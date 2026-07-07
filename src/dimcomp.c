#include "dimcomp.h"
#define _USE_MATH_DEFINES
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- Full-cycle PLC-signal flicker compensation. ----
 *
 * The host fires the dimmer once per FULL mains cycle (on the rising-edge
 * ZC) and synthesises the mid-cycle event from a timer compare. We see one
 * measurement T_full per call (~2*n_half ticks) and return a separate
 * corrected off-time for each of the two upcoming half-cycles.
 *
 * Model:
 *   f1(x) = A·sin(π·x/N),  f2(x) = B·sin(K_RATIO·π·x/N + π·a/8)
 *   δ_m   = -(-1)^m·C·sin(θ_m),   C = N·B/(πA)
 *   T_full_n = 2N + δ_{2n+2} − δ_{2n}
 *            ≈ 2N − K·cos(Θ_n + 0.3332π)
 * where Θ_n = θ_{2n}, K = 2·C·sin(0.3332π), per-full-cycle phase advance
 *   Θ_{n+1} − Θ_n = 0.6664π (mod 2π).
 *
 * Phase recovery from U_n and U_{n−1} (U = T_full − running mean):
 *   cos(Θ_n) = −(U_n + U_{n−1}) / (2K·cos(0.3332π))
 *   sin(Θ_n) = +(U_n − U_{n−1}) / (2K·sin(0.3332π))
 *
 * Each half-cycle's Δt is a linear mix of that phasor,
 *   Δt = sc_cos·cos(Θ) + sc_sin·sin(Θ),
 * with coefficients precomputed in dimcomp_set_brightness (see there for
 * the energy-balance derivation). Since (U_n,U_{n−1}) yields the phasor one
 * cycle in the past, the stored coefficients also absorb the rotation that
 * advances it to each upcoming crossing, so the hot path stays a plain
 * multiply-accumulate.
 */

#define K_RATIO     4.3332f
#define B_OVER_A    (1.0f / 41.07142f)

#define ALPHA_T_LOG2     5
#ifndef ALPHA_TH_LOG2
#define ALPHA_TH_LOG2    0   /* rotation-aware EMA on (sin,cos)(Θ); see hdr */
#endif
#ifndef VMAG_ON_NUM
#define VMAG_ON_NUM      3   /* on  threshold = NUM/10 of full tone amplitude */
#endif
#ifndef VMAG_OFF_NUM
#define VMAG_OFF_NUM     2   /* off threshold = NUM/10 (hysteresis)           */
#endif
#define DT_LIMIT_NUM     3
#define DT_LIMIT_DEN     8

/* Pair-magnitude detection thresholds: (fraction of full tone amplitude)²
 * in the Q30 units of mag² = sin²+cos² of the recovered phasor. */
#define MAG2_ON  ((uint32_t)((32768.0f * VMAG_ON_NUM  / 10) * \
                             (32768.0f * VMAG_ON_NUM  / 10)))
#define MAG2_OFF ((uint32_t)((32768.0f * VMAG_OFF_NUM / 10) * \
                             (32768.0f * VMAG_OFF_NUM / 10)))
#define S2_FLOOR         0.02f
#define WARMUP_CYCLES    4

/* Fractional bits of the stored sc_* coefficients. The hot-path mix
 *   sc * phasor(Q15)  must fit in int32, so the largest representable
 * correction is 2^31 >> (15+SC_SHIFT) ticks. Don't raise this: SC_SHIFT=8
 * caps |Δt| at 256 ticks, but near-full brightness (small sin² at the
 * window end => big mul) needs ~270, and the mix then WRAPS, flipping the
 * correction sign on the worst-phase cycles - seen as a +1% energy spike
 * every 3rd cycle. SC_SHIFT=2 keeps quarter-tick coefficient resolution
 * (quantisation ~0.1 tick on typical Δt, invisible) and allows |Δt| up to
 * 16384. */
#define SC_SHIFT         2

/* Round float to nearest int32. */
#define ROUND_I32(x) ((int32_t)((x) >= 0 ? (x) + 0.5f : (x) - 0.5f))

void dimcomp_init(dimcomp_t *d, uint16_t n_half_ticks)
{
    d->n_half = n_half_ticks;
    /* T_avg tracks the FULL-cycle period, ~ 2·n_half_ticks. */
    d->t_avg_q8 = (int32_t)((uint32_t)n_half_ticks << 9);  /* 2*N << 8 */

    float n   = (float)n_half_ticks;
    float pi  = (float)M_PI;
    float ha  = fmodf(K_RATIO * pi, 2.0f * pi);   /* per-half-cycle ≈ 0.3332π */
    float s33 = sinf(ha);
    float c33 = cosf(ha);
    float K   = 2.0f * n * B_OVER_A * s33 / pi;   /* peak |U_n|              */
    float Kc  = 2.0f * K * c33;                   /* denom for cos(Θ) recov  */
    float Ks  = 2.0f * K * s33;                   /* denom for sin(Θ) recov  */

    /* Q15 scales for phase recovery.
     *   cos_q15 = -(U_n + U_{n-1}) * scale_cos
     *   sin_q15 = +(U_n - U_{n-1}) * scale_sin */
    d->scale_cos = (int32_t)(32768.0f / Kc + 0.5f);
    d->scale_sin = (int32_t)(32768.0f / Ks + 0.5f);

    /* Per-cycle rotation of the (sin,cos)(Θ) phasor: 2·K_RATIO·π mod 2π. */
    float rot = fmodf(2.0f * K_RATIO * pi, 2.0f * pi);
    d->rot_cos_q15 = (int16_t)(cosf(rot) * 32767.0f);
    d->rot_sin_q15 = (int16_t)(sinf(rot) * 32767.0f);

    /* OC3 offset coefficients.
     * oc3 = δ_second − δ_first = 2·C·cos(α/2)·sin(Θ_n + 2.5α)
     *     = base·[cos(2.5α)·sinΘ + sin(2.5α)·cosΘ]
     * where α = K_RATIO·π mod 2π and C = N·(B/A)/π. */
    float C    = (float)n_half_ticks * B_OVER_A / pi;
    float base = 2.0f * C * cosf(0.5f * ha);          /* peak |oc3_offset| */
    float phs  = 2.5f * ha;
    d->oc3_sin_coef = ROUND_I32(base * cosf(phs) * 256.0f);
    d->oc3_cos_coef = ROUND_I32(base * sinf(phs) * 256.0f);

    d->t_prev_zc        = 0;
    d->v_prev           = 0;
    d->sin_th_q15       = 0;
    d->cos_th_q15       = 0;
    d->warmup           = WARMUP_CYCLES;
    d->smoothing_seeded = 0;
    d->signal_active    = 0;
    d->oc3_offset       = 0;

    d->t_dim_nominal = 0;
    d->t_dim_second  = 0;
    d->sc_sin        = 0;
    d->sc_cos        = 0;
    d->sc_sin2       = 0;
    d->sc_cos2       = 0;
}

void dimcomp_set_brightness(dimcomp_t *d, uint16_t t_dim_nominal,
                            uint16_t win_start_offset)
{
    if (t_dim_nominal == 0 || t_dim_nominal >= d->n_half) {
        d->sc_sin = 0;
        d->sc_cos = 0;
        d->sc_sin2 = 0;
        d->sc_cos2 = 0;
        d->t_dim_nominal = t_dim_nominal;
        return;
    }

    float pi     = (float)M_PI;
    float tau    = (float)t_dim_nominal / (float)d->n_half;
    float pi_tau = pi * tau;

    float km = K_RATIO - 1.0f;
    float kp = K_RATIO + 1.0f;

    float p_tilde = (cosf(km*pi_tau) - 1.0f)/(2.0f*km*pi) -
                    (cosf(kp*pi_tau) - 1.0f)/(2.0f*kp*pi);
    float q_tilde = sinf(km*pi_tau)/(2.0f*km*pi) -
                    sinf(kp*pi_tau)/(2.0f*kp*pi);
    /* Companion shape functions for the cos(πy/N)-weighted part of the
     * window (nonzero once the window starts W > 0 ticks past the ZC):
     *   P̃c = ∫cos·cos / N,  Q̃c = ∫cos·sin(K·) / N  over [0, τ]. */
    float pc_tilde = (sinf(kp*pi_tau)/kp + sinf(km*pi_tau)/km) / (2.0f*pi);
    float qc_tilde = ((1.0f - cosf(kp*pi_tau))/kp +
                      (1.0f - cosf(km*pi_tau))/km) / (2.0f*pi);

    /* Two physically distinct contributions per half-cycle window (even M),
     * for a window [W, W+T] starting W ticks after its ZC crossing
     * (W = base offset the host adds: detector-asymmetry centering + guard):
     *
     *   cross term : Δt_x = -2·mul·[(cosW̄·P̃ + sinW̄·P̃c)·sin(θ_win)
     *                            + (cosW̄·Q̃ + sinW̄·Q̃c)·cos(θ_win)]
     *                W̄ = π·W/N; θ_win = crossing phase advanced by
     *                γ = K_RATIO·π·W/N; mul = (B/A)·N / sin²(π(W+T)/N)
     *                (the energy/Δt slope is taken at the window END).
     *   boundary   : Δt_b = +C·(1 - sin²W̄/sin²(π(W+T)/N))·sin(θ_cross),
     *                C = N·(B/A)/π — cancels the jitter of the window
     *                start, which rides the CROSSING displacement δ
     *                regardless of the constant offset W, so it keeps the
     *                un-advanced crossing phase. The (1 - ...) factor is
     *                the f1²(start)/f1²(end) leak of a start shift.
     *
     * (For W = 0 the two would fold into the single shape function
     * P̃_eff = P̃ - sin²(πτ)/(2π).) For M odd both terms flip sign.
     *
     * Coefficients are stored ×SC_SCALE (Q-format with extra fractional
     * bits) so small Δt values aren't lost to int truncation. The hot path
     * shifts by (15 + SC_SHIFT) instead of 15 to recover ticks.
     */
    float wbar = pi * (float)win_start_offset / (float)d->n_half;
    float cw   = cosf(wbar);
    float sw   = sinf(wbar);
    float se   = sinf(pi_tau + wbar);     /* f1 envelope at window END   */
    float s2e  = se * se;
    if (s2e < S2_FLOOR) s2e = S2_FLOOR;

    float mul    = B_OVER_A * (float)d->n_half / s2e;
    float cr_sin = -2.0f * mul * (cw * p_tilde + sw * pc_tilde);
    float cr_cos = -2.0f * mul * (cw * q_tilde + sw * qc_tilde);
    float cb     = B_OVER_A * (float)d->n_half / pi
                 * (1.0f - sw * sw / s2e);

    /* We estimate (sin,cos)(Θ_{n-1}) where Θ = θ at the start of a FULL cycle.
     *   First  upcoming crossing:  M = 2n+2 → θ_M = Θ_{n-1} + 2α  (α = K_RATIO·π mod 2π)
     *   Second upcoming crossing:  M = 2n+3 → θ_M = Θ_{n-1} + 3α
     * Cross terms get the extra window-start advance γ; boundary terms don't.
     * Apply each rotation and the (-1)^M sign to bake everything into the
     * stored coefficients. */
    float alpha = fmodf(K_RATIO * pi, 2.0f * pi);       /* 0.3332π */
    float gamma = K_RATIO * pi * (float)win_start_offset / (float)d->n_half;
    float r1b = fmodf(2.0f * alpha, 2.0f * pi);
    float r2b = fmodf(3.0f * alpha, 2.0f * pi);
    float r1x = fmodf(r1b + gamma, 2.0f * pi);
    float r2x = fmodf(r2b + gamma, 2.0f * pi);

    /* For Δt = a·cos(θ+r) + b·sin(θ+r) expressed in terms of (sin θ, cos θ):
     *   new_cos = a·cos r + b·sin r
     *   new_sin = b·cos r − a·sin r            with (a,b) = (cos-coef, sin-coef).
     * Boundary contribution C·sin(θ+r) is the (a,b) = (0, C) case. */
    /* Round-to-nearest at the Q-scale. */
    #define SC_SCALE ((float)(1 << SC_SHIFT))
    int32_t new_cos1 = ROUND_I32(SC_SCALE * ( cr_cos * cosf(r1x) + cr_sin * sinf(r1x)
                                             + cb * sinf(r1b)));
    int32_t new_sin1 = ROUND_I32(SC_SCALE * ( cr_sin * cosf(r1x) - cr_cos * sinf(r1x)
                                             + cb * cosf(r1b)));
    int32_t new_cos2 = ROUND_I32(SC_SCALE * (-cr_cos * cosf(r2x) - cr_sin * sinf(r2x)
                                             - cb * sinf(r2b)));
    int32_t new_sin2 = ROUND_I32(SC_SCALE * (-cr_sin * cosf(r2x) + cr_cos * sinf(r2x)
                                             - cb * cosf(r2b)));
    #undef SC_SCALE

    d->sc_sin  = new_sin1;
    d->sc_cos  = new_cos1;
    d->sc_sin2 = new_sin2;
    d->sc_cos2 = new_cos2;
    d->t_dim_nominal = t_dim_nominal;
}

static inline int16_t clamp_q15(int32_t x)
{
    if (x >  32767) return  32767;
    if (x < -32768) return -32768;
    return (int16_t)x;
}

uint16_t dimcomp_on_zc(dimcomp_t *d, uint32_t timestamp_ticks)
{
    /* 1. Full-cycle period since last call. */
    int32_t T_n = (int32_t)(timestamp_ticks - d->t_prev_zc);
    d->t_prev_zc = timestamp_ticks;

    /* 2. Sanity gate, centred on ~2N. */
    int32_t T_avg = d->t_avg_q8 >> 8;
    if (T_n < (T_avg >> 1) || T_n > T_avg + (T_avg >> 1)) {
        /* Gap or glitch in the call stream — e.g. the host stopped calling
         * between signal bursts. v_prev and the phasor describe the signal
         * from before the gap, and correcting from that stale phase blinks
         * the light at the start of every burst, so fall back to fresh
         * detection. warmup=1 spends one call re-seeding v_prev before the
         * phase recovery is trusted again. */
        d->warmup           = 1;
        d->signal_active    = 0;
        d->smoothing_seeded = 0;
        d->t_dim_second = d->t_dim_nominal;
        d->oc3_offset = 0;
        return d->t_dim_nominal;
    }

    /* 3. U_n = T_n − T_avg, using T_avg from BEFORE the EMA update:
     *    updating first would fold α·T_n into the mean and shrink the
     *    recovered amplitude by 1-α = 1-1/32 ≈ 3%. */
    int32_t U_n = T_n - T_avg;

    /* Update the running mean. */
    d->t_avg_q8 += ((T_n << 8) - d->t_avg_q8) >> ALPHA_T_LOG2;

    /* 4. Warmup. */
    if (d->warmup) {
        d->warmup--;
        d->v_prev = U_n;
        d->t_dim_second = d->t_dim_nominal;
        d->oc3_offset = 0;
        return d->t_dim_nominal;
    }

    /* 5. Recover cos(Θ_{M-1}), sin(Θ_{M-1}) in Q15. Runs before detection:
     *    two consecutive samples of the 120°-per-cycle beat cosine determine
     *    amplitude AND phase, and detection keys on the amplitude. */
    int32_t sum_u  = U_n + d->v_prev;
    int32_t diff_u = U_n - d->v_prev;
    d->v_prev = U_n;

    int16_t cos_new = clamp_q15(-sum_u  * d->scale_cos);
    int16_t sin_new = clamp_q15( diff_u * d->scale_sin);

    /* 6. Pair-magnitude hysteretic signal-presence detector.
     *    mag² = sin²+cos² ≈ (amplitude/K)² in Q30, phase-INDEPENDENT for a
     *    steady tone. (Thresholding |U_n| alone would sample cos Θ, which
     *    can sit near a beat zero for the first 1-2 burst cycles, delaying
     *    detection by up to 2 full cycles — a visible flash at HDO burst
     *    onset.) A pair straddling the onset only shrinks mag² (partial
     *    sample), and the unnormalized phasor then also shrinks the applied
     *    correction by the same factor, so engaging early degrades
     *    gracefully. Release: two clean samples drop mag² below MAG2_OFF
     *    within 2 cycles of burst end. */
    uint32_t mag2 = (uint32_t)((int32_t)cos_new * cos_new)
                  + (uint32_t)((int32_t)sin_new * sin_new);
    if (d->signal_active) {
        if (mag2 < MAG2_OFF) d->signal_active = 0;
    } else {
        if (mag2 > MAG2_ON) {
            d->signal_active    = 1;
            d->smoothing_seeded = 0;
        }
    }

    if (!d->signal_active) {
        d->t_dim_second = d->t_dim_nominal;
        d->oc3_offset = 0;
        return d->t_dim_nominal;
    }

    /* 7. Rotation-aware smoothing (rotation by 2·K_RATIO·π per cycle).
     *    With ALPHA_TH_LOG2=0 this degenerates to "use latest sample". */
    if (!d->smoothing_seeded) {
        d->sin_th_q15      = sin_new;
        d->cos_th_q15      = cos_new;
        d->smoothing_seeded = 1;
    } else {
        int32_t sin_pred = ((int32_t)d->rot_cos_q15 * d->sin_th_q15
                          + (int32_t)d->rot_sin_q15 * d->cos_th_q15) >> 15;
        int32_t cos_pred = ((int32_t)d->rot_cos_q15 * d->cos_th_q15
                          - (int32_t)d->rot_sin_q15 * d->sin_th_q15) >> 15;
        d->sin_th_q15 = (int16_t)(sin_pred + (((int32_t)sin_new - sin_pred) >> ALPHA_TH_LOG2));
        d->cos_th_q15 = (int16_t)(cos_pred + (((int32_t)cos_new - cos_pred) >> ALPHA_TH_LOG2));
    }

    /* 8. Compute Δt for BOTH upcoming half-cycles. The coefficients carry the
     *    rotation that maps the estimated phasor at Θ_{n-1} forward to the
     *    correct θ_M (2α for the first half, 3α for the second) plus the
     *    (-1)^M sign for the second.
     *
     * Coefficients are stored ×2^SC_SHIFT, phasor in Q15, so shift back by
     * (15+SC_SHIFT) with round-to-nearest. Asymmetric truncation would bias
     * negative Δt values one tick more negative than positive ones. The mix
     * must fit int32 - see the SC_SHIFT note for the resulting |Δt| bound. */
    int32_t mix1 = d->sc_cos  * (int32_t)d->cos_th_q15
                 + d->sc_sin  * (int32_t)d->sin_th_q15;
    int32_t mix2 = d->sc_cos2 * (int32_t)d->cos_th_q15
                 + d->sc_sin2 * (int32_t)d->sin_th_q15;
    int32_t dt1 = (mix1 + (1 << (14 + SC_SHIFT))) >> (15 + SC_SHIFT);
    int32_t dt2 = (mix2 + (1 << (14 + SC_SHIFT))) >> (15 + SC_SHIFT);

    /* 9. OC3 offset: shift the mid-cycle event from N to N+(δ₂−δ₁) so the
     *    dimmer fires at the predicted ACTUAL second-half ZC, not just N
     *    ticks after the first ZC. Same Q8 coefficient/Q15 phasor scheme. */
    int32_t mix3 = d->oc3_sin_coef * (int32_t)d->sin_th_q15
                 + d->oc3_cos_coef * (int32_t)d->cos_th_q15;
    int32_t oc3 = (mix3 + (1 << 22)) >> 23;
    if (oc3 >  300) oc3 =  300;   /* sanity clamp — peak |offset| = 2C·cos(α/2) ≈ 134 */
    if (oc3 < -300) oc3 = -300;
    d->oc3_offset = (int16_t)oc3;

    /* 10. Clamp and apply (independently for each half-cycle). */
    int32_t lim = ((int32_t)d->t_dim_nominal * DT_LIMIT_NUM) / DT_LIMIT_DEN;
    if (dt1 >  lim) dt1 =  lim;
    if (dt1 < -lim) dt1 = -lim;
    if (dt2 >  lim) dt2 =  lim;
    if (dt2 < -lim) dt2 = -lim;

    int32_t t1 = (int32_t)d->t_dim_nominal + dt1;
    int32_t t2 = (int32_t)d->t_dim_nominal + dt2;
    if (t1 < 0)                  t1 = 0;
    if (t1 > (int32_t)d->n_half) t1 = (int32_t)d->n_half;
    if (t2 < 0)                  t2 = 0;
    if (t2 > (int32_t)d->n_half) t2 = (int32_t)d->n_half;

    d->t_dim_second = (uint16_t)t2;
    return (uint16_t)t1;
}
