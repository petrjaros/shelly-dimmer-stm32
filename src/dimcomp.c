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
 * measurement T_full per call (~2*n_half ticks). Same brightness_adj/Δt is
 * applied to both half-cycles within the upcoming full cycle.
 *
 * Model:
 *   f1(x) = A·sin(π·x/N),  f2(x) = B·sin(K_RATIO·π·x/N + π·a/8)
 *   δ_m   = -(-1)^m·C·sin(θ_m),   C = N·B/(πA)
 *   T_full_n = 2N + δ_{2n+2} − δ_{2n}
 *            ≈ 2N − K·cos(Θ_n + 0.3332π)
 * where Θ_n = θ_{2n}, K = 2·C·sin(0.3332π), per-full-cycle phase advance
 *   Θ_{n+1} − Θ_n = 0.6664π (mod 2π).
 *
 * Phase recovery from U_n and U_{n−1}:
 *   cos(Θ_n) = −(U_n + U_{n−1}) / (2K·cos(0.3332π))
 *   sin(Θ_n) = +(U_n − U_{n−1}) / (2K·sin(0.3332π))
 *
 * Δt for the upcoming cycle (sum of two half-cycle contributions, common Δt):
 *   Δt = (B/A)·N/sin²(πτ) · [P̃_eff·cos(Θ_M + π/6) − Q̃·sin(Θ_M + π/6)]
 *   P̃_eff = P̃ − sin²(πτ)/(2π)   (same as half-cycle)
 *   P̃, Q̃ are the standard cross-term shape functions.
 *
 * Since (U_n,U_{n−1}) yields Θ_{M−1} but we need Θ_M, the stored
 * coefficients absorb a rotation by γ = 2.5·(K_RATIO·π mod 2π) so the hot
 * path multiplies by cos(Θ_{M−1}), sin(Θ_{M−1}) as-is.
 */

#define K_RATIO     4.3332f
#define B_OVER_A    (1.0f / 41.07142f)

#define ALPHA_T_LOG2     5
#ifndef ALPHA_TH_LOG2
#define ALPHA_TH_LOG2    0   /* rotation-aware EMA on (sin,cos)(Θ); see hdr */
#endif
#ifndef VMAG_ON_NUM
#define VMAG_ON_NUM      5   /* on  threshold = NUM/10 · K                  */
#endif
#ifndef VMAG_OFF_NUM
#define VMAG_OFF_NUM     3   /* off threshold = NUM/10 · K  (hysteresis)    */
#endif
#define DT_LIMIT_NUM     3
#define DT_LIMIT_DEN     8
#define S2_FLOOR         0.02f
#define WARMUP_CYCLES    4

void dimcomp_init(dimcomp_t *d, uint16_t n_half_ticks)
{
    d->n_half = n_half_ticks;
    /* T_avg now tracks the FULL-cycle period, ~ 2·n_half_ticks. */
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

    /* Detection thresholds scaled to peak |U_n| = K. */
    d->v_thresh_on  = (int16_t)(K * (float)VMAG_ON_NUM  / 10.0f + 0.5f);
    d->v_thresh_off = (int16_t)(K * (float)VMAG_OFF_NUM / 10.0f + 0.5f);

    /* OC3 offset coefficients.
     * oc3 = δ_second − δ_first = 2·C·cos(α/2)·sin(Θ_n + 2.5α)
     *     = base·[cos(2.5α)·sinΘ + sin(2.5α)·cosΘ]
     * where α = K_RATIO·π mod 2π and C = N·(B/A)/π. */
    float C    = (float)n_half_ticks * B_OVER_A / pi;
    float base = 2.0f * C * cosf(0.5f * ha);          /* peak |oc3_offset| */
    float phs  = 2.5f * ha;
    d->oc3_sin_coef = (int32_t)(base * cosf(phs) * 256.0f + (base * cosf(phs) >= 0 ? 0.5f : -0.5f));
    d->oc3_cos_coef = (int32_t)(base * sinf(phs) * 256.0f + (base * sinf(phs) >= 0 ? 0.5f : -0.5f));

    d->t_prev_zc        = 0;
    d->v_prev           = 0;
    d->sin_th_q15       = 0;
    d->cos_th_q15       = 0;
    d->parity           = 0;  /* unused in full-cycle */
    d->warmup           = WARMUP_CYCLES;
    d->smoothing_seeded = 0;
    d->signal_active    = 0;
    d->u_recent[0]      = 0;
    d->u_recent[1]      = 0;
    d->u_recent[2]      = 0;
    d->u_idx            = 0;
    d->oc3_offset       = 0;

    d->t_dim_nominal = 0;
    d->sc_sin        = 0;
    d->sc_cos        = 0;
}

void dimcomp_set_brightness(dimcomp_t *d, uint16_t t_dim_nominal)
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

    float s  = sinf(pi_tau);
    float s2 = s * s;
    if (s2 < S2_FLOOR) s2 = S2_FLOOR;

    float km = K_RATIO - 1.0f;
    float kp = K_RATIO + 1.0f;

    float p_tilde = (cosf(km*pi_tau) - 1.0f)/(2.0f*km*pi) -
                    (cosf(kp*pi_tau) - 1.0f)/(2.0f*kp*pi);
    float q_tilde = sinf(km*pi_tau)/(2.0f*km*pi) -
                    sinf(kp*pi_tau)/(2.0f*kp*pi);
    float p_eff   = p_tilde - s2 / (2.0f * pi);

    /* Half-cycle "physical" coefficients (per half-cycle Δt, for M even):
     *   Δt_even = sc_sin_h · sin(θ_M) + sc_cos_h · cos(θ_M)
     * with sc_sin_h = -2·mul·P̃_eff,  sc_cos_h = -2·mul·Q̃,  mul = (B/A)·N/sin²(πτ).
     * For M odd, Δt carries an extra (-1)^M = -1.
     *
     * Coefficients are stored ×SC_SCALE (Q-format with extra fractional bits)
     * so small Δt values aren't lost to int truncation. The hot path shifts
     * by (15 + SC_SHIFT) instead of 15 to recover ticks.
     */
    float mul      = B_OVER_A * (float)d->n_half / s2;
    float sc_sin_h = -2.0f * mul * p_eff;
    float sc_cos_h = -2.0f * mul * q_tilde;

    /* We estimate (sin,cos)(Θ_{n-1}) where Θ = θ at the start of a FULL cycle.
     *   First  upcoming half-cycle:  M = 2n+2 → θ_M = Θ_{n-1} + 2α  (α = K_RATIO·π mod 2π)
     *   Second upcoming half-cycle:  M = 2n+3 → θ_M = Θ_{n-1} + 3α
     * Apply each rotation and the (-1)^M sign to bake everything into the
     * stored coefficients. */
    float alpha = fmodf(K_RATIO * pi, 2.0f * pi);       /* 0.3332π */
    float rot1  = fmodf(2.0f * alpha, 2.0f * pi);       /* 0.6664π */
    float rot2  = fmodf(3.0f * alpha, 2.0f * pi);       /* 0.9996π */
    float c1 = cosf(rot1), s1 = sinf(rot1);
    float c2 = cosf(rot2), s2r = sinf(rot2);

    /* For Δt = a·cos(θ+r) + b·sin(θ+r) expressed in terms of (sin θ, cos θ):
     *   new_cos = a·cos r + b·sin r
     *   new_sin = b·cos r − a·sin r           with (a,b) = (sc_cos_h, sc_sin_h) */
    /* Round-to-nearest at the Q-scale. */
    #define SC_SCALE 256.0f
    #define ROUND(x) ((int32_t)((x) >= 0 ? (x) + 0.5f : (x) - 0.5f))
    int32_t new_cos1 = ROUND(SC_SCALE * ( sc_cos_h * c1 + sc_sin_h * s1));
    int32_t new_sin1 = ROUND(SC_SCALE * ( sc_sin_h * c1 - sc_cos_h * s1));
    int32_t new_cos2 = ROUND(SC_SCALE * (-sc_cos_h * c2 - sc_sin_h * s2r));
    int32_t new_sin2 = ROUND(SC_SCALE * (-sc_sin_h * c2 + sc_cos_h * s2r));
    #undef ROUND
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

    /* 2. Sanity gate (now centred on ~2N). */
    int32_t T_avg = d->t_avg_q8 >> 8;
    if (T_n < (T_avg >> 1) || T_n > T_avg + (T_avg >> 1)) {
        d->t_dim_second = d->t_dim_nominal;
        d->oc3_offset = 0;
        return d->t_dim_nominal;
    }

    /* 3. U_n = T_n − T_avg, using T_avg BEFORE the EMA update so U_n isn't
     *    biased by α·T_n (which would shrink the recovered amplitude by
     *    1-α = 1-1/32 ≈ 3% under the previous "update first, subtract after"
     *    ordering). */
    int32_t U_n = T_n - T_avg;

    /* Update the running mean. */
    d->t_avg_q8 += ((T_n << 8) - d->t_avg_q8) >> ALPHA_T_LOG2;

    /* 4b. Update the 3-cycle |U_n| ring. The signal is a cosine sampled at
     *     120° per full mains cycle, so any 3 consecutive |U_n| values must
     *     include at least one ≥ K·√3/2 ≈ 0.866·K regardless of starting
     *     phase. A simple max-of-3 with hysteresis is therefore guaranteed
     *     to catch any real burst within 1-3 cycles AND release within 3
     *     cycles of burst end — no decay tuning needed. */
    int32_t u_abs = U_n < 0 ? -U_n : U_n;
    d->u_recent[d->u_idx] = (int16_t)(u_abs > 32767 ? 32767 : u_abs);
    d->u_idx = (uint8_t)((d->u_idx + 1) % 3);

    /* Warmup. */
    if (d->warmup) {
        d->warmup--;
        d->v_prev = U_n;
        d->t_dim_second = d->t_dim_nominal;
        d->oc3_offset = 0;
        return d->t_dim_nominal;
    }

    /* 4c. Max-of-3 hysteretic signal-presence detector. */
    int16_t u_max = d->u_recent[0];
    if (d->u_recent[1] > u_max) u_max = d->u_recent[1];
    if (d->u_recent[2] > u_max) u_max = d->u_recent[2];
    if (d->signal_active) {
        if (u_max < d->v_thresh_off) d->signal_active = 0;
    } else {
        if (u_max > d->v_thresh_on) {
            d->signal_active    = 1;
            d->smoothing_seeded = 0;
        }
    }

    if (!d->signal_active) {
        d->v_prev = U_n;
        d->t_dim_second = d->t_dim_nominal;
        d->oc3_offset = 0;
        return d->t_dim_nominal;
    }

    /* 5. Recover cos(Θ_{M-1}), sin(Θ_{M-1}) in Q15. */
    int32_t sum_u  = U_n + d->v_prev;
    int32_t diff_u = U_n - d->v_prev;
    d->v_prev = U_n;

    int16_t cos_new = clamp_q15(-sum_u  * d->scale_cos);
    int16_t sin_new = clamp_q15( diff_u * d->scale_sin);

    /* 6. Rotation-aware smoothing (rotation by 2·K_RATIO·π per cycle).
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

    /* 7. Compute Δt for BOTH upcoming half-cycles. The coefficients carry the
     *    rotation that maps the estimated phasor at Θ_{n-1} forward to the
     *    correct θ_M (2α for the first half, 3α for the second) plus the
     *    (-1)^M sign for the second.
     *
     * Coefficients are stored ×256 (Q8), phasor in Q15, so shift back by 23
     * with round-to-nearest. Asymmetric truncation would bias negative Δt
     * values one tick more negative than positive ones. */
    int32_t mix1 = d->sc_cos  * (int32_t)d->cos_th_q15
                 + d->sc_sin  * (int32_t)d->sin_th_q15;
    int32_t mix2 = d->sc_cos2 * (int32_t)d->cos_th_q15
                 + d->sc_sin2 * (int32_t)d->sin_th_q15;
    int32_t dt1 = (mix1 + (1 << 22)) >> 23;
    int32_t dt2 = (mix2 + (1 << 22)) >> 23;

    /* 7b. OC3 offset: shift the mid-cycle event from N to N+(δ₂−δ₁) so the
     *     dimmer fires at the predicted ACTUAL second-half ZC, not just N
     *     ticks after the first ZC. Same Q8 coefficient/Q15 phasor scheme. */
    int32_t mix3 = d->oc3_sin_coef * (int32_t)d->sin_th_q15
                 + d->oc3_cos_coef * (int32_t)d->cos_th_q15;
    int32_t oc3 = (mix3 + (1 << 22)) >> 23;
    if (oc3 >  127) oc3 =  127;   /* sanity clamp — peak |offset| ≈ K ≈ 15 */
    if (oc3 < -128) oc3 = -128;
    d->oc3_offset = (int16_t)oc3;

    /* 8. Clamp and apply (independently for each half-cycle). */
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
