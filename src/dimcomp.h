/*
 * dimcomp - Per-cycle compensation of a phase-cut dimmer against an HF data
 *           signal that rides on the mains.
 *
 * Model:
 *   mains : f1(x) = A * sin(pi * x / N)
 *   signal: f2(x) = B * sin(K_RATIO * pi * x / N  +  pi*a/8)
 *   v(x)  = f1(x) + f2(x)
 *
 * N    = timer ticks per mains HALF-cycle
 * B/A  = signal-to-mains amplitude ratio (~1/36.5 measured)
 * a    = slow phase drift we estimate from inter-ZC timing
 *
 * Sampling model is FULL-cycle: the host fires the dimmer once per mains
 * cycle on the rising-edge ZC and synthesises the mid-cycle event from a
 * timer compare. dimcomp_on_zc receives one timestamp per full cycle
 * (T_full ~ 2N) and returns *two* corrected off-times — one per half-cycle —
 * because the two half-cycles of the upcoming mains cycle see different
 * f2 phases (0.3332*pi apart).
 *
 * Usage:
 *     static dimcomp_t dc;
 *
 *     // Once at startup:
 *     dimcomp_init(&dc, 1000);              // ticks/half-cycle
 *
 *     // When user changes brightness (off the ISR — uses sin/cos):
 *     dimcomp_set_brightness(&dc, 100, 0);  // t_dim, window-start offset
 *
 *     // On every rising-edge ZC:
 *     uint16_t t1 = dimcomp_on_zc(&dc, now_ticks);   // first half-cycle
 *     uint16_t t2 = dc.t_dim_second;                 // second half-cycle
 *     schedule_off_at(t1);
 *     schedule_off_at(N_HALF + t2);                  // OC3-driven event
 *
 * When the data signal is absent, both t1 and t2 equal t_dim_nominal.
 *
 * Thread/IRQ safety: dimcomp_set_brightness writes the sc_* coefficients
 * and t_dim_nominal LAST, so a tolerance-check sync gate in the ISR can
 * compare against t_dim_nominal to know whether the coefficients are live.
 * The on_zc hot path is integer-only and ISR-safe.
 */

#ifndef DIMCOMP_H
#define DIMCOMP_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    /* --- configuration (set by init) --- */
    uint16_t n_half;          /* nominal timer ticks per mains half-cycle  */
    int32_t  scale_sin;       /* Q15-per-tick: sin_q15 = +(U_n-U_{n-1})*scale_sin */
    int32_t  scale_cos;       /* Q15-per-tick: cos_q15 = -(U_n+U_{n-1})*scale_cos */
    int16_t  rot_cos_q15;     /* cos/sin(2*K_RATIO*pi mod 2pi), Q15:       */
    int16_t  rot_sin_q15;     /*   per-FULL-cycle phasor rotation of theta */

    /* --- recomputed on brightness change --- */
    uint16_t t_dim_nominal;   /* what the user asked for                   */
    int32_t  sc_sin;          /* FIRST  half-cycle Δt coefficients         */
    int32_t  sc_cos;          /*   Δt = (sc_sin*sin + sc_cos*cos) >> 15    */
    int32_t  sc_sin2;         /* SECOND half-cycle Δt coefficients         */
    int32_t  sc_cos2;
    uint16_t t_dim_second;    /* corrected off-time for 2nd half-cycle     */
    int32_t  oc3_sin_coef;    /* coef of sin(Θ_n) for OC3 offset, Q8 ticks */
    int32_t  oc3_cos_coef;    /* coef of cos(Θ_n) for OC3 offset, Q8 ticks */
    int16_t  oc3_offset;      /* signed ticks to add to N for OC3 timing.  */
                              /* Predicts (δ_second − δ_first) so the mid-  */
                              /* cycle event fires at the ACTUAL 2nd ZC.    */

    /* --- runtime state --- */
    uint32_t t_prev_zc;       /* last zero-crossing timestamp (free-run)   */
    int32_t  t_avg_q8;        /* running mean FULL-cycle period, Q8 ticks  */
    int32_t  v_prev;          /* previous U_n                              */
    int16_t  sin_th_q15;      /* smoothed sin(theta_n), Q15                */
    int16_t  cos_th_q15;      /* smoothed cos(theta_n), Q15                */
    uint8_t  warmup;          /* cycles remaining before correction kicks in */
    uint8_t  smoothing_seeded;/* 0 until phasor EMA has been initialized   */
    uint8_t  signal_active;   /* 1 = signal detected, correction enabled.  */
                              /* Hysteresis on the recovered phasor's      */
                              /* magnitude (pair of U_n samples => tone    */
                              /* amplitude, phase-independent).            */
} dimcomp_t;

/* Initialize. n_half_ticks = nominal timer ticks per mains half-cycle. */
void     dimcomp_init(dimcomp_t *d, uint16_t n_half_ticks);

/* Call when the brightness setting OR the window-start offset changes (uses
 * floats; not for ISR). win_start_offset = constant ticks between a ZC
 * crossing and the start of its conduction window (detector-asymmetry
 * centering + guard offsets the host applies). The 216 Hz signal advances
 * ~0.08 rad per 100 ticks, so an unaccounted offset of a few hundred ticks
 * visibly misphases the correction. */
void     dimcomp_set_brightness(dimcomp_t *d, uint16_t t_dim_nominal,
                                uint16_t win_start_offset);

/* Call from the zero-crossing ISR. Returns the corrected t_dim (ticks) to
 * use for the half-cycle that just began. Integer-only, no transcendentals. */
uint16_t dimcomp_on_zc(dimcomp_t *d, uint32_t timestamp_ticks);

#ifdef __cplusplus
}
#endif

#endif /* DIMCOMP_H */
