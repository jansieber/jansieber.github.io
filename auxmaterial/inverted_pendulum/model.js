/* =====================================================================
   model.js -- mathematical core of the delayed cart-pendulum animation

   Written in a deliberately C-like subset of JavaScript, so that it
   reads like C:
     - all numbers are IEEE doubles (like C `double`);
     - `new Float64Array(n)` is a fixed-length array `double a[n]`,
       zero-initialised;
     - `var` declares a variable, `function f(a, b)` a function
       (no type annotations); arrays are passed by reference, as in C;
     - Math.sin, Math.cos, Math.floor, Math.sqrt, Math.abs, Math.atan2,
       Math.PI are the usual <math.h> functions.
   Nothing in this file touches the screen; the drawing code in the
   HTML file only calls sim_reset(), sim_step() and critical_delay().

   Model (x = pendulum angle from upright, y = cart position):
     (M+m) y'' + m l (x'' cos x - x'^2 sin x) = F(t)
     (4/3) l x'' + y'' cos x - g sin x       = 0,     l = L/2
     F(t) = Kp x(t-tau) + Kd x'(t-tau) + kp y(t-tau) + kd y'(t-tau)
   ===================================================================== */

/* ---- physical parameters (SI units) ---- */
var M = 1.0;          /* cart mass [kg]                           */
var m = 0.3;          /* pendulum (homogeneous rod) mass [kg]     */
var L = 1.0;          /* rod length [m]                           */
var l = 0.5 * L;      /* distance pivot -> centre of mass [m]     */
var g = 9.81;         /* gravity [m/s^2]                          */

/* ---- feedback gains: gain[0..3] = Kp, Kd, kp, kd ---- */
var gain = new Float64Array(4);
gain[0] = 20.0;  gain[1] = 7.0;  gain[2] = 0.3;  gain[3] = 1.5;

/* ---- delay [s] ---- */
var tau = 0.178;

/* ---- numerical parameters ---- */
var DT     = 0.0005;            /* RK4 time step [s]                     */
var TAUMAX = 0.30;              /* largest delay the history can hold    */
var NHIST  = 610;               /* history length, > TAUMAX/DT + 2       */
var NS     = 4;                 /* state dimension: z = (x, x', y, y')   */

/* ---- simulation state ---- */
var z     = new Float64Array(NS);          /* current state            */
var hist  = new Float64Array(NHIST * NS);  /* ring buffer of past z    */
var head  = 0;                             /* row of hist holding z(t) */
var t     = 0.0;                           /* current time             */
var force = 0.0;                           /* last applied force F     */

/* ---------------------------------------------------------------------
   rhs: right-hand side of the first-order system z' = f(z, F).
   The two equations of motion are linear in (x'', y''):
       [ M+m      m l cos x ] [y'']   [ F + m l x'^2 sin x ]
       [ cos x    (4/3) l   ] [x''] = [ g sin x            ]
   which is solved with Cramer's rule.
   --------------------------------------------------------------------- */
function rhs(zin, F, dz)
{
    var x  = zin[0], xd = zin[1];
    var s  = Math.sin(x), c = Math.cos(x);
    var a11 = M + m,  a12 = m * l * c;
    var a21 = c,      a22 = 4.0 * l / 3.0;
    var r1  = F + m * l * xd * xd * s;
    var r2  = g * s;
    var det = a11 * a22 - a12 * a21;
    var ydd = (r1 * a22 - a12 * r2) / det;
    var xdd = (a11 * r2 - a21 * r1) / det;
    dz[0] = xd;      /* x'  */
    dz[1] = xdd;     /* x'' */
    dz[2] = zin[3];  /* y'  */
    dz[3] = ydd;     /* y'' */
}

/* ---------------------------------------------------------------------
   past_state: state at time t - lag (0 <= lag <= TAUMAX), written to out.
   Linear interpolation between the stored grid values z(t - j*DT).
   --------------------------------------------------------------------- */
function past_state(lag, out)
{
    var i, j, theta, r0, r1;
    if (lag <= 0.0) {                 /* only happens if tau < DT */
        for (i = 0; i < NS; i++) out[i] = z[i];
        return;
    }
    j     = Math.floor(lag / DT);
    theta = lag / DT - j;
    r0 = ((head - j)     % NHIST + NHIST) % NHIST;   /* row of z(t - j DT)     */
    r1 = ((head - j - 1) % NHIST + NHIST) % NHIST;   /* row of z(t - (j+1) DT) */
    for (i = 0; i < NS; i++)
        out[i] = (1.0 - theta) * hist[r0 * NS + i] + theta * hist[r1 * NS + i];
}

/* feedback law: F = Kp x + Kd x' + kp y + kd y'  evaluated at a delayed state */
function control(zd)
{
    return gain[0] * zd[0] + gain[1] * zd[1] + gain[2] * zd[2] + gain[3] * zd[3];
}

/* work arrays for sim_step (allocated once, like static arrays in C) */
var k1 = new Float64Array(NS), k2 = new Float64Array(NS);
var k3 = new Float64Array(NS), k4 = new Float64Array(NS);
var ztmp = new Float64Array(NS), zdel = new Float64Array(NS);

/* ---------------------------------------------------------------------
   sim_step: advance z from t to t+DT with the classical RK4 method.
   The delayed force is evaluated at the stage times t, t+DT/2, t+DT,
   i.e. from the history at t-tau, t+DT/2-tau, t+DT-tau.
   Returns 1 if the pendulum has fallen over (|x| > pi/2), else 0.
   --------------------------------------------------------------------- */
function sim_step()
{
    var i, F0, Fh, F1;

    past_state(tau,            zdel);  F0 = control(zdel);
    past_state(tau - 0.5 * DT, zdel);  Fh = control(zdel);
    past_state(tau - DT,       zdel);  F1 = control(zdel);

    rhs(z, F0, k1);
    for (i = 0; i < NS; i++) ztmp[i] = z[i] + 0.5 * DT * k1[i];
    rhs(ztmp, Fh, k2);
    for (i = 0; i < NS; i++) ztmp[i] = z[i] + 0.5 * DT * k2[i];
    rhs(ztmp, Fh, k3);
    for (i = 0; i < NS; i++) ztmp[i] = z[i] + DT * k3[i];
    rhs(ztmp, F1, k4);
    for (i = 0; i < NS; i++)
        z[i] += DT / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);

    /* store new state in the ring buffer */
    head = (head + 1) % NHIST;
    for (i = 0; i < NS; i++) hist[head * NS + i] = z[i];
    t += DT;
    force = F0;

    return (Math.abs(z[0]) > 0.5 * Math.PI) ? 1 : 0;
}

/* sim_reset: constant initial history z(s) = (x0, 0, 0, 0) for s <= 0 */
function sim_reset(x0)
{
    var r, i;
    z[0] = x0;  z[1] = 0.0;  z[2] = 0.0;  z[3] = 0.0;
    for (r = 0; r < NHIST; r++)
        for (i = 0; i < NS; i++) hist[r * NS + i] = z[i];
    head = 0;  t = 0.0;  force = 0.0;
}

/* =====================================================================
   Linear stability of the upright equilibrium x = y = 0.
   Linearisation gives the characteristic equation
      P(lambda) = lambda^4 - a lambda^2 + exp(-lambda tau) Q(lambda) = 0,
      Q(lambda) = [ (K(lambda) - (4/3) l k(lambda)) lambda^2 + g k(lambda) ] / D,
   with K(lambda) = Kp + Kd lambda, k(lambda) = kp + kd lambda,
        D = l (4M+m)/3,  a = (M+m) g / D.
   A root lambda = i w lies on the imaginary axis iff
      |Q(i w)| = w^4 + a w^2         (modulus condition)
      exp(-i w tau) = -(w^4 + a w^2) / Q(i w)   (phase condition).
   ===================================================================== */
var D_lin = l * (4.0 * M + m) / 3.0;
var a_lin = (M + m) * g / D_lin;

/* real and imaginary part of Q(i w), written to q[0], q[1] */
function Q_of_iw(w, q)
{
    var Kp = gain[0], Kd = gain[1], kp = gain[2], kd = gain[3];
    var c  = 4.0 * l / 3.0, w2 = w * w;
    q[0] = (-w2 * (Kp - c * kp) + g * kp) / D_lin;
    q[1] = (-w2 * w * (Kd - c * kd) + g * kd * w) / D_lin;
}

/* h(w) = |Q(iw)|^2 - (w^4 + a w^2)^2 ; zeros are candidate Hopf frequencies */
var qtmp = new Float64Array(2);
function hfun(w)
{
    var p = w * w * w * w + a_lin * w * w;
    Q_of_iw(w, qtmp);
    return qtmp[0] * qtmp[0] + qtmp[1] * qtmp[1] - p * p;
}

/* Stability for tau = 0 (Routh-Hurwitz) of
   lambda^4 + c3 lambda^3 + c2 lambda^2 + c1 lambda + c0.
   Zero roots caused by kp = 0 (and kd = 0) are cart translations; they
   are factored out, so "stable" then means: stable up to cart drift.  */
function stable_without_delay()
{
    var c  = 4.0 * l / 3.0;
    var c3 = (gain[1] - c * gain[3]) / D_lin;
    var c2 = (gain[0] - c * gain[2]) / D_lin - a_lin;
    var c1 = g * gain[3] / D_lin;
    var c0 = g * gain[2] / D_lin;
    if (c0 == 0.0 && c1 == 0.0)             /* lambda^2 (lambda^2 + c3 lambda + c2) */
        return (c3 > 0.0 && c2 > 0.0) ? 1 : 0;
    if (c0 == 0.0)                          /* lambda (lambda^3 + c3 lambda^2 + c2 lambda + c1) */
        return (c3 > 0.0 && c1 > 0.0 && c3 * c2 > c1) ? 1 : 0;
    return (c3 > 0.0 && c2 > 0.0 && c1 > 0.0 && c0 > 0.0 &&
            c3 * c2 > c1 && c3 * c2 * c1 > c1 * c1 + c3 * c3 * c0) ? 1 : 0;
}

/* ---------------------------------------------------------------------
   critical_delay: smallest tau > 0 at which roots cross the imaginary
   axis.  Scans w in (0, 120] for sign changes of h, refines each by
   bisection, and gets tau from the phase condition.
   Writes out[0] = tau_c, out[1] = omega_c; returns 1 on success,
   0 if no crossing exists, -1 if the system is unstable already at tau=0.
   --------------------------------------------------------------------- */
function critical_delay(out)
{
    var w, wprev, hprev, hcur, lo, hi, flo, mid, fm, k, p, den, zr, zi, th, ta;
    var best = 1e30, wbest = 0.0;

    if (!stable_without_delay()) return -1;

    wprev = 1e-3;  hprev = hfun(wprev);
    for (w = 0.005; w < 120.0; w += 0.005) {
        hcur = hfun(w);
        if ((hcur > 0.0) != (hprev > 0.0)) {
            lo = wprev;  hi = w;  flo = hprev;
            for (k = 0; k < 50; k++) {                 /* bisection */
                mid = 0.5 * (lo + hi);  fm = hfun(mid);
                if ((fm > 0.0) == (flo > 0.0)) { lo = mid; flo = fm; } else hi = mid;
            }
            w = 0.5 * (lo + hi);
            Q_of_iw(w, qtmp);
            p   = w * w * w * w + a_lin * w * w;
            den = qtmp[0] * qtmp[0] + qtmp[1] * qtmp[1];
            zr  = -p * qtmp[0] / den;                  /* -p / Q = exp(-i w tau) */
            zi  =  p * qtmp[1] / den;
            th  = -Math.atan2(zi, zr);
            th  = ((th % (2.0 * Math.PI)) + 2.0 * Math.PI) % (2.0 * Math.PI);
            ta  = th / w;
            if (ta < best) { best = ta; wbest = w; }
            w = hi;
        }
        hprev = hcur;  wprev = w;
    }
    if (best > 1e29) return 0;
    out[0] = best;  out[1] = wbest;
    return 1;
}

/* allow `node` to load this file for testing (ignored in the browser) */
if (typeof module !== "undefined") module.exports = {
    gain: gain, z: z, sim_reset: sim_reset, sim_step: sim_step,
    critical_delay: critical_delay, setTau: function (v) { tau = v; },
    getT: function () { return t; }
};
