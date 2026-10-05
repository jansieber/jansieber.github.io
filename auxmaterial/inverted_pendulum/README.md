# Delayed cart–pendulum animation

This is a browser animation of an inverted pendulum balanced on a cart by
delayed PD feedback. Open `cart_pendulum_delay.html` in any modern browser.
It needs no server and no build step. `model.js` must sit in the same folder.

## Requests (history)

1. A simple animation of an inverted pendulum on a cart, balanced with
   delayed feedback.
   * Homogeneous rod of mass m, cart of mass M.
   * Angle x, cart position y. The cart moves horizontally.
   * PD feedback with delay τ: large gains on x and x', much smaller gains on
     y and y'.
   * Choose a regime in which the system starts to oscillate after a Hopf
     bifurcation as τ increases.
2. Make the page a portable standalone file.
3. Put the repository under git before making changes (initial commit = version 1).
4. Show the equations as static MathML, so no script library or network is
   needed for the mathematics.
5. Move the mathematics (right-hand side, RK4 loop, history and
   critical-delay computation) into its own file. Write it so that someone who
   knows C but not JavaScript can read it.
6. The falling over corresponds to a connecting orbit. Bifurcation analysis
   predicts that for larger K_d a symmetry breaking of the periodic orbit
   comes first. Find this regime by experimentation.
7. This README.

## Model

* x: angle from upright. y: cart position.
* l = L/2: distance from the pivot to the centre of mass.
* The rod's moment of inertia about its centre of mass is mL²/12.

```
(M+m) y'' + m l (x'' cos x − x'² sin x) = F(t)
(4/3) l x'' + y'' cos x − g sin x       = 0
F(t) = Kp x(t−τ) + Kd x'(t−τ) + kp y(t−τ) + kd y'(t−τ)
```

Parameters: M = 1 kg, m = 0.3 kg, L = 1 m, g = 9.81 m/s².

Default gains: Kp = 20, Kd = 7, kp = 0.3, kd = 1.5. The cart gains have the
same sign as the angle gains. With the opposite sign the system is unstable
even at τ = 0.

Linearising about the upright position gives the characteristic equation

```
λ⁴ − aλ² + e^(−λτ) [ (K(λ) − (4/3) l k(λ)) λ² + g k(λ) ] / D = 0,
K(λ) = Kp + Kd λ,  k(λ) = kp + kd λ,  D = l(4M+m)/3,  a = (M+m)g/D.
```

* With the default gains: τ_c ≈ 0.1715 s and ω_c ≈ 6.11 rad/s.
* The Hopf bifurcation is supercritical: the amplitude grows like √(τ−τ_c).
* The pendulum falls at about 1.1 τ_c.

## Files

| file | content |
|---|---|
| `cart_pendulum_delay.html` | Page, user interface, drawing and static MathML equations. Only calls the functions in `model.js`. |
| `model.js` | The mathematics, in C-like JavaScript: `rhs` (equations of motion, 2×2 solve), `past_state` (linear interpolation in the history ring buffer), `sim_step` (RK4, Δt = 0.5 ms, delayed force evaluated at the stage times), `critical_delay` (τ_c and ω_c from the modulus and phase conditions on the imaginary axis), `stable_without_delay` (Routh–Hurwitz at τ = 0). Can also be loaded with `node` for testing. |
| `experiments/sweep.c` | The same model in C, used for the parameter sweeps below. |
| `experiments/sweep.sh` | Sweeps τ = r·τ_c and stops at the first fall. Needs `node` and a compiled `sweep`. |
| `experiments/tauc.js` | Prints τ_c and ω_c for the given gains, using `model.js`. |

```
cd experiments
cc -O2 -o sweep sweep.c -lm
sh sweep.sh 20 13 0 0 1.01 1.4 0.02      # Kp Kd kp kd r0 r1 dr [T]
```

`sweep` reports statistics over the last 40 % of a 300 s run:

* `amp`: max |x|.
* `mean`: mean tilt x̄.
* `spread`: the spread of the local maxima of x. Zero means the motion is
  periodic.
* `yrange` and `vbar`: the cart's range of motion and its mean velocity.

The initial history is constant: (x, x', y, y') = (0.05, 0, 0.01, 0).

## Symmetry breaking (request 6)

The system is invariant under (x, y) ↦ (−x, −y). The Hopf orbit is symmetric,
so its mean tilt is x̄ = 0. Symmetry breaking shows up as a periodic orbit with
x̄ ≠ 0.

**Angle equation alone (kp = kd = 0).** In this case the x-equation does not
depend on y. The first point in each row where x̄ ≠ 0 is marked.

| Kd | τ_c [s] | symmetric orbit up to | asymmetric orbit (x̄ ≠ 0) | falls at |
|---|---|---|---|---|
| 4   | 0.164 | 1.12 τ_c | not seen | 1.14 τ_c |
| 5.5 | 0.171 | 1.16 τ_c | 1.18 τ_c (x̄ = −0.17 at 1.20) | 1.26 τ_c |
| 7   | 0.151 | 1.12 τ_c | 1.14 τ_c (x̄ = ±0.13) | 1.28 τ_c |
| 10  | 0.111 | 1.04 τ_c | 1.06 τ_c (x̄ = 0.10) | 1.28 τ_c |
| 13  | 0.086 | 1.01 τ_c | 1.03–1.05 τ_c (x̄ = 0.17 at 1.05) | 1.27 τ_c |

So for small Kd the pendulum falls from the symmetric orbit. As Kd grows, a
symmetry-breaking pitchfork of the periodic orbit appears before the fall and
moves towards τ_c. The asymmetric orbits stay periodic (spread = 0) until
just before the fall.

* On the asymmetric orbit the mean force is Kp x̄ ≠ 0.
* So the cart accelerates steadily towards the lean, and the cart velocity
  grows without bound.
* Physically, the controller keeps a leaning pendulum up by accelerating
  under it.

**With cart feedback (kp = 0.3, kd = 1.5).** Averaging the cart equation over
one period requires Kp x̄ + kp ȳ + kd ȳ' = 0. This cannot hold for long while
x̄ ≠ 0 and the cart gains are small. What you see instead:

* The cart runs off. Its own feedback pushes the mean tilt to the other side,
  and the cart runs back.
* The result is a fast pendulum oscillation with a slowly flipping mean tilt
  and large, slow cart sweeps.
* The sweeps have a period of about 5 s and an amplitude of 2–3 m.
* This looks like the symmetry-breaking pitchfork unfolded by the slow cart
  variables (a relaxation oscillation between the two asymmetric branches).

| Kd | τ_c [s] | symmetric up to | cart sweeps from | falls at |
|---|---|---|---|---|
| 7  | 0.1715 | 1.08 τ_c | (modulated at 1.10 τ_c) | 1.12 τ_c |
| 10 | 0.1244 | 1.04 τ_c | 1.06 τ_c | 1.18 τ_c |
| 13 | 0.0938 | 1.01 τ_c | 1.03 τ_c | 1.15 τ_c |

These are simulation results on a grid of step 0.02 τ_c (0.01 τ_c near τ_c),
not a continuation. Locating the pitchfork of periodic orbits exactly needs
DDE continuation, for example in DDE-BIFTOOL.

## User interface (minimal)

* **Scenarios:**
  * Hopf oscillation: Kd = 7, τ = 1.04 τ_c.
  * Symmetry breaking with cart feedback: Kd = 13, τ = 1.06 τ_c.
  * Symmetry breaking, pendulum only: Kd = 13, kp = kd = 0, τ = 1.08 τ_c. The
    view follows the cart, which keeps accelerating.
* **τ slider and buttons** for fixed multiples of τ_c.
* **Gain sliders** under "Feedback gains". τ_c is recomputed whenever a gain
  changes.
* **Scene:** the dashed amber outline is the delayed state that the controller
  sees.
* **Traces:** x(t), y(t), and the mean tilt x̄ over one Hopf period (amber).
  The status pill says "tilted" when |x̄| > 0.03 rad.
* **Internet use:** the page loads its fonts from Google Fonts. Offline, it
  falls back to system fonts and everything else still works.
