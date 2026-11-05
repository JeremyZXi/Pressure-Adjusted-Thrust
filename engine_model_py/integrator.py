"""Generic 4th-order Runge-Kutta step for any state vector.

The state s is a numpy array: [h, v] in the flight sim, [M_ox, U_ox, m_fuel] in the
static fire, and all five combined later. f(t, s) returns ds/dt, an array of the
same length. RK4 samples the slope four times across the step and takes a
weighted average (same 1-2-2-1 / 6 weights as Simpson's rule).
"""


def rk4_step(f, t, s, dt):
    k1 = f(t, s)
    k2 = f(t + 0.5 * dt, s + 0.5 * dt * k1)
    k3 = f(t + 0.5 * dt, s + 0.5 * dt * k2)
    k4 = f(t + dt, s + dt * k3)
    return s + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
