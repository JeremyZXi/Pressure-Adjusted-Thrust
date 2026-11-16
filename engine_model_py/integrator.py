"""One 4th-order Runge-Kutta step for any state vector s, with f(t, s) = ds/dt."""


def rk4_step(f, t, s, dt):
    k1 = f(t, s)
    k2 = f(t + 0.5 * dt, s + 0.5 * dt * k1)
    k3 = f(t + 0.5 * dt, s + 0.5 * dt * k2)
    k4 = f(t + dt, s + dt * k3)
    return s + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
