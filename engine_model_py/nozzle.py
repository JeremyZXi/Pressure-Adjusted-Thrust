"""Ideal converging-diverging nozzle: thrust coefficient C_F.

Thrust = C_F * Pc * A_throat. C_F packs together two contributions:
  * momentum: exhaust leaving fast, which depends only on gamma and the expansion
    ratio eps = A_exit / A_throat
  * pressure: (P_exit - P_ambient) * A_exit. This is the same term the C++ sim
    adds by hand as (P_ref - P_atm) * A_exit

Assumes 1D isentropic flow with constant gamma, no flow separation.
(A strongly over-expanded nozzle at low Pc would really separate, so the
pressure term is too pessimistic late in the burn.)
"""
from functools import lru_cache
import math


def area_ratio(M, gamma):
    """A / A_throat for Mach number M (the isentropic area-Mach relation)."""
    g = gamma
    return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** ((g + 1.0) / (2.0 * (g - 1.0)))


@lru_cache(maxsize=None)
def exit_mach(eps, gamma):
    """Supersonic exit Mach number for expansion ratio eps (bisection on M > 1)."""
    lo, hi = 1.0, 50.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if area_ratio(mid, gamma) < eps:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def exit_pressure_ratio(eps, gamma):
    """P_exit / Pc. Fixed by the nozzle geometry and gamma, not by Pc."""
    M = exit_mach(eps, gamma)
    return (1.0 + 0.5 * (gamma - 1.0) * M * M) ** (-gamma / (gamma - 1.0))


def thrust_coefficient(Pc, Pa, eps, gamma):
    if Pc <= 0.0:
        return 0.0
    g = gamma
    pe = exit_pressure_ratio(eps, g)
    momentum = math.sqrt(2.0 * g * g / (g - 1.0)
                         * (2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0))
                         * (1.0 - pe ** ((g - 1.0) / g)))
    pressure = (pe - Pa / Pc) * eps
    return momentum + pressure
