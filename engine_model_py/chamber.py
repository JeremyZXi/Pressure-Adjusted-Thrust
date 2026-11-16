"""Need to find the combustion chamber pressure that balances the flows!!

Flow in from both injectors gets smaller as chamber pressure Pc rises, because the
pressure drop across each injector shrinks:
    mdot_in(Pc) = spi(ox:   P_ox   -> Pc) + spi(fuel: P_fuel -> Pc)

Flow out gets larger as Pc rises:
    mdot_out(Pc) = Pc * A_throat / c_star

the characterisitic velocity (im calling it c_star) measures how much chamber pressure the
combustion makes per unit of mass flow. We treat it as a constant input. 


A possible improvement would be to make it depend on the oxidizer-to-fuel ratio (O/F) using NASA CEA.

Steady state is where IN = OUT. One curve falls and the other rises, so they
cross exactly once.
"""
from dataclasses import dataclass

from injector import spi_mass_flow


@dataclass
class ChamberSolution:
    Pc: float          # chamber pressure [Pa]
    mdot_ox: float     # oxidizer mass flow [kg/s]
    mdot_fuel: float   # fuel mass flow [kg/s]

    @property
    def mdot(self):
        return self.mdot_ox + self.mdot_fuel

    @property
    def OF(self):
        return self.mdot_ox / self.mdot_fuel if self.mdot_fuel > 0 else float('inf')


def solve_chamber(P_ox, rho_ox, P_fuel, rho_fuel, CdA_ox, CdA_fuel, A_throat, c_star):
    """Return the ChamberSolution where mdot_in(Pc) == mdot_out(Pc).

    Residual r(Pc) = mdot_in(Pc) - mdot_out(Pc) is >= 0 at Pc = 0 and < 0 at
    Pc = max(P_ox, P_fuel), and falls monotonically in between, so bisection on
    that bracket always finds the single root (to 1 Pa). If both feed pressures
    are zero, returns Pc = 0 with no flow.
    """
    hi = max(P_ox, P_fuel)
    if hi <= 0.0:
        return ChamberSolution(0.0, 0.0, 0.0)   # both sides empty: no flow, no pressure

    def flows(Pc):
        return (spi_mass_flow(CdA_ox, rho_ox, P_ox, Pc),
                spi_mass_flow(CdA_fuel, rho_fuel, P_fuel, Pc))

    def residual(Pc):
        # flow in through the injectors minus flow out through the choked throat
        m_ox, m_f = flows(Pc)
        return (m_ox + m_f) - Pc * A_throat / c_star

    # residual(0) >= 0 and residual(hi) < 0, and it falls monotonically in between,
    # so bisection always converges to the single crossing.
    lo = 0.0
    while hi - lo > 1.0:          # 1 Pa tolerance
        mid = 0.5 * (lo + hi)
        if residual(mid) > 0:     # more flowing in than out -> Pc must be higher
            lo = mid
        else:
            hi = mid
    Pc = 0.5 * (lo + hi)
    m_ox, m_f = flows(Pc)
    return ChamberSolution(Pc, m_ox, m_f)
