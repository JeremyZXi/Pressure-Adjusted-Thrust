import math


def spi_mass_flow(CdA, rho, P_up, P_down):
    
    """Mass flow through an injector.

    CdA    : discharge coefficient times total orifice area [m^2]
    rho    : propellant density upstream [kg/m^3]
    P_up   : upstream (tank/manifold) pressure [Pa]
    P_down : downstream (chamber) pressure [Pa]
    """
    dP = P_up - P_down

    if dP < 0:
        dP = 0
    return CdA * math.sqrt(2 * dP * rho)