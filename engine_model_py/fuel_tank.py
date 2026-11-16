"""Fuel tank models (IPA density treated as constant).

PistonFuelTank: N2O gas pushes a free piston, so P_fuel = P_tank - piston friction,
and the N2O volume grows as IPA leaves (the N2O does work P*dV on the piston).
RegulatedFuelTank: fuel held at a constant set pressure.
"""


class PistonFuelTank:
    feeds_from_ox_tank = True

    def __init__(self, fuel_mass, density, piston_dp=0.0):
        self.m0 = fuel_mass        # loaded fuel [kg]
        self.rho = density         # fuel density [kg/m^3]
        self.piston_dp = piston_dp # pressure lost to piston seal friction [Pa]

    def pressure(self, m_fuel, P_ox_tank):
        """Pressure feeding the fuel injector; zero once the piston bottoms out."""
        return max(P_ox_tank - self.piston_dp, 0.0) if m_fuel > 0.0 else 0.0

    def volume_freed(self, m_fuel):
        """Volume the piston has swept so far = extra volume now holding N2O [m^3]."""
        return (self.m0 - m_fuel) / self.rho


class RegulatedFuelTank:
    feeds_from_ox_tank = False

    def __init__(self, pressure, fuel_mass, density):
        self.P = pressure          # regulated tank pressure [Pa]
        self.m0 = fuel_mass
        self.rho = density

    def pressure(self, m_fuel, P_ox_tank=None):
        return self.P if m_fuel > 0.0 else 0.0

    def volume_freed(self, m_fuel):
        return 0.0                 # separate tank: N2O volume is unaffected
