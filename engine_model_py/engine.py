"""Connects the tank, injectors, chamber and nozzle.

State: s = [M_ox, U_ox, m_fuel] (N2O mass, N2O internal energy, fuel mass).
Each derivative call: tank (M, U, V) -> P_tank -> feed pressures -> Pc and flows -> ds/dt.
"""
import numpy as np

import config as C
from chamber import solve_chamber
from fuel_tank import PistonFuelTank, RegulatedFuelTank
from n2o_props import SaturationTable
from nozzle import thrust_coefficient
from ox_tank import EquilibriumOxTank

G0 = 9.80665


class LiquidEngine:
    def __init__(self, config, props=None):
        self.cfg = config
        self.props = props or SaturationTable()
        self.T0 = C.ox_initial_temperature(config, self.props)
        self.V_ox0 = C.ox_tank_volume(config, self.props)
        self.ox_tank = EquilibriumOxTank(self.V_ox0, config.ox_mass, self.T0, self.props)
        if config.fuel_feed == 'piston':
            self.fuel_tank = PistonFuelTank(config.fuel_mass, config.fuel_density, config.piston_dp)
        elif config.fuel_feed == 'regulated':
            self.fuel_tank = RegulatedFuelTank(config.fuel_pressure, config.fuel_mass, config.fuel_density)
        else:
            raise ValueError(f'unknown fuel_feed {config.fuel_feed!r}')
        self.CdA_ox, self.CdA_fuel = C.injector_CdA(config, self.props)

    def initial_state(self):
        return np.array([self.ox_tank.M0, self.ox_tank.U0, self.fuel_tank.m0])

    def ox_volume(self, m_fuel):
        """Volume holding N2O: the tank, plus whatever the piston has swept."""
        return self.V_ox0 + self.fuel_tank.volume_freed(m_fuel)

    def evaluate(self, s):
        """Solve tank and chamber for state s. Returns (tank_state, chamber_solution, P_fuel)."""
        M_ox, U_ox, m_fuel = s
        tank = self.ox_tank.solve(M_ox, U_ox, self.ox_volume(m_fuel))
        # Feed pressure is zero once the liquid is gone
        P_ox = 0.0 if tank.liquid_gone else tank.P
        P_fuel = self.fuel_tank.pressure(m_fuel, tank.P)
        chamber = solve_chamber(
            P_ox, tank.rho_liq, P_fuel, self.fuel_tank.rho,
            self.CdA_ox, self.CdA_fuel, self.cfg.throat_area, self.cfg.c_star)
        return tank, chamber, P_fuel

    def derivatives(self, t, s):
        tank, chamber, _ = self.evaluate(s)
        # With a piston, the N2O volume grows as fuel leaves: dV/dt = mdot_fuel / rho_fuel,
        # and the N2O spends energy P*dV/dt pushing the piston.
        dV_dt = chamber.mdot_fuel / self.fuel_tank.rho if self.fuel_tank.feeds_from_ox_tank else 0.0
        return np.array([
            -chamber.mdot_ox,
            -chamber.mdot_ox * tank.h_liq - tank.P * dV_dt,
            -chamber.mdot_fuel,
        ])

    def thrust(self, Pc, ambient_pressure):
        cf = thrust_coefficient(Pc, ambient_pressure, self.cfg.expansion_ratio, self.cfg.gamma)
        return max(0.0, self.cfg.cf_efficiency * cf * Pc * self.cfg.throat_area)
