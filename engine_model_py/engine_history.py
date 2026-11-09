"""Calculate what happens inside the engine during the burn.

Start with the tank and fuel inputs from config.py. At each step, calculate
how much propellant flows through the injectors and find the chamber pressure
where that flow matches the nozzle outflow. Then update the tank energy and
remaining fuel and nitrous. Stop when the engine can no longer feed liquid.

rocket.py uses the saved time, chamber pressure, remaining propellant mass,
and tank pressure to calculate thrust and rocket mass during flight. You
don't need to run this file directly.
"""

import numpy as np

from engine_model_py.config import injector_CdA, ox_initial_temperature, ox_tank_volume
from engine_model_py.fuel_tank import PistonFuelTank, RegulatedFuelTank
from engine_model_py.injector import spi_mass_flow
from engine_model_py.integrator import rk4_step
from engine_model_py.n2o_props import SaturationTable
from engine_model_py.ox_tank import EquilibriumOxTank


class EngineCutoff(ValueError):
    """The current liquid-feed model cannot continue burning."""


def build_engine(cfg):
    """Set up the tanks and injectors from the engine configuration.

    Return the starting state [nitrous mass, nitrous energy, fuel mass],
    along with functions that calculate pressures, flows, and state changes.
    """
    props = SaturationTable()
    fill_volume = ox_tank_volume(cfg, props)
    fill_temperature = ox_initial_temperature(cfg, props)
    tank = EquilibriumOxTank(fill_volume, cfg.ox_mass, fill_temperature, props)
    cda_ox, cda_fuel = injector_CdA(cfg, props)
    if cfg.fuel_feed == 'piston':
        fuel = PistonFuelTank(cfg.fuel_mass, cfg.fuel_density, cfg.piston_dp)
    elif cfg.fuel_feed == 'regulated':
        fuel = RegulatedFuelTank(cfg.fuel_pressure, cfg.fuel_mass, cfg.fuel_density)
    else:
        raise ValueError('Unknown fuel feed configuration')

    def engine_state(state):
        """Find the chamber pressure where injector inflow equals nozzle outflow.

        Return the nitrous tank state, chamber pressure [Pa], and the
        nitrous and fuel flow rates [kg/s].
        """
        ox_mass, ox_energy, fuel_mass = state
        if ox_mass <= 0 or fuel_mass <= 0:
            raise EngineCutoff('Propellant depleted')
        ox_volume = fill_volume + fuel.volume_freed(fuel_mass)
        ox = tank.solve(ox_mass, ox_energy, ox_volume)
        if ox.liquid_gone or ox.m_vap < 0:
            raise EngineCutoff('Oxidizer left the liquid-feed range')
        fuel_pressure = fuel.pressure(fuel_mass, ox.P)

        def flows(pc):
            """Return nitrous and fuel flow [kg/s] at a trial chamber pressure [Pa]."""
            return (spi_mass_flow(cda_ox, ox.rho_liq, ox.P, pc),
                    spi_mass_flow(cda_fuel, cfg.fuel_density, fuel_pressure, pc))

        # Close the chamber mass balance: injector flow = Pc * At / c*.
        lo, hi = 0.0, max(ox.P, fuel_pressure)
        for _ in range(45):
            pc = (lo + hi) / 2
            nozzle_flow = pc * cfg.throat_area / cfg.c_star
            if sum(flows(pc)) > nozzle_flow:
                lo = pc
            else:
                hi = pc
        pc = (lo + hi) / 2
        ox_flow, fuel_flow = flows(pc)
        if min(ox_flow, fuel_flow) <= 0:
            raise EngineCutoff('One propellant stopped flowing')
        return ox, pc, ox_flow, fuel_flow

    def derivatives(_time, state):
        """Return how fast nitrous mass, tank energy, and fuel mass change.

        The tank loses energy as nitrous exits and, for piston feed, as it
        pushes the piston. RK4 uses these rates to advance the state.
        """
        ox, _, ox_flow, fuel_flow = engine_state(state)
        volume_rate = fuel_flow / cfg.fuel_density if fuel.feeds_from_ox_tank else 0.0
        energy_rate = -ox_flow * ox.h_liq - ox.P * volume_rate
        return np.array([-ox_flow, energy_rate, -fuel_flow])

    state = np.array([tank.M0, tank.U0, cfg.fuel_mass])
    return state, engine_state, derivatives


def generate_engine_history(cfg, dt=0.01, max_time=60.0):
    """Step the engine forward until it can no longer feed liquid.

    Each row contains time [s], chamber pressure [Pa], total remaining
    propellant [kg], and nitrous tank pressure [Pa].
    """
    state, evaluate, derivatives = build_engine(cfg)
    rows = []

    for step in range(int(max_time / dt)):
        time = step * dt
        try:
            ox, pc, _, _ = evaluate(state)
            propellant_mass = state[0] + state[2]
            rows.append((time, pc, propellant_mass, ox.P))
            state = rk4_step(derivatives, time, state, dt)
        except EngineCutoff:
            # Stop at the last valid sample, including cutoff during an RK4 step.
            if not rows:
                raise
            return np.array(rows)

    raise RuntimeError("No engine cutoff within the time limit")
