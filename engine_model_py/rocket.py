"""Vertical 1D flight driven by the liquid-feed engine model.

Tracks altitude and velocity with RK4 and exports the C++-compatible CSV.
Engine inputs live in config.py. Residual propellant remains aboard after cutoff.
"""
from pathlib import Path
import sys

# Support both direct IDE execution and python -m engine_model_py.rocket.
project_root = Path(__file__).resolve().parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from engine_model_py.config import EngineConfig
from engine_model_py.engine_history import generate_engine_history
from engine_model_py.nozzle import thrust_coefficient
import numpy as np
import pandas as pd

# Engine-model input parameters live in config.py.
cfg = EngineConfig()
ignition_time = 0.05
engine_table = generate_engine_history(cfg)
engine_times = engine_table[:, 0] + ignition_time
chamber_pressures = engine_table[:, 1]
propellant_masses = engine_table[:, 2]
tank_pressures = engine_table[:, 3]
cutoff_time = engine_times[-1]

# --- aerodynamics ---
dry_mass = 46.7 # kg
cg = 2.33  # m, center of gravity, assumed the same for both dry mass and propellant
# Aerodynamic coefficients at zero angle of attack, indexed by Mach number.
aero = pd.read_csv(project_root / "Eureka3_Machnum_RASAero2_sim.CSV")
straight_flight = aero.loc[aero["Alpha"] == 0, ["Mach", "CD", "CN", "CP"]]
straight_flight = straight_flight.sort_values("Mach")
mach_table = straight_flight["Mach"].to_numpy()
cd_table = straight_flight["CD"].to_numpy()
cn_table = straight_flight["CN"].to_numpy()
cp_table = straight_flight["CP"].to_numpy()
# Aerodynamic reference area.
area = 0.03243 # m^2
# Flight integration timestep.
dt = 0.01 # s


def g(z):
    """Gravitational acceleration [m/s^2] at altitude z [m], using inverse-square gravity."""
    g0 = 9.80665
    earth_radius = 6_371_000.0
    return g0 * (earth_radius / (earth_radius + z)) ** 2

def mach_number(z, v):
    """Signed Mach number for velocity v [m/s] at altitude z [m]."""
    air_temperature = 288.15 - 0.0065 * np.clip(z, 0, 11_000)  # Air temperature [K].
    speed_of_sound = np.sqrt(1.4 * 287.05 * air_temperature)  # Speed of sound [m/s].
    return v / speed_of_sound

# drag coefficient
def cd(z, v):
    """Drag coefficient at zero angle of attack, interpolated from the Mach table.

    z : altitude [m]
    v : vertical velocity [m/s]
    Returns dimensionless C_D."""
    # interpolate from the lookup table
    return np.interp(mach_number(z, v), mach_table, cd_table)

def ambient_pressure(z):
    """Ambient air pressure [Pa] at altitude z [m].

    Uses a linear temperature lapse rate below 11 km and a constant
    temperature above it. Negative altitude uses sea-level pressure."""

    h = max(0.0, z)
    R_air = 287.05
    g0 = 9.80665
    lapse_rate = 0.0065
    exponent = g0 / (lapse_rate * R_air)

    if h <= 11_000:
        temperature = 288.15 - lapse_rate * h
        return 101325.0 * (temperature / 288.15) ** exponent

    temperature = 216.65
    pressure_11km = 101325.0 * (temperature / 288.15) ** exponent
    return pressure_11km * np.exp(
        -g0 * (h - 11_000) / (R_air * temperature)
    )


def rho(z):
    """Air density [kg/m^3] at altitude z [m], from rho = P / (R_air * T)."""
    h = max(0.0, z)
    temperature = 288.15 - 0.0065 * h if h <= 11_000 else 216.65
    return ambient_pressure(h) / (287.05 * temperature)


def engine_thrust(t, z):
    """Nozzle thrust [N] at simulation time t [s] and altitude z [m].

    Interpolates chamber pressure from the engine history and includes
    ambient pressure through C_F. Returns zero before ignition and
    at or after liquid-feed cutoff."""
    if t < ignition_time or t >= cutoff_time:
        return 0.0
    pc = np.interp(t, engine_times, chamber_pressures)
    cf = thrust_coefficient(pc, ambient_pressure(z), cfg.expansion_ratio, cfg.gamma)
    return max(0.0, cfg.cf_efficiency * cf * pc * cfg.throat_area)

# Flight equations and integration.
def derivatives(t, state):
    """Flight-state derivative [v, a] for RK4.

    t     : simulation time [s]
    state : numpy array [altitude (m), vertical velocity (m/s)]
    Returns [dz/dt (m/s), dv/dt (m/s^2)]. Includes thrust, upward-flight
    drag and weight, with pad support before ignition."""
    z, v = state
    thrust = engine_thrust(t, z)
    mass = dry_mass + np.interp(t, engine_times, propellant_masses)
    drag = 0.0
    if v > 0.0:
        dynamic_pressure = 0.5 * rho(z) * v**2
        drag = dynamic_pressure * area * cd(z, v)
    weight = mass * g(z)
    a = (thrust - drag - weight) / mass
    if t < ignition_time and z <= 0.0 and v <= 0.0:
        a = 0.0  # Pad supports the fully loaded rocket until ignition.
    return np.array([v, a])  # [dz/dt, dv/dt]

def rk4_step(t, state, dt):
    """Advance [altitude, velocity] by one fourth-order Runge-Kutta step.

    t     : current time [s]
    state : numpy array [altitude (m), vertical velocity (m/s)]
    dt    : timestep [s]
    Returns the new state after dt."""
    k1 = derivatives(t, state)
    k2 = derivatives(t + dt / 2, state + dt * k1 / 2)
    k3 = derivatives(t + dt / 2, state + dt * k2 / 2)
    k4 = derivatives(t + dt, state + dt * k3)
    return state + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


t = 0.0
state = np.array([0.0, 0.0])  # [altitude, velocity], used to record current state during the loop
history = []  # One flight snapshot per timestep.

time_limit = 100
# Stop at apogee or the time limit.
while t < time_limit:
    z, v = state
    density = rho(z)
    mach = mach_number(z, v)
    mass = dry_mass + np.interp(t, engine_times, propellant_masses)
    thrust = engine_thrust(t, z)
    history.append((
        t, z, v, derivatives(t, state)[1], mass, cg, mach,
        0.5 * density * v**2, density,
        np.interp(mach, mach_table, cd_table),
        np.interp(mach, mach_table, cn_table),
        np.interp(mach, mach_table, cp_table),
        thrust, np.interp(t, engine_times, tank_pressures) / tank_pressures[0],
    ))
    if t > 0.1 and (v < 0.0 or z < 0.0):
        break
    if t < ignition_time:
        t = min(t + dt, ignition_time)
        continue
    state = rk4_step(t, state, dt)
    t += dt

flight_data = pd.DataFrame(history, columns=[
    "time", "altitude", "velocity", "acceleration", "mass", "CG",
    "mach", "q", "density", "CD", "CN", "CP", "thrust", "tankPressureRatio",
])
# Keep the Python output separate so it doesn't overwrite the C++ run.
output = Path(__file__).with_name("flight_data_python_engine.csv")
flight_data.to_csv(output, index=False, float_format="%.6g")


# Console summary.
print("\nSimulation Results (engine model; illustrative inputs):")
print(f"  Liquid-feed cutoff: {cutoff_time:.3f} s")
print(f"  Residual propellant: {propellant_masses[-1]:.3f} kg")
print(f"  Apogee: {flight_data['altitude'].max():.6g} m")
print(f"  Max Velocity: {flight_data['velocity'].max():.6g} m/s")
print(f"  Max Mach: {flight_data['mach'].max():.6g}")
print(f"  Flight data exported to: {output.name}")
print("\nExample completed successfully!")
