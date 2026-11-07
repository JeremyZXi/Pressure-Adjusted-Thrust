from pathlib import Path

import numpy as np
import pandas as pd

"""These are all the inputs required for this simulation.
I left them hardcoded here for now but essentially all of them are input variables
"""
# ---thrust---
times = [0.05, 8.50, 8.51, 9.50] # time t corresponding the thrust
thrust_points = [8133, 7603, 1200, 1000] # thrust in N
# nozzle data
nozzle_exit_diameter = 0.1241
nozzle_exit_area = np.pi * nozzle_exit_diameter**2 / 4.0
propellant_points = [44, 5, 5, 0] # remaining propellant in kg

# --- aerodynamics ---
dry_mass = 46.7 # kg
cg = 2.33  # m, center of gravity, assumed the same for both dry mass and propellant
# create a lookup table to get corresponding Cd under certain velocity in vertical motion
aero = pd.read_csv(Path(__file__).with_name("Eureka3_Machnum_RASAero2_sim.CSV"))
straight_flight = aero.loc[aero["Alpha"] == 0, ["Mach", "CD", "CN", "CP"]]
straight_flight = straight_flight.sort_values("Mach")
mach_table = straight_flight["Mach"].to_numpy()
cd_table = straight_flight["CD"].to_numpy()
cn_table = straight_flight["CN"].to_numpy()
cp_table = straight_flight["CP"].to_numpy()
# reference area
area = 0.03243 # m^2
# time step
dt = 0.01 # s


def g(z):
    """Calculate the gravitational acceleration at altitude z(m)
    Parameters
    ----------
    z : float
        Current altitude of the rocket in meter
    Returns
    -------
    float
        the gravitational acceleration at altitude z(m) in m/s^2
    """
    g0 = 9.80665
    earth_radius = 6_371_000.0
    return g0 * (earth_radius / (earth_radius + z)) ** 2
def mach_number(z, v):
    """Calculate current mach number
    Parameters
    ----------
    z : float
        Current altitude of the rocket in meter
    v : float
        Current velocity of the rocket in m/s
    Returns
    -------
    float
        speed v in mach at some height z
    """
    air_temperature = 288.15 - 0.0065 * np.clip(z, 0, 11_000) # calculate the current air temperature
    speed_of_sound = np.sqrt(1.4 * 287.05 * air_temperature) # calculate the sound of speed at current height
    return v / speed_of_sound

# drag coefficient
def cd(z, v):
    """Look up the drag coefficient
    Parameters
    ----------
    z : float
        Current altitude of the rocket in meter
    v : float
        Current velocity of the rocket in m/s
    Returns
    -------
    float
        Drag coefficient at altitude z and velocity v with 0 attack angle
    """
    # interpolate from the lookup table
    return np.interp(mach_number(z, v), mach_table, cd_table)

def ambient_pressure(z):
    """Ambient air pressure (Pa) at altitude z (m)."""
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
    """Air density (kg/m³) at altitude z (m)."""
    h = max(0.0, z)
    temperature = 288.15 - 0.0065 * h if h <= 11_000 else 216.65
    return ambient_pressure(h) / (287.05 * temperature)


def ambient_thrust_correction(thrust, z, A_exit, p_ref=101325.0):
    if thrust <= 0.0:
        return 0.0  # No thrust correction after burnout.

    return max(0.0, thrust + (p_ref - ambient_pressure(z)) * A_exit)

# ---solver for RK4 method---
def derivatives(t, state):
    z, v = state
    thrust = np.interp(t, times, thrust_points, left=0.0, right=0.0)
    thrust = ambient_thrust_correction(thrust, z, nozzle_exit_area)
    mass = dry_mass + np.interp(t, times, propellant_points, left=0.0, right=0.0)
    drag = 0.0
    if v > 0.0:
        dynamic_pressure = 0.5 * rho(z) * v**2
        drag = dynamic_pressure * area * cd(z, v)
    weight = mass * g(z)
    a = (thrust - drag - weight) / mass
    return np.array([v, a])  # [dz/dt, dv/dt]

def rk4_step(t, state, dt):
    k1 = derivatives(t, state)
    k2 = derivatives(t + dt / 2, state + dt * k1 / 2)
    k3 = derivatives(t + dt / 2, state + dt * k2 / 2)
    k4 = derivatives(t + dt, state + dt * k3)
    return state + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


t = 0.0
state = np.array([0.0, 0.0])  # [altitude, velocity], used to record current state during the loop
history = [] # keep

time_limit = 100
# run the program for time_limit seconds
while t < time_limit:
    z, v = state
    density = rho(z)
    mach = mach_number(z, v)
    mass = dry_mass + np.interp(t, times, propellant_points, left=0.0, right=0.0)
    thrust = np.interp(t, times, thrust_points, left=0.0, right=0.0)
    thrust = ambient_thrust_correction(thrust, z, nozzle_exit_area)
    history.append((
        t, z, v, derivatives(t, state)[1], mass, cg, mach,
        0.5 * density * v**2, density,
        np.interp(mach, mach_table, cd_table),
        np.interp(mach, mach_table, cn_table),
        np.interp(mach, mach_table, cp_table),
        thrust, 1.0,  # Tank correction disabled.
    ))
    if t > 0.1 and (v < 0.0 or z < 0.0):
        break
    state = rk4_step(t, state, dt)
    t += dt

flight_data = pd.DataFrame(history, columns=[
    "time", "altitude", "velocity", "acceleration", "mass", "CG",
    "mach", "q", "density", "CD", "CN", "CP", "thrust", "tankPressureRatio",
])
# Keep the Python output separate so it doesn't overwrite the C++ run.
output = Path(__file__).with_name("flight_data_python_ambient_pressure.csv")
flight_data.to_csv(output, index=False, float_format="%.6g")


# console output
print("\nSimulation Results:")
print(f"  Apogee: {flight_data['altitude'].max():.6g} m")
print(f"  Max Velocity: {flight_data['velocity'].max():.6g} m/s")
print(f"  Max Mach: {flight_data['mach'].max():.6g}")
print(f"  Flight data exported to: {output.name}")
print("\nExample completed successfully!")
