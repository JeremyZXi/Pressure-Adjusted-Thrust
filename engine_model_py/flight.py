"""1D vertical flight with pressure-adjusted thrust (tank blowdown + altitude).

    python3 flight.py [--no-show]      saves flight.png / flight.csv

Aero and atmosphere come from aero.py. Vehicle dry mass is a placeholder.
"""
import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from aero import AeroData, Atmosphere
from config import EngineConfig
from engine import LiquidEngine
from integrator import rk4_step

HERE = Path(__file__).resolve().parent

G0 = 9.80665
R_EARTH = 6371000.0
SEA_LEVEL = 101325.0


@dataclass
class VehicleConfig:
    dry_mass: float = 46.7                 # kg  TBD: value from 1DSimExample.cpp, NOT the Spaceshot vehicle
    aero_file: str = str(HERE / 'Eureka3_Machnum_RASAero2_sim.CSV')  # RASAero table (Eureka-3; replace for Spaceshot)
    reference_area: float = 0.03243        # m^2 (8 in diameter body), same default as the 1D framework


class PressureAdjustedFlight:
    def __init__(self, engine_cfg=None, vehicle_cfg=None, altitude_correction=True):
        self.engine = LiquidEngine(engine_cfg or EngineConfig())
        self.veh = vehicle_cfg or VehicleConfig()
        self.aero = AeroData(self.veh.aero_file, self.veh.reference_area)
        self.altitude_correction = altitude_correction  # False -> thrust always at sea-level C_F
        self.burning = True

    def mass(self, s):
        return self.veh.dry_mass + s[2] + max(s[4], 0.0)   # all N2O (liquid + vapor) stays on board

    def thrust(self, s, P_amb):
        if not self.burning:
            return 0.0, None
        tank, ch, _ = self.engine.evaluate(s[2:])
        return self.engine.thrust(ch.Pc, P_amb if self.altitude_correction else SEA_LEVEL), (tank, ch)

    def derivatives(self, t, s):
        h, v = s[0], s[1]
        atm = Atmosphere.get_atmosphere(h)
        F, _ = self.thrust(s, atm.pressure)
        cd = self.aero.get_coeffs(0.0, abs(v) / atm.speed_of_sound).CD
        drag = 0.5 * atm.density * v * abs(v) * self.aero.reference_area * cd
        m = self.mass(s)
        weight = m * G0 * (R_EARTH / (R_EARTH + max(h, 0.0))) ** 2
        engine_rates = self.engine.derivatives(t, s[2:]) if self.burning else np.zeros(3)
        return np.concatenate([[v, (F - drag - weight) / m], engine_rates])

    def run(self, dt=0.01, t_max=600.0):
        s = np.concatenate([[0.0, 0.0], self.engine.initial_state()])
        t, rows, burnout = 0.0, [], None
        while t <= t_max:
            atm = Atmosphere.get_atmosphere(s[0])
            F, eng = self.thrust(s, atm.pressure)
            if self.burning and (eng[0].liquid_gone or s[4] <= 0.0):
                self.burning, burnout = False, t           # engine shuts down; coast from here
                F, eng = 0.0, None
            rows.append({'t': t, 'h': s[0], 'v': s[1], 'mach': s[1] / atm.speed_of_sound,
                         'mass': self.mass(s), 'thrust': F, 'P_amb': atm.pressure,
                         'P_tank': eng[0].P if eng else np.nan, 'Pc': eng[1].Pc if eng else np.nan})
            if t > 1.0 and s[1] < 0.0:                      # apogee
                break
            if t < 1.0 and s[0] <= 0.0 and F < self.mass(s) * G0:
                s[0], s[1] = 0.0, max(s[1], 0.0)            # on the pad until thrust exceeds weight
            s = rk4_step(self.derivatives, t, s, dt)
            t += dt
        data = {k: np.array([r[k] for r in rows]) for k in rows[0]}
        return data, burnout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-show', action='store_true')
    args = parser.parse_args()

    adj, t_bo = PressureAdjustedFlight(altitude_correction=True).run()
    sl, _ = PressureAdjustedFlight(altitude_correction=False).run()
    i_bo = int(np.searchsorted(adj['t'], t_bo))
    print(f'Burnout at t = {t_bo:.2f} s, h = {adj["h"][i_bo]/1e3:.2f} km, v = {adj["v"][i_bo]:.0f} m/s')
    print(f'Thrust: {adj["thrust"][0]/1e3:.2f} kN at liftoff -> {adj["thrust"][i_bo-1]/1e3:.2f} kN at burnout '
          f'(sea-level-only model: {sl["thrust"][i_bo-1]/1e3:.2f} kN)')
    print(f'Apogee: {adj["h"].max()/1e3:.2f} km with altitude correction, '
          f'{sl["h"].max()/1e3:.2f} km with sea-level thrust only '
          f'({(adj["h"].max()-sl["h"].max())/1e3:+.2f} km)')
    print(f'Max velocity {adj["v"].max():.0f} m/s, max Mach {adj["mach"].max():.2f}')

    with open(HERE / 'flight.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(adj.keys()); w.writerows(zip(*adj.values()))

    import matplotlib
    if args.no_show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    burn = adj['t'] <= t_bo
    ax[0].plot(adj['t'][burn], adj['thrust'][burn] / 1e3, label='altitude-adjusted')
    ax[0].plot(sl['t'][sl['t'] <= t_bo], sl['thrust'][sl['t'] <= t_bo] / 1e3, ls='--', label='sea-level C_F only')
    ax[0].set_ylabel('Thrust (kN)'); ax[0].legend()
    ax[1].plot(adj['t'][burn], adj['P_tank'][burn] / 6894.757, label='N2O tank')
    ax[1].plot(adj['t'][burn], adj['Pc'][burn] / 6894.757, label='chamber')
    ax[1].set_ylabel('Pressure (psi)'); ax[1].legend()
    ax[2].plot(adj['t'], adj['h'] / 1e3, label='altitude-adjusted')
    ax[2].plot(sl['t'], sl['h'] / 1e3, ls='--', label='sea-level C_F only')
    ax[2].set_ylabel('Altitude (km)'); ax[2].legend()
    for a in ax:
        a.set_xlabel('Time (s)'); a.grid(alpha=.25)
    fig.suptitle('1D flight with pressure-adjusted thrust. Vehicle dry mass and propellant loads are placeholders')
    fig.savefig(HERE / 'flight.png', dpi=150)
    print(f'Saved {HERE / "flight.png"} and flight.csv')
    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
