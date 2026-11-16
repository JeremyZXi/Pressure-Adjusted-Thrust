"""Sea-level static fire.

    python3 static_fire.py [--no-show] [--eng]

Saves static_fire.png / static_fire.csv; --eng also writes thrust_curve.eng for OpenRocket.
The burn ends when the N2O liquid or the fuel runs out.
"""
import argparse
import csv
import math
from pathlib import Path

import numpy as np

from config import PSI, EngineConfig
from engine import G0, LiquidEngine
from integrator import rk4_step

HERE = Path(__file__).resolve().parent
SEA_LEVEL = 101325.0


def run(config=None, dt=0.01, t_max=60.0, ambient_pressure=SEA_LEVEL):
    engine = LiquidEngine(config or EngineConfig())
    s = engine.initial_state()
    t = 0.0
    rows = []
    while t <= t_max:
        tank, ch, P_fuel = engine.evaluate(s)
        if tank.liquid_gone or s[2] <= 0.0:
            reason = 'N2O liquid ran out' if tank.liquid_gone else 'fuel ran out'
            break
        rows.append({
            't': t, 'P_tank': tank.P, 'T_tank': tank.T, 'm_liq': tank.m_liq,
            'm_vap': tank.m_vap, 'm_fuel': s[2], 'P_fuel': P_fuel, 'Pc': ch.Pc,
            'mdot_ox': ch.mdot_ox, 'mdot_fuel': ch.mdot_fuel, 'OF': ch.OF,
            'thrust': engine.thrust(ch.Pc, ambient_pressure),
        })
        s = rk4_step(engine.derivatives, t, s, dt)
        t += dt
    else:
        reason = f'hit t_max = {t_max} s'
    data = {key: np.array([r[key] for r in rows]) for key in rows[0]}
    return data, reason, engine


def write_eng(d, cfg, path, name='SEB_CF_N2O_IPA', diameter_mm=200, length_mm=3000):
    """RASP .eng thrust curve (sea level) for OpenRocket. Header: name, diameter, length,
    delays, propellant mass, total mass, maker. Total mass = propellant only (hardware mass
    belongs in the rocket's own model). Points are thinned to ~0.1 s spacing."""
    t, F = d['t'], d['thrust']
    prop = (cfg.ox_mass - d['m_liq'][-1] - d['m_vap'][-1]) + (cfg.fuel_mass - d['m_fuel'][-1])
    keep = np.unique(np.concatenate([np.searchsorted(t, np.arange(0.0, t[-1], 0.1)), [len(t) - 1]]))
    with open(path, 'w') as f:
        f.write('; Pressure-adjusted N2O/IPA engine model (engine_model_py/static_fire.py)\n')
        f.write('; Sea-level thrust. Propellant loads in config.py may be placeholders.\n')
        f.write(f'{name} {diameter_mm} {length_mm} 0 {prop:.3f} {prop:.3f} SEB\n')
        for i in keep:
            f.write(f'   {max(t[i], 0.001):.3f} {F[i]:.1f}\n')
        f.write(f'   {t[-1] + 0.05:.3f} 0.0\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-show', action='store_true')
    parser.add_argument('--eng', action='store_true', help='also write thrust_curve.eng')
    args = parser.parse_args()

    cfg = EngineConfig()
    d, reason, eng = run(cfg)
    t = d['t']
    trapezoid = getattr(np, 'trapezoid', None) or np.trapz  # name changed in numpy 2
    impulse = trapezoid(d['thrust'], t)
    prop_used = (cfg.ox_mass - d['m_liq'][-1] - d['m_vap'][-1]) + (cfg.fuel_mass - d['m_fuel'][-1])
    area = math.pi * cfg.tank_inner_diameter ** 2 / 4

    print('--- Setup ---')
    print(f'Fuel feed: {cfg.fuel_feed}; N2O fill {cfg.ox_fill_pressure/PSI:.0f} psi -> '
          f'{eng.T0-273.15:.1f} C saturated')
    print(f'N2O tank {eng.V_ox0*1e3:.1f} L ({eng.V_ox0/area:.2f} m long at ID {cfg.tank_inner_diameter*1e3:.1f} mm); '
          f'IPA {cfg.fuel_mass/cfg.fuel_density*1e3:.1f} L ({cfg.fuel_mass/cfg.fuel_density/area:.2f} m)'
          '   <- check against real hardware')
    print(f'Injector CdA: ox {eng.CdA_ox*1e6:.1f} mm^2, fuel {eng.CdA_fuel*1e6:.1f} mm^2'
          + ('  (calibrated to RPA design point)' if cfg.CdA_ox is None else ''))
    print('--- At ignition vs RPA design point ---')
    print(f'Pc {d["Pc"][0]/PSI:.0f} psi (RPA 350) | mdot {d["mdot_ox"][0]+d["mdot_fuel"][0]:.2f} kg/s (RPA 6.6) | '
          f'O/F {d["OF"][0]:.2f} (RPA 4.5) | sea-level thrust {d["thrust"][0]/1e3:.2f} kN (RPA 13.56)')
    print('--- Burn ---')
    print(f'Ended at t = {t[-1]:.2f} s ({reason})')
    print(f'Tank pressure: {d["P_tank"][0]/PSI:.0f} -> {d["P_tank"][-1]/PSI:.0f} psi  '
          f'(ratio {d["P_tank"][-1]/d["P_tank"][0]:.3f})')
    print(f'Tank temperature: {d["T_tank"][0]-273.15:.1f} -> {d["T_tank"][-1]-273.15:.1f} C')
    print(f'Chamber pressure: {d["Pc"][0]/PSI:.0f} -> {d["Pc"][-1]/PSI:.0f} psi')
    print(f'Thrust: {d["thrust"][0]/1e3:.2f} -> {d["thrust"][-1]/1e3:.2f} kN')
    print(f'O/F: {d["OF"][0]:.2f} -> {d["OF"][-1]:.2f}')
    print(f'Total impulse: {impulse/1e3:.1f} kN*s, average sea-level Isp: {impulse/(prop_used*G0):.0f} s')
    print(f'Left over (last full step): {d["m_liq"][-1]:.2f} kg N2O liquid, {d["m_vap"][-1]:.2f} kg N2O vapor, '
          f'{d["m_fuel"][-1]:.2f} kg fuel')

    with open(HERE / 'static_fire.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(d.keys())
        writer.writerows(zip(*d.values()))

    if args.eng:
        write_eng(d, cfg, HERE / 'thrust_curve.eng')
        print(f'Wrote {HERE / "thrust_curve.eng"}')

    import matplotlib
    if args.no_show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    a = ax[0, 0]
    a.plot(t, d['P_tank'] / PSI, label='N2O tank')
    a.plot(t, d['P_fuel'] / PSI, ls=':', label=f'Fuel feed ({cfg.fuel_feed})')
    a.plot(t, d['Pc'] / PSI, label='Chamber')
    a.set_ylabel('Pressure (psi)'); a.legend()
    a = ax[0, 1]
    a.plot(t, d['mdot_ox'], label='N2O')
    a.plot(t, d['mdot_fuel'], label='IPA')
    a.set_ylabel('Mass flow (kg/s)'); a.legend(loc='center left')
    b = a.twinx(); b.plot(t, d['OF'], color='gray', ls='--'); b.set_ylabel('O/F (dashed)')
    a = ax[1, 0]
    a.plot(t, d['T_tank'] - 273.15, color='tab:red')
    a.set_ylabel('Tank temperature (C)', color='tab:red')
    b = a.twinx(); b.plot(t, d['m_liq'], color='tab:blue'); b.set_ylabel('Liquid N2O (kg)', color='tab:blue')
    a = ax[1, 1]
    a.plot(t, d['thrust'] / 1e3, label='Equilibrium + injector model')
    a.axhline(13.56, ls=':', color='gray', label='RPA design point (13.56 kN SL)')
    a.set_ylabel('Sea-level thrust (kN)'); a.legend()
    for a in ax.flat:
        a.set_xlabel('Time (s)'); a.grid(alpha=0.25)
    fig.suptitle('Static fire: Eureka CF engine (N2O / IPA, piston-fed). '
                 'Propellant loads and CdA are placeholders, see config.py')
    fig.savefig(HERE / 'static_fire.png', dpi=150)
    print(f'Saved {HERE / "static_fire.png"} and static_fire.csv')
    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
