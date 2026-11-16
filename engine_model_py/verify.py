"""Numerical checks: mass balance, timestep convergence, and the property
table vs CoolProp (skipped if CoolProp isn't installed).

    python3 verify.py
"""
import numpy as np

from config import EngineConfig
from engine import LiquidEngine
from static_fire import run

trap = getattr(np, 'trapezoid', None) or np.trapz


def main():
    cfg = EngineConfig()
    d, _, _ = run(cfg, dt=0.01)
    lost_ox = cfg.ox_mass - (d['m_liq'][-1] + d['m_vap'][-1])
    flowed_ox = trap(d['mdot_ox'], d['t'])
    lost_f = cfg.fuel_mass - d['m_fuel'][-1]
    flowed_f = trap(d['mdot_fuel'], d['t'])
    print(f'1. Mass balance  N2O: tank lost {lost_ox:.3f} kg, injector passed {flowed_ox:.3f} kg '
          f'({100*(flowed_ox/lost_ox-1):+.2f}%)')
    print(f'                 IPA: tank lost {lost_f:.3f} kg, injector passed {flowed_f:.3f} kg '
          f'({100*(flowed_f/lost_f-1):+.2f}%)')

    d2, _, _ = run(cfg, dt=0.005)
    I1, I2 = trap(d['thrust'], d['t']), trap(d2['thrust'], d2['t'])
    print(f'2. Timestep      dt=0.01 s: burn {d["t"][-1]:.3f} s, impulse {I1/1e3:.3f} kN*s | '
          f'dt=0.005 s: burn {d2["t"][-1]:.3f} s, impulse {I2/1e3:.3f} kN*s ({100*(I2/I1-1):+.3f}%)')

    try:
        from CoolProp.CoolProp import PropsSI
    except ImportError:
        print('3. Property table: skipped (pip install CoolProp to run)')
        return
    props = LiquidEngine(cfg).props
    worst = {}
    for T in np.linspace(235.0, 305.0, 141) + 0.0537:   # offset so points fall between table rows
        for name, out, q in [('P', 'P', 0), ('rho_l', 'D', 0), ('rho_v', 'D', 1), ('u_l', 'U', 0),
                             ('u_v', 'U', 1), ('h_l', 'H', 0)]:
            exact = PropsSI(out, 'T', T, 'Q', q, 'NitrousOxide')
            worst[name] = max(worst.get(name, 0.0), abs(props(name, T) / exact - 1))
    print('3. Property table vs CoolProp, worst relative error 235-305 K: '
          + ', '.join(f'{k} {100*v:.4f}%' for k, v in worst.items()))


if __name__ == '__main__':
    main()
