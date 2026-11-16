r"""Compare the model against measured test data (CSV, e.g. a Grafana export).

    python3 validate.py TEST.csv --time-col Time \
        --col "P_tank_psi=NOS Tank PT" --col "Pc_psi=Chamber PT" \
        --col "thrust_N=Load Cell Sum" --thrust-scale 9.80665 \
        --t-ignition "2026-11-08 13:41:06.200" --set ox_mass=30

Channels: P_tank_psi, P_fuel_psi, Pc_psi (psi) and thrust_N (N). --set overrides any
EngineConfig field (SI units, or a value ending in "psi"). Saves validation.png.
"""
import argparse
import csv
from dataclasses import fields
from datetime import datetime
from pathlib import Path

import numpy as np

from config import PSI, EngineConfig
from static_fire import run

HERE = Path(__file__).resolve().parent
CHANNELS = {  # standard name -> (model output key, model->display scale, axis label)
    'P_tank_psi': ('P_tank', 1 / PSI, 'N2O tank pressure (psi)'),
    'P_fuel_psi': ('P_fuel', 1 / PSI, 'Fuel feed pressure (psi)'),
    'Pc_psi': ('Pc', 1 / PSI, 'Chamber pressure (psi)'),
    'thrust_N': ('thrust', 1.0, 'Thrust (N)'),
}


def parse_time(values):
    try:
        return np.array([float(v) for v in values])
    except ValueError:
        stamps = [datetime.fromisoformat(v.strip().replace('Z', '+00:00')) for v in values]
        return np.array([(s - stamps[0]).total_seconds() for s in stamps]), stamps[0]


def load_test(path, time_col, colmap, thrust_scale):
    with open(path, newline='') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit(f'{path} has no data rows')
    if time_col not in rows[0]:
        raise SystemExit(f'time column {time_col!r} not found; columns are {list(rows[0])}')
    t = parse_time([r[time_col] for r in rows])
    t0_stamp = None
    if isinstance(t, tuple):
        t, t0_stamp = t
    data = {}
    for std in CHANNELS:
        src = colmap.get(std, std)
        if src in rows[0]:
            vals = []
            for r in rows:
                try:
                    vals.append(float(r[src]))
                except (TypeError, ValueError):
                    vals.append(np.nan)      # blank cells in sparse Grafana exports
            data[std] = np.array(vals) * (thrust_scale if std == 'thrust_N' else 1.0)
    if not data:
        raise SystemExit('no recognised channels; use --col NAME=COLUMN (see --help)')
    return t, t0_stamp, data


def apply_overrides(cfg, sets):
    names = {f.name for f in fields(cfg)}
    for item in sets:
        key, _, raw = item.partition('=')
        if key not in names:
            raise SystemExit(f'--set {key}: not an EngineConfig field ({sorted(names)})')
        raw = raw.strip()
        if raw.lower() in ('none', ''):
            val = None
        elif raw.lower().endswith('psi'):
            val = float(raw[:-3]) * PSI
        else:
            try:
                val = float(raw)
            except ValueError:
                val = raw            # e.g. fuel_feed=regulated
        setattr(cfg, key, val)
    return cfg


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('csv')
    p.add_argument('--time-col', default='time_s')
    p.add_argument('--col', action='append', default=[], help='STANDARD=COLUMN, e.g. "Pc_psi=Chamber PT"')
    p.add_argument('--t-ignition', default=None, help='main-valve open time in the file time base')
    p.add_argument('--thrust-scale', type=float, default=1.0, help='multiply thrust column (9.80665 for kgf)')
    p.add_argument('--ambient-psi', type=float, default=14.696)
    p.add_argument('--set', action='append', default=[], help='EngineConfig override, e.g. ox_mass=30')
    p.add_argument('--no-show', action='store_true')
    a = p.parse_args()

    colmap = dict(c.split('=', 1) for c in a.col)
    t_test, t0_stamp, test = load_test(a.csv, a.time_col, colmap, a.thrust_scale)
    if a.t_ignition is not None:
        try:
            t_ign = float(a.t_ignition)
        except ValueError:
            if t0_stamp is None:
                raise SystemExit('--t-ignition is a timestamp but the time column is numeric')
            t_ign = (datetime.fromisoformat(a.t_ignition.replace('Z', '+00:00')) - t0_stamp).total_seconds()
    else:
        t_ign = t_test[0]
    t_test = t_test - t_ign

    cfg = apply_overrides(EngineConfig(), a.set)
    sim, reason, _ = run(cfg, ambient_pressure=a.ambient_psi * PSI)
    t_end = sim['t'][-1]
    print(f'Model burn: 0 -> {t_end:.2f} s ({reason}). Comparing over the overlap with the test data.\n')
    print(f'{"channel":<12}{"points":>8}{"mean test":>12}{"mean model":>12}{"bias":>10}{"RMS err":>10}{"RMS %":>8}')

    import matplotlib
    if a.no_show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(test), 1, figsize=(9, 2.8 * len(test)), squeeze=False, constrained_layout=True)
    for ax, (std, meas) in zip(axes[:, 0], test.items()):
        key, scale, label = CHANNELS[std]
        model = sim[key] * scale
        win = (t_test >= 0) & (t_test <= t_end) & np.isfinite(meas)
        if win.sum() >= 2:
            pred = np.interp(t_test[win], sim['t'], model)
            err = pred - meas[win]
            rms = float(np.sqrt(np.mean(err ** 2)))
            m_test = float(np.mean(meas[win]))
            print(f'{std:<12}{win.sum():>8d}{m_test:>12.1f}{np.mean(pred):>12.1f}{np.mean(err):>10.1f}'
                  f'{rms:>10.1f}{100 * rms / max(abs(m_test), 1e-9):>7.1f}%')
        else:
            print(f'{std:<12}{int(win.sum()):>8d}   (no test points inside the model burn window)')
        ax.plot(t_test, meas, '.', ms=2, alpha=.5, label='test')
        ax.plot(sim['t'], model, lw=2, label='model')
        ax.set_ylabel(label); ax.grid(alpha=.25); ax.legend(loc='best')
    axes[-1, 0].set_xlabel('Time since main valves opened (s)')
    fig.suptitle(f'Model vs test: {Path(a.csv).name}')
    fig.savefig(HERE / 'validation.png', dpi=150)
    print(f'\nSaved {HERE / "validation.png"}')
    if not a.no_show:
        plt.show()


if __name__ == '__main__':
    main()
