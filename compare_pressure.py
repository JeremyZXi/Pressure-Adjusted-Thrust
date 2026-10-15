"""Compare the actual C++ simulator with and without ambient- and tank-pressure corrections.

Run: python3 compare_pressure.py
Requires matplotlib, numpy, and a C++14 compiler (clang++ or g++).
Inputs are read from 1DSimExample.cpp; the baseline disables both pressure corrections.
"""
import argparse
import csv
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def simulate(directory, source, compiler, corrected):
    """Build temporary copies; leave the original source and flight_data.csv alone."""
    if not corrected:
        source, count = re.subn(
            r'(const double nozzleExitArea\s*=\s*)[^;]+;',
            r'\g<1>0.0;', source,
        )
        if count != 1:
            raise ValueError('Expected one nozzleExitArea input in main()')
        source, count = re.subn(
            r'(const double tankPressureDropFraction\s*=\s*)[^;]+;',
            r'\g<1>0.0;', source,
        )
        if count != 1:
            raise ValueError('Expected one tankPressureDropFraction input in main()')
    with tempfile.TemporaryDirectory(prefix='rocket-comparison-') as work:
        work = Path(work)
        cpp = work / 'main.cpp'
        cpp.write_text(source)
        executable = work / 'flight'
        subprocess.run([compiler, '-std=c++14', '-O2', '-I', str(directory),
                        str(cpp), '-o', str(executable)], check=True)
        aero = 'Eureka3_Machnum_RASAero2_sim.CSV'
        (work / aero).symlink_to(directory / aero)
        result = subprocess.run([str(executable)], cwd=work, check=True,
                                capture_output=True, text=True)
        print(('WITH' if corrected else 'WITHOUT') + ' correction')
        print(result.stdout)
        with (work / 'flight_data.csv').open(newline='') as file:
            rows = list(csv.DictReader(file))
        return {key: [float(row[key]) for row in rows] for key in rows[0]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-show', action='store_true', help='Save without opening a window')
    args = parser.parse_args()
    import matplotlib
    if args.no_show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    directory = Path(__file__).resolve().parent
    compiler = shutil.which('clang++') or shutil.which('g++')
    if compiler is None:
        raise RuntimeError('Install clang++ or g++ to run the C++ simulator')
    source = (directory / '1DSimExample.cpp').read_text()
    baseline = simulate(directory, source, compiler, corrected=False)
    corrected = simulate(directory, source, compiler, corrected=True)

    # Compare at matching times, only while both trajectories are available.
    times = np.asarray(baseline['time'])
    times = times[times <= corrected['time'][-1]]
    fig, axes = plt.subplots(3, 2, figsize=(12, 10), constrained_layout=True)
    for row, (key, label) in enumerate([
        ('thrust', 'Thrust (N)'), ('altitude', 'Altitude (m)'),
        ('velocity', 'Velocity (m/s)'),
    ]):
        left, right = axes[row]
        left.plot(baseline['time'], baseline[key], label='Without correction')
        left.plot(corrected['time'], corrected[key], '--', label='With correction')
        left.set_ylabel(label)
        left.legend()
        delta = (np.interp(times, corrected['time'], corrected[key])
                 - np.interp(times, baseline['time'], baseline[key]))
        right.plot(times, delta, color='tab:green')
        right.axhline(0, color='gray', linewidth=.8)
        right.set_ylabel('Difference in ' + label.lower())
        for axis in (left, right):
            axis.set_xlabel('Time (s)')
            axis.grid(alpha=.25)
        if key == 'thrust':
            last_burn = max(t for t, thrust in zip(corrected['time'], corrected['thrust'])
                            if thrust > 0)
            left.set_xlim(0, last_burn + .5)
            right.set_xlim(0, last_burn + .5)
    axes[0, 0].set_title('Flight histories')
    axes[0, 1].set_title('With correction minus without correction')
    diameter = float(re.search(r'const double nozzleExitDiameter\s*=\s*([^;]+);', source).group(1))
    area = np.pi * diameter**2 / 4.0
    pressure = re.search(r'const double referencePressure\s*=\s*([^;]+);', source).group(1)
    fig.suptitle(f'Atmospheric + tank-pressure corrections | exit area {area:.7f} m² | reference {pressure} Pa\n'
                 'Exploratory linear tank model; example engine inputs', fontsize=13)
    output = directory / 'pressure_comparison.png'
    fig.savefig(output, dpi=170)
    for key, label in [('altitude', 'Apogee (m)'), ('velocity', 'Maximum velocity (m/s)')]:
        before, after = max(baseline[key]), max(corrected[key])
        print(f'{label}: {before:.3f} → {after:.3f} (difference {after-before:+.3f})')
    print(f'Plot saved to {output}')
    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
