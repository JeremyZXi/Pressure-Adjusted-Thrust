"""Saturated nitrous oxide (N2O) properties: generate a table once, then interpolate.

This file uses a table to interpolate because calling CoolProp for every time step
is super slow

Regenerate the table(requires CoolProp tho):
    python3 n2o_props.py
Everything else only needs numpy.

The CoolProps library has all of the fluid property info necessary

Units (SI): T [K], P [Pa], rho [kg/m^3], u and h [J/kg]
  u = specific internal energy, h = specific enthalpy = u + P/rho
  suffix _l = saturated liquid, _v = saturated vapor
"""
from pathlib import Path

import numpy as np

TABLE_PATH = Path(__file__).resolve().parent / 'n2o_saturation.csv'
COLUMNS = ['T', 'P', 'rho_l', 'rho_v', 'u_l', 'u_v', 'h_l', 'h_v']


def generate(path=TABLE_PATH, T_min=220.0, T_max=308.0, n=881):
    """Write the saturation table. T_max stays below the critical point (309.5 K)."""
    from CoolProp.CoolProp import PropsSI
    fluid = 'NitrousOxide'
    rows = []
    for T in np.linspace(T_min, T_max, n): #gives an array of temps points
        rows.append([
            T,
            PropsSI('P', 'T', T, 'Q', 0, fluid), #vapor pressure
            PropsSI('D', 'T', T, 'Q', 0, fluid), #liquid density 
            PropsSI('D', 'T', T, 'Q', 1, fluid), # vapor density
            PropsSI('U', 'T', T, 'Q', 0, fluid), #liquid internal energy
            PropsSI('U', 'T', T, 'Q', 1, fluid),#vapor internal energy
            PropsSI('H', 'T', T, 'Q', 0, fluid),#liquid enthalpy
            PropsSI('H', 'T', T, 'Q', 1, fluid), #vapor enthalpy
        ])

    #create a csv file with all the properties
    np.savetxt(path, np.array(rows), delimiter=',', header=','.join(COLUMNS),
               comments='', fmt='%.10g')
    print(f'Wrote {len(rows)} rows to {path}')


class SaturationTable:
    """Look up saturated N2O properties at temperature T by linear interpolation.

    Usage:
        props = SaturationTable()
        props('P', 293.15)       # vapor pressure at 20 C, Pa
        props('rho_l', 280.0)    # liquid density at 280 K
    """

    def __init__(self, path=TABLE_PATH):
        data = np.loadtxt(path, delimiter=',', skiprows=1)
        self._columns = {name: data[:, i] for i, name in enumerate(COLUMNS)}
        self.T_min = float(data[0, 0])
        self.T_max = float(data[-1, 0])

    def __call__(self, name, T):
        if not self.T_min <= T <= self.T_max:
            raise ValueError(f'T = {T:.2f} K is outside the table range '
                             f'[{self.T_min:.1f}, {self.T_max:.1f}] K')
        return float(np.interp(T, self._columns['T'], self._columns[name]))

    def T_from_P(self, P):
        """Saturation temperature for a given vapor pressure (inverse lookup)."""
        return float(np.interp(P, self._columns['P'], self._columns['T']))


if __name__ == '__main__':
    generate()
