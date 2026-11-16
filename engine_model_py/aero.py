"""Drag table (AeroData) and 1976 standard atmosphere, used by flight.py.
Copied from the team's 1D framework (flight_sim.py) so this folder runs on its own.
"""
from bisect import bisect_left
from dataclasses import dataclass
import csv
import math

class Interpolator2D:
    def __init__(self, x, y, z):
        if len(x) < 2 or len(y) < 2:
            raise ValueError('At least a 2 by 2 grid required')
        if len(z) != len(x) or any(len(row) != len(y) for row in z):
            raise ValueError('Grid dimensions do not match axes')
        self.x, self.y, self.z = x, y, z

    def interpolate(self, x, y):
        x = max(self.x[0], min(x, self.x[-1]))
        y = max(self.y[0], min(y, self.y[-1]))
        i, j = max(1, bisect_left(self.x, x)), max(1, bisect_left(self.y, y))
        wx = (x - self.x[i - 1]) / (self.x[i] - self.x[i - 1])
        wy = (y - self.y[j - 1]) / (self.y[j] - self.y[j - 1])
        r0 = self.z[i - 1][j - 1] * (1 - wx) + self.z[i][j - 1] * wx
        r1 = self.z[i - 1][j] * (1 - wx) + self.z[i][j] * wx
        return r0 * (1 - wy) + r1 * wy


@dataclass
class AeroCoeffs:
    CD: float
    CN: float
    CP: float


class AeroData:
    def __init__(self, filename, reference_area=0.03243):
        self.reference_area = reference_area
        with open(filename, newline='') as file:
            rows = list(csv.DictReader(file))
        if not rows:
            raise ValueError('No aerodynamic data found')
        aoas = sorted({float(row['Alpha']) for row in rows})
        machs = sorted({float(row['Mach']) for row in rows})
        aoa_index = {value: i for i, value in enumerate(aoas)}
        mach_index = {value: i for i, value in enumerate(machs)}
        self.curves = {}
        for name in ('CD', 'CN', 'CP'):
            grid = [[0.0 for _ in machs] for _ in aoas]
            for row in rows:
                grid[aoa_index[float(row['Alpha'])]][mach_index[float(row['Mach'])]] = float(row[name])
            self.curves[name] = Interpolator2D(aoas, machs, grid)

    def get_coeffs(self, aoa, mach):
        return AeroCoeffs(*(self.curves[name].interpolate(aoa, mach) for name in ('CD', 'CN', 'CP')))


@dataclass
class AtmosphericData:
    temperature: float
    pressure: float
    density: float
    speed_of_sound: float


class Atmosphere:
    @staticmethod
    def get_atmosphere(altitude):
        altitude = max(0.0, altitude)
        if altitude <= 11000.0:
            temperature = 288.15 - 0.0065 * altitude
            pressure = 101325.0 * (temperature / 288.15) ** (9.80665 / (0.0065 * 287.05))
        else:
            temperature = 216.65
            p_trop = 101325.0 * (216.65 / 288.15) ** (9.80665 / (0.0065 * 287.05))
            pressure = p_trop * math.exp(-9.80665 * (altitude - 11000.0) / (287.05 * temperature))
        return AtmosphericData(temperature, pressure, pressure / (287.05 * temperature), math.sqrt(1.4 * 287.05 * temperature))
