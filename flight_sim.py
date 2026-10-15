"""Vertical rocket simulator ported from the C++ framework (SI units).

Run: python3 flight_sim.py
Input and output paths are relative to this script, independent of IDE settings.
"""
from bisect import bisect_left
from dataclasses import dataclass, asdict
from pathlib import Path
import csv
import math


class Interpolator1D:
    def __init__(self, x, y):
        if len(x) != len(y) or len(x) < 2:
            raise ValueError('At least two matching x/y points required')
        if any(b <= a for a, b in zip(x, x[1:])):
            raise ValueError('x values must be strictly increasing')
        self.x, self.y = list(x), list(y)

    def interpolate(self, x):
        # Preserve the C++ behavior, including zero propellant before ignition.
        if x < self.x[0] or x > self.x[-1]:
            return 0.0
        i = max(1, bisect_left(self.x, x))
        fraction = (x - self.x[i - 1]) / (self.x[i] - self.x[i - 1])
        return self.y[i - 1] + fraction * (self.y[i] - self.y[i - 1])


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


class Engine:
    def __init__(self, time, thrust, propellant_mass, propellant_cg):
        self.thrust = Interpolator1D(time, thrust)
        self.propellant_mass = Interpolator1D(time, propellant_mass)
        self.propellant_cg = Interpolator1D(time, propellant_cg)


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


class Rocket:
    def __init__(self, dry_mass, dry_cg, aero_data, engine):
        self.dry_mass, self.dry_cg = dry_mass, dry_cg
        self.aero_data, self.engine = aero_data, engine

    def get_mass(self, time):
        return self.dry_mass + self.engine.propellant_mass.interpolate(time)

    def get_cg(self, time):
        prop_mass = self.engine.propellant_mass.interpolate(time)
        if prop_mass < 1e-6:
            return self.dry_cg
        return (self.dry_mass * self.dry_cg + prop_mass * self.engine.propellant_cg.interpolate(time)) / self.get_mass(time)


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


@dataclass
class FlightSnapshot:
    time: float
    altitude: float
    velocity: float
    acceleration: float
    mass: float
    CG: float
    mach: float
    q: float
    density: float
    CD: float
    CN: float
    CP: float


class FlightData:
    def __init__(self):
        self.snapshots = []

    def export_to_csv(self, filename):
        with open(filename, 'w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=list(FlightSnapshot.__dataclass_fields__))
            writer.writeheader()
            writer.writerows(asdict(snapshot) for snapshot in self.snapshots)

    def get_apogee(self):
        return max((s.altitude for s in self.snapshots), default=0.0)

    def get_max_velocity(self):
        return max((s.velocity for s in self.snapshots), default=0.0)

    def get_max_mach(self):
        return max((s.mach for s in self.snapshots), default=0.0)


class FlightSim:
    def __init__(self, rocket, dt=0.01, rail_length=10.0):
        if dt <= 0:
            raise ValueError('Timestep must be positive')
        self.rocket, self.dt, self.rail_length = rocket, dt, rail_length
        self.flight_data = FlightData()

    def evaluate_acceleration(self, time, altitude, velocity):
        rocket = self.rocket
        mass = rocket.get_mass(time)
        thrust = rocket.engine.thrust.interpolate(time)
        atm = Atmosphere.get_atmosphere(altitude)
        coeffs = rocket.aero_data.get_coeffs(0.0, velocity / atm.speed_of_sound)
        drag = 0.5 * atm.density * velocity ** 2 * rocket.aero_data.reference_area * coeffs.CD if velocity > 0 else 0.0
        weight = mass * 9.80665 * (6371000.0 / (6371000.0 + altitude)) ** 2
        return (thrust - drag - weight) / mass

    def rk4_step(self, time, altitude, velocity):
        dt = self.dt
        a1 = self.evaluate_acceleration(time, altitude, velocity)
        v2 = velocity + a1 * dt * 0.5
        a2 = self.evaluate_acceleration(time + dt * 0.5, altitude + velocity * dt * 0.5, v2)
        v3 = velocity + a2 * dt * 0.5
        a3 = self.evaluate_acceleration(time + dt * 0.5, altitude + v2 * dt * 0.5, v3)
        v4 = velocity + a3 * dt
        a4 = self.evaluate_acceleration(time + dt, altitude + v3 * dt, v4)
        return (altitude + (velocity + 2 * v2 + 2 * v3 + v4) * dt / 6,
                velocity + (a1 + 2 * a2 + 2 * a3 + a4) * dt / 6)

    def run(self, max_time=100.0):
        self.flight_data = FlightData()
        time, altitude, velocity = 0.0, 0.0, 0.0
        on_rail = True  # As in C++, rail exit is tracked without changing forces.
        while time < max_time:
            acceleration = self.evaluate_acceleration(time, altitude, velocity)
            atm = Atmosphere.get_atmosphere(altitude)
            mach = velocity / atm.speed_of_sound
            coeffs = self.rocket.aero_data.get_coeffs(0.0, mach)
            self.flight_data.snapshots.append(FlightSnapshot(
                time, altitude, velocity, acceleration, self.rocket.get_mass(time),
                self.rocket.get_cg(time), mach, 0.5 * atm.density * velocity ** 2,
                atm.density, coeffs.CD, coeffs.CN, coeffs.CP))
            if time > 0.1 and (velocity < 0 or altitude < 0):
                break
            if on_rail and altitude >= self.rail_length:
                on_rail = False
            altitude, velocity = self.rk4_step(time, altitude, velocity)
            time += self.dt
        return self.flight_data


def main():
    directory = Path(__file__).resolve().parent
    print('1D Flight Simulator - Python Example\n================================')
    time = [0.05, 8.5, 8.51, 9.5]
    engine = Engine(time, [8133.0, 7603.0, 1200.0, 1000.0],
                    [44.0, 5.0, 5.0, 0.0], [2.33] * 4)
    aero = AeroData(directory / 'Eureka3_Machnum_RASAero2_sim.CSV')
    rocket = Rocket(46.7, 2.33, aero, engine)
    data = FlightSim(rocket, dt=0.01, rail_length=15.0).run()
    output = directory / 'flight_data_python.csv'
    data.export_to_csv(output)
    print(f'Apogee: {data.get_apogee():.6f} m')
    print(f'Max Velocity: {data.get_max_velocity():.6f} m/s')
    print(f'Max Mach: {data.get_max_mach():.6f}')
    print(f'Flight data exported to: {output}')


if __name__ == '__main__':
    main()
