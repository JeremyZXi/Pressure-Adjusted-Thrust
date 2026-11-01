// 1D Flight Simulator
// Author: Marlow Nedelchev

#ifndef FLIGHTSIM_H
#define FLIGHTSIM_H

#include "Rocket.h"
#include "Atmosphere.h"
#include "FlightData.h"
#include "tank.h"
#include <string>
#include <algorithm>
#include <cmath>
#include <stdexcept>

class FlightSim {
private:
    Rocket rocket;
    double dt;
    double railLength;
    double nozzleExitArea;     // Nozzle exit area (m^2), not aerodynamic reference area
    double referencePressure;  // Ambient pressure of the supplied thrust curve (Pa)

    std::vector<double> tankHistoryTimes;
    std::vector<double> tankHistoryRatios;
    Interpolator1D tankHistory;

    double getTankPressureRatio(double time) const {
        if (!tankHistoryTimes.empty()) {
            if (time <= tankHistoryTimes.front()) return 1.0;
            if (time >= tankHistoryTimes.back()) return tankHistoryRatios.back();
            return tankHistory.interpolate(time);
        }
        // No tank correction unless an equilibrium history is configured.
        return 1.0;
    }

    double getAdjustedThrust(double time, double ambientPressure) const {
        const double baseThrust = rocket.getThrust(time);
        // The supplied curve uses positive thrust to identify the active burn.
        if (baseThrust <= 0.0) return 0.0;
        if (!tankHistoryTimes.empty() && time > tankHistoryTimes.back()) {
            throw std::domain_error("Active thrust extends beyond valid tank history");
        }
        // Exploratory assumption: reference-pressure thrust scales with tank pressure.
        // Assumes the baseline curve does not already include this pressure loss.
        const double pressureRatio = getTankPressureRatio(time);
        if (pressureRatio <= 0.0) return 0.0;
        const double thrust = baseThrust * pressureRatio
            + (referencePressure - ambientPressure) * nozzleExitArea;
        return thrust > 0.0 ? thrust : 0.0;
    }
    FlightData flightData;

    static constexpr double g0 = 9.80665;   // Sea level gravitational acceleration, m/s^2
    static constexpr double re = 6371000.0; // Earth radius, m

    struct State {
        double altitude;
        double velocity;
        double acceleration;
    };

    // Calculate acceleration from forces at given state
    double evaluateAcceleration(double time, double altitude, double velocity) const {
        double mass = rocket.getMass(time);

        // Get atmospheric properties
        auto atm = Atmosphere::getAtmosphere(altitude);
        double thrust = getAdjustedThrust(time, atm.pressure);
        double rho = atm.density;
        double speedOfSound = atm.speedOfSound;

        // Calculate Mach number
        double mach = 0.0;
        if (speedOfSound > 0) {
            mach = velocity / speedOfSound;
        }

        // Get aero coefficients (AoA = 0 for 1D vertical flight)
        auto aeroCoeffs = rocket.getAeroCoeffs(0.0, mach);
        double CD = aeroCoeffs.CD;

        // Calculate drag force: Fd = 0.5 * rho * v^2 * A * CD
        double drag = 0.0;
        if (velocity > 0) {
            double q = 0.5 * rho * velocity * velocity;
            drag = q * rocket.getReferenceArea() * CD;
        }

        // Calculate weight
        double weight = mass * g0 * std::pow((re / (re + altitude)), 2);

        // Net force and acceleration
        double netForce = thrust - drag - weight;
        return netForce / mass;
    }

    // RK4 integration step
    State rk4Step(double time, const State& state, double dt_step) const {
        // RK4 coefficients for velocity (k) and altitude (m)
        // k1: acceleration at current state
        double k1_v = state.velocity;
        double k1_a = evaluateAcceleration(time, state.altitude, state.velocity);

        // k2: acceleration at midpoint using k1
        double h2 = state.altitude + k1_v * dt_step * 0.5;
        double v2 = state.velocity + k1_a * dt_step * 0.5;
        double k2_v = v2;
        double k2_a = evaluateAcceleration(time + dt_step * 0.5, h2, v2);

        // k3: acceleration at midpoint using k2
        double h3 = state.altitude + k2_v * dt_step * 0.5;
        double v3 = state.velocity + k2_a * dt_step * 0.5;
        double k3_v = v3;
        double k3_a = evaluateAcceleration(time + dt_step * 0.5, h3, v3);

        // k4: acceleration at endpoint using k3
        double h4 = state.altitude + k3_v * dt_step;
        double v4 = state.velocity + k3_a * dt_step;
        double k4_v = v4;
        double k4_a = evaluateAcceleration(time + dt_step, h4, v4);

        // Integrate: y_next = y + (k1 + 2*k2 + 2*k3 + k4) * dt / 6
        State nextState;
        nextState.altitude = state.altitude + (k1_v + 2*k2_v + 2*k3_v + k4_v) * dt_step / 6.0;
        nextState.velocity = state.velocity + (k1_a + 2*k2_a + 2*k3_a + k4_a) * dt_step / 6.0;
        nextState.acceleration = evaluateAcceleration(time + dt_step, nextState.altitude, nextState.velocity);

        return nextState;
    }

    // Create a flight snapshot from current state
    FlightSnapshot createSnapshot(double time, const State& state,
                                   double thrust, double drag, double mass, double mach,
                                   double q, double rho, const AeroData::AeroCoeffs& coeffs) const {
        FlightSnapshot snap;
        snap.thrust = thrust;
        snap.tankPressureRatio = getTankPressureRatio(time);
        snap.time = time;
        snap.altitude = state.altitude;
        snap.velocity = state.velocity;
        snap.acceleration = state.acceleration;
        snap.mass = mass;
        snap.CG = rocket.getCG(time);
        snap.mach = mach;
        snap.q = q;
        snap.density = rho;
        snap.CD = coeffs.CD;
        snap.CN = coeffs.CN;
        snap.CP = coeffs.CP;
        return snap;
    }

public:
    FlightSim(const Rocket& rocket_, double dt_, double railLength_,
              double nozzleExitArea_, double referencePressure_)
        : rocket(rocket_), dt(dt_), railLength(railLength_),
          nozzleExitArea(nozzleExitArea_), referencePressure(referencePressure_) {
        if (!std::isfinite(dt) || dt <= 0.0 ||
            !std::isfinite(nozzleExitArea) || nozzleExitArea < 0.0 ||
            !std::isfinite(referencePressure) || referencePressure < 0.0) {
            throw std::invalid_argument("Require positive timestep and nonnegative finite area/pressure");
        }
    }

    // Precompute from the supplied N2O mass curve; force evaluations only look up
    // ratios, so RK4 never advances tank state multiple times per timestep.
    // This mode retains the existing approximate pressure-to-thrust scaling.
    void configureAdiabaticTank(double volume, double initialTemperature,
                               const std::vector<double>& times,
                               const std::vector<double>& masses,
                               double gasExponent = 1.3) {
        if (times.size() < 2 || times.size() != masses.size()) {
            throw std::invalid_argument("Tank mass curve requires matching time/mass points");
        }
        for (std::size_t i = 0; i < times.size(); i++) {
            if (!std::isfinite(times[i]) || times[i] < 0.0 ||
                !std::isfinite(masses[i]) || masses[i] < 0.0) {
                throw std::invalid_argument("Invalid tank mass curve value");
            }
            if (i > 0 && (times[i] <= times[i - 1] || masses[i] > masses[i - 1])) {
                throw std::invalid_argument("Tank times must increase and masses must not increase");
            }
            const double rocketMass = rocket.getMass(times[i]) - rocket.getDryMass();
            if (std::abs(rocketMass - masses[i]) > 1e-8 * std::max(1.0, masses[i])) {
                throw std::invalid_argument("Tank mass curve must match engine propellant mass curve");
            }
        }
        tank oxidizerTank(volume, masses.front(), initialTemperature, gasExponent);
        std::vector<double> historyTimes;
        std::vector<double> historyRatios;
        historyTimes.push_back(times.front());
        historyRatios.push_back(1.0);
        // Resolve each original mass-curve segment, including short plateaus.
        for (std::size_t segment = 1; segment < times.size(); segment++) {
            const double duration = times[segment] - times[segment - 1];
            const double count = std::ceil(duration / dt);
            if (!std::isfinite(count) || count > 1000000.0) {
                throw std::invalid_argument("Tank history exceeds one million steps per segment");
            }
            const int steps = static_cast<int>(count);
            for (int step = 1; step <= steps; step++) {
                const double fraction = static_cast<double>(step) / steps;
                const double time = times[segment - 1] + duration * fraction;
                double mass = masses[segment - 1]
                    + fraction * (masses[segment] - masses[segment - 1]);
                if (step == steps) {
                    mass = masses[segment];
                }
                try {
                    oxidizerTank.advanceToMass(mass);
                } catch (const std::exception& error) {
                    throw std::runtime_error("Tank history failed at t=" + std::to_string(time)
                        + " s, remaining N2O=" + std::to_string(mass) + " kg: " + error.what()
                        + ". Check tank inputs, saturation-property bounds and gas exponent.");
                }
                historyTimes.push_back(time);
                historyRatios.push_back(oxidizerTank.getPressureRatio());
            }
        }
        // Commit only after the whole history has been computed successfully.
        tankHistory = Interpolator1D(historyTimes, historyRatios);
        tankHistoryTimes = historyTimes;
        tankHistoryRatios = historyRatios;
    }

    // Run the simulation
    FlightData run(double maxTime = 100.0) {
        flightData.clear();

        State state;
        state.altitude = 0.0;
        state.velocity = 0.0;
        state.acceleration = 0.0;

        double time = 0.0;
        bool onRail = true;

        while (time < maxTime) {
            // Calculate acceleration for current state
            state.acceleration = evaluateAcceleration(time, state.altitude, state.velocity);

            // Get properties for snapshot
            double mass = rocket.getMass(time);
            auto atm = Atmosphere::getAtmosphere(state.altitude);
            double thrust = getAdjustedThrust(time, atm.pressure);
            double rho = atm.density;
            double speedOfSound = atm.speedOfSound;
            double mach = state.velocity / speedOfSound;
            double q = 0.5 * rho * state.velocity * state.velocity;
            auto aeroCoeffs = rocket.getAeroCoeffs(0.0, mach);
            double drag = q * rocket.getReferenceArea() * aeroCoeffs.CD;

            // Create and store snapshot
            FlightSnapshot snap = createSnapshot(time, state, thrust, drag, mass,
                                                 mach, q, rho, aeroCoeffs);
            flightData.addSnapshot(snap);

            // Apogee
            if (state.velocity < 0 && time > 0.1) {
                break;
            }

            // Just to be safe
            if (state.altitude < 0 && time > 0.1) {
                break;
            }

            // Check rail exit
            if (onRail && state.altitude >= railLength) {
                onRail = false;
            }

            // Integrate to next timestep
            state = rk4Step(time, state, dt);
            time += dt;
        }

        return flightData;
    }

    // Get the flight data
    const FlightData& getFlightData() const {
        return flightData;
    }
};

#endif // FLIGHTSIM_H
