// Author: Jeremy Zhang
// Rigid adiabatic N2O tank: equilibrium liquid withdrawal, followed by
// an approximate constant-exponent gas blowdown anchored at liquid depletion.
#ifndef INC_1D_FRAMEWORK_TANK_H
#define INC_1D_FRAMEWORK_TANK_H

#include <cmath>
#include <stdexcept>

#include "nitrous_properties.h"

class tank {
public:
    // Supply consistent N2O saturation data over [minimumT, maximumT], in K.
    // The interval must lie strictly below the critical temperature.

    // gasExponent is an assumed polytropic exponent, not real-gas N2O data.
    tank(double volume, double initialMass, double initialTemperature,
         double gasExponent = 1.3) {
        initialize(volume, initialMass, initialTemperature,
                   nitrousProperties, 183.0, 305.0, gasExponent);
    }

    // Optional custom saturation data provider.
    tank(double volume, double initialMass, double initialTemperature,
         NitrousProperties (*properties)(double), double minimumT, double maximumT,
         double gasExponent = 1.3) {
        initialize(volume, initialMass, initialTemperature,
                   properties, minimumT, maximumT, gasExponent);
    }

private:
    void initialize(double volume, double initialMass, double initialTemperature,
                    NitrousProperties (*properties)(double),
                    double minimumT, double maximumT, double gasExponent) {
        volume_ = volume;
        mass_ = initialMass;
        temperature_ = initialTemperature;
        properties_ = properties;
        minimumT_ = minimumT;
        maximumT_ = maximumT;
        gasExponent_ = gasExponent;
        gasOnly_ = false;
        empty_ = false;
        gasStartMass_ = 0.0;
        gasStartTemperature_ = 0.0;
        gasStartPressure_ = 0.0;
        gasCv_ = 0.0;
        gasEnergyOffset_ = 0.0;
        if (!positive(volume_) || !positive(mass_) || !properties_ ||
            !std::isfinite(gasExponent_) || gasExponent_ <= 1.0 ||
            !positive(minimumT_) || !positive(maximumT_) || minimumT_ >= maximumT_ ||
            !std::isfinite(temperature_) || temperature_ < minimumT_ ||
            temperature_ > maximumT_)
            throw std::invalid_argument("Invalid tank inputs, property temperature bounds or gas exponent (must exceed 1)");
        const NitrousProperties p = propertiesAt(temperature_);
        quality_ = qualityAt(mass_, p);
        requireTwoPhase(quality_);
        energy_ = equilibriumEnergy(mass_, p);
        if (!std::isfinite(energy_))
            throw std::invalid_argument("Nonfinite initial tank energy");
        initialPressure_ = p.pressure;
        pressure_ = initialPressure_;
        maxMassStep_ = initialMass * 0.001;
        if (quality_ == 1.0) {
            startGasDischarge();
        }
    }

public:
    // Drive this with the unchanged N2O mass curve. Zero mass is allowed.
    // Work on a copy so a failed update leaves the original tank untouched.
    void advanceToMass(double remainingMass) {
        if (!std::isfinite(remainingMass) || remainingMass < 0.0 || remainingMass > mass_) {
            throw std::invalid_argument("Remaining N2O mass must be finite, nonnegative and nonincreasing");
        }
        if (remainingMass == mass_) {
            return;
        }
        tank next = *this;
        next.advanceState(remainingMass);
        *this = next;
    }

    //accessor method
    double getPressure() const { return pressure_; }
    double pressureRatio() const { return pressure_ / initialPressure_; }
    double getPressureRatio() const { return pressureRatio(); }
    double getTemperature() const {
        if (empty_) {
            throw std::domain_error("An empty tank has no fluid temperature");
        }
        return temperature_;
    }
    double getMass() const { return mass_; }
    double getInternalEnergy() const { return energy_; }
    double getVaporFraction() const {
        if (empty_) {
            throw std::domain_error("An empty tank has no vapor fraction");
        }
        return quality_;
    }
    double getLiquidMass() const { return mass_ * (1.0 - quality_); }
    double getVaporMass() const { return mass_ * quality_; }
    bool isGasOnly() const { return gasOnly_ && !empty_; }
    bool isEmpty() const { return empty_; }
    double getGasExponent() const { return gasExponent_; }
    double getGasStartMass() const { return gasStartMass_; }
    double getGasStartTemperature() const { return gasStartTemperature_; }
    double getGasStartPressure() const { return gasStartPressure_; }
    // Outlet enthalpy for conservation diagnostics, with the common energy reference.
    double getOutletEnthalpy() const {
        if (empty_) {
            throw std::domain_error("An empty tank has no outlet enthalpy");
        }
        if (gasOnly_) {
            return gasEnergyOffset_ + gasExponent_ * gasCv_ * temperature_;
        }
        return propertiesAt(temperature_).liquidEnthalpy;
    }

private:
    double volume_, mass_, temperature_;
    // Pointer to an ordinary function: NitrousProperties functionName(double).
    NitrousProperties (*properties_)(double);
    double minimumT_, maximumT_;
    double quality_, energy_, pressure_, initialPressure_, maxMassStep_;
    bool gasOnly_, empty_;
    double gasExponent_;
    double gasStartMass_, gasStartTemperature_, gasStartPressure_;
    double gasCv_, gasEnergyOffset_;

    void advanceState(double remainingMass) {
        if (gasOnly_) {
            advanceGas(remainingMass);
            return;
        }
        const double startMass = mass_;
        const int steps = static_cast<int>(std::ceil((mass_ - remainingMass) / maxMassStep_));
        for (int i = 1; i <= steps; i++) {
            double nextMass;
            if (i == steps) {
                nextMass = remainingMass;
            } else {
                nextMass = startMass - (startMass - remainingMass) * i / steps;
            }
            double newTemperature = 0.0;
            double newQuality = 2.0;
            const double oldEnthalpy = propertiesAt(temperature_).liquidEnthalpy;
            if (nextMass > 0.0) {
                newTemperature = solveTemperature(nextMass, energy_, mass_ - nextMass, oldEnthalpy);
                newQuality = qualityAt(nextMass, propertiesAt(newTemperature));
            }
            if (newQuality >= 1.0) {
                // Resolve the last-liquid event inside this substep, then discharge
                // the remaining requested mass as gas. Do not clamp vapor quality.
                finishLiquidDischarge(nextMass, oldEnthalpy);
                advanceGas(remainingMass);
                return;
            }
            requireTwoPhase(newQuality);
            mass_ = nextMass;
            temperature_ = newTemperature;
            quality_ = newQuality;
            const NitrousProperties p = propertiesAt(temperature_);
            energy_ = equilibriumEnergy(mass_, p);
            pressure_ = p.pressure;
        }
    }

    void finishLiquidDischarge(double lowerMass, double oldEnthalpy) {
        // Bisection in remaining mass: locate x=1 using the same liquid energy
        // balance as the preceding steps. The starting state remains fixed.
        double lo = lowerMass;
        double hi = mass_;
        double eventMass = mass_;
        double eventTemperature = temperature_;
        bool converged = false;
        for (int i = 0; i < 100; i++) {
            eventMass = lo + 0.5 * (hi - lo);
            eventTemperature = solveTemperature(eventMass, energy_, mass_ - eventMass, oldEnthalpy);
            const double x = qualityAt(eventMass, propertiesAt(eventTemperature));
            if (std::abs(x - 1.0) <= 1e-11) {
                converged = true;
                break;
            }
            if (x >= 1.0) {
                lo = eventMass;
            } else {
                hi = eventMass;
            }
        }
        if (!converged) {
            throw std::runtime_error("Liquid-depletion event did not converge");
        }
        const NitrousProperties p = propertiesAt(eventTemperature);
        mass_ = eventMass;
        temperature_ = eventTemperature;
        quality_ = 1.0;
        pressure_ = p.pressure;
        energy_ = mass_ * p.vaporEnergy;
        startGasDischarge();
    }

    void startGasDischarge() {
        gasOnly_ = true;
        gasStartMass_ = mass_;
        gasStartTemperature_ = temperature_;
        gasStartPressure_ = pressure_;
        // Anchor a calorically perfect gas surrogate to the real saturation
        // state. R_eff=P*V/(m*T) includes a fixed depletion-state Z factor.
        // cv and the energy offset make pressure, temperature and U continuous.
        const double effectiveR = pressure_ * volume_ / (mass_ * temperature_);
        gasCv_ = effectiveR / (gasExponent_ - 1.0);
        gasEnergyOffset_ = energy_ / mass_ - gasCv_ * temperature_;
        if (!positive(gasCv_) || !std::isfinite(gasEnergyOffset_)) {
            throw std::runtime_error("Invalid gas-discharge reference state");
        }
    }

    void advanceGas(double remainingMass) {
        // Exact solution of the surrogate rigid-tank balance dU=h_out*dm:
        // T/T*=(m/m*)^(gamma-1), P/P*=(m/m*)^gamma.
        if (remainingMass == 0.0) {
            mass_ = 0.0;
            energy_ = 0.0;
            pressure_ = 0.0;
            quality_ = 0.0; // No mass remains; vapor fraction is undefined.
            empty_ = true;
            return; // Retain last gas T for diagnostics; no empty-tank temperature.
        }
        const double ratio = remainingMass / gasStartMass_;
        temperature_ = gasStartTemperature_ * std::pow(ratio, gasExponent_ - 1.0);
        pressure_ = gasStartPressure_ * std::pow(ratio, gasExponent_);
        mass_ = remainingMass;
        quality_ = 1.0;
        energy_ = mass_ * (gasEnergyOffset_ + gasCv_ * temperature_);
        if (!positive(temperature_) || !positive(pressure_) || !std::isfinite(energy_)) {
            throw std::runtime_error("Nonfinite gas-discharge state");
        }
    }

    static bool positive(double value) {
        return std::isfinite(value) && value > 0.0;
    }

    NitrousProperties propertiesAt(double t) const {
        const NitrousProperties p = properties_(t);
        if (!positive(p.pressure) || !positive(p.liquidDensity) ||
            !positive(p.vaporDensity) || p.liquidDensity <= p.vaporDensity ||
            !std::isfinite(p.liquidEnergy) || !std::isfinite(p.vaporEnergy) ||
            p.vaporEnergy <= p.liquidEnergy || !std::isfinite(p.liquidEnthalpy))
            throw std::invalid_argument("Invalid N2O saturation properties");
        return p;
    }

    double qualityAt(double m, const NitrousProperties& p) const {
        return (volume_ / m - 1.0 / p.liquidDensity) / (1.0 / p.vaporDensity - 1.0 / p.liquidDensity);
    }

    double equilibriumEnergy(double m, const NitrousProperties& p) const {
        const double x = qualityAt(m, p);
        return m * ((1.0 - x) * p.liquidEnergy + x * p.vaporEnergy);
    }

    static void requireTwoPhase(double x) {
        if (!std::isfinite(x) || x < 0.0 || x > 1.0)
            throw std::domain_error("Invalid initial or liquid-phase equilibrium state; "
            "check nitrous mass, tank volume and temperature");
    }

    double energyResidual(double trialTemperature, double newMass,
                          double oldEnergy, double removedMass,
                          double oldEnthalpy) const {
        const NitrousProperties p = propertiesAt(trialTemperature);
        const double result = equilibriumEnergy(newMass, p) - oldEnergy
            + removedMass * 0.5 * (oldEnthalpy + p.liquidEnthalpy);
        if (!std::isfinite(result)) {
            throw std::runtime_error("Nonfinite tank energy residual");
        }
        return result;
    }

    double solveTemperature(double newMass, double oldEnergy,
                            double removedMass, double oldEnthalpy) const {
        // Trial quality may leave [0,1] during root finding; the final state must not.
        double lo = minimumT_, hi = maximumT_;
        double flo = energyResidual(lo, newMass, oldEnergy, removedMass, oldEnthalpy), fhi = energyResidual(hi, newMass, oldEnergy, removedMass, oldEnthalpy);
        if (flo == 0.0) return lo;
        if (fhi == 0.0) return hi;
        if ((flo < 0.0) == (fhi < 0.0))
            throw std::domain_error("No bracketed equilibrium temperature in property range");
        for (int i = 0; i < 100; ++i) {
            const double mid = lo + (hi - lo) * 0.5;
            const double fm = energyResidual(mid, newMass, oldEnergy, removedMass, oldEnthalpy);
            if (fm == 0.0 || hi - lo <= 1e-9) return mid;
            if ((flo < 0.0) == (fm < 0.0)) {
                lo = mid;
                flo = fm;
            } else {
                hi = mid;
            }
        }
        throw std::runtime_error("Tank equilibrium temperature did not converge");
    }
};

#endif // INC_1D_FRAMEWORK_TANK_H
