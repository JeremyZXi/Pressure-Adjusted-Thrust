// Engine.h
// Engine module for rocket propulsion
// SI Units: m, kg, s, N, Pa, K

#ifndef ENGINE_H
#define ENGINE_H

#include "Interpolator.h"
#include <vector>

// Engine module - simplified for now, will be expanded later
class Engine {
private:
    Interpolator1D thrustCurve;      // time -> thrust (N)
    Interpolator1D propMassCurve;    // time -> propellant mass (kg)
    Interpolator1D propCGCurve;      // time -> propellant CG (m from nose)

public:
    Engine() {}

    // Set thrust curve from data points
    void setThrustCurve(const std::vector<double>& time, const std::vector<double>& thrust) {
        thrustCurve = Interpolator1D(time, thrust);
    }

    // Set thrust curve from an existing Interpolator1D
    void setThrustCurve(const Interpolator1D& interpolator) {
        thrustCurve = interpolator;
    }

    // Set propellant mass curve
    void setPropellantMassCurve(const std::vector<double>& time, const std::vector<double>& mass) {
        propMassCurve = Interpolator1D(time, mass);
    }

    // Set propellant mass curve from an existing Interpolator1D
    void setPropellantMassCurve(const Interpolator1D& interpolator) {
        propMassCurve = interpolator;
    }

    // Set propellant CG curve
    void setPropellantCGCurve(const std::vector<double>& time, const std::vector<double>& cg) {
        propCGCurve = Interpolator1D(time, cg);
    }

    // Set propellant CG curve from an existing Interpolator1D
    void setPropellantCGCurve(const Interpolator1D& interpolator) {
        propCGCurve = interpolator;
    }

    // Get thrust at given time
    double getThrust(double time) const {
        return thrustCurve.interpolate(time);
    }

    // Get propellant mass at given time
    double getPropellantMass(double time) const {
        return propMassCurve.interpolate(time);
    }

    // Get propellant CG at given time
    double getPropellantCG(double time) const {
        return propCGCurve.interpolate(time);
    }
};

#endif // ENGINE_H
