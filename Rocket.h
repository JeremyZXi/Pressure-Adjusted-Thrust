// Rocket class containing all vehicle properties
// Author: Marlow Nedelchev

#ifndef ROCKET_H
#define ROCKET_H

#include "AeroData.h"
#include "Engine.h"

// Rocket class containing all vehicle properties
class Rocket {
private:
    double dryMass;
    double dryCG;
    AeroData aeroData;
    Engine engine;

public:
    Rocket() : dryMass(0), dryCG(0) {}

    Rocket(double dryMass_, double dryCG_, const AeroData& aeroData_, const Engine& engine_)
        : dryMass(dryMass_), dryCG(dryCG_), aeroData(aeroData_), engine(engine_) {}


    double getMass(double time) const {
        return dryMass + engine.getPropellantMass(time);
    }

    double getCG(double time) const {
        double propMass = engine.getPropellantMass(time);
        double propCG = engine.getPropellantCG(time);

        if (propMass < 1e-6) {
            return dryCG;
        }

        double totalMass = dryMass + propMass;
        return (dryMass * dryCG + propMass * propCG) / totalMass;
    }

    double getThrust(double time) const {
        return engine.getThrust(time);
    }

    AeroData::AeroCoeffs getAeroCoeffs(double aoa_deg, double mach) const {
        return aeroData.getCoeffs(aoa_deg, mach);
    }

    double getReferenceArea() const {
        return aeroData.getReferenceArea();
    }

    double getDryMass() const { return dryMass; }
    double getDryCG() const { return dryCG; }
};

#endif 
