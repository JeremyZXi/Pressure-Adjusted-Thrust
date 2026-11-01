#ifndef NITROUS_PROPERTIES_H
#define NITROUS_PROPERTIES_H
#include <cmath>
#include <stdexcept>
#include "nitrous_saturation_table.h"

struct NitrousProperties {
    double pressure;       // Pa
    double liquidDensity;  // kg/m^3
    double vaporDensity;   // kg/m^3
    double liquidEnergy;   // J/kg
    double vaporEnergy;    // J/kg
    double liquidEnthalpy; // J/kg, same reference as internal energy
};

// Linear interpolation of CoolProp 8.0.0 HEOS saturation data.
// Deliberately restricted to 183-305 K, away from the critical point.
inline NitrousProperties nitrousProperties(double temperature) {
    if (!std::isfinite(temperature) || temperature < 183.0 || temperature > 305.0) {
        throw std::out_of_range("N2O saturation table supports 183 to 305 K");
    }
    int index = static_cast<int>(std::floor((temperature - 183.0) / 0.1));
    if (index >= 1220) {
        index = 1219;
    }
    double fraction = (temperature - nitrousSaturationTable[index][0])
        / (nitrousSaturationTable[index + 1][0] - nitrousSaturationTable[index][0]);
    double values[5];
    for (int column = 1; column < 6; column++) {
        values[column - 1] = nitrousSaturationTable[index][column]
            + fraction * (nitrousSaturationTable[index + 1][column]
                          - nitrousSaturationTable[index][column]);
    }
    NitrousProperties p;
    p.pressure = values[0];
    p.liquidDensity = values[1];
    p.vaporDensity = values[2];
    p.liquidEnergy = values[3];
    p.vaporEnergy = values[4];
    p.liquidEnthalpy = p.liquidEnergy + p.pressure / p.liquidDensity;
    return p;
}
#endif
