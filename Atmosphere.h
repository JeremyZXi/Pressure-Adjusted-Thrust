// Atmosphere.h
// 1976 US Standard Atmosphere model
// SI Units: m, kg, s, N, Pa, K

#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include <cmath>

// 1976 US Standard Atmosphere model
class Atmosphere {
private:
    // Sea level constants
    static constexpr double T0 = 288.15;      // Sea level temperature (K)
    static constexpr double P0 = 101325.0;    // Sea level pressure (Pa)
    static constexpr double rho0 = 1.225;     // Sea level density (kg/m^3)
    static constexpr double g0 = 9.80665;     // Gravitational acceleration (m/s^2)
    static constexpr double R = 287.05;       // Specific gas constant for air (J/(kg*K))
    static constexpr double gamma = 1.4;      // Ratio of specific heats for air

    // Troposphere constants
    static constexpr double L = 0.0065;       // Temperature lapse rate (K/m)
    static constexpr double h_trop = 11000.0; // Tropopause altitude (m)
    static constexpr double T_trop = 216.65;  // Tropopause temperature (K)

public:
    struct AtmosphericData {
        double temperature;    // K
        double pressure;       // Pa
        double density;        // kg/m^3
        double speedOfSound;   // m/s
    };

    // Get atmospheric properties at given altitude (m)
    static AtmosphericData getAtmosphere(double altitude) {
        AtmosphericData data;

        if (altitude < 0.0) altitude = 0.0;  // Clamp to ground level

        if (altitude <= h_trop) {
            // Troposphere (0 - 11 km)
            data.temperature = T0 - L * altitude;
            data.pressure = P0 * std::pow(data.temperature / T0, g0 / (L * R));
            data.density = data.pressure / (R * data.temperature);
        } else {
            // Isothermal layer (11 - 25 km) - simplified for stratosphere
            data.temperature = T_trop;
            double P_trop = P0 * std::pow(T_trop / T0, g0 / (L * R));
            data.pressure = P_trop * std::exp(-g0 * (altitude - h_trop) / (R * T_trop));
            data.density = data.pressure / (R * data.temperature);
        }

        // Speed of sound: a = sqrt(gamma * R * T)
        data.speedOfSound = std::sqrt(gamma * R * data.temperature);

        return data;
    }

    // Get density at altitude (convenience method)
    static double getDensity(double altitude) {
        return getAtmosphere(altitude).density;
    }

    // Get speed of sound at altitude (convenience method)
    static double getSpeedOfSound(double altitude) {
        return getAtmosphere(altitude).speedOfSound;
    }
};

#endif // ATMOSPHERE_H
