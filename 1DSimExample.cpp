// 1Dsim.cpp
// 1D flight sim framework - Example implementation
// Author: Marlow Nedelchev
// SI Units: m, kg, s, N, Pa, K

#include <iostream>
#include <fstream>
#include "FlightSim.h"
//#include "TradeStudyExample.h"

int main() {
    std::cout << "1D Flight Simulator - Example" << std::endl;
    std::cout << "================================\n" << std::endl;

    try {
        // atmosphere pressure-correction inputs.
        // Adopted E3 sizing-model exit diameter (not verified).
        const double nozzleExitDiameter = 0.1241;  // m
        const double nozzleExitArea = 3.141592653589793
            * nozzleExitDiameter * nozzleExitDiameter / 4.0;  // m^2
        // Assumes the original thrust curve is referenced to standard sea level.
        const double referencePressure = 101325.0;  // Pa




        const bool useAdiabaticTank =true; //indicator to disable/enable pressure tank
        //Jeremy Zhang: I totally made these up
        const double tankVolume = 65.0 / 1000.0;    // m^3
        const double initialTankTemperature = 20.0 + 273.15; // K
        const double gasExponent = 1.3; // Assumed constant gas-tail exponent (>1)

        // E3 example flight
        std::cout << "Setting up example rocket..." << std::endl;

        Engine engine;

        std::vector<double> thrustTime = {0.05, 8.5, 8.51, 9.5};
        std::vector<double> thrustValues = {8133.0, 7603.0, 1200.0, 1000.0};
        engine.setThrustCurve(thrustTime, thrustValues);

        std::vector<double> propMass = {44.0, 5.0, 5.0, 0.0};
        engine.setPropellantMassCurve(thrustTime, propMass);

        std::vector<double> propCG = {2.33, 2.33, 2.33, 2.33};
        engine.setPropellantCGCurve(thrustTime, propCG);

        AeroData aeroData("Eureka3_Machnum_RASAero2_sim.CSV");

        // Create rocket: dry mass, CG pos
        Rocket rocket(46.7, 2.33, aeroData, engine);

        std::cout << "Running simulation..." << std::endl;
        FlightSim sim(rocket, 0.01, 15.0, nozzleExitArea, referencePressure);  // 0.01s timestep, 15m rail
        if (useAdiabaticTank) {
            sim.configureAdiabaticTank(tankVolume, initialTankTemperature,
                                       thrustTime, propMass, gasExponent);
        }
        FlightData data = sim.run();

        data.exportToCSV("flight_data.csv");

        // Print summary
        std::cout << "\nSimulation Results:" << std::endl;
        std::cout << "  Apogee: " << data.getApogee() << " m" << std::endl;
        std::cout << "  Max Velocity: " << data.getMaxVelocity() << " m/s" << std::endl;
        std::cout << "  Max Mach: " << data.getMaxMach() << std::endl;
        std::cout << "  Flight data exported to: flight_data.csv" << std::endl;

        std::cout << "\nExample completed successfully!" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
