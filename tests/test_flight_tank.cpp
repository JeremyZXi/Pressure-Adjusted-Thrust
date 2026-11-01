#include "FlightSim.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
int main() {
    Engine engine;
    const std::vector<double> times = {0.05, 8.5, 8.51, 9.5};
    const std::vector<double> masses = {44.0, 5.0, 5.0, 0.0};
    engine.setThrustCurve(times, {8133.0,7603.0,1200.0,1000.0});
    engine.setPropellantMassCurve(times, masses);
    engine.setPropellantCGCurve(times, {2.33,2.33,2.33,2.33});
    AeroData aero("Eureka3_Machnum_RASAero2_sim.CSV");
    Rocket rocket(46.7,2.33,aero,engine);
    FlightSim sim(rocket,0.01,15.0,0.0121,101325.0);
    sim.configureAdiabaticTank(0.065,293.15,times,masses,1.3);
    FlightData result = sim.run(12.0);
    bool sawGasTail = false;
    bool sawEmpty = false;
    double previousRatio = 1.0;
    const std::vector<FlightSnapshot>& rows = result.getSnapshots();
    for (std::size_t i=0; i<rows.size(); i++) {
        const FlightSnapshot& row = rows[i];
        if (!std::isfinite(row.altitude) || !std::isfinite(row.tankPressureRatio)
            || row.tankPressureRatio < 0 || row.tankPressureRatio > previousRatio + 1e-10) {
            throw std::runtime_error("Invalid flight pressure history");
        }
        previousRatio = row.tankPressureRatio;
        if (std::abs(row.mass-rocket.getMass(row.time)) > 1e-10) {
            throw std::runtime_error("Flight mass curve changed");
        }
        if (row.time>8.4 && row.time<9.4 && row.tankPressureRatio>0 && row.thrust>0) {
            sawGasTail = true;
        }
        if (row.time>9.5) {
            sawEmpty = true;
            if (row.tankPressureRatio!=0 || row.thrust!=0) {
                throw std::runtime_error("Empty tank produces thrust");
            }
        }
    }
    if (!sawGasTail || !sawEmpty) throw std::runtime_error("Incomplete gas-tail coverage");
    // Positive supplied thrust at an empty endpoint must not retain ambient thrust.
    FlightSim endpoint(rocket,0.05,15.0,0.0121,101325.0);
    endpoint.configureAdiabaticTank(0.065,293.15,times,masses);
    FlightData endpointResult = endpoint.run(10.0);
    const std::vector<FlightSnapshot>& endpointRows = endpointResult.getSnapshots();
    for (std::size_t i=0; i<endpointRows.size(); i++) {
        if (endpointRows[i].tankPressureRatio==0 && endpointRows[i].thrust!=0) {
            throw std::runtime_error("Ambient correction survives empty endpoint");
        }
    }
    std::cout << "Original mass curve, gas-tail flight, finite monotone pressure and burnout checks passed\n";
}
