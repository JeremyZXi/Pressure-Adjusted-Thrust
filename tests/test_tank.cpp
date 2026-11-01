#include "tank.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

// Explicit checks remain active in Release builds (unlike assert).
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
NitrousProperties synthetic(double t) {
    double pressure = 1e6 * std::exp((t - 280.0) / 30.0);
    NitrousProperties p = {pressure, 800.0, pressure / (180.0 * t),
                          2000.0 * t, 2000.0 * t + 200000.0,
                          2000.0 * t + pressure / 800.0};
    return p;
}
int main() {
    tank fixture(0.01, 5, 280, synthetic, 200, 300);
    double oldEnergy = fixture.getInternalEnergy();
    double oldH = fixture.getOutletEnthalpy();
    fixture.advanceToMass(4.999);
    double newH = fixture.getOutletEnthalpy();
    check(std::abs(fixture.getInternalEnergy() - oldEnergy + 0.001*(oldH+newH)/2) < 0.01,
          "Liquid energy balance failed");
    tank real(0.065, 44, 293.15); // Test fixture only; not measured vehicle inputs.
    check(std::abs(real.getPressure() - 5052509.28) < 5.0, "Initial property mismatch");
    real.advanceToMass(30);
    check(!real.isGasOnly() && !real.isEmpty(), "Early gas transition");
    NitrousProperties p = nitrousProperties(real.getTemperature());
    check(std::abs(real.getLiquidMass()/p.liquidDensity
        + real.getVaporMass()/p.vaporDensity - 0.065) < 1e-12, "Volume balance failed");
    real.advanceToMass(5.6);
    check(real.isGasOnly(), "No gas transition");
    check(real.getLiquidMass() == 0.0 && real.getVaporMass() == 5.6, "Gas mass lost");
    const double eventMass = real.getGasStartMass();
    const double eventT = real.getGasStartTemperature();
    const double eventP = real.getGasStartPressure();
    p = nitrousProperties(eventT);
    check(std::abs(eventMass / p.vaporDensity - 0.065) < 1e-10, "Transition volume discontinuity");
    check(eventP == p.pressure, "Transition pressure discontinuity");
    const double gamma = real.getGasExponent();
    const double cv = eventP * 0.065 / (eventMass * eventT * (gamma - 1.0));
    const double offset = p.vaporEnergy - cv * eventT;
    real.advanceToMass(2.0);
    check(std::abs(real.getPressure()/eventP - std::pow(2.0/eventMass,gamma)) < 1e-12,
          "Gas pressure law failed");
    check(std::abs(real.getTemperature()/eventT - std::pow(2.0/eventMass,gamma-1)) < 1e-12,
          "Gas temperature law failed");
    check(std::abs(real.getInternalEnergy()-2.0*(offset+cv*real.getTemperature())) < 1e-7,
          "Gas energy reference discontinuity");
    // Finite-difference check of dU/dm=h_out for the gas model.
    tank plus = real;
    tank minus = real;
    plus.advanceToMass(1.9 + 1e-5);
    minus.advanceToMass(1.9 - 1e-5);
    tank middle = real;
    middle.advanceToMass(1.9);
    const double numericalH = (plus.getInternalEnergy()-minus.getInternalEnergy()) / 2e-5;
    check(std::abs(numericalH-middle.getOutletEnthalpy()) < 0.01, "Gas energy balance failed");
    // Refining the caller's mass increments must preserve the transition.
    tank fine(0.065,44,293.15);
    for (int i=1; i<=3840; i++) fine.advanceToMass(44.0-i*0.01);
    check(std::abs(fine.getGasStartMass()-eventMass)/eventMass < 1e-5,
          "Transition failed step refinement");
    check(std::abs(fine.getPressure()-eventP*std::pow(5.6/eventMass,gamma))/eventP < 1e-5,
          "Gas pressure failed step refinement");
    bool rejected = false;
    const double savedEnergy = real.getInternalEnergy();
    try { real.advanceToMass(3.0); } catch (const std::invalid_argument&) { rejected=true; }
    check(rejected && real.getMass()==2.0 && real.getInternalEnergy()==savedEnergy,
          "Failed update modified state");
    rejected = false;
    try { real.advanceToMass(-0.1); } catch (const std::invalid_argument&) { rejected=true; }
    check(rejected, "Negative mass accepted");
    real.advanceToMass(0.0);
    check(real.isEmpty() && !real.isGasOnly(), "Empty state missing");
    check(real.getMass()==0 && real.getPressure()==0 && real.getPressureRatio()==0
        && real.getInternalEnergy()==0 && real.getLiquidMass()==0 && real.getVaporMass()==0,
        "Empty state not zero");
    real.advanceToMass(0.0); // Repeated empty update is valid.
    rejected=false;
    try { real.getTemperature(); } catch (const std::domain_error&) { rejected=true; }
    check(rejected, "Empty tank reported a physical temperature");
    tank direct(0.065,44,293.15);
    direct.advanceToMass(0.0);
    check(direct.isEmpty(), "Single update to empty failed");
    rejected=false;
    try { tank invalid(0.065,44,293.15,1.0); }
    catch (const std::invalid_argument&) { rejected=true; }
    check(rejected, "Invalid gas exponent accepted");
    rejected=false;
    try { nitrousProperties(std::numeric_limits<double>::quiet_NaN()); }
    catch (const std::out_of_range&) { rejected=true; }
    check(rejected, "NaN property input accepted");
    nitrousProperties(183.0);
    nitrousProperties(305.0);
    std::cout << "Liquid/gas energy, transition continuity, refinement, mass, empty state and input checks passed\n";
    std::cout << "Liquid depletion: " << eventMass << " kg at " << eventT << " K\n";
}
