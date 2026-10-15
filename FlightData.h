// FlightData.h
// Data structures for storing and exporting flight simulation data
// SI Units: m, kg, s, N, Pa, K

#ifndef FLIGHTDATA_H
#define FLIGHTDATA_H

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

// Snapshot of flight state at a single timestep
struct FlightSnapshot {
    double tankPressureRatio;  // P_tank / P_initial (dimensionless)
    double thrust;  // Pressure-adjusted thrust (N)
    double time;
    double altitude;
    double velocity;
    double acceleration;

    double mass;
    double CG;

    double mach;
    double q;
    double density;

    double CD;
    double CN;
    double CP;

    FlightSnapshot() : tankPressureRatio(1), thrust(0), time(0), altitude(0), velocity(0), acceleration(0),
                       mass(0), CG(0), mach(0), q(0), density(0),
                       CD(0), CN(0), CP(0) {}
};


class FlightData {
private:
    std::vector<FlightSnapshot> snapshots;

public:
    void addSnapshot(const FlightSnapshot& snapshot) {
        snapshots.push_back(snapshot);
    }

    const std::vector<FlightSnapshot>& getSnapshots() const {
        return snapshots;
    }

    void clear() {
        snapshots.clear();
    }

    bool exportToCSV(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
            return false;
        }

        file << "time,altitude,velocity,acceleration,mass,CG,mach,q,density,CD,CN,CP,thrust,tankPressureRatio\n";

        for (const auto& snap : snapshots) {
            file << snap.time << ","
                 << snap.altitude << ","
                 << snap.velocity << ","
                 << snap.acceleration << ","
                 << snap.mass << ","
                 << snap.CG << ","
                 << snap.mach << ","
                 << snap.q << ","
                 << snap.density << ","
                 << snap.CD << ","
                 << snap.CN << ","
                 << snap.CP << ","
                 << snap.thrust << ","
                 << snap.tankPressureRatio << "\n";
        }

        file.close();
        return true;
    }

    double getApogee() const {
        double maxAlt = 0.0;
        for (const auto& snap : snapshots) {
            if (snap.altitude > maxAlt) maxAlt = snap.altitude;
        }
        return maxAlt;
    }

    double getMaxVelocity() const {
        double maxVel = 0.0;
        for (const auto& snap : snapshots) {
            if (snap.velocity > maxVel) maxVel = snap.velocity;
        }
        return maxVel;
    }

    double getMaxMach() const {
        double maxMach = 0.0;
        for (const auto& snap : snapshots) {
            if (snap.mach > maxMach) maxMach = snap.mach;
        }
        return maxMach;
    }
};

#endif 
