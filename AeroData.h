// Aerodynamic coefficient table with 2D interpolation
// Author: Marlow Nedelchev

#ifndef AERODATA_H
#define AERODATA_H

#include "Interpolator.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>

// Aerodynamic coefficient table with 2D interpolation
// CSV requires: AOA, Mach Num, CD, CN, CP
class AeroData {
private:
    Interpolator2D CD_interp;
    Interpolator2D CP_interp;
    Interpolator2D CN_interp;
    double referenceArea;

    // Parse CSV and create interpolators, returns true on success
    bool loadFromCSV(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open aero file " << filename << std::endl;
            return false;
        }

        std::vector<double> aoa_raw, mach_raw, cd_raw, cp_raw, cn_raw;
        std::string line;

        // Skip header
        std::getline(file, line);

        // Read all data
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string cell;
            std::vector<double> row;

            while (std::getline(ss, cell, ',')) {
                row.push_back(std::stod(cell));
            }

            // Set for RasAero table currently
            if (row.size() >= 12) {
                mach_raw.push_back(row[0]);
                aoa_raw.push_back(row[1]);
                cd_raw.push_back(row[2]);
                cn_raw.push_back(row[8]);
                cp_raw.push_back(row[12]);
            }
        }

        file.close();

        if (aoa_raw.empty()) {
            std::cerr << "Error: No data read from aero file" << std::endl;
            return false;
        }

        // Extract unique AoA and Mach values
        std::vector<double> aoa_unique, mach_unique;
        for (double val : aoa_raw) {
            if (std::find(aoa_unique.begin(), aoa_unique.end(), val) == aoa_unique.end()) {
                aoa_unique.push_back(val);
            }
        }
        for (double val : mach_raw) {
            if (std::find(mach_unique.begin(), mach_unique.end(), val) == mach_unique.end()) {
                mach_unique.push_back(val);
            }
        }

        std::sort(aoa_unique.begin(), aoa_unique.end());
        std::sort(mach_unique.begin(), mach_unique.end());

        // Organize data into 2D grids [AoA][Mach]
        size_t n_aoa = aoa_unique.size();
        size_t n_mach = mach_unique.size();

        std::vector<std::vector<double>> cd_grid(n_aoa, std::vector<double>(n_mach, 0.0));
        std::vector<std::vector<double>> cp_grid(n_aoa, std::vector<double>(n_mach, 0.0));
        std::vector<std::vector<double>> cn_grid(n_aoa, std::vector<double>(n_mach, 0.0));

        for (size_t i = 0; i < aoa_raw.size(); ++i) {
            auto aoa_it = std::find(aoa_unique.begin(), aoa_unique.end(), aoa_raw[i]);
            auto mach_it = std::find(mach_unique.begin(), mach_unique.end(), mach_raw[i]);

            size_t aoa_idx = aoa_it - aoa_unique.begin();
            size_t mach_idx = mach_it - mach_unique.begin();

            cd_grid[aoa_idx][mach_idx] = cd_raw[i];
            cp_grid[aoa_idx][mach_idx] = cp_raw[i];
            cn_grid[aoa_idx][mach_idx] = cn_raw[i];
        }

        // Create interpolators
        CD_interp = Interpolator2D(aoa_unique, mach_unique, cd_grid);
        CP_interp = Interpolator2D(aoa_unique, mach_unique, cp_grid);
        CN_interp = Interpolator2D(aoa_unique, mach_unique, cn_grid);

        return true;
    }

public:
    AeroData() : referenceArea(0.0) {}

    AeroData(const std::string& csvFilename, double refArea = 0.03243) {
        // Default reference area for 8" diameter
        referenceArea = refArea;
        if (!loadFromCSV(csvFilename)) {
            throw std::runtime_error("Failed to load aero data from " + csvFilename);
        }
    }

    // Get coefficients at given AoA and Mach number
    struct AeroCoeffs {
        double CD;
        double CN;
        double CP; 
    };

    AeroCoeffs getCoeffs(double aoa_deg, double mach) const {
        AeroCoeffs coeffs;
        coeffs.CD = CD_interp.interpolate(aoa_deg, mach);
        coeffs.CN = CN_interp.interpolate(aoa_deg, mach);
        coeffs.CP = CP_interp.interpolate(aoa_deg, mach);
        return coeffs;
    }

    double getReferenceArea() const { return referenceArea; }
};

#endif 
