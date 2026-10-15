// Custom 1D and 2D interpolators
// Author: Marlow Nedelchev

#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include <vector>
#include <algorithm>
#include <stdexcept>

// 1D linear interpolator (thrust curves, shifting cg, etc)
class Interpolator1D {
private:
    std::vector<double> x_data;
    std::vector<double> y_data;

public:
    Interpolator1D() {}

    Interpolator1D(const std::vector<double>& x, const std::vector<double>& y) {
        if (x.size() != y.size() || x.size() < 2) {
            throw std::invalid_argument("Interpolator1D requires at least 2 points with matching x and y sizes");
        }
        x_data = x;
        y_data = y;
    }

    // Linear interpolation returning 0 outside bounds
    // Useful for thrust curves that naturally go to 0 outside the defined range
    double interpolate(double x) const {
        if (x_data.empty()) return 0.0;

        // Return 0 outside bounds
        if (x < x_data.front() || x > x_data.back()) return 0.0;

        // Binary search for the interval
        auto it = std::lower_bound(x_data.begin(), x_data.end(), x);
        size_t idx = it - x_data.begin();

        // Handle exact match
        if (idx > 0 && x_data[idx] == x) return y_data[idx];

        // Linear interpolation
        if (idx == 0) idx = 1;
        double x0 = x_data[idx - 1];
        double x1 = x_data[idx];
        double y0 = y_data[idx - 1];
        double y1 = y_data[idx];

        return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
    }

    bool isEmpty() const { return x_data.empty(); }
};


// 2D bilinear interpolator (aero tables)
class Interpolator2D {
private:
    std::vector<double> x_data; // 1D
    std::vector<double> y_data; // 1D
    std::vector<std::vector<double>> z_data; // 2D 

public:
    Interpolator2D() {}

    Interpolator2D(const std::vector<double>& x, const std::vector<double>& y,
                   const std::vector<std::vector<double>>& z) {
        if (x.size() < 2 || y.size() < 2) {
            throw std::invalid_argument("Interpolator2D requires at least 2x2 grid");
        }
        if (z.size() != x.size()) {
            throw std::invalid_argument("z_data rows must match x_data size");
        }
        for (const auto& row : z) {
            if (row.size() != y.size()) {
                throw std::invalid_argument("All z_data rows must match y_data size");
            }
        }

        x_data = x;
        y_data = y;
        z_data = z;
    }

    // Bilinear interpolation with clamping
    double interpolate(double x, double y) const {
        if (x_data.empty() || y_data.empty()) return 0.0;

        // Clamp x to bounds
        double x_clamped = std::max(x_data.front(), std::min(x, x_data.back()));
        double y_clamped = std::max(y_data.front(), std::min(y, y_data.back()));

        // Find x interval
        auto x_it = std::lower_bound(x_data.begin(), x_data.end(), x_clamped);
        size_t i1 = x_it - x_data.begin();
        if (i1 == 0) i1 = 1;
        size_t i0 = i1 - 1;

        // Find y interval
        auto y_it = std::lower_bound(y_data.begin(), y_data.end(), y_clamped);
        size_t j1 = y_it - y_data.begin();
        if (j1 == 0) j1 = 1;
        size_t j0 = j1 - 1;

        // Bilinear interpolation
        double x0 = x_data[i0], x1 = x_data[i1];
        double y0 = y_data[j0], y1 = y_data[j1];

        double Q00 = z_data[i0][j0];
        double Q10 = z_data[i1][j0];
        double Q01 = z_data[i0][j1];
        double Q11 = z_data[i1][j1];

        double wx = (x_clamped - x0) / (x1 - x0);
        double wy = (y_clamped - y0) / (y1 - y0);

        double R0 = Q00 * (1 - wx) + Q10 * wx;
        double R1 = Q01 * (1 - wx) + Q11 * wx;

        return R0 * (1 - wy) + R1 * wy;
    }

    bool isEmpty() const { return x_data.empty() || y_data.empty(); }
};

#endif 
