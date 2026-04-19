#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <cmath>
#include <string>
#include <numeric>

namespace py = pybind11;

inline double dist_sq(const double* a, const double* b) {
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    double dz = a[2] - b[2];
    return dx*dx + dy*dy + dz*dz;
}

int func1(py::array_t<double> atoms_py, 
          py::array_t<double> water_coords_py, 
          py::array_t<double> radiuses_py, 
          double r) {
    auto atoms = atoms_py.unchecked<2>();
    auto waters = water_coords_py.unchecked<2>();
    auto radii = radiuses_py.unchecked<1>();
    size_t n_atoms = atoms.shape(0);
    size_t n_waters = waters.shape(0);
    int total_number = 0;
    #pragma omp parallel for reduction(+:total_number)
    for (int i = 0; i < (int)n_atoms; ++i) {
        const double* pos_i = &atoms(i, 0);
        double rad_i = radii(i);
        std::vector<size_t> valid_waters;
        std::vector<double> dist_i_w;
        for (size_t w = 0; w < n_waters; ++w) {
            double d2 = dist_sq(pos_i, &waters(w, 0));
            double limit = rad_i + r;
            if (d2 < limit * limit) {
                valid_waters.push_back(w);
                dist_i_w.push_back(std::sqrt(d2));
            }
        }
        if (valid_waters.empty()) continue;
        std::vector<size_t> other_indices;
        for (size_t j = 0; j < n_atoms; ++j) {
            if (i == (int)j) continue;
            double d2 = dist_sq(pos_i, &atoms(j, 0));
            double limit = radii(i) + 2 * r + radii(j);
            if (d2 <= limit * limit) {
                other_indices.push_back(j);
            }
        }
        if (other_indices.empty()) continue;
        for (size_t idx = 0; idx < valid_waters.size(); ++idx) {
            size_t w_idx = valid_waters[idx];
            double d_iw = dist_i_w[idx];
            const double* pos_w = &waters(w_idx, 0);
            
            bool all_condition = true;
            for (size_t j : other_indices) {
                double d_wj = std::sqrt(dist_sq(pos_w, &atoms(j, 0)));
                if (!(d_iw * radii(j) < d_wj * rad_i)) {
                    all_condition = false;
                    break; 
                }
            }
            if (all_condition) {
                total_number++;
            }
        }
    }
    return total_number;
}


double func2(
    const std::vector<std::array<double, 3>>& atoms,
    const std::vector<double>& radiuses,
    const std::vector<int>& sel,
    double r = 0,
    int n_points = 1000
) {
    int n = atoms.size();
    std::vector<std::vector<double>> distances(n, std::vector<double>(n));
    std::vector<std::vector<bool>> masks(n, std::vector<bool>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double dx = atoms[i][0] - atoms[j][0];
            double dy = atoms[i][1] - atoms[j][1];
            double dz = atoms[i][2] - atoms[j][2];
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            distances[i][j] = dist;
            double threshold = radiuses[i] + 2*r + radiuses[j];
            if (dist <= threshold) {
                masks[i][j] = true;
            }
        }
    }
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist01(0.0, 1.0);
    auto random_points_in_sphere = [&](int n_points, double radius) {
        std::vector<std::array<double, 3>> points(n_points);
        std::vector<double> norms(n_points);
        for (int i = 0; i < n_points; ++i) {
            double rr = radius * std::cbrt(dist01(rng));
            double theta = 2 * M_PI * dist01(rng);
            double phi = std::acos(1 - 2 * dist01(rng));
            double x = rr * std::sin(phi) * std::cos(theta);
            double y = rr * std::sin(phi) * std::sin(theta);
            double z = rr * std::cos(phi);
            points[i] = {x, y, z};
            norms[i] = rr;
        }
        return std::make_pair(points, norms);
    };
    double total_volume = 0.0;
    for (int idx : sel) {
        const auto& atom = atoms[idx];
        double radius = radiuses[idx];
        const auto& mask_row = masks[idx];
        auto[points, point_norms] = random_points_in_sphere(n_points, radius + r);
        for (auto& p : points) {
            p[0] += atom[0];
            p[1] += atom[1];
            p[2] += atom[2];
        }
        std::vector<std::array<double, 3>> neighbors;
        std::vector<double> neighbor_radiuses;
        for (int j = 0; j < n; ++j) {
            if (mask_row[j]) {
                neighbors.push_back(atoms[j]);
                neighbor_radiuses.push_back(radiuses[j]);
            }
        }
        int count = 0;
        for (int i = 0; i < n_points; ++i) {
            bool all_good = true;
            for (size_t j = 0; j < neighbors.size(); ++j) {
                double dx = points[i][0] - neighbors[j][0];
                double dy = points[i][1] - neighbors[j][1];
                double dz = points[i][2] - neighbors[j][2];
                double d = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (point_norms[i] * neighbor_radiuses[j] >= d * radius) {
                    all_good = false;
                    break;
                }
            }
            if (all_good) ++count;
        }
        double sphere_volume = (4.0 / 3.0) * M_PI * std::pow(radius + r, 3);
        total_volume += sphere_volume * (static_cast<double>(count) / n_points);
    }

PYBIND11_MODULE(module_monte, m) {
    m.doc() = "Monte Carlo volume calculation module";

    m.def("func1", &func1, "Fast water distance calculation",
          py::arg("atoms_py"), 
          py::arg("water_coords_py"), 
          py::arg("radiuses_py"), 
          py::arg("r"));

    m.def("func2", &func2, "Monte Carlo volume calculation",
          py::arg("atoms"),
          py::arg("radiuses"),
          py::arg("sel"),
          py::arg("r") = 0.0,
          py::arg("n_points") = 1000);
}
    return total_volume;
}
