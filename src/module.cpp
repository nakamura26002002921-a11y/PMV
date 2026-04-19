#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <string>
#include <numeric>
#include <omp.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> AABB_tree;

namespace py = pybind11;

using namespace std;
using namespace Eigen;
using Point3D = std::array<double, 3>;
using Region = std::vector<int>;


double complex_volume_estimate(const Matrix<double, 3, 3>& vectors, double radius) {
    auto cone_volume = [&](const vector<Vector3d>& p, double radius) -> double {
        vector<Vector3d> cross(3);
        for (int i = 0; i < 3; i++) cross[i] = p[(i+2)%3].cross(p[i]);
        vector<double> angle(3);
        for (int i = 0; i < 3; i++) angle[i] = acos((-cross[(i+2)%3]).dot(cross[i]) / (cross[(i+2)%3].norm() * cross[i].norm()));
        double sum_angle = angle[0] + angle[1] + angle[2];
        return fabs(((sum_angle - M_PI) * pow(radius, 3)) / 3.0);
    };

    auto pyramid_volume = [&](const Vector3d& p1, const Vector3d& p2, const Vector3d& center,
                              const Vector3d& p1_n, const Vector3d& p2_n, const Vector3d& center_n,
                              double radius) -> double {
        double pyramid = fabs(center.dot(p1.cross(p2))) / 6.0;
        vector<Vector3d> arr = {p1_n, p2_n, center_n};
        double cone = cone_volume(arr, radius);
        return ((p1 - center).cross(p2 - center)).dot(center) > 0 ? cone - pyramid : -cone + pyramid;
    };

    auto sphere_volume = [&](const Vector3d& p1, const Vector3d& p2, const Vector3d& center, double volume) -> double {
        Vector3d v1 = p1 - center;
        Vector3d v2 = p2 - center;
        Vector3d cross_product = v1.cross(v2);
        double cos_theta = v1.dot(v2) / (v1.norm() * v2.norm());
        cos_theta = max(-1.0, min(1.0, cos_theta));
        double theta = acos(cos_theta);
        return cross_product.dot(center) > 0 ? volume * theta / (2 * M_PI) : -volume * theta / (2 * M_PI);
    };

    auto line_sphere_intersection = [&](const Vector3d& p1, const Vector3d& p2, double radius) -> vector<Vector3d> {
        Vector3d d = p2 - p1;
        double a = d.dot(d), b = 2*p1.dot(d), c = p1.dot(p1)-radius*radius;
        double disc = b*b - 4*a*c;
        vector<Vector3d> intersections;
        if (disc < 0) return intersections;
        double sqrt_disc = sqrt(disc);
        for (double t : {(-b - sqrt_disc)/(2*a), (-b + sqrt_disc)/(2*a)}) 
            if (t >= 0 && t <= 1) intersections.push_back(p1 + t*d);
        return intersections;
    };

    Eigen::Matrix<double, 3, 3> vectors_n;
    for (int i = 0; i < 3; i++)
        vectors_n.row(i) = radius * (vectors.row(i) / vectors.row(i).norm());

    Eigen::Vector3d normal_vector = (vectors.row(1) - vectors.row(0)).cross(vectors.row(2) - vectors.row(0));
    double d = -normal_vector.dot(vectors.row(0));
    Eigen::Vector3d center = (-d / normal_vector.dot(normal_vector)) * normal_vector;
    double center_norm = center.norm();

    std::vector<Eigen::Vector3d> vn = {
        vectors_n.row(0),
        vectors_n.row(1),
        vectors_n.row(2)
    };
    double total_volume = cone_volume(vn, radius);

    if (center_norm > radius) return total_volume;
    if ((vectors.rowwise().norm().array() < radius).all()) {
        return std::fabs(vectors.row(0).dot(vectors.row(1).cross(vectors.row(2)))) / 6.0;
    }
    double center_radius = std::sqrt(radius * radius - center_norm * center_norm);
    Eigen::Vector3d center_n = radius * (center / center_norm);
    std::vector<std::vector<Eigen::Vector3d>> intersections_list;
    for (int i = 0; i < 3; i++) {
        intersections_list.push_back(
            line_sphere_intersection(vectors.row((i + 2) % 3), vectors.row(i), radius)
        );
    }

    double ex_volume = 0;
    double ex_sphere = 2*M_PI*(1-center_norm/radius)*pow(radius,3)/3 - center_norm*M_PI*center_radius*center_radius/3;

    for (int i = 0; i < 3; i++) {
        Eigen::Vector3d v1 = vectors.row((i + 2) % 3);
        Eigen::Vector3d v2 = vectors.row(i);
        auto arr = intersections_list[i];
        if (v1.norm() > radius) {
            if (v2.norm() > radius) {
                if (arr.size() == 2) {
                    ex_volume += sphere_volume(v1, arr[0], center, ex_sphere);
                    ex_volume += pyramid_volume(arr[0], arr[1], center, arr[0], arr[1], center_n, radius);
                    ex_volume += sphere_volume(arr[1], v2, center, ex_sphere);
                } else {
                    ex_volume += sphere_volume(v1, v2, center, ex_sphere);
                }
            } else {
                ex_volume += sphere_volume(v1, arr[0], center, ex_sphere);
                ex_volume += pyramid_volume(arr[0], v2, center, arr[0], vn[i], center_n, radius);
            }
        } else {
            if (v2.norm() > radius) {
                ex_volume += pyramid_volume(v1, arr[0], center, vn[(i + 2) % 3], arr[0], center_n, radius);
                ex_volume += sphere_volume(arr[0], v2, center, ex_sphere);
            } else {
                ex_volume += pyramid_volume(v1, v2, center, vn[(i + 2) % 3], vn[i], center_n, radius);
            }
        }
    }

    total_volume -= fabs(ex_volume);
    return total_volume;
}

double compute_volume(const std::vector<Region>& REGIONS,
    const std::vector<Point3D>& VERTICES,
    const std::vector<Point3D>& POINTS,
    double RADIUS) {
    double total_volume = 0.0;

    for (const auto& s : REGIONS) {
        if (s.empty() || std::find(s.begin(), s.end(), -1) != s.end()) continue;
        std::vector<Point> v;
        for (int idx : s) {
            if (idx < 0 || idx >= (int)VERTICES.size()) continue;
            const auto& pt = VERTICES[idx];
            v.emplace_back(pt[0], pt[1], pt[2]);
        }
        Polyhedron mesh;
        CGAL::convex_hull_3(v.begin(), v.end(), mesh);
        AABB_tree tree(faces(mesh).first, faces(mesh).second, mesh);
        tree.accelerate_distance_queries();
        CGAL::Side_of_triangle_mesh<Polyhedron, K> point_inside(mesh);
        for (const auto& p : POINTS) {
            Point cg_p(p[0], p[1], p[2]);
            if (point_inside(cg_p) == CGAL::ON_BOUNDED_SIDE) {
		double vlm = 0.0;
                for (auto f = mesh.facets_begin(); f != mesh.facets_end(); ++f) {
                    auto h = f->halfedge();
		    Eigen::Matrix3d v_tri;
                    int i = 0;
                    for (auto v_o : CGAL::vertices_around_face(h, mesh)) {
                        const auto& v_p = v_o->point();
                        v_tri.col(i) << v_p.x() - cg_p.x(), v_p.y() - cg_p.y(), v_p.z() - cg_p.z();
			++i;
                    }
                    total_volume += complex_volume_estimate(v_tri.transpose(), RADIUS);
                }
            }
        }
    }
    return total_volume;
}

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
        const double* atom_i = &atoms(i, 0);
        double radius_i = radii(i);
        std::vector<int> neighbor_atoms;
        for (size_t j = 0; j < n_atoms; ++j) {
            if (i == (int)j) continue;
            double d2 = dist_sq(atom_i, &atoms(j, 0));
            double threshold = radius_i + 2*r + radii(j);
            if (d2 <= threshold * threshold) {
                neighbor_atoms.push_back(j);
            }
        }
        if (neighbor_atoms.empty()) continue;
        std::vector<int> valid_waters;
        std::vector<double> dist_iw;
        for (size_t w = 0; w < n_waters; ++w) {
            double d2 = dist_sq(atom_i, &waters(w, 0));
            double limit = radius_i + r;
            if (d2 < limit * limit) {
                valid_waters.push_back(w);
                dist_iw.push_back(std::sqrt(d2));
            }
        }
        if (valid_waters.empty()) continue;
        int local_count = 0;
        for (size_t wi = 0; wi < valid_waters.size(); ++wi) {
            const double* pos_w = &waters(valid_waters[wi], 0);
            double d_iw_val = dist_iw[wi];
            bool ok = true;
            for (size_t j : neighbor_atoms) {
                double d_wj = std::sqrt(dist_sq(pos_w, &atoms(j, 0)));
                if (!(d_iw_val * radii(j) < d_wj * radius_i)) {
                    ok = false;
                    break;
                }
            }
            if (ok) local_count++;
        }
        total_number += local_count;
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
    return total_volume;
}

PYBIND11_MODULE(pmv_calc, m) {
    m.doc() = "pmv module";

    m.def("compute_volume", &compute_volume,
          py::arg("REGIONS"),
          py::arg("VERTICES"),
          py::arg("POINTS"),
          py::arg("RADIUS"));

    m.def("func1", &func1,
          py::arg("atoms_py"),
          py::arg("water_coords_py"),
          py::arg("radiuses_py"),
          py::arg("r"));

    m.def("func2", &func2,
          py::arg("atoms"),
          py::arg("radiuses"),
          py::arg("sel"),
          py::arg("r") = 0.0,
          py::arg("n_points") = 1000);
}
