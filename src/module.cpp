#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>


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

PYBIND11_MODULE(pmv_calc, m) {
    m.doc() = "compute_volume binding";
    m.def("compute_volume", &compute_volume,
          py::arg("REGIONS"),
          py::arg("VERTICES"),
          py::arg("POINTS"),
          py::arg("RADIUS"));
}
