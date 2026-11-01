// Self tests for the transfer and the cage. They run on the host and in the guest.
#pragma once

#include <cmath>
#include <iostream>
#include <set>
#include <utility>

#include "cage.hpp"
#include "transfer.hpp"

namespace rwt {

#define RWT_CHECK(cond)                                                                      \
	do {                                                                                     \
		if (!(cond)) {                                                                       \
			std::cerr << __func__ << ": check failed: " #cond " (line " << __LINE__ << ")" << std::endl; \
			return false;                                                                    \
		}                                                                                    \
	} while (0)

inline bool test_find_closest_point_on_surface() {
	Eigen::MatrixXd vertices(3, 3);
	vertices << -1, -1, -1, 1, -1, -1, 1, 1, -1;
	Eigen::MatrixXi triangles(1, 3);
	triangles << 0, 1, 2;
	Eigen::MatrixXd points(3, 3);
	points << 0, 0, 0, 2, 2, 2, 0, 0, -2;
	Eigen::VectorXd expected_distances(3);
	expected_distances << 1.0, 11.0, 1.0;
	Eigen::MatrixXd expected_points(3, 3);
	expected_points << 0, 0, -1, 1, 1, -1, 0, 0, -1;
	Eigen::MatrixXd expected_barycentric(3, 3);
	expected_barycentric << 0.5, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0, 0.5;

	Eigen::VectorXd sqrD;
	Eigen::VectorXi I;
	Eigen::MatrixXd C, B;
	find_closest_point_on_surface(points, vertices, triangles, sqrD, I, C, B);
	RWT_CHECK(sqrD.isApprox(expected_distances, 1e-6));
	RWT_CHECK((I.array() == 0).all());
	RWT_CHECK(C.isApprox(expected_points, 1e-6));
	RWT_CHECK(B.isApprox(expected_barycentric, 1e-6));

	// The same triangle translated by (1, 2, 3).
	Eigen::MatrixXd moved_vertices = vertices;
	moved_vertices.rowwise() += Eigen::RowVector3d(1, 2, 3);
	Eigen::MatrixXd moved_points(3, 3);
	moved_points << 1, 2, 2, 3, 4, 5, 1, 2, 1;
	Eigen::VectorXd moved_distances(3);
	moved_distances << 0.0, 11.0, 1.0;
	Eigen::MatrixXd moved_expected(3, 3);
	moved_expected << 1, 2, 2, 2, 3, 2, 1, 2, 2;
	find_closest_point_on_surface(moved_points, moved_vertices, triangles, sqrD, I, C, B);
	RWT_CHECK(sqrD.isApprox(moved_distances, 1e-6));
	RWT_CHECK((I.array() == 0).all());
	RWT_CHECK(C.isApprox(moved_expected, 1e-6));
	RWT_CHECK(B.isApprox(expected_barycentric, 1e-6));
	return true;
}

inline bool test_interpolate_attribute_from_bary() {
	Eigen::MatrixXd attributes(5, 2);
	attributes << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
	Eigen::MatrixXd bary(2, 3);
	bary << 0.2, 0.5, 0.3, 0.6, 0.3, 0.1;
	Eigen::VectorXi primitives(2);
	primitives << 0, 1;
	Eigen::MatrixXi triangles(2, 3);
	triangles << 0, 1, 2, 2, 3, 4;
	Eigen::MatrixXd expected(2, 2);
	expected << 1 * 0.2 + 3 * 0.5 + 5 * 0.3, 2 * 0.2 + 4 * 0.5 + 6 * 0.3,
			5 * 0.6 + 7 * 0.3 + 9 * 0.1, 6 * 0.6 + 8 * 0.3 + 10 * 0.1;
	Eigen::MatrixXd result;
	interpolate_attribute_from_bary(attributes, bary, primitives, triangles, result);
	RWT_CHECK(result.isApprox(expected, 1e-6));
	return true;
}

inline bool test_normalize_vector() {
	Eigen::VectorXd vector(3);
	vector << 3, 4, 0;
	Eigen::VectorXd expected(3);
	expected << 0.6, 0.8, 0;
	RWT_CHECK(normalize_vector(vector).isApprox(expected, 1e-6));
	return true;
}

inline bool test_find_matches_closest_surface() {
	Eigen::MatrixXd source_vertices(3, 3);
	source_vertices << 0, 0, 0, 1, 0, 0, 0, 1, 0;
	Eigen::MatrixXi source_triangles(1, 3);
	source_triangles << 0, 1, 2;
	Eigen::MatrixXd source_normals(3, 3);
	source_normals << 0, 0, 1, 0, 0, 1, 0, 0, 1;
	Eigen::MatrixXd source_weights(3, 2);
	source_weights << 1, 0, 0, 1, 0.5, 0.5;

	Eigen::MatrixXd target_vertices(2, 3);
	target_vertices << 0.1, 0.1, 0, 2, 2, 2; // The second is too far.
	Eigen::MatrixXi target_triangles(1, 2);
	target_triangles << 0, 1;
	Eigen::MatrixXd target_normals(2, 3);
	target_normals << 0, 0, 1, 1, 0, 0; // The second is at the wrong angle.

	Eigen::MatrixXd expected_weights(2, 2);
	expected_weights << 0.85, 0.15, 0.25, 0.75;

	Eigen::MatrixXd target_weights;
	BoolArray matched;
	find_matches_closest_surface(source_vertices, source_triangles, source_normals, target_vertices, target_triangles,
			target_normals, source_weights, 0.5, 10, target_weights, matched);
	RWT_CHECK(matched(0) && !matched(1));
	RWT_CHECK(target_weights.isApprox(expected_weights, 1e-6));
	return true;
}

inline bool test_find_matches_ignores_flipped_normals() {
	Eigen::MatrixXd source_vertices(3, 3);
	source_vertices << 0, 0, 0, 1, 0, 0, 0, 1, 0;
	Eigen::MatrixXi source_triangles(1, 3);
	source_triangles << 0, 1, 2;
	Eigen::MatrixXd source_normals(3, 3);
	source_normals << 0, 0, 1, 0, 0, 1, 0, 0, 1;
	Eigen::MatrixXd source_weights(3, 1);
	source_weights << 1, 1, 1;
	Eigen::MatrixXd target_vertices(1, 3);
	target_vertices << 0.2, 0.2, 0.01;
	Eigen::MatrixXi target_triangles(1, 3);
	target_triangles << 0, 0, 0;
	Eigen::MatrixXd target_normals(1, 3);
	target_normals << 0, 0, -1; // Flipped, still parallel.

	Eigen::MatrixXd target_weights;
	BoolArray matched;
	find_matches_closest_surface(source_vertices, source_triangles, source_normals, target_vertices, target_triangles,
			target_normals, source_weights, 0.5, 10, target_weights, matched);
	RWT_CHECK(matched(0));
	return true;
}

inline bool test_is_valid_array() {
	Eigen::MatrixXd valid(2, 2);
	valid << 1, 2, 3, 4;
	Eigen::MatrixXd invalid(2, 2);
	invalid << std::numeric_limits<double>::quiet_NaN(), 2, std::numeric_limits<double>::infinity(), 4;
	RWT_CHECK(is_valid_array(valid));
	RWT_CHECK(!is_valid_array(invalid));
	return true;
}

inline bool test_inpaint() {
	Eigen::MatrixXd V2(4, 3);
	V2 << 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0;
	Eigen::MatrixXi F2(2, 3);
	F2 << 0, 1, 2, 1, 2, 3;
	Eigen::MatrixXd W2(4, 2);
	W2 << 1, 0, 0, 1, 0.5, 0.5, 0, 0;
	BoolArray matched(4);
	matched << true, true, true, false;
	Eigen::MatrixXd expected(4, 2);
	expected << 1, 0, 0, 1, 0.5, 0.5, 0.117647, 0.882353;

	Eigen::MatrixXd inpainted;
	RWT_CHECK(inpaint(V2, F2, W2, matched, inpainted));
	RWT_CHECK(inpainted.isApprox(expected, 1e-5));

	// Bad input fails instead of crashing.
	Eigen::MatrixXi bad_faces(1, 3);
	bad_faces << 0, 1, 9;
	RWT_CHECK(!inpaint(V2, bad_faces, W2, matched, inpainted));
	RWT_CHECK(!inpaint(V2, F2, W2, BoolArray::Constant(4, false), inpainted));
	return true;
}

inline bool test_smooth() {
	Eigen::MatrixXd vertices(5, 3);
	vertices << 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 0; // The fifth is in no face.
	Eigen::MatrixXi faces(2, 3);
	faces << 0, 1, 2, 1, 2, 3;
	Eigen::MatrixXd weights(5, 2);
	weights << 1, 0, 0, 1, 0.5, 0.5, 0.25, 0.75, 0.1, 0.9;
	BoolArray matched(5);
	matched << true, true, true, false, false;
	Eigen::MatrixXd expected(5, 2);
	expected << 0.85, 0.15, 0.10666667, 0.89333333, 0.48044444, 0.51955556, 0.25871111, 0.74128889, 0.1, 0.9;

	Eigen::MatrixXd smoothed;
	BoolArray to_smooth;
	smooth(smoothed, to_smooth, vertices, faces, weights, matched, 1.5, 1, 0.2);
	RWT_CHECK(smoothed.isApprox(expected, 1e-6));
	RWT_CHECK(to_smooth.all());
	return true;
}

// Source weights vary along x; the right half of the target is lifted out of reach and must be inpainted.
inline bool test_robust_weight_transfer() {
	const int n = 5;
	Eigen::MatrixXd V1(n * n, 3), V2(n * n, 3), N(n * n, 3), W1(n * n, 2);
	Eigen::MatrixXi F(2 * (n - 1) * (n - 1), 3);
	for (int y = 0; y < n; ++y) {
		for (int x = 0; x < n; ++x) {
			const int i = y * n + x;
			const double fx = x / double(n - 1), fy = y / double(n - 1);
			V1.row(i) << fx, fy, 0;
			V2.row(i) << fx, fy, fx > 0.5 ? 1.0 : 0.0;
			N.row(i) << 0, 0, 1;
			W1.row(i) << fx, 1.0 - fx;
		}
	}
	int f = 0;
	for (int y = 0; y + 1 < n; ++y) {
		for (int x = 0; x + 1 < n; ++x) {
			const int a = y * n + x, b = a + 1, c = a + n, d = c + 1;
			F.row(f++) << a, b, c;
			F.row(f++) << b, d, c;
		}
	}
	const double distance_threshold = 0.05 * (V2.colwise().maxCoeff() - V2.colwise().minCoeff()).norm();

	Eigen::MatrixXd W2;
	BoolArray matched;
	find_matches_closest_surface(V1, F, N, V2, F, N, W1, distance_threshold * distance_threshold, 30, W2, matched);
	int expected_matched = 0;
	for (int i = 0; i < n * n; ++i) {
		const bool near = V2(i, 2) == 0.0;
		RWT_CHECK(matched(i) == near);
		expected_matched += near;
	}
	RWT_CHECK(matched.count() == expected_matched);

	Eigen::MatrixXd inpainted;
	RWT_CHECK(inpaint(V2, F, W2, matched, inpainted));
	RWT_CHECK(is_valid_array(inpainted));
	RWT_CHECK(inpainted.rows() == n * n && inpainted.cols() == 2);
	for (int i = 0; i < n * n; ++i) {
		RWT_CHECK(std::abs(inpainted.row(i).sum() - 1.0) < 1e-6); // Weights stay a partition of unity.
		if (matched(i)) {
			RWT_CHECK(inpainted.row(i).isApprox(W2.row(i), 1e-9)); // Matched vertices are untouched.
		}
	}

	Eigen::MatrixXd smoothed;
	BoolArray to_smooth;
	smooth(smoothed, to_smooth, V2, F, inpainted, matched, distance_threshold);
	RWT_CHECK(is_valid_array(smoothed));
	RWT_CHECK(smoothed.rows() == n * n);
	return true;
}

inline bool test_create_diamond() {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	create_diamond(V, F);
	RWT_CHECK(V.rows() == 6 && F.rows() == 8);
	RWT_CHECK(std::abs(V.col(2).maxCoeff() - 0.005) < 1e-12);
	RWT_CHECK(std::abs(V.col(2).minCoeff() + 0.005) < 1e-12);
	return true;
}

// Signed volume by the divergence theorem, and a check that every edge is shared by two faces in opposite order.
inline bool hull_is_closed_and_outward(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double &volume) {
	volume = 0.0;
	std::set<std::pair<int, int>> edges;
	for (int f = 0; f < F.rows(); ++f) {
		const Eigen::Vector3d a = V.row(F(f, 0)), b = V.row(F(f, 1)), c = V.row(F(f, 2));
		volume += a.dot(b.cross(c)) / 6.0;
		for (int k = 0; k < 3; ++k) {
			if (!edges.insert({ F(f, k), F(f, (k + 1) % 3) }).second) {
				return false; // A directed edge used twice.
			}
		}
	}
	for (const auto &e : edges) {
		if (edges.find({ e.second, e.first }) == edges.end()) {
			return false;
		}
	}
	return true;
}

inline bool test_convex_hull() {
	Eigen::MatrixXd P(10, 3);
	P << 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
			0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1,
			0.5, 0.5, 0.5, // Inside.
			0.5, 0.5, 0; // On a face.
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	RWT_CHECK(convex_hull(P, V, F));
	double volume = 0.0;
	RWT_CHECK(V.rows() == 8);
	RWT_CHECK(hull_is_closed_and_outward(V, F, volume));
	RWT_CHECK(std::abs(volume - 1.0) < 1e-9);

	Eigen::MatrixXd flat(4, 3);
	flat << 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0;
	RWT_CHECK(!convex_hull(flat, V, F)); // No volume.
	return true;
}

inline bool test_create_cage() {
	Eigen::MatrixXd target(8, 3);
	target << 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1;
	Eigen::MatrixXd V1, V2;
	Eigen::MatrixXi F1, F2;
	RWT_CHECK(create_cage(target, 200, 7, V1, F1));
	RWT_CHECK(create_cage(target, 200, 7, V2, F2));
	RWT_CHECK(V1.isApprox(V2) && F1 == F2); // Repeats with the same seed.
	double volume = 0.0;
	RWT_CHECK(hull_is_closed_and_outward(V1, F1, volume));
	RWT_CHECK(volume > 1.0); // The diamonds grow the cube.
	// Every target vertex lies inside the cage.
	for (int i = 0; i < target.rows(); ++i) {
		for (int f = 0; f < F1.rows(); ++f) {
			const Eigen::Vector3d a = V1.row(F1(f, 0)), b = V1.row(F1(f, 1)), c = V1.row(F1(f, 2));
			const Eigen::Vector3d n = (b - a).cross(c - a).normalized();
			RWT_CHECK(n.dot(Eigen::Vector3d(target.row(i)) - a) <= 1e-9);
		}
	}
	return true;
}

// Runs every self test; returns the number that failed.
inline int run_self_tests() {
	struct Entry {
		const char *name;
		bool (*run)();
	};
	const Entry tests[] = {
		{ "find_closest_point_on_surface", test_find_closest_point_on_surface },
		{ "interpolate_attribute_from_bary", test_interpolate_attribute_from_bary },
		{ "normalize_vector", test_normalize_vector },
		{ "find_matches_closest_surface", test_find_matches_closest_surface },
		{ "find_matches_ignores_flipped_normals", test_find_matches_ignores_flipped_normals },
		{ "is_valid_array", test_is_valid_array },
		{ "inpaint", test_inpaint },
		{ "smooth", test_smooth },
		{ "robust_weight_transfer", test_robust_weight_transfer },
		{ "create_diamond", test_create_diamond },
		{ "convex_hull", test_convex_hull },
		{ "create_cage", test_create_cage },
	};
	int failed = 0;
	for (const Entry &t : tests) {
		const bool ok = t.run();
		std::cout << (ok ? "PASS " : "FAIL ") << t.name << std::endl;
		failed += ok ? 0 : 1;
	}
	return failed;
}

} // namespace rwt
