// Property tests on witness-cpp: each property is checked over a ladder of generated inputs,
// and each is paired with a planted-false control that must be caught.
#pragma once

#include <array>
#include <functional>
#include <string>
#include <vector>

#include "self_test.hpp"
#include "witness/ladder.h"

namespace rwt {

// A small ladder: the guest runs these in an emulator.
inline constexpr std::array<witness::Level, 3> PROPERTY_LADDER = { {
		{ 0, 64, 256, 25 },
		{ 1, 512, 1024, 100 },
		{ 2, 4000, 4096, 300 },
} };

// A jittered n by n grid with random skin weights and a random matched flag per vertex.
struct GridCase {
	int n = 3;
	int bones = 2;
	std::vector<double> z; // height of each vertex
	std::vector<double> weights; // n*n*bones, each row sums to 1
	std::vector<char> matched; // vertex 0 is always matched
};

inline GridCase gen_grid_case(witness::RNG &rng, const witness::Level &) {
	GridCase c;
	c.n = rng.int_range(3, 6);
	c.bones = rng.int_range(2, 4);
	const int count = c.n * c.n;
	for (int i = 0; i < count; ++i) {
		c.z.push_back(rng.float_range(-0.05, 0.05));
		double sum = 0.0;
		std::vector<double> row;
		for (int b = 0; b < c.bones; ++b) {
			row.push_back(rng.float_range(0.01, 1.0));
			sum += row.back();
		}
		for (double v : row) {
			c.weights.push_back(v / sum);
		}
		c.matched.push_back(i == 0 || rng.float_range(0.0, 1.0) < 0.6);
	}
	return c;
}

inline void grid_mesh(const GridCase &c, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	const int n = c.n;
	V.resize(n * n, 3);
	F.resize(2 * (n - 1) * (n - 1), 3);
	for (int y = 0; y < n; ++y) {
		for (int x = 0; x < n; ++x) {
			V.row(y * n + x) << x / double(n - 1), y / double(n - 1), c.z[y * n + x];
		}
	}
	int f = 0;
	for (int y = 0; y + 1 < n; ++y) {
		for (int x = 0; x + 1 < n; ++x) {
			const int a = y * n + x, b = a + 1, d = a + n, e = d + 1;
			F.row(f++) << a, b, d;
			F.row(f++) << b, e, d;
		}
	}
}

inline Eigen::MatrixXd grid_weights(const GridCase &c) {
	Eigen::MatrixXd W(c.n * c.n, c.bones);
	for (int i = 0; i < c.n * c.n; ++i) {
		for (int b = 0; b < c.bones; ++b) {
			W(i, b) = c.weights[i * c.bones + b];
		}
	}
	return W;
}

inline BoolArray grid_matched(const GridCase &c) {
	BoolArray m(c.n * c.n);
	for (int i = 0; i < c.n * c.n; ++i) {
		m(i) = c.matched[i] != 0;
	}
	return m;
}

inline Eigen::MatrixXd gen_points_matrix(witness::RNG &rng, int count) {
	Eigen::MatrixXd P(count, 3);
	for (int i = 0; i < count; ++i) {
		P.row(i) << rng.float_range(-1, 1), rng.float_range(-1, 1), rng.float_range(-1, 1);
	}
	return P;
}

struct PointsCase {
	int count = 4;
	uint32_t seed = 0;
	uint64_t points_seed = 0;
};

inline PointsCase gen_points_case(witness::RNG &rng, const witness::Level &) {
	PointsCase c;
	c.count = rng.int_range(4, 30);
	c.seed = (uint32_t)rng.next_u64();
	c.points_seed = rng.next_u64();
	return c;
}

inline Eigen::MatrixXd points_of(const PointsCase &c) {
	witness::RNG rng(c.points_seed ? c.points_seed : 1);
	return gen_points_matrix(rng, c.count);
}

inline bool hull_contains(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &points) {
	for (int i = 0; i < points.rows(); ++i) {
		for (int f = 0; f < F.rows(); ++f) {
			const Eigen::Vector3d a = V.row(F(f, 0)), b = V.row(F(f, 1)), c = V.row(F(f, 2));
			const Eigen::Vector3d n = (b - a).cross(c - a).normalized();
			if (n.dot(Eigen::Vector3d(points.row(i)) - a) > 1e-9) {
				return false;
			}
		}
	}
	return true;
}

template <class T>
inline bool holds(const char *name, witness::Generator<T> gen, std::function<bool(const T &)> predicate) {
	const witness::Trial t = witness::resolve_with_ladder<T>(name, PROPERTY_LADDER, std::move(gen), std::move(predicate),
			&witness::no_shrink<T>, std::function<void(std::ostream &, const T &)>());
	if (t.outcome == witness::Outcome::FOUND) {
		std::cerr << t.message << std::endl;
		return false;
	}
	return true;
}

// A planted-false property: the ladder must report a counterexample, or the check above proves nothing.
template <class T>
inline bool is_caught(const char *name, witness::Generator<T> gen, std::function<bool(const T &)> predicate) {
	const witness::Trial t = witness::resolve_with_ladder<T>(name, PROPERTY_LADDER, std::move(gen), std::move(predicate),
			&witness::no_shrink<T>, std::function<void(std::ostream &, const T &)>());
	return t.outcome == witness::Outcome::FOUND;
}

inline bool inpaint_keeps_partition_of_unity(const GridCase &c) {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	grid_mesh(c, V, F);
	const Eigen::MatrixXd W = grid_weights(c);
	const BoolArray matched = grid_matched(c);
	Eigen::MatrixXd out;
	if (!inpaint(V, F, W, matched, out) || !is_valid_array(out)) {
		return false;
	}
	for (int i = 0; i < out.rows(); ++i) {
		if (std::abs(out.row(i).sum() - 1.0) > 1e-6) {
			return false;
		}
		if (matched(i) && !out.row(i).isApprox(W.row(i), 1e-9)) {
			return false;
		}
	}
	return true;
}

inline bool inpaint_changes_nothing(const GridCase &c) {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	grid_mesh(c, V, F);
	const Eigen::MatrixXd W = grid_weights(c);
	Eigen::MatrixXd out;
	return inpaint(V, F, W, grid_matched(c), out) && out.isApprox(W, 1e-9);
}

inline bool smooth_is_a_convex_average(const GridCase &c) {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	grid_mesh(c, V, F);
	const Eigen::MatrixXd W = grid_weights(c);
	const BoolArray matched = grid_matched(c);
	Eigen::MatrixXd out;
	BoolArray to_smooth;
	smooth(out, to_smooth, V, F, W, matched, 0.4);
	for (int i = 0; i < out.rows(); ++i) {
		if (!matched(i) && !to_smooth(i)) {
			return false; // Every unmatched vertex is smoothed.
		}
		if (!to_smooth(i) && !out.row(i).isApprox(W.row(i), 1e-12)) {
			return false; // Others stay put.
		}
		if (std::abs(out.row(i).sum() - 1.0) > 1e-9) {
			return false;
		}
	}
	for (int b = 0; b < W.cols(); ++b) {
		if (out.col(b).minCoeff() < W.col(b).minCoeff() - 1e-9 || out.col(b).maxCoeff() > W.col(b).maxCoeff() + 1e-9) {
			return false; // A convex average never leaves the input range.
		}
	}
	return true;
}

inline bool smooth_changes_nothing(const GridCase &c) {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	grid_mesh(c, V, F);
	const Eigen::MatrixXd W = grid_weights(c);
	Eigen::MatrixXd out;
	BoolArray to_smooth;
	smooth(out, to_smooth, V, F, W, grid_matched(c), 0.4);
	return out.isApprox(W, 1e-12);
}

inline bool identical_meshes_match_everywhere(const GridCase &c) {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	grid_mesh(c, V, F);
	const Eigen::MatrixXd W = grid_weights(c);
	const Eigen::MatrixXd N = Eigen::MatrixXd::Ones(V.rows(), 1) * Eigen::RowVector3d(0, 0, 1);
	Eigen::MatrixXd W2;
	BoolArray matched;
	find_matches_closest_surface(V, F, N, V, F, N, W, 1e-12, 1.0, W2, matched);
	return matched.all() && W2.isApprox(W, 1e-9);
}

inline bool far_target_matches_nowhere(const GridCase &c) {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	grid_mesh(c, V, F);
	const Eigen::MatrixXd W = grid_weights(c);
	const Eigen::MatrixXd N = Eigen::MatrixXd::Ones(V.rows(), 1) * Eigen::RowVector3d(0, 0, 1);
	Eigen::MatrixXd lifted = V;
	lifted.col(2).array() += 10.0;
	Eigen::MatrixXd W2;
	BoolArray matched;
	find_matches_closest_surface(V, F, N, lifted, F, N, W, 1.0, 30.0, W2, matched);
	return !matched.any();
}

inline bool hull_is_valid(const PointsCase &c) {
	const Eigen::MatrixXd P = points_of(c);
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	double volume = 0.0;
	return convex_hull(P, V, F) && hull_is_closed_and_outward(V, F, volume) && volume > 0.0 && V.rows() <= P.rows() &&
			// Euler: V - E + F = 2 with E = 3F/2 gives F = 2V - 4.
			F.rows() == 2 * V.rows() - 4 && hull_contains(V, F, P);
}

inline bool hull_is_empty(const PointsCase &c) {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	return convex_hull(points_of(c), V, F) && F.rows() == 0;
}

inline bool cage_contains_target_and_repeats(const PointsCase &c) {
	const Eigen::MatrixXd P = points_of(c);
	Eigen::MatrixXd V1, V2;
	Eigen::MatrixXi F1, F2;
	if (!create_cage(P, 10, c.seed, V1, F1) || !create_cage(P, 10, c.seed, V2, F2)) {
		return false;
	}
	double volume = 0.0;
	return V1.isApprox(V2) && F1 == F2 && hull_is_closed_and_outward(V1, F1, volume) && volume > 0.0 && hull_contains(V1, F1, P);
}

inline int run_property_tests() {
	struct Entry {
		const char *name;
		bool (*run)();
	};
	const Entry properties[] = {
		{ "property inpaint keeps weights a partition of unity", [] {
			 return holds<GridCase>("inpaint partition of unity", &gen_grid_case, &inpaint_keeps_partition_of_unity);
		 } },
		{ "property smooth is a convex average", [] {
			 return holds<GridCase>("smooth convex average", &gen_grid_case, &smooth_is_a_convex_average);
		 } },
		{ "property identical meshes match everywhere", [] {
			 return holds<GridCase>("identical meshes match", &gen_grid_case, &identical_meshes_match_everywhere);
		 } },
		{ "property a far target matches nowhere", [] {
			 return holds<GridCase>("far target matches nowhere", &gen_grid_case, &far_target_matches_nowhere);
		 } },
		{ "property convex hull is closed, outward and contains its points", [] {
			 return holds<PointsCase>("convex hull valid", &gen_points_case, &hull_is_valid);
		 } },
		{ "property cage contains its mesh and repeats", [] {
			 return holds<PointsCase>("cage contains target", &gen_points_case, &cage_contains_target_and_repeats);
		 } },
		// Planted-false controls: each must be caught.
		{ "control inpaint changing nothing is caught", [] {
			 return is_caught<GridCase>("inpaint changes nothing", &gen_grid_case, &inpaint_changes_nothing);
		 } },
		{ "control smooth changing nothing is caught", [] {
			 return is_caught<GridCase>("smooth changes nothing", &gen_grid_case, &smooth_changes_nothing);
		 } },
		{ "control an empty hull is caught", [] {
			 return is_caught<PointsCase>("hull is empty", &gen_points_case, &hull_is_empty);
		 } },
	};
	int failed = 0;
	for (const Entry &p : properties) {
		const bool ok = p.run();
		std::cout << (ok ? "PASS " : "FAIL ") << p.name << std::endl;
		failed += ok ? 0 : 1;
	}
	return failed;
}

// Every example test and every property test; returns the number that failed.
inline int run_all_tests() {
	return run_self_tests() + run_property_tests();
}

} // namespace rwt
