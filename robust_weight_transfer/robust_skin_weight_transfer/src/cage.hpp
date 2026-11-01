// Cage building: a convex hull of a mesh's vertices, each grown by a small diamond.
// Pure Eigen, so a host build can test it.
#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cstdint>
#include <map>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

namespace rwt {

/**
 * A small diamond (octahedron), scaled down from meters.
 *
 *  V: 6 by 3 vertices
 *  F: 8 by 3 triangles
 */
inline void create_diamond(Eigen::MatrixXd &V, Eigen::MatrixXi &F, double scale_factor = 0.005) {
	V.resize(6, 3);
	V << 0, 0, 1, // Top vertex
			1, 0, 0, // Side vertices
			-1, 0, 0,
			0, 1, 0,
			0, -1, 0,
			0, 0, -1; // Bottom vertex
	V *= scale_factor;
	F.resize(8, 3);
	F << 0, 1, 3,
			0, 3, 2,
			0, 2, 4,
			0, 4, 1,
			5, 3, 1,
			5, 2, 3,
			5, 4, 2,
			5, 1, 4;
}

/**
 * Convex hull of a point set by incremental insertion.
 *
 *  P: #P by 3 points
 *  V: #V by 3 hull vertices, a subset of P
 *  F: #F by 3 outward-facing (counter-clockwise from outside) triangles
 *  Returns false when the points are fewer than four or do not span a volume.
 */
inline bool convex_hull(const Eigen::MatrixXd &P, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	const int n = (int)P.rows();
	auto pt = [&](int i) -> Eigen::Vector3d { return P.row(i).transpose(); };
	if (n < 4) {
		return false;
	}
	const double diagonal = (P.colwise().maxCoeff() - P.colwise().minCoeff()).norm();
	const double eps = 1e-9 * diagonal;
	if (!(diagonal > 0.0)) {
		return false;
	}

	// Seed with a tetrahedron: extreme point, farthest point, farthest from their line, farthest from their plane.
	int i0 = 0;
	for (int i = 1; i < n; ++i) {
		if (P(i, 0) < P(i0, 0)) {
			i0 = i;
		}
	}
	int i1 = -1;
	double best = eps;
	for (int i = 0; i < n; ++i) {
		const double d = (pt(i) - pt(i0)).norm();
		if (d > best) {
			best = d;
			i1 = i;
		}
	}
	if (i1 < 0) {
		return false;
	}
	const Eigen::Vector3d e = pt(i1) - pt(i0);
	int i2 = -1;
	best = eps * diagonal;
	for (int i = 0; i < n; ++i) {
		const double d = e.cross(pt(i) - pt(i0)).norm();
		if (d > best) {
			best = d;
			i2 = i;
		}
	}
	if (i2 < 0) {
		return false;
	}
	const Eigen::Vector3d plane_normal = e.cross(pt(i2) - pt(i0)).normalized();
	int i3 = -1;
	best = eps;
	for (int i = 0; i < n; ++i) {
		const double d = std::abs(plane_normal.dot(pt(i) - pt(i0)));
		if (d > best) {
			best = d;
			i3 = i;
		}
	}
	if (i3 < 0) {
		return false;
	}

	struct Face {
		int a, b, c;
		Eigen::Vector3d normal;
		double offset;
	};
	auto make_face = [&](int a, int b, int c, const Eigen::Vector3d &inside) {
		Face f{ a, b, c, Eigen::Vector3d::Zero(), 0.0 };
		const Eigen::Vector3d pa = pt(a);
		Eigen::Vector3d nrm = (pt(b) - pa).cross(pt(c) - pa).normalized();
		if (nrm.dot(inside - pa) > 0.0) { // Face the interior point away.
			std::swap(f.b, f.c);
			nrm = -nrm;
		}
		f.normal = nrm;
		f.offset = nrm.dot(pa);
		return f;
	};

	const Eigen::Vector3d centroid = (pt(i0) + pt(i1) + pt(i2) + pt(i3)) / 4.0;
	std::vector<Face> faces;
	faces.push_back(make_face(i0, i1, i2, centroid));
	faces.push_back(make_face(i0, i1, i3, centroid));
	faces.push_back(make_face(i0, i2, i3, centroid));
	faces.push_back(make_face(i1, i2, i3, centroid));

	for (int p = 0; p < n; ++p) {
		if (p == i0 || p == i1 || p == i2 || p == i3) {
			continue;
		}
		const Eigen::Vector3d point = pt(p);
		std::vector<bool> visible(faces.size(), false);
		bool any_visible = false;
		for (size_t f = 0; f < faces.size(); ++f) {
			if (faces[f].normal.dot(point) - faces[f].offset > eps) {
				visible[f] = true;
				any_visible = true;
			}
		}
		if (!any_visible) {
			continue;
		}
		// A directed edge is on the horizon when the face across it is not visible.
		std::map<std::pair<int, int>, size_t> edge_owner;
		for (size_t f = 0; f < faces.size(); ++f) {
			edge_owner[{ faces[f].a, faces[f].b }] = f;
			edge_owner[{ faces[f].b, faces[f].c }] = f;
			edge_owner[{ faces[f].c, faces[f].a }] = f;
		}
		std::vector<std::pair<int, int>> horizon;
		for (size_t f = 0; f < faces.size(); ++f) {
			if (!visible[f]) {
				continue;
			}
			const int edges[3][2] = { { faces[f].a, faces[f].b }, { faces[f].b, faces[f].c }, { faces[f].c, faces[f].a } };
			for (const auto &edge : edges) {
				auto across = edge_owner.find({ edge[1], edge[0] });
				if (across != edge_owner.end() && !visible[across->second]) {
					horizon.push_back({ edge[0], edge[1] });
				}
			}
		}
		std::vector<Face> kept;
		for (size_t f = 0; f < faces.size(); ++f) {
			if (!visible[f]) {
				kept.push_back(faces[f]);
			}
		}
		faces.swap(kept);
		for (const auto &edge : horizon) {
			faces.push_back(make_face(edge.first, edge.second, p, centroid));
		}
	}

	// Keep only the vertices the hull uses.
	std::vector<int> remap(n, -1);
	std::vector<int> used;
	F.resize(faces.size(), 3);
	for (size_t f = 0; f < faces.size(); ++f) {
		const int ids[3] = { faces[f].a, faces[f].b, faces[f].c };
		for (int k = 0; k < 3; ++k) {
			if (remap[ids[k]] < 0) {
				remap[ids[k]] = (int)used.size();
				used.push_back(ids[k]);
			}
			F(f, k) = remap[ids[k]];
		}
	}
	V.resize(used.size(), 3);
	for (size_t i = 0; i < used.size(); ++i) {
		V.row(i) = P.row(used[i]);
	}
	return true;
}

/**
 * A cage around a mesh: the convex hull of its vertices and of a diamond placed on
 * each of num_points randomly chosen vertices.
 *
 *  target_vertices: #V by 3
 *  num_points: how many vertices get a diamond; clamped to #V
 *  seed: seed of the random choice, so a run repeats
 *  V, F: the cage hull
 */
inline bool create_cage(const Eigen::MatrixXd &target_vertices, int num_points, uint32_t seed,
		Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	Eigen::MatrixXd diamond_vertices;
	Eigen::MatrixXi diamond_faces;
	create_diamond(diamond_vertices, diamond_faces);

	std::vector<int> indices(target_vertices.rows());
	std::iota(indices.begin(), indices.end(), 0);
	std::mt19937 rng(seed);
	std::shuffle(indices.begin(), indices.end(), rng);
	const int count = std::max(0, std::min(num_points, (int)indices.size()));

	Eigen::MatrixXd points(target_vertices.rows() + count * diamond_vertices.rows(), 3);
	points.topRows(target_vertices.rows()) = target_vertices;
	int row = (int)target_vertices.rows();
	for (int k = 0; k < count; ++k) {
		for (int d = 0; d < diamond_vertices.rows(); ++d) {
			points.row(row++) = diamond_vertices.row(d) + target_vertices.row(indices[k]);
		}
	}
	return convex_hull(points, V, F);
}

} // namespace rwt
