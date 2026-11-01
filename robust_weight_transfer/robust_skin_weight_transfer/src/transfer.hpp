// Robust skin weights transfer via weight inpainting (SIGGRAPH Asia 2023).
// Pure Eigen, libigl and geometry-central; no engine types, so a host build can test it.
#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <queue>
#include <set>
#include <tuple>
#include <vector>

#include <igl/adjacency_list.h>
#include <igl/barycentric_coordinates.h>
#include <igl/invert_diag.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/slice_mask.h>

#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/tufted_laplacian.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace rwt {

namespace gc = geometrycentral;

using BoolArray = Eigen::Array<bool, Eigen::Dynamic, 1>;

/**
 * Find the closest points on the surface of the V,F mesh for a number of points.
 *
 *  P: #P by 3 points
 *  V: #V by 3 mesh vertices
 *  F: #F by 3 mesh triangle indices
 *  sqrD: #P smallest squared distances
 *  I: #P primitive indices of the closest triangles
 *  C: #P by 3 closest points
 *  B: #P by 3 barycentric coordinates of the closest points
 */
inline void find_closest_point_on_surface(const Eigen::MatrixXd &P, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
		Eigen::VectorXd &sqrD, Eigen::VectorXi &I, Eigen::MatrixXd &C, Eigen::MatrixXd &B) {
	igl::point_mesh_squared_distance(P, V, F, sqrD, I, C);

	Eigen::MatrixXi F_closest = F(I, Eigen::indexing::all);
	Eigen::MatrixXd V1 = V(F_closest(Eigen::indexing::all, 0), Eigen::indexing::all);
	Eigen::MatrixXd V2 = V(F_closest(Eigen::indexing::all, 1), Eigen::indexing::all);
	Eigen::MatrixXd V3 = V(F_closest(Eigen::indexing::all, 2), Eigen::indexing::all);

	igl::barycentric_coordinates(C, V1, V2, V3, B);
}

/**
 * Interpolate per-vertex attributes A via barycentric coordinates B of the F[I,:] vertices.
 *
 *  A: #V by N per-vertex attributes
 *  B: #B by 3 barycentric coordinates
 *  I: #B primitive indices containing the closest point
 *  F: #F by 3 mesh triangle indices
 *  A_out: #B by N interpolated attributes
 */
inline void interpolate_attribute_from_bary(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
		const Eigen::VectorXi &I, const Eigen::MatrixXi &F, Eigen::MatrixXd &A_out) {
	Eigen::MatrixXi F_closest = F(I, Eigen::indexing::all);

	Eigen::MatrixXd a1 = A(F_closest.col(0), Eigen::indexing::all);
	Eigen::MatrixXd a2 = A(F_closest.col(1), Eigen::indexing::all);
	Eigen::MatrixXd a3 = A(F_closest.col(2), Eigen::indexing::all);

	a1 = a1.array().colwise() * B.col(0).array();
	a2 = a2.array().colwise() * B.col(1).array();
	a3 = a3.array().colwise() * B.col(2).array();

	A_out = a1 + a2 + a3;
}

inline Eigen::VectorXd normalize_vector(const Eigen::VectorXd &p_vector) {
	return p_vector / p_vector.norm();
}

/**
 * For each vertex on the target mesh find a match on the source mesh.
 *
 *  V1, F1, N1: source mesh vertices, triangles, per-vertex normals
 *  V2, F2, N2: target mesh vertices, triangles, per-vertex normals
 *  W1: #V1 by num_bones source skin weights
 *  distance_threshold_squared: squared distance threshold
 *  angle_threshold_degrees: normal threshold
 *  W2: #V2 by num_bones weights copied from the source by the closest point method
 *  Matched: #V2 flags; true where the closest point passes both thresholds
 */
inline void find_matches_closest_surface(const Eigen::MatrixXd &V1, const Eigen::MatrixXi &F1, const Eigen::MatrixXd &N1,
		const Eigen::MatrixXd &V2, const Eigen::MatrixXi &F2, const Eigen::MatrixXd &N2,
		const Eigen::MatrixXd &W1, double distance_threshold_squared, double angle_threshold_degrees,
		Eigen::MatrixXd &W2, BoolArray &Matched) {
	Matched = BoolArray::Constant(V2.rows(), false);
	Eigen::VectorXd sqrD;
	Eigen::VectorXi I;
	Eigen::MatrixXd C, B;
	find_closest_point_on_surface(V2, V1, F1, sqrD, I, C, B);

	// Interpolate the skin weights and normals at each closest point.
	if (W1.size() > 0) {
		interpolate_attribute_from_bary(W1, B, I, F1, W2);
	} else {
		W2 = Eigen::MatrixXd::Zero(V2.rows(), 0);
	}
	Eigen::MatrixXd N1_interpolated;
	interpolate_attribute_from_bary(N1, B, I, F1, N1_interpolated);

	for (int row = 0; row < V2.rows(); ++row) {
		const double n1_norm = N1_interpolated.row(row).norm();
		const double n2_norm = N2.row(row).norm();
		if (n1_norm == 0.0 || n2_norm == 0.0) {
			continue; // No usable normal, so no match.
		}
		// The absolute value of the dot product ignores flipped normals.
		double cos_angle = std::abs(N1_interpolated.row(row).dot(N2.row(row))) / (n1_norm * n2_norm);
		cos_angle = std::min(1.0, cos_angle);
		const double deg_angle = std::acos(cos_angle) * (180.0 / M_PI);

		if (sqrD(row) <= distance_threshold_squared && deg_angle <= angle_threshold_degrees) {
			Matched(row) = true;
		}
	}
}

inline bool is_valid_array(const Eigen::MatrixXd &p_matrix) {
	return p_matrix.allFinite();
}

inline bool is_valid_array(const Eigen::SparseMatrix<double> &p_matrix) {
	for (int k = 0; k < p_matrix.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(p_matrix, k); it; ++it) {
			if (!std::isfinite(it.value())) {
				return false;
			}
		}
	}
	return true;
}

/**
 * Robust mesh Laplacian, after the robust_laplacian package: an intrinsically
 * mollified tufted Laplacian (positive semi-definite) and a lumped mass matrix.
 */
inline bool mesh_laplacian(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
		Eigen::SparseMatrix<double> &L, Eigen::SparseMatrix<double> &M, double mollify_factor = 1e-5) {
	std::vector<std::vector<size_t>> polygons(F.rows());
	for (int i = 0; i < F.rows(); ++i) {
		polygons[i] = { (size_t)F(i, 0), (size_t)F(i, 1), (size_t)F(i, 2) };
	}
	gc::surface::SurfaceMesh mesh(polygons);
	if ((int)mesh.nVertices() != V.rows()) {
		return false; // Unreferenced trailing vertices; the matrices would not match V.
	}
	gc::surface::VertexData<gc::Vector3> positions(mesh);
	for (int i = 0; i < V.rows(); ++i) {
		positions[i] = gc::Vector3{ V(i, 0), V(i, 1), V(i, 2) };
	}
	gc::surface::VertexPositionGeometry geometry(mesh, positions);
	std::tie(L, M) = gc::surface::buildTuftedLaplacian(mesh, geometry, mollify_factor);
	return true;
}

/**
 * Inpaint weights for the target vertices that found no good match.
 *
 *  V2: #V2 by 3 target vertices
 *  F2: #F2 by 3 target triangles
 *  W2: #V2 by num_bones weights copied from the source
 *  Matched: #V2 flags
 *  W_inpainted: #V2 by num_bones final weights
 *  Returns true on success.
 */
inline bool inpaint(const Eigen::MatrixXd &V2, const Eigen::MatrixXi &F2, const Eigen::MatrixXd &W2,
		const BoolArray &Matched, Eigen::MatrixXd &W_inpainted) {
	if (V2.cols() != 3 || F2.cols() != 3 || F2.rows() == 0 || F2.maxCoeff() >= V2.rows() || F2.minCoeff() < 0) {
		return false;
	}
	if (W2.rows() != V2.rows() || Matched.size() != V2.rows()) {
		return false;
	}
	if (!Matched.any()) {
		return false; // Nothing to anchor the solve to.
	}

	Eigen::SparseMatrix<double> L, M;
	if (!mesh_laplacian(V2, F2, L, M)) {
		return false;
	}
	L = -L; // Flip the sign of the Laplacian.

	Eigen::SparseMatrix<double> Minv;
	igl::invert_diag(M, Minv);
	Eigen::SparseMatrix<double> Q = -L + L * Minv * L;
	if (!is_valid_array(Q)) {
		return false;
	}

	Eigen::SparseMatrix<double> Aeq;
	Eigen::VectorXd Beq;
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(L.rows(), W2.cols());

	Eigen::VectorXi b_all = Eigen::VectorXi::LinSpaced(V2.rows(), 0, V2.rows() - 1);
	Eigen::VectorXi b;
	igl::slice_mask(b_all, Matched, 1, b);
	Eigen::MatrixXd bc;
	igl::slice_mask(W2, Matched, 1, bc);

	igl::min_quad_with_fixed_data<double> mqwf;
	if (!igl::min_quad_with_fixed_precompute(Q, b, Aeq, true, mqwf)) {
		return false;
	}
	return igl::min_quad_with_fixed_solve(mqwf, B, bc, Beq, W_inpainted);
}

/**
 * Smooth weights in the inpainted areas and their close neighbours.
 *
 *  W2_smoothed: #V2 by num_bones new weights
 *  VIDs_to_smooth: #V2 flags for the vertices that were smoothed
 *  V2, F2: target mesh
 *  W2: #V2 by num_bones weights to smooth
 *  Matched: #V2 flags
 *  dDISTANCE_THRESHOLD: neighbours closer than this to an unmatched vertex are smoothed too
 */
inline void smooth(Eigen::MatrixXd &W2_smoothed, BoolArray &VIDs_to_smooth, const Eigen::MatrixXd &V2,
		const Eigen::MatrixXi &F2, const Eigen::MatrixXd &W2, const BoolArray &Matched,
		const double dDISTANCE_THRESHOLD, const int num_smooth_iter_steps = 10, const double smooth_alpha = 0.2) {
	const BoolArray NotMatched = !Matched;
	VIDs_to_smooth = NotMatched;

	std::vector<std::vector<int>> adj_list;
	igl::adjacency_list(F2, adj_list);
	const int adj_size = (int)adj_list.size();

	// Mark every vertex reachable from VID through edges whose far end is within the distance of VID.
	auto mark_points_within_distance = [&](const int VID, const double distance) {
		if (VID < 0 || VID >= adj_size) {
			return;
		}
		std::queue<int> queue;
		queue.push(VID);
		std::set<int> visited{ VID };
		while (!queue.empty()) {
			const int current = queue.front();
			queue.pop();
			if (current >= adj_size) {
				continue;
			}
			for (const int neighbor : adj_list[current]) {
				if (visited.find(neighbor) == visited.end() && (V2.row(VID) - V2.row(neighbor)).norm() < distance) {
					visited.insert(neighbor);
					queue.push(neighbor);
				}
			}
		}
		for (const int v : visited) {
			VIDs_to_smooth[v] = true;
		}
	};

	for (int i = 0; i < V2.rows(); ++i) {
		if (NotMatched[i]) {
			mark_points_within_distance(i, dDISTANCE_THRESHOLD);
		}
	}

	W2_smoothed = W2;
	for (int step = 0; step < num_smooth_iter_steps; ++step) {
		for (int i = 0; i < V2.rows(); ++i) {
			if (!VIDs_to_smooth[i] || i >= adj_size || adj_list[i].empty()) {
				continue;
			}
			const std::vector<int> &neighbors = adj_list[i];
			const double count = (double)neighbors.size();
			Eigen::VectorXd new_weight = (1.0 - smooth_alpha) * W2_smoothed.row(i).transpose();
			for (const int n : neighbors) {
				new_weight += (W2_smoothed.row(n).transpose() / count) * smooth_alpha;
			}
			W2_smoothed.row(i) = new_weight.transpose();
		}
	}
}

} // namespace rwt
