#include <api.hpp>

#include <cstdint>
#include <iostream>
#include <vector>

#include "generated_api.hpp"
#include "src/cage.hpp"
#include "src/property_test.hpp"
#include "src/transfer.hpp"

using namespace rwt;

extern "C" Variant robust_weight_transfer(Mesh source_mesh, Mesh target_mesh, Dictionary arguments, Array matched_array, Array interpolated_weights_array, Array inpainted_weights_array, Array smoothed_weights_array) {
    if (!arguments.has("verbose") || !arguments.has("angle_threshold_degrees") || !arguments.has("distance_threshold") || !arguments.has("source_mesh_surface") || !arguments.has("target_mesh_surface")) {
        std::cerr << "Missing required arguments" << std::endl;
        return false;
    }

    bool verbose = arguments["verbose"].value();
    double angle_threshold_degrees = arguments["angle_threshold_degrees"].value();
    double distance_threshold = arguments["distance_threshold"].value();
    int64_t source_mesh_surface = arguments["source_mesh_surface"].value();
    int64_t target_mesh_surface = arguments["target_mesh_surface"].value();

    if (verbose) { std::cout << "arguments:\n" << std::endl; }

    Array source_mesh_arrays = source_mesh.surface_get_arrays(source_mesh_surface);
    if (source_mesh_arrays.size() <= Mesh::ARRAY_VERTEX || source_mesh_arrays.size() <= Mesh::ARRAY_INDEX || source_mesh_arrays.size() <= Mesh::ARRAY_NORMAL || source_mesh_arrays.size() <= Mesh::ARRAY_WEIGHTS) {
        std::cerr << "Source mesh arrays are incomplete" << std::endl;
        return false;
    }

    Variant vertices_1_variant = source_mesh_arrays[Mesh::ARRAY_VERTEX].get();
    if (vertices_1_variant == Nil) {
        std::cerr << "vertices_1_variant is null" << std::endl;
        return false;
    }
    PackedArray<Vector3> vertices_1_ref = vertices_1_variant;

    Variant faces_1_variant = source_mesh_arrays[Mesh::ARRAY_INDEX].get();
    if (faces_1_variant == Nil) {
        std::cerr << "faces_1_variant is null" << std::endl;
        return false;
    }
    PackedArray<int32_t> faces_1_ref = faces_1_variant;

    Variant normals_1_variant = source_mesh_arrays[Mesh::ARRAY_NORMAL].get();
    if (normals_1_variant == Nil) {
        std::cerr << "normals_1_variant is null" << std::endl;
        return false;
    }
    PackedArray<Vector3> normals_1_ref = normals_1_variant;

    Variant skin_weights_variant = source_mesh_arrays[Mesh::ARRAY_WEIGHTS].get();
    if (skin_weights_variant == Nil) {
        std::cerr << "skin_weights_variant is null" << std::endl;
        return false;
    }
    PackedArray<float> skin_weights_ref = skin_weights_variant;

    if (verbose) { std::cout << "source_mesh_arrays:\n" << std::endl; }

    if (vertices_1_ref.is_empty() || faces_1_ref.is_empty() || normals_1_ref.is_empty() || skin_weights_ref.is_empty()) {
        std::cerr << "One or more source mesh arrays are empty" << std::endl;
        return false;
    }

    std::vector<Vector3> vertices_1 = vertices_1_ref.fetch();
    std::vector<int32_t> faces_1 = faces_1_ref.fetch();
    std::vector<Vector3> normals_1 = normals_1_ref.fetch();
    std::vector<float> skin_weights_1 = skin_weights_ref.fetch();

    Array target_mesh_arrays = target_mesh.surface_get_arrays(target_mesh_surface);
    if (target_mesh_arrays.size() <= Mesh::ARRAY_VERTEX || target_mesh_arrays.size() <= Mesh::ARRAY_INDEX || target_mesh_arrays.size() <= Mesh::ARRAY_NORMAL) {
        std::cerr << "Target mesh arrays are incomplete" << std::endl;
        return false;
    }
    Variant vertices_2_variant = target_mesh_arrays[Mesh::ARRAY_VERTEX].get();
    if (vertices_2_variant == Nil) {
        std::cerr << "vertices_2_variant is null" << std::endl;
        return false;
    }
    PackedArray<Vector3> vertices_2_ref = vertices_2_variant;

    Variant faces_2_variant = target_mesh_arrays[Mesh::ARRAY_INDEX].get();
    if (faces_2_variant == Nil) {
        std::cerr << "faces_2_variant is null" << std::endl;
        return false;
    }
    PackedArray<int32_t> faces_2_ref = faces_2_variant;

    Variant normals_2_variant = target_mesh_arrays[Mesh::ARRAY_NORMAL].get();
    if (normals_2_variant == Nil) {
        std::cerr << "normals_2_variant is null" << std::endl;
        return false;
    }
    PackedArray<Vector3> normals_2_ref = normals_2_variant;

    if (verbose) { std::cout << "target_mesh_arrays:\n" << std::endl; }

    if (vertices_2_ref.is_empty() || faces_2_ref.is_empty() || normals_2_ref.is_empty()) {
        std::cerr << "One or more target mesh arrays are empty" << std::endl;
        return false;
    }

    std::vector<Vector3> vertices_2 = vertices_2_ref.fetch();
    std::vector<int32_t> faces_2 = faces_2_ref.fetch();
    std::vector<Vector3> normals_2 = normals_2_ref.fetch();

    Eigen::MatrixXd vertices_1_eigen(vertices_1.size(), 3);
    for (int i = 0; i < vertices_1.size(); ++i) {
        Vector3 v = vertices_1[i];
        vertices_1_eigen(i, 0) = v.x;
        vertices_1_eigen(i, 1) = v.y;
        vertices_1_eigen(i, 2) = v.z;
    }
    if (verbose) { std::cout << "vertices_1_eigen:\n" << vertices_1_eigen << std::endl; }

    Eigen::MatrixXi faces_1_eigen(faces_1.size() / 3, 3);
    for (int i = 0; i < faces_1.size() / 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            faces_1_eigen(i, j) = faces_1[i * 3 + j];
        }
    }
    if (verbose) { std::cout << "faces_1_eigen:\n" << faces_1_eigen << std::endl; }

    Eigen::MatrixXd normals_1_eigen(normals_1.size(), 3);
    for (int i = 0; i < normals_1.size(); ++i) {
        Vector3 n = normals_1[i];
        normals_1_eigen(i, 0) = n.x;
        normals_1_eigen(i, 1) = n.y;
        normals_1_eigen(i, 2) = n.z;
    }
    if (verbose) { std::cout << "normals_1_eigen:\n" << normals_1_eigen << std::endl; }

    if (skin_weights_1.empty()) {
        std::cerr << "skin_weights array is empty" << std::endl;
        return false;
    }
    int32_t num_bones = skin_weights_1.size() / vertices_1.size();
    Eigen::MatrixXd skin_weights_eigen(vertices_1.size(), num_bones);
    for (int i = 0; i < vertices_1.size(); ++i) {
        for (int j = 0; j < num_bones; ++j) {
            skin_weights_eigen(i, j) = skin_weights_1[i * num_bones + j];
        }
    }

    if (verbose) { std::cout << "skin_weights_eigen:\n" << skin_weights_eigen << std::endl; }

    Eigen::MatrixXd vertices_2_eigen(vertices_2.size(), 3);
    for (int i = 0; i < vertices_2.size(); ++i) {
        Vector3 v = vertices_2[i];
        vertices_2_eigen(i, 0) = v.x;
        vertices_2_eigen(i, 1) = v.y;
        vertices_2_eigen(i, 2) = v.z;
    }
    if (verbose) { std::cout << "vertices_2_eigen:\n" << vertices_2_eigen << std::endl; }

    Eigen::MatrixXi faces_2_eigen(faces_2.size() / 3, 3);
    for (int i = 0; i < faces_2.size() / 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            faces_2_eigen(i, j) = faces_2[i * 3 + j];
        }
    }
    if (verbose) { std::cout << "faces_2_eigen:\n" << faces_2_eigen << std::endl; }

    Eigen::MatrixXd normals_2_eigen(normals_2.size(), 3);
    for (int i = 0; i < normals_2.size(); ++i) {
        Vector3 n = normals_2[i];
        normals_2_eigen(i, 0) = n.x;
        normals_2_eigen(i, 1) = n.y;
        normals_2_eigen(i, 2) = n.z;
    }
    if (verbose) { std::cout << "normals_2_eigen:\n" << normals_2_eigen << std::endl; }

    // Section 3.1 Closest Point Matching
    if (verbose) { std::cout << "Distance threshold: " << distance_threshold << std::endl; }

    Eigen::MatrixXd W2_eigen = Eigen::MatrixXd::Zero(vertices_2_eigen.rows(), skin_weights_eigen.cols());
    Eigen::Array<bool, Eigen::Dynamic, 1> Matched_eigen;

    find_matches_closest_surface(vertices_1_eigen, faces_1_eigen, normals_1_eigen, vertices_2_eigen, faces_2_eigen, normals_2_eigen, skin_weights_eigen, distance_threshold * distance_threshold, angle_threshold_degrees, W2_eigen, Matched_eigen);
    if (verbose) { std::cout << "Matched_eigen:\n" << Matched_eigen << std::endl; }
    if (verbose) { std::cout << "W2_eigen:\n" << W2_eigen << std::endl; }

    // Section 3.2 Skinning Weights Inpainting
    Eigen::MatrixXd W_inpainted;
    bool success = inpaint(vertices_2_eigen, faces_2_eigen, W2_eigen, Matched_eigen, W_inpainted);
    if (verbose) { std::cout << "Inpainting success: " << success << std::endl; }
    if (verbose) { std::cout << "W_inpainted:\n" << W_inpainted << std::endl; }

    // Optional smoothing
    Eigen::MatrixXd W2_smoothed;
    Eigen::Array<bool, Eigen::Dynamic, 1> VIDs_to_smooth;
    smooth(W2_smoothed, VIDs_to_smooth, vertices_2_eigen, faces_2_eigen, W_inpainted, Matched_eigen, distance_threshold, 10, 0.2);
    if (verbose) { std::cout << "W2_smoothed:\n" << W2_smoothed << std::endl; }
    if (verbose) { std::cout << "VIDs_to_smooth:\n" << VIDs_to_smooth << std::endl; }
    
    if (verbose) { std::cout << "Matched: " << Matched_eigen.transpose() << std::endl; }
    if (verbose) { std::cout << "Interpolated Skin Weights: " << W2_eigen << std::endl; }
    if (verbose) { std::cout << "Inpainted Weights: " << W_inpainted << std::endl; }
    if (verbose) { std::cout << "Smoothed Inpainted Weights: " << W2_smoothed << std::endl; }

    for (int i = 0; i < Matched_eigen.size(); ++i) {
        matched_array.append(Matched_eigen(i));
    }
    if (verbose) { std::cout << "Matched array appended." << std::endl; }
    for (int i = 0; i < W2_eigen.size(); ++i) {
        interpolated_weights_array.append(Variant(W2_eigen(i)));
    }
    if (verbose) { std::cout << "Interpolated weights array appended." << std::endl; }

    for (int i = 0; i < W_inpainted.size(); ++i) {
        inpainted_weights_array.append(Variant(W_inpainted(i)));
    }
    if (verbose) { std::cout << "Inpainted weights array appended." << std::endl; }

    for (int i = 0; i < W2_smoothed.size(); ++i) {
        smoothed_weights_array.append(Variant(W2_smoothed(i)));
    }
    if (verbose) { std::cout << "Smoothed weights array appended." << std::endl; }

    return true;
}

/**
 * Build a cage around a mesh: the convex hull of its vertices and of a small diamond on
 * each of num_points randomly chosen vertices.
 *
 * arguments: target_mesh_surface (int), num_points (int, default 200), seed (int, default 0)
 * out_vertices: flat x, y, z values of the cage vertices
 * out_faces: flat vertex indices, three per triangle
 */
extern "C" Variant create_cage(Mesh target_mesh, Dictionary arguments, Array out_vertices, Array out_faces) {
	if (!arguments.has("target_mesh_surface")) {
		std::cerr << "Missing required arguments" << std::endl;
		return false;
	}
	int64_t target_mesh_surface = arguments["target_mesh_surface"].value();
	int64_t num_points = arguments.has("num_points") ? (int64_t)arguments["num_points"].value() : 200;
	int64_t seed = arguments.has("seed") ? (int64_t)arguments["seed"].value() : 0;

	Array mesh_arrays = target_mesh.surface_get_arrays(target_mesh_surface);
	if (mesh_arrays.size() <= Mesh::ARRAY_VERTEX) {
		std::cerr << "Target mesh arrays are incomplete" << std::endl;
		return false;
	}
	Variant vertices_variant = mesh_arrays[Mesh::ARRAY_VERTEX].get();
	if (vertices_variant == Nil) {
		std::cerr << "vertices_variant is null" << std::endl;
		return false;
	}
	PackedArray<Vector3> vertices_ref = vertices_variant;
	if (vertices_ref.is_empty()) {
		std::cerr << "Target mesh vertices are empty" << std::endl;
		return false;
	}
	std::vector<Vector3> vertices = vertices_ref.fetch();

	Eigen::MatrixXd vertices_eigen(vertices.size(), 3);
	for (size_t i = 0; i < vertices.size(); ++i) {
		vertices_eigen(i, 0) = vertices[i].x;
		vertices_eigen(i, 1) = vertices[i].y;
		vertices_eigen(i, 2) = vertices[i].z;
	}

	Eigen::MatrixXd cage_vertices;
	Eigen::MatrixXi cage_faces;
	if (!rwt::create_cage(vertices_eigen, (int)num_points, (uint32_t)seed, cage_vertices, cage_faces)) {
		std::cerr << "The vertices do not span a volume" << std::endl;
		return false;
	}
	for (int i = 0; i < cage_vertices.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			out_vertices.append(Variant(cage_vertices(i, j)));
		}
	}
	for (int i = 0; i < cage_faces.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			out_faces.append(Variant((int64_t)cage_faces(i, j)));
		}
	}
	return true;
}

// The example tests: fast, well inside the per-call instruction cap.
extern "C" Variant run_tests() {
	const int failed = run_self_tests();
	if (failed == 0) {
		std::cout << "All tests passed!" << std::endl;
		return Variant(0);
	}
	std::cerr << failed << " tests failed." << std::endl;
	return Variant(1);
}

// The property tests, in their own call so each stays inside the instruction cap.
extern "C" Variant run_property_tests() {
	const int failed = rwt::run_property_tests();
	if (failed == 0) {
		std::cout << "All property tests passed!" << std::endl;
		return Variant(0);
	}
	std::cerr << failed << " property tests failed." << std::endl;
	return Variant(1);
}
