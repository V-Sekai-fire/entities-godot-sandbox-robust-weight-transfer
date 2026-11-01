#!/usr/bin/env -S godot --headless -s
# Loads the guest ELF into a Sandbox and checks it: its own tests, then the entry points on a known case.
#
#   godot --headless --path <project with the godot_sandbox addon> -s <this file> -- <robust_weight_transfer.elf>
#
# Ends with "RESULT: PASS" or "RESULT: FAIL". The weight arrays the guest returns are column-major.
extends SceneTree

var failed := 0


func _initialize() -> void:
	var args := OS.get_cmdline_user_args()
	if args.is_empty():
		print("FAIL: no ELF given")
		quit(2)
		return
	var sb := Sandbox.new()
	sb.set("memory_max", 4096)
	sb.set("references_max", 4096)
	sb.set("allocations_max", 1000000)
	sb.set_redirect_stdout(func(line): print(line))
	sb.load_buffer(FileAccess.get_file_as_bytes(args[0]))

	for fn in ["run_tests", "run_property_tests"]:
		_expect(fn, sb.has_function(fn) and sb.vmcall(fn) == 0)

	var call_args := {"verbose": false, "angle_threshold_degrees": 30.0, "distance_threshold": 0.05 * sqrt(3.0),
		"source_mesh_surface": 0, "target_mesh_surface": 0}
	# The right two columns of the target are out of reach of the source: 15 of 25 vertices match.
	var matched := []
	var interpolated := []
	var inpainted := []
	var smoothed := []
	var ok = sb.vmcall("robust_weight_transfer", _grid(5, 0.0), _grid(5, 1.0), call_args, matched, interpolated, inpainted, smoothed)
	_expect("robust_weight_transfer returns true", ok == true)
	_expect("matched 15 of 25", _count(matched) == 15 and matched.size() == 25)
	_expect("100 inpainted weights", inpainted.size() == 100)
	var sums_ok := true
	for i in 25:
		var sum := 0.0
		for bone in 4:
			sum += inpainted[i + bone * 25]
		if absf(sum - 1.0) > 1e-4:
			sums_ok = false
	_expect("each vertex's weights sum to 1", sums_ok)

	# Control: a target equal to the source matches everywhere.
	matched = []
	sb.vmcall("robust_weight_transfer", _grid(5, 0.0), _grid(5, 0.0), call_args, matched, [], [], [])
	_expect("identical meshes match 25 of 25", _count(matched) == 25)

	# Cage around a unit cube's corners: a closed hull has faces = 2 * vertices - 4.
	var cube := PackedVector3Array()
	for x in 2:
		for y in 2:
			for z in 2:
				cube.append(Vector3(x, y, z))
	var cube_arrays := []
	cube_arrays.resize(Mesh.ARRAY_MAX)
	cube_arrays[Mesh.ARRAY_VERTEX] = cube
	var cube_mesh := ArrayMesh.new()
	cube_mesh.add_surface_from_arrays(Mesh.PRIMITIVE_POINTS, cube_arrays)
	var cage_vertices := []
	var cage_faces := []
	ok = sb.vmcall("create_cage", cube_mesh, {"target_mesh_surface": 0, "num_points": 8, "seed": 3}, cage_vertices, cage_faces)
	var vertex_count: int = cage_vertices.size() / 3
	_expect("create_cage closed hull", ok == true and vertex_count >= 8 and cage_faces.size() / 3 == 2 * vertex_count - 4)

	print("RESULT: %s" % ("PASS" if failed == 0 else "FAIL"))
	quit(0 if failed == 0 else 1)


func _expect(what: String, condition: bool) -> void:
	print("%s %s" % ["PASS" if condition else "FAIL", what])
	if not condition:
		failed += 1


func _count(values: Array) -> int:
	var n := 0
	for v in values:
		if v:
			n += 1
	return n


# An n by n grid in the plane z = 0, the vertices right of x = 0.5 lifted by lift_right, weights (x, 1 - x, 0, 0).
func _grid(n: int, lift_right: float) -> ArrayMesh:
	var vertices := PackedVector3Array()
	var normals := PackedVector3Array()
	var weights := PackedFloat32Array()
	var bones := PackedInt32Array()
	var indices := PackedInt32Array()
	for y in n:
		for x in n:
			var fx := x / float(n - 1)
			vertices.append(Vector3(fx, y / float(n - 1), lift_right if fx > 0.5 else 0.0))
			normals.append(Vector3(0, 0, 1))
			weights.append_array(PackedFloat32Array([fx, 1.0 - fx, 0.0, 0.0]))
			bones.append_array(PackedInt32Array([0, 1, 2, 3]))
	for y in n - 1:
		for x in n - 1:
			var a := y * n + x
			indices.append_array(PackedInt32Array([a, a + 1, a + n, a + 1, a + n + 1, a + n]))
	var arrays := []
	arrays.resize(Mesh.ARRAY_MAX)
	arrays[Mesh.ARRAY_VERTEX] = vertices
	arrays[Mesh.ARRAY_NORMAL] = normals
	arrays[Mesh.ARRAY_WEIGHTS] = weights
	arrays[Mesh.ARRAY_BONES] = bones
	arrays[Mesh.ARRAY_INDEX] = indices
	var mesh := ArrayMesh.new()
	mesh.add_surface_from_arrays(Mesh.PRIMITIVE_TRIANGLES, arrays)
	return mesh
