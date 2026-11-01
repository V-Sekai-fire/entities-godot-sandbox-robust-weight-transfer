# Robust Skin Weights Transfer via Weight Inpainting

![Teaser](https://www.dgp.toronto.edu/~rinat/projects/RobustSkinWeightsTransfer/teaser.jpg)
Sample code for the Siggraph Asia 2023 Technical Communications paper - [Robust Skin Weights Transfer via Weight Inpainting](https://www.dgp.toronto.edu/~rinat/projects/RobustSkinWeightsTransfer/index.html)

## Building the guest
The guest builds with the clang and lld in the `riscv64-sysroot` checkout (pixi) against `guest-runtime`'s sandbox API,
both beside this repository, and is built at double precision (`.sandbox_variant` reads `SBXV 28000000`):
```sh
(cd ../../../riscv64-sysroot && pixi install)
./compile.sh            # WS=<folder holding riscv64-sysroot and guest-runtime> if they are elsewhere
```
The ELF lands next to `robust_weight_transfer/`, named `robust_weight_transfer`; the project loads it as `robust_weight_transfer.elf`.

## Dependencies
[Eigen](https://eigen.tuxfamily.org), [libigl](https://libigl.github.io) and
[geometry-central](https://geometry-central.net) (for the robust tufted Laplacian), all under `thirdparty/`.
The tests use [witness-cpp](https://github.com/V-Sekai-fire/repository-witness-cpp).

## Layout
- `src/transfer.hpp`: closest point matching, weight inpainting and smoothing.
- `src/cage.hpp`: the diamond, a convex hull, and the cage built from them.
- `src/self_test.hpp`, `src/property_test.hpp`: example tests and property tests.
- `robust_weight_transfer.cpp`: the guest entry points `robust_weight_transfer`, `create_cage`, `run_tests` and `run_property_tests`.

## Running the tests on the host
No RISC-V toolchain or engine is needed:
```sh
cmake -S host_test -B build-host -GNinja -DCMAKE_BUILD_TYPE=Release
cmake --build build-host
./build-host/host_test
```
The same tests run in the guest through the `run_tests` and `run_property_tests` calls, each inside the default instruction cap. The guest uses a shorter ladder of trials than the host. The weight arrays the guest returns are column-major: all vertices for bone 0, then bone 1, and so on.

## Checking the built guest
`tests/run_guest.gd` loads the ELF into a Sandbox in a headless double Godot editor with the godot_sandbox addon and runs the guest's tests and its entry points on a known case:
```sh
godot --headless --path <project with the addon> -s <path to>/tests/run_guest.gd -- <robust_weight_transfer.elf>
```
It ends with `RESULT: PASS` or `RESULT: FAIL`.

## Body to garment transfer
`create_cage` builds a cage around a garment mesh, and `robust_weight_transfer` moves the skin weights of
a body onto it.

## Other 3rd party implementations

Blender addon (by Sent From Space) - https://sentfromspacevr.gumroad.com/l/robust-weight-transfer

## Cite

If you use this code for an academic publication, cite it as:

```bib
@inproceedings{abdrashitov2023robust,
author = {Abdrashitov, Rinat and Raichstat, Kim and Monsen, Jared and Hill, David},
title = {Robust Skin Weights Transfer via Weight Inpainting},
year = {2023},
isbn = {9798400703140},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3610543.3626180},
doi = {10.1145/3610543.3626180},
booktitle = {SIGGRAPH Asia 2023 Technical Communications},
articleno = {25},
numpages = {4},
location = {<conf-loc>, <city>Sydney</city>, <state>NSW</state>, <country>Australia</country>, </conf-loc>},
series = {SA '23}
}
```
