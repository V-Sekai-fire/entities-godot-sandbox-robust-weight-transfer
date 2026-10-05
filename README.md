# entities-godot-sandbox-robust-weight-transfer

An engine project that runs robust skin weight transfer, by weight inpainting, inside a sandboxed RISC-V guest.

## What it is for

The main scene hands a skinned avatar mesh and a target mesh to the guest, which matches vertices, interpolates their weights, inpaints the weights the match missed, and smooths the result. The guest's C++ is a subrepo of `godot-robust-skin-weights-transfer`, after the paper "Robust Skin Weights Transfer via Weight Inpainting".

## Build and run

The guest compiles with a riscv64 GCC toolchain:

```sh
cd robust_weight_transfer/robust_skin_weight_transfer && ./compile.sh
```

Then open the project in the editor and run its main scene.

## Licence

MIT; see LICENSE.md. The transfer code carries its own MIT LICENSE.
