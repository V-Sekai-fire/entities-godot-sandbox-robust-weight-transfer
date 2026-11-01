#!/bin/sh
# Builds robust_weight_transfer.elf with the clang in riscv64-sysroot (pixi) against guest-runtime's sandbox API.
# Expects the manifest layout: riscv64-sysroot and guest-runtime beside this repository.
set -e
here=$(cd "$(dirname "$0")" && pwd)
ws=$(cd "$here/../../.." && pwd)
PATH="$ws/riscv64-sysroot/.pixi/envs/default/bin:$PATH"
cmake -S "$here" -B "$here/build" -G Ninja -DCMAKE_BUILD_TYPE=Release -DUSE_DOUBLE_PRECISION=ON \
  -DCMAKE_TOOLCHAIN_FILE="$ws/riscv64-sysroot/toolchain.cmake"
cmake --build "$here/build"
