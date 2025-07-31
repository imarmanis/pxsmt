# PxSMT

An old SMT-based Model Checking project that handles shared-memory persistency and consistency semantics via a custom theory solver,
implemented as Z3 theory propagator.

Persistency is verified via invariants of a read-only `rec()` function that can observe any possible post-crash state.

The expected inputs are (limited fragments of) C programs that have been compiled into LLVM-IR.

# Requirements

- LLVM/Clang toolchain version 11
- Z3 version 4.8.17
- Boost

# Build Instructions

Running `make` produces the `main` executable

# Example Input
```
#include "verify.h"

unsigned i = 0;
int y = 0;
int a[3];

void t() {
    faa(&i, 1);
}

void rec() {
    assert((a[i] == 1) || (y!= 1));
}

int main() {
    i = 0;
    y = 0;
    a[0] = 1;
    a[1] = 0;
    a[2] = 0;

    flush(&a[0]);
    flush(&a[1]);
    y = 1;
    par(t, t);
}
```
and `clang -S -O0 -Xclang -disable-O0-optnone -emit-llvm -mclflushopt -Iinclude input.c`
