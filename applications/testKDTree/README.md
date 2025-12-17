# testKDTree

## Purpose
`testKDTree` is a validation and benchmarking application for the OpenFOAM-native
`KDTree` library.

It is designed to **verify correctness, robustness, and performance**
of nearest-neighbour queries under conditions representative of
CFD and ML-coupled workflows.

The application intentionally separates **accuracy validation**
from **performance measurement** to remain tractable at large scales.

---

## What Is Tested

### 1. Correctness
Nearest-neighbour results are validated against a brute-force search:

- Full brute-force validation for small datasets
- Subsampled brute-force validation for large datasets (up to millions of points)
- Exact squared-distance matching within numerical tolerance

### 2. Robustness
The following pathological and non-ideal cases are explicitly tested:

- Very small datasets
- Duplicate points
- Uniform spatial distributions
- Strongly clustered / wake-like distributions

### 3. Performance & Scaling
The benchmark suite evaluates:

- Tree build time
- Average query time
- Speedup relative to brute-force search
- Scaling behaviour up to **5 million points (3D)**

---

## Build

From this directory:

```bash
wmake
