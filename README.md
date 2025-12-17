# algorithms

This repository contains **standalone, OpenFOAM-native algorithm implementations**
intended for reuse in CFD workflows and ML–CFD coupling.

# Important

All code in this repository is intended for research and development use.
Users are responsible for validating suitability, performance,
and resource usage in their own environments.



All algorithms are implemented in **pure OpenFOAM C++ style** and are designed to:
- Avoid external dependencies
- Scale to large datasets
- Be suitable for integration into solvers, utilities, or ML pipelines

Each algorithm lives in `src/` and may optionally provide a corresponding
validation or benchmark application in `applications/`.

---

## Current Contents

- **KDTree**
  - Exact nearest-neighbour search
  - Optimized for large point clouds (millions of points)
  - Validated against brute-force search
  - Benchmark application included

---

## Build

Each module follows standard OpenFOAM conventions.

Example:

```bash
wmake libso src/KDTree
wmake applications/testKDTree
