# KDTree (OpenFOAM Native)

## Description
High-performance, memory-efficient KDTree implementation designed for
CFD-scale datasets and ML coupling inside OpenFOAM.

Features:
- O(N) average build using QuickSelect
- O(log N) average nearest-neighbour queries (low dimensionality)
- Leaf-bucket optimization
- No STL containers in public API
- Fully OpenFOAM-native types

## Build Instructions

From this directory:

```bash
wmake libso
