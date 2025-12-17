```markdown
# testKDTree

## Purpose
Validation and benchmarking application for the KDTree library.

The test suite covers:
- Correctness vs brute-force search
- Stability on degenerate datasets
- Performance scaling up to 5M points
- Distribution sensitivity (uniform vs clustered)

## Build

```bash
wmake
