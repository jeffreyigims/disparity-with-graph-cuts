# KZ2: Stereo Disparity Estimation via Graph Cuts

## Overview

This repository implements the **Kolmogorov-Zabih (KZ2) stereo matching algorithm**, which computes a dense disparity map from a rectified stereo image pair. Disparity maps encode per-pixel depth information: given two horizontally-aligned images of the same scene (left and right), the algorithm finds how far each pixel in the left image must shift horizontally to match the corresponding pixel in the right image. Larger disparity = closer to the camera.

The algorithm frames stereo matching as an **energy minimization** problem and solves it exactly (for each expansion step) using **graph cuts** (min-cut/max-flow). It handles **occlusions** explicitly -- pixels visible in one image but not the other are labeled as occluded rather than forced into a bad match.

### Academic References

| # | Paper | Authors | Venue |
|---|-------|---------|-------|
| 1 | "Kolmogorov and Zabih's Graph Cuts Stereo Matching Algorithm" | Kolmogorov, Monasse, Tan (2014) | IPOL |
| 2 | "Computing Visual Correspondence with Occlusions using Graph Cuts" | Kolmogorov, Zabih (2001) | ICCV/IJCV |
| 3 | "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision" | Boykov, Kolmogorov (2004) | IEEE PAMI |

---

## Repository Structure

```
disparity-with-graph-cuts/
|-- CMakeLists.txt              # Top-level CMake build (project: KZ2)
|-- README.txt                  # Original documentation
|-- LICENSE.txt                 # GPLv3
|-- images/
|   |-- scene_l.png             # Sample left stereo image
|   |-- scene_r.png             # Sample right stereo image
|   `-- disp.png                # Reference disparity output
`-- src/
    |-- CMakeLists.txt          # Build config, dependency resolution
    |-- main.cpp                # Entry point, CLI parsing, parameter setup
    |-- match.h / match.cpp     # Core Match class: disparity storage, I/O, orchestration
    |-- kz2.cpp                 # Alpha-expansion loop + graph construction
    |-- data.cpp                # Data cost (Birchfield-Tomasi) + smoothness penalties
    |-- statistics.cpp          # Automatic computation of occlusion cost K
    |-- image.h / image.cpp     # Custom image data structure + PNG/TIFF/PGM I/O
    |-- cmdLine.h               # Lightweight command-line parser
    |-- nan.h                   # NaN constant for marking occluded pixels
    |-- io_png.h / io_png.c     # libPNG wrapper (C)
    |-- io_tiff.h / io_tiff.c   # libTIFF wrapper (C)
    |-- energy/
    |   |-- energy.h            # Energy minimization API (thin wrapper over graph)
    |   `-- test_energy.cpp     # Unit test for the energy module
    |-- maxflow/
    |   |-- graph.h             # Templated graph structure (nodes, arcs, BFS trees)
    |   |-- graph.cpp           # Graph construction (add_node, add_edge, add_tweights)
    |   `-- maxflow.cpp         # Boykov-Kolmogorov max-flow/min-cut algorithm
    `-- third_party/            # Vendored zlib, libPNG, libTIFF, libjpeg sources
```

---

## Algorithm Deep Dive

### 1. Energy Formulation

The algorithm minimizes a global energy function over all pixels:

```
E(d) = sum_p DataCost(p, d_p) + sum_{p,q neighbors} SmoothnessCost(p, q) + sum_p OcclusionCost(p)
```

- **Data cost** (`data.cpp`): How well does pixel `p` in the left image match pixel `p + d` in the right image? Uses the **Birchfield-Tomasi** dissimilarity measure (robust to sub-pixel sampling). Supports both **L1** (absolute difference) and **L2** (squared difference) norms. A hard cutoff (`CUTOFF=30`) caps outlier penalties.

- **Smoothness cost** (`data.cpp`): Penalizes neighboring pixels that have different disparities. Uses a **Potts model** with two penalty levels:
  - `lambda1` (higher): applied when both left and right image intensities are similar across the pixel boundary (no edge detected)
  - `lambda2` (lower): applied when an intensity edge is detected (diff >= `edgeThresh`), allowing disparity discontinuities at object boundaries

- **Occlusion cost** (`K`): Fixed penalty for labeling a pixel as occluded (no match). Balances between forcing bad matches and allowing too many occlusions.

### 2. Alpha-Expansion (kz2.cpp)

The energy is minimized via **alpha-expansion moves**, an iterative approach:

1. Start with all pixels labeled as occluded
2. For each disparity label `alpha` in `[dMin, dMax]`:
   - Build a graph encoding the energy change if any subset of pixels switches to disparity `alpha`
   - Solve the **min-cut** (via max-flow) to find the optimal subset
   - If energy decreases, accept the move
3. Repeat until no label expansion reduces energy, or `maxIter` is reached
4. Optionally **randomize** the label order each iteration (`-r` flag)

The convergence detector uses a `done[]` array: once an expansion of label `alpha` fails to reduce energy, it is skipped until some other expansion succeeds (which resets all labels as candidates again).

### 3. Graph Construction (kz2.cpp)

For each alpha-expansion, a binary-labeled graph is built:

- **`build_nodes`**: Creates variables for each pixel. For pixel `p` with current disparity `d != alpha`:
  - A variable for keeping `(p, p+d)` active (in A^0)
  - A variable for activating `(p, p+alpha)` (in A^alpha)
  - If `d == alpha`, the assignment is constant (no variable needed)

- **`build_smoothness`**: Adds edges between neighboring pixels to enforce the smoothness term. Handles all combinations of current and candidate disparities.

- **`build_uniqueness`**: Adds infinite-capacity edges (via `forbid01`) to enforce the constraint that each pixel has at most one match, and each right-image pixel is matched by at most one left-image pixel.

### 4. Max-Flow / Min-Cut (maxflow/)

Implements the **Boykov-Kolmogorov algorithm** with two search trees growing from SOURCE and SINK:

1. **Grow**: Extend trees via BFS until they meet (augmenting path found)
2. **Augment**: Push flow along the path (bottleneck capacity)
3. **Adopt**: Reconnect orphaned nodes after saturated edges break tree paths

The graph is templated (`Graph<captype, tcaptype, flowtype>`) and instantiated as `Graph<short, short, int>` for this application.

---

## Explicit Optimizations and Efficiencies

### Memory Pre-Allocation
**File**: `kz2.cpp:183`
```cpp
Energy e(2*imSizeL.x*imSizeL.y, 12*imSizeL.x*imSizeL.y);
```
Pre-allocates `2n` nodes and `12n` edges (where `n` = pixel count) via `std::vector::reserve()`. These are the theoretical maximums, avoiding costly reallocations during graph construction. The README notes this trades elegance for speed compared to on-demand allocation.

### Integer Arithmetic via Fractional Denominators
**File**: `main.cpp:28-78`

All energy terms use **short integer** arithmetic internally (the graph capacities are `short`). Floating-point parameters (`K`, `lambda1`, `lambda2`) are approximated as fractions with a shared denominator (up to `MAX_DENOM = 16`). The denominator is chosen to minimize the sum of relative approximation errors. This avoids floating-point overhead in the inner loop and keeps values below `2^15` to prevent `short` overflow.

### Birchfield-Tomasi Sub-Pixel Preprocessing
**File**: `data.cpp:93-169`

Rather than computing intensity intervals on-the-fly for each pixel match, the algorithm **precomputes** per-pixel min/max intensity ranges (`imLeftMin`, `imLeftMax`, etc.) from 4-connected neighbors. This `InitSubPixel()` call happens once, making every subsequent `data_penalty_gray/color()` call a simple table lookup + comparison instead of recomputing neighbor interpolations.

### Intensity Cutoff
**File**: `data.cpp:42`
```cpp
static int CUTOFF = 30;
```
Clamps the maximum data penalty, preventing outlier pixels from dominating the energy function. Also keeps values small enough that `L2` cost (`CUTOFF^2 = 900 < 2^10`) fits comfortably within the `short` representation.

### Convergence Tracking with `done[]` Array
**File**: `kz2.cpp:235-258`

Instead of blindly iterating over all disparities each round, the algorithm tracks which labels have already been tried without improvement. When any expansion succeeds, all labels are re-enabled (since the configuration changed). This can skip many redundant max-flow computations in later iterations.

### Auto-Detection of Grayscale Images
**File**: `main.cpp:33-52`

Color images where all channels are identical are automatically converted to grayscale, reducing the data cost computation from 3-channel to 1-channel (3x fewer operations per pixel match).

### Inline Graph Operations
**Files**: `energy/energy.h`, `maxflow/graph.h`

All `Energy` methods and many `Graph` methods are declared `inline` in headers. The templated graph implementation is `#include`-ed directly into the header (not separately compiled), enabling full inlining and compiler optimization of the hot max-flow inner loops.

### Efficient Active Node Queue
**File**: `maxflow/maxflow.cpp:27-58`

The active node list uses an intrusive linked list (the `node::next` pointer doubles as both "next in queue" and "is in queue" flag). Nodes with null parents are lazily skipped rather than explicitly removed, avoiding costly list removal operations.

### Minimal Neighborhood Definition
**File**: `kz2.cpp:29-30`
```cpp
const struct Coord NEIGHBORS[] = { Coord(-1,0), Coord(0,1) };
```
Only half the 4-connected neighborhood is stored; the other half comes from reversed iteration. This halves the smoothness edge construction work while maintaining full connectivity.

### Heuristic K Computation
**File**: `statistics.cpp:25-61`

Instead of requiring manual tuning, `K` is automatically estimated by sampling the k-th smallest data penalty across all disparities for each pixel, then averaging. This adapts the occlusion cost to the actual image pair without user intervention.

---

## How to Build and Run

### Prerequisites

- **CMake** >= 2.6
- **C/C++ compiler** (GCC, Clang, or MSVC)
- Optional: system-installed **libPNG** and **libTIFF** (with `-dev` headers). If not found, CMake builds them from the vendored `third_party/` sources automatically.

### Build (Linux / macOS)

```bash
cd /path/to/disparity-with-graph-cuts
mkdir Build && cd Build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

### Build (Windows / MSVC)

1. Open CMake GUI, set source to the repo root, binaries to a new folder
2. Click "Configure" then "Generate", selecting your MSVC version
3. Open the generated `KZ2.sln` in Visual Studio
4. Set to "Release" mode, then Build Solution (F7)

### Run

```bash
bin/KZ2 [options] im1.png im2.png dMin dMax [dispMap.tif]
```

**Required arguments:**
| Argument | Description |
|----------|-------------|
| `im1.png` | Left stereo image (PNG or PGM) |
| `im2.png` | Right stereo image (PNG or PGM) |
| `dMin` | Minimum disparity (integer, can be negative) |
| `dMax` | Maximum disparity (integer) |
| `dispMap.tif` | (Optional) Output disparity map as float TIFF |

**Options:**
| Flag | Description | Default |
|------|-------------|---------|
| `-o disp.png` | Save viewable disparity map (gray = disparity, cyan = occluded) | none |
| `-i N` | Max iterations | 4 |
| `-r` | Randomize label order each iteration | false |
| `-c L1\|L2` | Data cost norm | L2 |
| `-l lambda` | Smoothness weight | K/5 |
| `--lambda1 l1` | Smoothness cost (no edge) | 3*lambda |
| `--lambda2 l2` | Smoothness cost (across edge) | lambda |
| `-t threshold` | Intensity diff threshold for edge detection | 8 |
| `-k K` | Occlusion penalty | auto-computed |

### Example

```bash
# From the Build directory:
bin/KZ2 ../images/scene_l.png ../images/scene_r.png -15 0 disp.tif -o disp.png
```

**Output:**
- `disp.tif` -- Float TIFF with raw disparity values (NaN for occluded pixels)
- `disp.png` -- Viewable 8-bit image: gray levels map to disparity range (64--255), cyan marks occluded pixels

### Parameter-Only Mode

If no output file is specified (no `dispMap.tif` and no `-o`), the program just prints the auto-computed `K` and `lambda` values without running the optimization:

```bash
bin/KZ2 ../images/scene_l.png ../images/scene_r.png -15 0
# Output: K=... lambda=...
```

---

## Input Requirements

- Images must be **epipolar rectified** (corresponding points lie on the same scanline)
- Both images should have the same height (cropped to the shorter if different)
- Supported formats: PNG (if libPNG available), PGM/PPM (always supported)
- The disparity range `[dMin, dMax]` is scene-dependent and must be specified manually

## Output Format

- **Float TIFF** (`.tif`): Precise disparity values, occluded pixels as NaN. Requires a specialized viewer.
- **Scaled PNG** (`.png` via `-o`): Human-viewable, gray levels proportional to disparity, occluded pixels in cyan (R=0, G=B=255).

---

## Technical Notes

- The code compiles under **C++98** (`-std=c++98`) and C89 (`-std=c89`) for maximum portability
- The graph uses `std::vector` for node/arc storage with index-based references, which is ~10-20% slower than the original raw-pointer implementation (Kolmogorov's "Match 3.4") due to vector indexing overhead vs direct pointer chasing
- The random alpha ordering may produce slightly different results across runs due to `srand(time(NULL))`
- Licensed under **GPLv3**
