# mtDNA_mt2_XL

A simulation study estimating mitochondrial DNA (mtDNA) variance components using structural equation modeling via OpenMx. Simulated pedigree data are generated across combinations of variance components and pedigree structures, and models are fit to recover those components.

## Dependencies

- R packages: `OpenMx`, `mvtnorm`

---

## Project Structure

```
mtDNA_mt2_XL/
├── InitialData/          # Pedigree structures, variance combinations, prep scripts
├── Functions/            # Core simulation and helper functions
├── Analysis/             # Power analysis plots and figures
├── Result_*/             # Simulation output (one folder per run configuration)
│   ├── cXpY/             # Results for variance combo X, pedigree Y
│   └── Analysis_*/       # Summary analyses for that result set
├── RunThis*.R            # Top-level scripts to launch simulations
├── compress_results.ps1  # Compress results for storage/sharing (maintainers)
└── decompress_results.R  # Extract results after cloning (users)
```

---

## Naming Convention

Result subfolders follow the pattern `cXpY`:

- **X** = row index in `InitialData/VarianceComb.csv` (the variance component combination)
- **Y** = row index in `InitialData/Pedigrees.csv` (the pedigree structure)

For example, `c3p7` contains results for variance combination 3 run on pedigree structure 7.

---

## Workflow

### 1. Prepare initial data

Before running simulations, generate the fixed pedigrees and load variance combinations:

```r
source("InitialData/InitialDataPrep.R")
```

This reads `Pedigrees.csv` and `VarianceComb.csv`, simulates pedigree structures, and saves everything to `InitialData/FixedPedVar.RData`.

### 2. Run simulations

Use the appropriate `RunThis*.R` script depending on the simulation configuration:

| Script | Description | Output |
|---|---|---|
| `RunThis.R` | Base simulation (n=10,000) | `Result2/` |
| `RunThis_r2.R` | 20k simulation | `Result_20k/` |
| `RunThis_a2mt2.R` | A2/Mt2 model only (n=20,000) | `ResultA2Mt2_20k_2/` |
| `RunThis_rdPed.R` | Random pedigree, full model | `Result_rdped_full/` |
| `RunThis_PedCheck.R` | Pedigree structure checks | `Result_Check/` |

Each script loops over variance combinations (rows of `VarianceComb.csv`) and pedigree structures (rows of `Pedigrees.csv`), writing one `cXpY/` subfolder of results per combination.

### 3. Analyse results

Analysis scripts are in the `Analysis/` folder and generate power plots grouped by chain/partition combinations (e.g., `PowPlotc1p1234.R`).

---

## Input Data

### `InitialData/Pedigrees.csv`

Defines pedigree structures used in simulations.

| Column | Description |
|---|---|
| `k` | Kids per couple |
| `G` | Number of generations |
| `p` | Sex ratio |
| `r` | Marriage rate |

### `InitialData/VarianceComb.csv`

Defines variance component combinations to simulate.

| Column | Description |
|---|---|
| `a2` | Additive genetic variance |
| `cN2` | Nuclear common environment |
| `cE2` | Extended common environment |
| `mt2` | Mitochondrial variance |
| `d2` | Dominance variance |
| `j2` | Additive × mitochondrial interaction |
| `e2` | Residual variance |

---

## File Compression

Result subfolders are stored as `.tar.gz` archives to keep the repository size manageable. `Analysis*` subfolders are excluded.

### After cloning — decompress results

```r
source("decompress_results.R")
```

Or from a terminal at the repo root:

```bash
Rscript decompress_results.R
```

This extracts all `.tar.gz` archives across every `Result*` directory.

### For maintainers — compress results

Run from the repo root in PowerShell:

```powershell
& ".\compress_results.ps1"
```

---

## Seeds

| Run | Seed |
|---|---|
| rdped_v1 | 11 |
| rdped_v2 | 12 |
| rdped_v3 | 13 |
| rdped_20k | 20 |
| rdped_full | 100 |
