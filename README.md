# Asymmetric Archimedean BNP reproducibility stack

This repository reproduces manuscript Sections 4.1--4.4. It provides exact
Khoudraji simulation, a Poisson--Dirichlet sampler, deterministic study jobs,
compact summaries, clustering, contour comparisons, and the tables and figures
used in the revised manuscript.

## Setup and checks

Run commands from the repository root:

```sh
Rscript scripts/install_dependencies.R
Rscript scripts/make_manifest.R
Rscript -e 'testthat::test_dir("tests/testthat")'
```

`salso` is optional. When available, the code reports SALSO partitions under
both Binder and variation-of-information (VI) loss. Otherwise it selects the
retained partition minimizing squared Binder loss to the posterior similarity
matrix (the Dahl least-squares fallback).

## Study manifest and deterministic data

`config/job_manifest.csv` is explicit and resumable. Models sharing a
scenario/truth/sample-size setting have the same `data_seed`, so they fit
exactly one deterministic data set. Each fit has its own `fit_seed`.

- Simulation 1: all four alpha truths, `n=200,500`, free and exchangeable
  Clayton/Frank/Gumbel BNP kernels.
- Simulation 2: free and exchangeable Clayton kernels on 20 seeded data sets
  (`mirror_rep01`--`mirror_rep20`, `n=500`) from the mirror-asymmetric mixture
  `0.5 CLA(50; 1, .65) + 0.5 CLA(50; .65, 1)`. The components share `theta`
  and Kendall tau and differ only in the direction of asymmetry, so no
  exchangeable kernel can separate them. For this experiment only,
  `theta ~ Gamma(4, .08)` and `b+a ~ Gamma(2,10)` replace the defaults.
  Each job stores recovery metrics against the generating labels (ARI of the
  Binder, Binder `K<=2`, VI and PSM-cut partitions, per-draw ARI, and
  within/between posterior similarity). Replicate 1 is the featured data set.
- Simulation 3: free and exchangeable Gumbel kernels.
- Simulation 4: single and two-component Gaussian truths at `n=500`, with
  free and exchangeable Clayton/Frank/Gumbel kernels.

Smoke mode uses 24 iterations (8 burn-in, thin 4). Full mode uses exactly
15,000 iterations (5,000 burn-in, thin 5). To smoke one representative job
from every study:

```sh
Rscript scripts/run_manifest.R --mode smoke --subset studies --workers 4
```

To run the complete production study:

```sh
Rscript scripts/run_manifest.R --mode full --subset all --workers 4
Rscript scripts/aggregate_results.R full
Rscript scripts/manuscript_results.R
```

Each worker is an independent `Rscript` process. This avoids inheriting a
`sourceCpp` DLL through a fork. Existing completed summary RDS files are
skipped; use `--force` to rerun.

## Outputs

Outputs are under `results/{smoke,full}/`. The repository includes the
completed compact `full` results; full fitted chains are deliberately excluded.

- `raw/`: full fitted chains (ignored by git);
- `summaries/`: compact RDS and CSV files for LPML, alpha/theta, pairwise
  posterior-predictive tau, and clusters;
- `diagnostics/`: traces and acceptance information;
- `figures/`: observed/predictive contour PDFs for bivariate jobs;
- `latex/`: per-job and aggregate table fragments.
- `manuscript/`: the four final tables and comparison figures used in Sections
  4.1--4.4.

The retained allocation vector `z` is stored at every saved iteration. The
posterior similarity matrix and selected partition are included in each
compact RDS.

## Statistical implementation notes

The exact product construction is

`Z_j = max(X_j^(1/alpha_j), Y_j^(1/(1-alpha_j)))`,

where `X` is sampled from the inner Archimedean copula using `copula` and the
`Y_j` are independent uniforms. Alpha zero returns `Y_j` exactly and alpha one
returns `X_j` exactly.

Parameter support is dimension-aware. The fitted Clayton kernel is restricted
to positive theta in every dimension because this repository uses the
absolutely continuous strict-generator density; negative bivariate Clayton
has a singular boundary component. Frank must be positive above dimension
two but permits either sign in dimension two; Gumbel and Joe require theta at
least one. Theta and alpha have separate MH decisions and acceptance diagnostics.
The theta random walk includes the truncated-uniform Hastings correction.
Batch adaptation is diminishing and stops after burn-in. The transformed
update for `b+a` enforces `b > -a`.

LPML uses a full leave-one-out Poisson--Dirichlet predictive density at each
retained state:

`[sum_k (n_{-i,k}-a)c_k(u_i) + (b+a*m_{-i}) m0(u_i)] / (b+n-1)`.

The base density `m0` is Monte Carlo integrated under `G0`; the number of
draws is recorded in `cpo_method`. CPO uses the harmonic mean of those
full-mixture ordinates. This is more expensive but avoids evaluating only the
atom occupied by observation `i`.

Posterior-predictive tau is computed from draws from the complete predictive
mixture. Consequently it includes the cross-component terms (`k != j`) in the
mixture Kendall-tau formula. `fit_single_khoud()` supplies the optional
single-component maximum-likelihood comparison.

## Main simulation findings

- The free model recovers substantial bivariate asymmetry at `n=500`; under
  exchangeability its alpha coordinates concentrate near one.
- In the revised two-component asymmetric Clayton experiment, free alpha gains
  76.36 LPML units and the two-group Binder summary attains ARI 0.672, close to
  the known-parameter oracle benchmark 0.712. The unconstrained VI summary is
  more conservative and merges the free fit into one group.
- In dimension three, free alpha recovers `(1, 0.7, 0.2)`, reproduces all three
  unequal pairwise taus, and improves LPML by 470.35.
- For Gaussian truths, nested Frank is best. The Khoudraji layer cannot replace
  a generator capable of representing both positive and negative dependence.

## Scope and caveats

The compiled full likelihood is implemented for dimensions two and three,
which covers the four studies. Numerical tau integration can be slower near
parameter boundaries; Monte Carlo tau is available as a validated
alternative. The base-measure part of CPO and posterior-predictive tau are
Monte Carlo approximations with deterministic seeds. Clayton is restricted to
positive theta because its negative bivariate form contains a singular boundary
component that is outside the absolutely continuous likelihood implemented
here.

## Legacy Simulation 2 (manuscript v2)

The earlier single-data-set Simulation 2 (`CLA(30; 1, .95)` and
`CLA(1; .05, 1)`) is archived, with its summaries and figures, under
`results/full/legacy_v2_sim2/` (fitted chains excluded). It is no longer in
the job manifest and no longer feeds the manuscript tables.

## Real-data refits (occupancy and red wine)

```sh
Rscript scripts/run_real_data.R --mode full --cores 9 --seeds 1,2,3
Rscript scripts/real_data_tables.R
```

Uses the same sampler as Simulations 1--4. Data are in `data/real/`
(Candanedo and Feldheim occupancy; Cortez et al. red wine). Exchangeable
(`alpha = 1`) fits are run for three seeds; Khoudraji fits for seed 1.
`results/real/full/real_data_summary.csv` holds the LPML summaries.

## Sampler conventions

- **Finite adaptation.** The MH tuning constants are updated in batches of 50
  iterations only while `it <= burn` (step `min(.05, j^-0.6)`, target
  acceptance .35); they are frozen afterwards. `fit$tuning_trace` and
  `fit$adaptation` record this, and `tests/testthat` checks it.
- **LPML.** `fit$lpml` is the harmonic mean of leave-one-out PD predictive
  densities (used in Simulations 1--4); `fit$lpml_comp` is the
  component-specific harmonic mean of Pan et al. (2025) (used for the
  real-data tables). `lpml_component()` recomputes the latter from stored
  draws; `results/full/summaries/lpml_both_estimators.csv` compares both.
- Plot labels use "Exchangeable" and "Khoudraji"; `scripts/regenerate_job_figures.R`
  rebuilds the per-job PDFs from stored raw fits.
