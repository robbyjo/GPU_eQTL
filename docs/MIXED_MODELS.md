# Fitted Gaussian mixed models and pedigree analyses

GPU eQTL now supports two fitted mixed-model modes. Both estimate variance
components by REML using the pinned JLinAlg 0.3.0 library and JDistlib 0.10.2.
The supplied-covariance `--cohort-model` GLS path remains a separate option.

| Mode | Fitting and test | Intended use |
| --- | --- | --- |
| `--mixed-model score` (alias `genesis`) | Fit one REML null per trait; reuse its covariance for all predictor score tests; chi-square(1) inference | Association scans |
| `--mixed-model reml` | Fit the null, then refit REML with each predictor in the fixed-effect design; Wald inference | Slower per-association fitting and sensitivity checks |

Supplying `--pedigree` or `--random-covariates` defaults to `score`. Without
either, the default remains the existing fixed-effect analysis. Explicit
`--mixed-model none` rejects random/pedigree terms.

## GENESIS relationship

GENESIS's Gaussian `fitNullModel` estimates variance components by average-information
REML. Its usual score scan reuses the fitted null covariance. This implementation
uses that statistical strategy, with JLinAlg doing the REML optimization. It does
not invoke R, reproduce every GENESIS option, or use GENESIS's separate
`fast.Score.SE` diagonal/variance-ratio approximation. See the
[GENESIS null-model documentation](https://github.com/UW-GAC/GENESIS/blob/devel/man/fitNullModel.Rd).

The checked-out baseline had no fitted REML or pedigree execution path:
random terms were ignored and pedigree input was not consumed by model fitting.
The new modes implement those features explicitly.

For each trait, with fixed design X and fitted covariance V, the score path forms
`P = V^-1 - V^-1 X (X' V^-1 X)^-1 X' V^-1`.
For predictor g, `U = g' P y`, `I = g' P g`, and `Q = U²/I`.
The reported effect and standard error are `U/I` and `1/sqrt(I)`; the p-value
uses chi-square(1). Variance components are estimated under the null, so this
is asymptotic inference. It is not the fixed-effect t test.

The refit path estimates V with `[X,g]` in the mean model. The default
`--reml-df residual` reports `N - rank([X,g])` degrees of freedom. That denominator
is an approximation even though variance estimation is REML.
`--reml-df satterthwaite` and `--reml-df kenward-roger` select JLinAlg's
coefficient-specific inference; the latter also adjusts the inference covariance.
`--reml-max-iterations` defaults to 200 per fit. Nonconvergence stops execution.

## Inputs and commands

Pedigrees are headered, tab-separated files with exactly these columns:

```text
ID	SIRE	DAM
father	0	0
mother	0	0
child	father	mother
```

`0` and `.` mean unknown parent; blank values are errors. Parents must have
their own rows. Input order may place children before parents. Duplicate IDs,
cycles, self-parenting, missing parents, and absent analysis individuals are
fatal. Unknown parents are assumed unrelated and noninbred. JLinAlg constructs
the numerator relationship matrix A, equal to twice the kinship matrix.
The fitted covariance includes `sigma_A² A + sigma_e² I`.

By default, pedigree IDs match aligned genotype sample IDs after any explicitly
requested literal prefix removal. For repeated observations, use
`--pedigree-id-column individual` to map each unique observation ID through the
covariate table to a pedigree individual. Multiple observations can share that
individual. Matrix sample IDs must still be unique.

`--random-covariates family,batch` names random-intercept **grouping columns**
in the covariate file. Numeric-looking labels are group labels too; this option
does not implement numeric random slopes. Each grouping contributes a fitted
covariance basis with entries one for observations in the same group and zero
otherwise. Grouping terms may accompany a pedigree. Redundant or confounded
variance components fail the fitter's identifiability checks.

```sh
java -Xmx8g -jar target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar \
  --genotype genotypes.csv --expression traits.csv --output score.csv \
  --covariates covariates.csv --fixed-covariates Age,Sex \
  --pedigree pedigree.tsv --mixed-model score --backend auto

java -Xmx8g -jar target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar \
  --genotype genotypes.csv --expression traits.csv --output reml.csv \
  --covariates covariates.csv --fixed-covariates Age,Sex \
  --pedigree pedigree.tsv --mixed-model reml --reml-df satterthwaite
```

Legacy INI equivalents are `mixed_model`, `pedigree_file`, `pedigree_id_column`,
`covariate_random`, `reml_df`, and `reml_max_iterations`. INI covariate lists are
whitespace-separated. Mixed analyses require modern headered CSV traits and
headered CSV or VCF/BCF predictors; legacy TPED/headerless inputs are rejected.

Existing alignment, categorical fixed-effect encoding, aligned VCF QC and
explicit predictor filling policies run before fitting. Missing grouping or
individual-mapping entries participate in the selected covariate-missing policy.
No missingness policy is inferred from pedigree relationships.

## Output, restart, and resource use

Mixed output uses a distinct CSV schema:

`SNP,Probe,N,DF,Effect,SE,Statistic,log10P,Test,Status,VarianceComponents,REMLLogLik,Iterations,Singular`

Rows follow trait source order, then predictor source order. `DF` is one for
score chi-square tests or the Wald denominator DF for refits. `Test` identifies
the distinction. Unestimable predictors are reported with explicit status,
NaN effect/SE and p=1; a constant or fixed-effect-explained trait is fatal.
The variance column lists named estimates for the relevant null or full fit.
`Singular` marks a component at or below 1e-6 of the largest fitted variance.
The optimizer variance bounds are 1e-10 and 1e10 times the trait's initial
fixed-effect residual variance, so changing response units does not impose an
absolute variance floor. Near-zero estimates are not
evidence that a random effect is well identified.

`--threshold pval P` filters the actual mixed-model p-value, with `0 < P <= 1`.
R-squared thresholds/output, legacy rounding, and `--df-offset` are rejected.

One atomic checkpoint part contains one complete trait scan. Resume is bound to
input SHA-256 hashes, configuration, sample permutations, dependency versions,
and compute-device/runtime identity. Source identity is checked before committing
parts and hashes are checked again before final assembly. An interrupted trait
is recomputed. `--checkpoint-dir`, `--resume`, and `--keep-checkpoints` apply;
the default directory is `OUTPUT.mixed-checkpoint`.

One trait's dense null state is resident at a time. The fitter needs several
N-by-N workspaces and refuses an estimated allocation above 65% of the JVM heap.
`--genotype-block-rows` bounds predictor tiles (default 256 here). The score
path sends P-times-predictor tiles through the existing backend interfaces,
with FP64, rounded active tiles and exclusive pooled contexts. It currently
submits those tiles sequentially; multiple devices do not imply concurrent
multi-GPU speedup. REML fitting and per-predictor refitting use JLinAlg's
portable CPU backend. Refitting skips GPU initialization.

This is a dense implementation, not a sparse pedigree solver. It rereads the
predictor source for each trait. Very large samples or many traits can therefore
be expensive despite null reuse. General `--block-size`, `--threads`,
`--expression-block-rows`, `--trait-cache` and prepared fixed-effect caches do
not alter this one-trait schedule. Full private-cohort and multi-GPU scaling
have not been established.

Current exclusions: fitted mixed-model burden/SKAT/SKAT-O, exact trait-pattern
deletion, supplied cohort GLS combined with fitted random effects, GLMM/binary
outcomes, random slopes, and FP32. For traits with missing entries, choose an
explicit filling/exclusion policy or a complete input. A default `pattern`
policy is accepted when there are no missing trait values.

## Dependency packaging and verification

Upstream releases are not available from Maven Central. During `initialize`,
the normal Maven wrapper build downloads pinned release artifacts to ignored
`.dependencies/`, verifies SHA-256 on every build, and incorporates only
`org/jlinalg/**` and `jdistlib/**` Java entries into the application. It does not
register upstream GPU/native service providers or bundle their vendor runtimes.
The existing CUDA/OpenCL/CPU backends continue to own association execution.
All jars retain the upstream GPL notices under `META-INF/licenses` alongside
the existing application and oneMKL notices. See `pom.xml` for hashes and URLs.

The legacy `net.sourceforge.jdistlib` namespace remains for established
fixed-effect paths; modern `jdistlib` supplies the new inference and tail code.
CPU QR/rank policy and the existing GPU GEMM engine have not been replaced.
After changing dependencies, run `./mvnw test` (Windows: `.\mvnw.cmd test`)
before refreshing an Eclipse/m2e project. First build requires network access;
subsequent clean builds reuse the verified `.dependencies` artifacts.

`scripts/generate-mixed-reference.R` independently fits Gaussian REML using
base R Cholesky/optimization and checks the random-intercept refit against
`nlme`. Committed references cover grouping, pedigree and combined covariance,
both null score and full refit. Hardware tests compare the actual CUDA, OpenCL
and CPU projection paths to the scalar calculation. These comparisons are not
a claim of byte-for-byte parity with every GENESIS or R small-sample DF mode.
