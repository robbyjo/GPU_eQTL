# CPU eigenvalues and overlapping set-test tails

Set tests keep the selected fixed-effect or pedigree model. In particular,
`--set-null-fitter spectral` still fits a separate two-component REML null for
each response and exact observed mask. It does not replace the eigendecomposition
of each window's full score covariance or its weighted chi-square tail.

For SKAT, the production runner uses CPU LAPACK `dsyevd` in values-only mode.
SKAT-O requests eigenvectors as well. The normalized symmetric covariance,
positive-semidefinite check and numerical rank tolerance are unchanged.
Bundled OpenBLAS provides LAPACK, including when association products use CUDA,
OpenCL or oneMKL. `-Deqtl.cpu.blas=java` selects the portable reference. No GPU
eigensolver is claimed. The existing exclusive spectral algebra scope owns and
restores the full OpenBLAS thread setting; matrix products and eigen calls stay
on the producer thread, without nested native workers.

The positive gamma series retains its relative remainder check and 131,072-term
hard ceiling. A preferred 2,048-term budget now hands expensive spectra to the
existing tilted characteristic inversion earlier. Both compute the same
weighted chi-square law. Rank-one, equal-eigenvalue and rank-two special cases,
the inversion's relative-error and evaluation limits, and explicit numerical
failures are preserved. This is a choice between checked numerical methods;
moment matching, eigenvalue truncation and relaxed tail accuracy are not used.
If early inversion fails, the series can still use its original full budget;
the coefficient cache grows only for that recovery. Failure of both methods
remains fatal, with both diagnostics preserved. Cancellation never triggers a
fallback calculation.

`--threads` controls bounded CPU result workers in the production set runner.
Zero/omitted selects at most eight; the count is capped by available processors,
sets per block, 256 workers and one worker per 32 MiB of JVM heap. That last cap
is an admission limit, not a full workload memory estimate. One producer retains
the genotype rows and dense covariance for the current window. Only compact
scores, eigenvalues, optional burden projections and audit text enter the CPU
queue, with at most one pending result per worker. While those workers compute
p-values, the producer can prepare the next window's CPU/GPU products. Output
is consumed in source order and committed only after all tails and source
checks in the block succeed. Failed tasks cancel pending work; contexts close
after workers exit. Device contexts remain exclusive. The existing product
tile schedule is unchanged; this does not introduce concurrent multi-GPU tiles.

Budget the **sum** of the producer's native BLAS team and CPU result workers,
across all concurrent JVMs. For example, a process allocated eight cores can use
`-Deqtl.cpu.threads=4` before `-jar`, and `--threads 4` among its analysis options.
More workers or a larger heap does not guarantee better throughput. Leave heap
headroom for input tiles, pedigree geometry and native workspaces. Benchmarks
must separate component timing from complete analysis throughput.

Use fresh result/checkpoint paths after this change: the set runner and native
eigen identities are versioned. Compatible raw dosage caches and fitted-null
identities are unchanged. Changing only the CPU result worker count does not
change the result signature; compatible resumes remain byte-identical.
Existing running jobs and installed releases are unaffected until a later,
separately chosen deployment.

The synthetic component benchmark is
[`SetNumericsBenchmark.java`](../scripts/benchmarks/SetNumericsBenchmark.java).
From a checkout containing this change, compile it against the candidate, then
run the same class against the previous release and the candidate in separate
JVMs. This PowerShell example uses three replicates and reports median times:

```powershell
$candidate = 'target/gpu-eqtl-1.9.2.8-all.jar'
$baseline = 'previous/gpu-eqtl-1.9.2.8-all.jar'
javac --release 17 -cp $candidate -d .dependencies/bench-classes scripts/benchmarks/SetNumericsBenchmark.java
$classes = '.dependencies/bench-classes'
java -Xmx2g -cp "$classes;$baseline" gov.nih.eqtl.settest.SetNumericsBenchmark eigen-java 3
java -Xmx2g -Deqtl.cpu.blas=openblas -Deqtl.cpu.threads=4 -cp "$classes;$candidate" gov.nih.eqtl.settest.SetNumericsBenchmark eigen-native 3
java -Xmx2g -cp "$classes;$baseline" gov.nih.eqtl.settest.SetNumericsBenchmark tail 3
java -Xmx2g -cp "$classes;$candidate" gov.nih.eqtl.settest.SetNumericsBenchmark tail 3
```

Place the unchanged previous release at the baseline path and build the candidate
with the Maven wrapper first. Use `:` as the classpath separator on POSIX.
`eigen-vectors` benchmarks native full eigenvectors for SKAT-O. Record JAR hashes,
JDK, hardware, thread counts and complete output, and run timing comparisons
without competing tests or analyses on the benchmark host. These fixtures
measure individual numerical stages, not full-cohort completion time.
