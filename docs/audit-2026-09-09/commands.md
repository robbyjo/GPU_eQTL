# Audit commands

Run from the repository root. Candidate JARs and optional runtime libraries
are deliberately kept under ignored `target/`; checksums are in the report.
Maven may need access to its normal local dependency cache and Maven Central.

```powershell
# Independent references (ordinary unit tests use checked-in fixtures).
New-Item -ItemType Directory -Force target\audit-r-lib
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' -e 'install.packages("CompQuadForm",lib="target/audit-r-lib",repos="https://cloud.r-project.org")'
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' scripts\generate-audit-reference.R

# Optional, when compatible CUDA 12 libraries are missing from the host PATH.
python -m pip install --target target\cuda12 --no-deps nvidia-cublas-cu12==12.9.2.10 nvidia-cuda-runtime-cu12==12.9.79 nvidia-cuda-nvrtc-cu12==12.9.86
$env:PATH=(Resolve-Path target\cuda12\nvidia\cublas\bin).Path+';'+(Resolve-Path target\cuda12\nvidia\cuda_runtime\bin).Path+';'+(Resolve-Path target\cuda12\nvidia\cuda_nvrtc\bin).Path+';'+$env:PATH

.\mvnw.cmd -q package
.\mvnw.cmd -q -Pcpu-mkl-windows-x86_64 '-Dtest=CpuBackendTest,GpuKernelIntegrationTest' test

# CPU module and supplied-covariance benchmarks.
java -cp 'target\test-classes;target\classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar' gov.nih.eqtl.QModuleAuditBenchmark
java -cp 'target\test-classes;target\classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar' gov.nih.eqtl.QCohortAuditBenchmark

# Identical harness on baseline classes saved before this audit's edits.
java -cp 'target\test-classes;target\baseline-classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar' gov.nih.eqtl.QModuleAuditBenchmark
java -cp 'target\test-classes;target\baseline-classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar' gov.nih.eqtl.QCohortAuditBenchmark

# Current transfer-inclusive GPU/CPU association products.
java --enable-native-access=ALL-UNNAMED -cp 'target\test-classes;target\classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar' gov.nih.gpu.GpuBackendBenchmark

# Standalone comparison against the checksum-verified published candidate JARs.
javac -cp 'target\classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar;target\dependency-audit\JLinAlg-0.3.0-library.jar;target\dependency-audit\jdistlib-all-0.10.2.jar' -d target\dependency-audit scripts\QDependencyAudit.java
java -cp 'target\dependency-audit;target\classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar;target\dependency-audit\JLinAlg-0.3.0-library.jar;target\dependency-audit\jdistlib-all-0.10.2.jar' org.jlinalg.settest.QDependencyAudit
javac -cp target\dependency-audit\jdistlib-all-0.10.2.jar -d target\dependency-audit scripts\QDependencyGpuAudit.java
java --enable-native-access=ALL-UNNAMED -cp 'target\dependency-audit;target\dependency-audit\jdistlib-all-0.10.2.jar' QDependencyGpuAudit

# Packaged smoke; --pedigree must fail with an explicit unsupported-model error.
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo
java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --pedigree unimplemented.tsv
```

In a fresh shell without the optional CUDA PATH additions, packaged automatic
discovery also completed successfully on this host using its NVIDIA OpenCL
fallback. Baseline and optimized CPU benchmarks ran sequentially, not
concurrently with other audit workloads.
