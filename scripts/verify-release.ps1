param(
    [ValidateSet("universal", "windows-x86_64", "linux-x86_64", "linux-arm64", "macosx-x86_64", "macosx-arm64")]
    [string]$Platform = "universal",
    [switch]$SkipTests,
    [switch]$SkipBackendSmoke,
    [switch]$AllowDirty
)

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $projectRoot
$gitCommit = (git rev-parse HEAD).Trim()
$outputTimestamp = (git show -s --format=%cI HEAD).Trim()
$dirty = -not [string]::IsNullOrWhiteSpace((git status --porcelain --untracked-files=normal | Out-String))
if ($dirty -and -not $AllowDirty) {
    throw "Release verification requires a clean worktree; commit/stash changes or use -AllowDirty for development-only validation"
}

if ($Platform -eq "universal") {
    $profile = "cpu-openblas-universal"
    $buildDirectory = Join-Path $projectRoot "target"
    $classifier = "all"
    $platformArgument = @()
} else {
    $profile = "cpu-openblas-platform"
    $buildDirectory = Join-Path $projectRoot ("target-openblas-" + $Platform)
    $classifier = "openblas-$Platform-all"
    $platformArgument = @("-Dopenblas.platform=$Platform")
}

$mavenArguments = @("-P$profile", "-Dproject.build.outputTimestamp=$outputTimestamp")
$mavenArguments += $platformArgument
if ($SkipTests) { $mavenArguments += "-DskipTests" }
$mavenArguments += @("clean", "package")
& .\mvnw.cmd @mavenArguments
if ($LASTEXITCODE -ne 0) { throw "Maven release build failed" }

$sourceJar = Join-Path $buildDirectory "gpu-eqtl-2.0.0-SNAPSHOT-$classifier.jar"
if (-not (Test-Path -LiteralPath $sourceJar -PathType Leaf)) {
    throw "Expected shaded jar was not produced: $sourceJar"
}
$releaseDirectory = Join-Path $projectRoot ("target\release\" + $Platform)
New-Item -ItemType Directory -Force -Path $releaseDirectory | Out-Null
$releaseJar = Join-Path $releaseDirectory ("gpu-eqtl-2.0.0-SNAPSHOT-$Platform.jar")
Copy-Item -LiteralPath $sourceJar -Destination $releaseJar -Force

$entries = & jar tf $releaseJar
if ($LASTEXITCODE -ne 0) { throw "The packaged jar cannot be listed" }
$requiredEntries = @("META-INF/LICENSE", "META-INF/LICENSE_EXCEPTION")
foreach ($entry in $requiredEntries) {
    if ($entries -notcontains $entry) { throw "Packaged jar omits $entry" }
}
if (-not ($entries | Where-Object { $_ -like "META-INF/licenses/*" })) {
    throw "Packaged jar omits third-party license notices"
}

if (-not $SkipBackendSmoke) {
    & java -jar $releaseJar --printbackendinfo
    if ($LASTEXITCODE -ne 0) { throw "Packaged backend-information smoke test failed" }
}

$hash = (Get-FileHash -LiteralPath $releaseJar -Algorithm SHA256).Hash.ToLowerInvariant()
$checksumFile = Join-Path $releaseDirectory "SHA256SUMS"
Set-Content -LiteralPath $checksumFile -Encoding ascii -NoNewline -Value "$hash  $(Split-Path -Leaf $releaseJar)`n"
$javaVersion = (& java -version 2>&1 | Select-Object -First 1).ToString()
$manifest = @(
    "artifact=$(Split-Path -Leaf $releaseJar)",
    "sha256=$hash",
    "git_commit=$gitCommit",
    "maven_profile=$profile",
    "platform=$Platform",
    "build_os=$([System.Runtime.InteropServices.RuntimeInformation]::OSDescription)",
    "java=$javaVersion",
    "output_timestamp=$outputTimestamp",
    "dirty_worktree=$dirty",
    "backend_smoke=$(-not $SkipBackendSmoke)",
    "tests=$(-not $SkipTests)"
)
Set-Content -LiteralPath (Join-Path $releaseDirectory "RELEASE-MANIFEST.txt") -Encoding utf8 -Value $manifest
Write-Host "Verified release artifact: $releaseJar"
Write-Host "SHA-256: $hash"
