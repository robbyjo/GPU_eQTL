#!/usr/bin/env sh
set -eu

platform="${1:-universal}"
project_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$project_root"
git_commit=$(git rev-parse HEAD)
output_timestamp=$(git show -s --format=%cI HEAD)
dirty=false
if test -n "$(git status --porcelain --untracked-files=normal)"; then
  dirty=true
  if test "${ALLOW_DIRTY:-0}" != 1; then
    echo "Release verification requires a clean worktree; commit/stash changes or set ALLOW_DIRTY=1 for development-only validation" >&2
    exit 2
  fi
fi

case "$platform" in
  universal)
    profile=cpu-openblas-universal
    build_directory=target
    classifier=all
    platform_argument=
    ;;
  windows-x86_64|linux-x86_64|linux-arm64|macosx-x86_64|macosx-arm64)
    profile=cpu-openblas-platform
    build_directory="target-openblas-$platform"
    classifier="openblas-$platform-all"
    platform_argument="-Dopenblas.platform=$platform"
    ;;
  *)
    echo "Unsupported platform: $platform" >&2
    exit 2
    ;;
esac

./mvnw "-P$profile" "-Dproject.build.outputTimestamp=$output_timestamp" \
  ${platform_argument:+"$platform_argument"} clean package
source_jar="$build_directory/gpu-eqtl-2.0.0-SNAPSHOT-$classifier.jar"
test -f "$source_jar"
release_directory="target/release/$platform"
mkdir -p "$release_directory"
release_jar="$release_directory/gpu-eqtl-2.0.0-SNAPSHOT-$platform.jar"
cp "$source_jar" "$release_jar"

entries=$(jar tf "$release_jar")
printf '%s\n' "$entries" | grep -Fx 'META-INF/LICENSE' >/dev/null
printf '%s\n' "$entries" | grep -Fx 'META-INF/LICENSE_EXCEPTION' >/dev/null
printf '%s\n' "$entries" | grep -E '^META-INF/licenses/.+' >/dev/null
java -jar "$release_jar" --printbackendinfo

(
  cd "$release_directory"
  sha256sum "$(basename "$release_jar")" > SHA256SUMS
)
java_version=$(java -version 2>&1 | sed -n '1p')
hash=$(awk '{print $1}' "$release_directory/SHA256SUMS")
{
  echo "artifact=$(basename "$release_jar")"
  echo "sha256=$hash"
  echo "git_commit=$git_commit"
  echo "maven_profile=$profile"
  echo "platform=$platform"
  echo "build_os=$(uname -srm)"
  echo "java=$java_version"
  echo "output_timestamp=$output_timestamp"
  echo "dirty_worktree=$dirty"
  echo "backend_smoke=true"
  echo "tests=true"
} > "$release_directory/RELEASE-MANIFEST.txt"
echo "Verified release artifact: $release_jar"
echo "SHA-256: $hash"
