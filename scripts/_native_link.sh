# Sourced by scripts/run_pybind11_tests.sh: the checks around the compiled module it
# links into traminex/ for the length of a run (docs/DEV_BUILD.md).

# The runner builds in a fresh `mktemp -d "${NATIVE_BUILD_PREFIX}XXXXXX"` directory and
# writes its process ID there under NATIVE_PID_FILE.
NATIVE_BUILD_PREFIX=/tmp/traminex-r1-pybind11.
NATIVE_PID_FILE=runner.pid

# Whether the runner build directory holding the link target $1 belongs to a run that is
# over: the directory sits directly under the prefix, its NATIVE_PID_FILE names a process
# that no longer exists, and no process has the target mapped (a test process the killed
# runner left behind would). A directory without the file (an older runner) is not over.
native_run_is_over() {
  local target=$1 build=${1%/*} pid
  [[ "$build" == "$NATIVE_BUILD_PREFIX"* && "${build#"$NATIVE_BUILD_PREFIX"}" != */* ]] || return 1
  [[ -f "$build/$NATIVE_PID_FILE" ]] || return 1
  pid=$(<"$build/$NATIVE_PID_FILE")
  [[ "$pid" =~ ^[0-9]+$ && ! -e "/proc/$pid" ]] || return 1
  ! grep -qsF -- "$target" /proc/[0-9]*/maps
}

# Remove each traminex/_pb_kernels*.so link a killed run left behind (TX-AUD-20261010-088):
# a run killed before its EXIT trap keeps its link, live while its build directory exists
# and dead once /tmp is cleared, and the next run would refuse it. A link into one of the
# runner's build directories is removed when its target is gone, or when its run is over
# (above), together with that build directory. A link of a run in progress, a regular
# file, or a link pointing anywhere else is left alone.
clear_dead_native_links() {
  local package=$1 module target
  for module in "$package"/_pb_kernels*.so; do
    [[ -L "$module" ]] || continue
    target=$(readlink -- "$module")
    [[ "$target" == "$NATIVE_BUILD_PREFIX"*/* ]] || continue
    if [[ ! -e "$module" ]]; then
      echo "Removing a dead link left by an interrupted run: $module -> $target" >&2
      rm -f -- "$module"
    elif native_run_is_over "$target"; then
      echo "Removing the link and build directory of a killed run: $module -> $target" >&2
      rm -f -- "$module"
      rm -rf -- "${target%/*}"
    fi
  done
}

# Refuse to replace an existing module or link at the path the runner is about to link.
refuse_existing_module() {
  local candidate=$1
  if [[ -e "$candidate" || -L "$candidate" ]]; then
    echo "Refusing to replace existing module: $candidate" >&2
    return 1
  fi
}
