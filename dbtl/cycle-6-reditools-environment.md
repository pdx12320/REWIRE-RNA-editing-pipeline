# Cycle 6 — Make REDItools2 and MPI reproducible

## Design

### Question

Could REDItools2 be executed consistently on all six libraries despite its Python 2 and MPI dependencies?

REDItools2 was the central substitution-calling tool, but its runtime environment was less straightforward than the Python 3 wrappers used elsewhere in the project.

The environment had to satisfy three conditions:

1. Python 2.7-compatible dependencies were installed;
2. `mpi4py` was compiled against the same MPI implementation used at runtime;
3. the number of requested workers fit the available allocation.

## Build

A dedicated environment definition was retained in:

```text
pipeline/env/reditools2_py2.yml
```

Core setup:

```bash
conda env create -f pipeline/env/reditools2_py2.yml
conda activate reditools2_py2
```

When the old solver could not install every package, compatible versions were installed explicitly:

```bash
curl -L https://bootstrap.pypa.io/pip/2.7/get-pip.py -o /tmp/get-pip-py2.py
python /tmp/get-pip-py2.py
python -m pip install "setuptools<45" "wheel<0.35"
python -m pip install \
  "pysam==0.15.4" \
  "sortedcontainers==2.2.2" \
  "psutil==5.6.7" \
  "netifaces==0.10.9"

export MPICC="$(command -v mpicc)"
python -m pip install --no-cache-dir --no-binary=mpi4py "mpi4py==3.0.3"
```

The REDItools2 launcher receives the environment interpreter explicitly:

```bash
"$CONDA_PREFIX/bin/python"
```

rather than relying on the generic `python` command.

## Test

### Dependency failures

Observed errors showed that the Python 2 environment lacked one or more required modules:

```text
mpi4py
pysam
sortedcontainers
psutil
netifaces
```

The issue was not in the RNA data or command-line parameters. It was an incomplete runtime environment.

### MPI allocation failure

A separate test produced an insufficient-slots error when the requested process count exceeded the scheduler allocation.

Allocation inspection:

```bash
mpirun --display-allocation -np 1 hostname
```

The final run used a process count that left capacity for the controller and operating system rather than occupying every advertised core.

### Multi-layer execution

The final launcher separated three forms of parallelism:

```text
MPI REDItools2 processes
concurrent coverage jobs
compression threads
```

Example:

```bash
nohup bash pipeline/scripts/rna/run_reditools_all_samples.sh \
  "$PROJECT" "$REF" "$REDITOOLS" "$CONDA_PREFIX/bin/python" \
  30 8 8 \
  > "$PROJECT/logs/reditools.log" 2>&1 &
```

The values `30 8 8` correspond to MPI processes, concurrent coverage jobs and compression threads in the frozen command pattern.

## Learn

### Lesson 1 — Interpreter identity must be explicit

A successful package installation is not enough when the workflow later invokes another interpreter from `PATH`. Passing the full environment interpreter removed this ambiguity.

### Lesson 2 — MPI must be consistent end to end

The following commands must refer to the same MPI implementation:

```text
mpirun
mpicc
mpi4py build
```

Mixing implementations can produce import or runtime failures even when each component is installed.

### Lesson 3 — Parallelism is a resource contract

More workers do not automatically make the run safer or faster. The revised workflow records the requested process count and checks it against the actual allocation.

### Lesson 4 — Environment setup belongs in reproducibility evidence

Because dependency versions determined whether the core caller could run, the environment file and fallback installation commands were retained in the repository rather than treated as temporary server notes.

## Final role in the pipeline

This cycle converted REDItools2 from a manually repaired tool into a repeatable pipeline component:

```text
dedicated Python 2 environment
→ explicit interpreter
→ consistent MPI build
→ allocation-aware worker count
→ per-sample logged execution
```
