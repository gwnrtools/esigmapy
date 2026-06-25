#!/usr/bin/env python3
"""Compare surrogate behavior and line-level profiling across two git commits.

This harness is intentionally standalone:

- It resolves two commits into detached git worktrees.
- Each worker subprocess runs with an explicit PYTHONPATH that points at the
  chosen worktree root.
- Each worker starts in a temporary cwd outside any worktree so cwd-based import
  shadowing cannot mask a PYTHONPATH bug.
- Each worker prints and asserts the imported esigmapy path before any surrogate
  work begins.

Outputs:

- one line_profiler text report per commit / parameter-set pair
- one raw JSON artifact containing worker measurements and comparison metrics
- one flattened CSV artifact for easy diffing in later optimization passes
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


DEFAULT_PARAM_SETS = [
    {
        "M": 10.0,
        "q": 2.3,
        "reference_eccentricity": 0.43,
        "reference_mean_anomaly": 1.3,
        "delta_t": 0.000244140625,
        "t_start": -136.0,
    },
    {
        "M": 10.0,
        "q": 2.3,
        "reference_eccentricity": 0.43,
        "reference_mean_anomaly": 1.3,
        "delta_t": 0.000244140625,
        "t_start": -1.0,
    },
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _default_config_path() -> Path:
    return Path(__file__).resolve().with_name("surrogate_compare_params.json")


def _default_output_dir(ref_commit: str, head_commit: str) -> Path:
    stamp = time.strftime("%Y%m%d-%H%M%S", time.gmtime())
    safe_ref = ref_commit.replace("/", "_")
    safe_head = head_commit.replace("/", "_")
    return _repo_root() / "tests" / "perf" / "results" / f"esigmapy-surrogate-compare-{safe_ref}-vs-{safe_head}-{stamp}"


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if dataclasses.is_dataclass(value):
        return dataclasses.asdict(value)
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    raise TypeError(f"Object of type {type(value)!r} is not JSON serializable")


def _run(cmd: list[str], *, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd is not None else None,
        text=True,
        capture_output=True,
        check=True,
    )


def _run_git(repo_root: Path, *args: str) -> str:
    result = _run(["git", "-C", str(repo_root), *args], cwd=repo_root)
    return result.stdout.strip()


def _ensure_param_config(config_path: Path) -> list[dict[str, Any]]:
    if config_path.exists():
        with config_path.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    with config_path.open("w", encoding="utf-8") as fh:
        json.dump(DEFAULT_PARAM_SETS, fh, indent=2)
        fh.write("\n")
    return DEFAULT_PARAM_SETS


def _ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _sha256_array(array: Any) -> str:
    import numpy as np

    arr = np.asarray(array)
    return hashlib.sha256(np.ascontiguousarray(arr).view(np.uint8)).hexdigest()


def _stats_total_seconds(stats: Any, code_obj: Any) -> float:
    key = (code_obj.co_filename, code_obj.co_firstlineno, code_obj.co_name)
    timings = stats.timings.get(key)
    if timings is None:
        return 0.0
    total = sum(line_time for _, _, line_time in timings)
    return total * stats.unit


def _path_is_within(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def _relative_to_root(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _safe_mean(values: Iterable[float]) -> float:
    values = list(values)
    return sum(values) / len(values) if values else 0.0


def _relative_diff(a: Any, b: Any, *, floor: float = 1e-15) -> float:
    import numpy as np

    lhs = np.asarray(a)
    rhs = np.asarray(b)
    diff = np.abs(lhs - rhs)
    denom = np.maximum(np.maximum(np.abs(lhs), np.abs(rhs)), floor)
    return float(np.max(diff / denom)) if diff.size else 0.0


def _absolute_diff(a: Any, b: Any) -> float:
    import numpy as np

    diff = np.abs(np.asarray(a) - np.asarray(b))
    return float(np.max(diff)) if diff.size else 0.0


def _normalized_overlap(a: Any, b: Any) -> float:
    import numpy as np

    lhs = np.asarray(a).ravel()
    rhs = np.asarray(b).ravel()
    lhs_norm = float(np.vdot(lhs, lhs).real)
    rhs_norm = float(np.vdot(rhs, rhs).real)
    if lhs_norm <= 0.0 or rhs_norm <= 0.0:
        return float("nan")
    overlap = abs(np.vdot(lhs, rhs)) / (lhs_norm * rhs_norm) ** 0.5
    return float(overlap)


def _mismatch(a: Any, b: Any) -> float:
    overlap = _normalized_overlap(a, b)
    if overlap != overlap:
        return float("nan")
    return float(1.0 - overlap)


def _max_diff_metrics(ref: Any, head: Any) -> dict[str, float]:
    return {
        "max_abs_diff": _absolute_diff(ref, head),
        "max_rel_diff": _relative_diff(ref, head),
        "mismatch": _mismatch(ref, head),
    }


def _serializable_array_summary(array: Any) -> dict[str, Any]:
    import numpy as np

    arr = np.asarray(array)
    return {
        "shape": list(arr.shape),
        "dtype": str(arr.dtype),
        "sha256": _sha256_array(arr),
        "l2_norm": float(np.linalg.norm(arr.ravel())),
        "mean_abs": float(np.mean(np.abs(arr))) if arr.size else 0.0,
        "max_abs": float(np.max(np.abs(arr))) if arr.size else 0.0,
    }


@dataclass
class WorktreeInfo:
    label: str
    requested_commit: str
    resolved_commit: str
    root: Path


@dataclass
class ScenarioOutputs:
    time_grid: Any | None
    waveform: Any | None
    amp: Any | None
    phase: Any | None
    avg_wall_time_sec: float


def _make_worktree(repo_root: Path, label: str, commit: str, parent_dir: Path) -> WorktreeInfo:
    resolved_commit = _run_git(repo_root, "rev-parse", "--verify", f"{commit}^{{commit}}")
    worktree_root = parent_dir / f"{label}-{resolved_commit[:12]}"
    if worktree_root.exists():
        shutil.rmtree(worktree_root)
    _run_git(repo_root, "worktree", "add", "--detach", str(worktree_root), resolved_commit)
    return WorktreeInfo(
        label=label,
        requested_commit=commit,
        resolved_commit=resolved_commit,
        root=worktree_root,
    )


def _remove_worktree(repo_root: Path, worktree_root: Path) -> None:
    try:
        _run_git(repo_root, "worktree", "remove", "--force", str(worktree_root))
    except subprocess.CalledProcessError:
        pass
    shutil.rmtree(worktree_root, ignore_errors=True)


def _call_count_for_profiler(repetitions: int) -> dict[str, int]:
    return {
        "CircularSurrogate.__call__": repetitions * 3,
        "EccentricSurrogate.__call__": repetitions,
    }


def _run_worker(
    *,
    repo_root: Path,
    worktree: WorktreeInfo,
    param_index: int,
    param_set: dict[str, Any],
    repetitions: int,
    output_dir: Path,
) -> Path:
    worker_root = _ensure_dir(output_dir / worktree.label / f"param_{param_index:02d}")
    worker_json = worker_root / "worker_results.json"
    profile_txt = worker_root / "line_profiler.txt"

    env = os.environ.copy()
    env["PYTHONPATH"] = (
        f"{worktree.root}{os.pathsep}{env['PYTHONPATH']}"
        if env.get("PYTHONPATH")
        else str(worktree.root)
    )
    env["ESIGMASUR_DATA_PATH"] = str(worktree.root / "esigmapy" / "surrogate" / "data")

    child_script = r"""
import json
import os
import re
import time
from pathlib import Path

import numpy as np
from line_profiler import LineProfiler

import esigmapy

print(esigmapy.__file__, flush=True)

expected_root = Path(os.environ["ESIGMAPY_WORKTREE_ROOT"]).resolve()
package_file = Path(esigmapy.__file__).resolve()
if expected_root not in package_file.parents:
    raise RuntimeError(f"esigmapy imported from {package_file}, expected worktree root {expected_root}")

import esigmapy.surrogate.surrogate as surrogate_module

print(surrogate_module.__file__, flush=True)
module_file = Path(surrogate_module.__file__).resolve()
if expected_root not in module_file.parents:
    raise RuntimeError(f"surrogate module imported from {module_file}, expected worktree root {expected_root}")

out_dir = Path(os.environ["ESIGMASUR_OUTPUT_DIR"]).resolve()
out_dir.mkdir(parents=True, exist_ok=True)
profile_txt = Path(os.environ["ESIGMASUR_PROFILE_TXT"]).resolve()
worker_json = out_dir / "worker_results.json"

param = json.loads(os.environ["ESIGMASUR_PARAM_JSON"])
repetitions = int(os.environ["ESIGMASUR_REPETITIONS"])

circ = surrogate_module.CircularSurrogate(
    data_dir=str(expected_root / "esigmapy" / "surrogate" / "data" / "circ_sur_data")
)
ecc = surrogate_module.EccentricSurrogate(
    ecc_data_dir=str(expected_root / "esigmapy" / "surrogate" / "data" / "ecc_sur_data"),
    circ_data_dir=str(expected_root / "esigmapy" / "surrogate" / "data" / "circ_sur_data"),
)

lp = LineProfiler()
lp.add_function(surrogate_module.CircularSurrogate.__call__)
lp.add_function(surrogate_module.EccentricSurrogate.__call__)

M = param["M"]
q = param["q"]
reference_mean_anomaly = param["reference_mean_anomaly"]
reference_eccentricity = param["reference_eccentricity"]
delta_t = param["delta_t"]
t_start = param["t_start"]

def run_batch(label, func):
    first = None
    wall_start = time.perf_counter()
    for idx in range(repetitions):
        result = func()
        if idx == 0:
            first = result
    elapsed = time.perf_counter() - wall_start
    return first, elapsed / repetitions, elapsed

lp.enable_by_count()
try:
    circ_full, circ_full_avg_wall, circ_full_total_wall = run_batch(
        "circular_full",
        lambda: circ(
            M=M,
            q=q,
            delta_t=delta_t,
            t_start=t_start,
            remove_initial_phase=False,
            return_amp_phase_only=False,
            return_orbital_variables=False,
        ),
    )
    circ_amp_phase, circ_amp_phase_avg_wall, circ_amp_phase_total_wall = run_batch(
        "circular_amp_phase",
        lambda: circ(
            M=M,
            q=q,
            delta_t=delta_t,
            t_start=t_start,
            remove_initial_phase=False,
            return_amp_phase_only=True,
            return_orbital_variables=False,
        ),
    )
    ecc_full, ecc_full_avg_wall, ecc_full_total_wall = run_batch(
        "eccentric_full",
        lambda: ecc(
            M=M,
            params=(
                q,
                reference_eccentricity,
                reference_mean_anomaly,
            ),
            delta_t=delta_t,
            t_start=t_start,
            return_orbital_variables=False,
        ),
    )
finally:
    lp.disable_by_count()

with profile_txt.open("w", encoding="utf-8") as fh:
    lp.print_stats(stream=fh)

profile_text = profile_txt.read_text(encoding="utf-8")
profile_totals = [float(value) for value in re.findall(r"^Total time:\s+([0-9.]+)\s+s$", profile_text, flags=re.MULTILINE)]
if len(profile_totals) < 2:
    raise RuntimeError(f"Could not parse profiler totals from {profile_txt}")

profiler_totals = {
    "CircularSurrogate.__call__": profile_totals[0],
    "EccentricSurrogate.__call__": profile_totals[1],
}

def arr_summary(arr):
    arr = np.asarray(arr)
    return {
        "shape": list(arr.shape),
        "dtype": str(arr.dtype),
        "sha256": __import__("hashlib").sha256(np.ascontiguousarray(arr).view(np.uint8)).hexdigest(),
        "l2_norm": float(np.linalg.norm(arr.ravel())),
        "mean_abs": float(np.mean(np.abs(arr))) if arr.size else 0.0,
        "max_abs": float(np.max(np.abs(arr))) if arr.size else 0.0,
    }

def pack_full(name, result, avg_wall, total_wall):
    t_grid, waveform = result
    artifact_path = out_dir / f"{name}.npz"
    np.savez_compressed(artifact_path, time_grid=t_grid, waveform=waveform)
    return {
        "artifact_path": str(artifact_path),
        "time_grid": arr_summary(t_grid),
        "waveform": arr_summary(waveform),
        "avg_wall_time_sec": avg_wall,
        "total_wall_time_sec": total_wall,
    }

def pack_amp_phase(name, time_grid, result, avg_wall, total_wall):
    amp, phase = result
    artifact_path = out_dir / f"{name}.npz"
    np.savez_compressed(artifact_path, time_grid=time_grid, amp=amp, phase=phase)
    return {
        "artifact_path": str(artifact_path),
        "time_grid": arr_summary(time_grid),
        "amp": arr_summary(amp),
        "phase": arr_summary(phase),
        "avg_wall_time_sec": avg_wall,
        "total_wall_time_sec": total_wall,
    }

payload = {
    "import_guard": {
        "esigmapy_file": str(package_file),
        "surrogate_module_file": str(module_file),
        "worktree_root": str(expected_root),
    },
    "param_set": param,
    "repetitions": repetitions,
    "profiling": {
        "report_path": str(profile_txt),
        "call_counts": {
            "CircularSurrogate.__call__": repetitions * 3,
            "EccentricSurrogate.__call__": repetitions,
        },
        "function_totals_sec": profiler_totals,
    },
    "scenarios": {
        "circular_full": pack_full(
            "circular_full", circ_full, circ_full_avg_wall, circ_full_total_wall
        ),
        "circular_amp_phase": pack_amp_phase(
            "circular_amp_phase",
            circ_full[0],
            circ_amp_phase,
            circ_amp_phase_avg_wall,
            circ_amp_phase_total_wall,
        ),
        "eccentric_full": pack_full(
            "eccentric_full", ecc_full, ecc_full_avg_wall, ecc_full_total_wall
        ),
    },
}

with worker_json.open("w", encoding="utf-8") as fh:
    json.dump(payload, fh, indent=2, sort_keys=True, default=str)
    fh.write("\n")

print(worker_json, flush=True)
"""

    child_env = env.copy()
    child_env["ESIGMAPY_WORKTREE_ROOT"] = str(worktree.root)
    child_env["ESIGMASUR_WORKTREE_ROOT"] = str(worktree.root)
    child_env["ESIGMASUR_OUTPUT_DIR"] = str(worker_root)
    child_env["ESIGMASUR_PROFILE_TXT"] = str(profile_txt)
    child_env["ESIGMASUR_PARAM_JSON"] = json.dumps(param_set)
    child_env["ESIGMASUR_REPETITIONS"] = str(repetitions)

    child_cwd = Path(tempfile.mkdtemp(prefix="esigmapy-surrogate-worker-"))
    try:
        result = subprocess.run(
            [sys.executable, "-c", child_script],
            cwd=str(child_cwd),
            env=child_env,
            text=True,
            capture_output=True,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        print(exc.stdout)
        print(exc.stderr, file=sys.stderr)
        raise
    finally:
        shutil.rmtree(child_cwd, ignore_errors=True)

    stdout_lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if not stdout_lines:
        raise RuntimeError("Worker did not emit the expected import-path guard output.")
    expected_prefix = str(worktree.root.resolve())
    path_lines = [line for line in stdout_lines if line.startswith(expected_prefix)]
    if len(path_lines) < 2:
        raise RuntimeError(
            "Import guard failed: expected esigmapy and surrogate module paths "
            f"under {expected_prefix}, got stdout={stdout_lines!r}"
        )

    return worker_json


def _comparison_metrics(ref: dict[str, Any], head: dict[str, Any]) -> dict[str, Any]:
    import numpy as np

    def load_arrays(payload: dict[str, Any]) -> dict[str, np.ndarray]:
        artifact = payload["artifact_path"]
        with np.load(artifact) as data:
            return {name: data[name] for name in data.files}

    ref_arrays = load_arrays(ref)
    head_arrays = load_arrays(head)

    ref_time = np.asarray(ref_arrays["time_grid"])
    head_time = np.asarray(head_arrays["time_grid"])
    time_grid_match = ref_time.shape == head_time.shape and np.array_equal(ref_time, head_time)

    full_ref = np.asarray(ref_arrays.get("waveform", []))
    full_head = np.asarray(head_arrays.get("waveform", []))
    amp_ref = np.asarray(ref_arrays.get("amp", []))
    amp_head = np.asarray(head_arrays.get("amp", []))
    phase_ref = np.asarray(ref_arrays.get("phase", []))
    phase_head = np.asarray(head_arrays.get("phase", []))

    amp_phase_ref = amp_ref * np.exp(-1j * phase_ref) if amp_ref.size else amp_ref
    amp_phase_head = amp_head * np.exp(-1j * phase_head) if amp_head.size else amp_head

    return {
        "time_grid_match": time_grid_match,
        "full_waveform": _max_diff_metrics(full_ref, full_head),
        "amp": _max_diff_metrics(amp_ref, amp_head),
        "phase": _max_diff_metrics(phase_ref, phase_head),
    }


def _load_worker_result(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _flatten_rows(
    *,
    ref_commit: str,
    head_commit: str,
    param_index: int,
    param_set: dict[str, Any],
    ref_payload: dict[str, Any],
    head_payload: dict[str, Any],
    comparison: dict[str, Any],
) -> list[dict[str, Any]]:
    ref_circ_full = ref_payload["scenarios"]["circular_full"]
    ref_circ_amp = ref_payload["scenarios"]["circular_amp_phase"]
    ref_ecc_full = ref_payload["scenarios"]["eccentric_full"]
    head_circ_full = head_payload["scenarios"]["circular_full"]
    head_circ_amp = head_payload["scenarios"]["circular_amp_phase"]
    head_ecc_full = head_payload["scenarios"]["eccentric_full"]

    def profile_avg(payload: dict[str, Any], name: str) -> float:
        total = payload["profiling"]["function_totals_sec"][name]
        count = payload["profiling"]["call_counts"][name]
        return total / count if count else float("nan")

    rows = [
        {
            "param_index": param_index,
            "scenario": "circular_full",
            "commit_role": "ref",
            "commit": ref_commit,
            "avg_wall_time_sec": ref_circ_full["avg_wall_time_sec"],
            "profile_avg_time_sec": profile_avg(ref_payload, "CircularSurrogate.__call__"),
            "time_grid_match": comparison["time_grid_match"],
            "max_abs_diff": comparison["full_waveform"]["max_abs_diff"],
            "max_rel_diff": comparison["full_waveform"]["max_rel_diff"],
            "mismatch": comparison["full_waveform"]["mismatch"],
        },
        {
            "param_index": param_index,
            "scenario": "circular_full",
            "commit_role": "head",
            "commit": head_commit,
            "avg_wall_time_sec": head_circ_full["avg_wall_time_sec"],
            "profile_avg_time_sec": profile_avg(head_payload, "CircularSurrogate.__call__"),
            "time_grid_match": comparison["time_grid_match"],
            "max_abs_diff": comparison["full_waveform"]["max_abs_diff"],
            "max_rel_diff": comparison["full_waveform"]["max_rel_diff"],
            "mismatch": comparison["full_waveform"]["mismatch"],
        },
        {
            "param_index": param_index,
            "scenario": "circular_amp_phase",
            "commit_role": "ref",
            "commit": ref_commit,
            "avg_wall_time_sec": ref_circ_amp["avg_wall_time_sec"],
            "profile_avg_time_sec": profile_avg(ref_payload, "CircularSurrogate.__call__"),
            "time_grid_match": comparison["time_grid_match"],
            "max_abs_diff_amp": comparison["amp"]["max_abs_diff"],
            "max_rel_diff_amp": comparison["amp"]["max_rel_diff"],
                "max_abs_diff_phase": comparison["phase"]["max_abs_diff"],
                "max_rel_diff_phase": comparison["phase"]["max_rel_diff"],
            },


        {
            "param_index": param_index,
            "scenario": "circular_amp_phase",
            "commit_role": "head",
            "commit": head_commit,
            "avg_wall_time_sec": head_circ_amp["avg_wall_time_sec"],
            "profile_avg_time_sec": profile_avg(head_payload, "CircularSurrogate.__call__"),
            "time_grid_match": comparison["time_grid_match"],
            "max_abs_diff_amp": comparison["amp"]["max_abs_diff"],
            "max_rel_diff_amp": comparison["amp"]["max_rel_diff"],
                "max_abs_diff_phase": comparison["phase"]["max_abs_diff"],
                "max_rel_diff_phase": comparison["phase"]["max_rel_diff"],
            },


        {
            "param_index": param_index,
            "scenario": "eccentric_full",
            "commit_role": "ref",
            "commit": ref_commit,
            "avg_wall_time_sec": ref_ecc_full["avg_wall_time_sec"],
            "profile_avg_time_sec": ref_payload["profiling"]["function_totals_sec"]["EccentricSurrogate.__call__"]
            / ref_payload["profiling"]["call_counts"]["EccentricSurrogate.__call__"],
            "time_grid_match": comparison["time_grid_match"],
            "max_abs_diff": comparison["full_waveform"]["max_abs_diff"],
            "max_rel_diff": comparison["full_waveform"]["max_rel_diff"],
            "mismatch": comparison["full_waveform"]["mismatch"],
        },
        {
            "param_index": param_index,
            "scenario": "eccentric_full",
            "commit_role": "head",
            "commit": head_commit,
            "avg_wall_time_sec": head_ecc_full["avg_wall_time_sec"],
            "profile_avg_time_sec": head_payload["profiling"]["function_totals_sec"]["EccentricSurrogate.__call__"]
            / head_payload["profiling"]["call_counts"]["EccentricSurrogate.__call__"],
            "time_grid_match": comparison["time_grid_match"],
            "max_abs_diff": comparison["full_waveform"]["max_abs_diff"],
            "max_rel_diff": comparison["full_waveform"]["max_rel_diff"],
            "mismatch": comparison["full_waveform"]["mismatch"],
        },
    ]
    return rows


def _print_summary(
    *,
    ref_commit: str,
    head_commit: str,
    param_index: int,
    param_set: dict[str, Any],
    ref_payload: dict[str, Any],
    head_payload: dict[str, Any],
    comparison: dict[str, Any],
) -> None:
    print(f"Parameter set {param_index}: {json.dumps(param_set, sort_keys=True)}")
    for scenario in ("circular_full", "circular_amp_phase", "eccentric_full"):
        ref_s = ref_payload["scenarios"][scenario]
        head_s = head_payload["scenarios"][scenario]
        if scenario == "circular_amp_phase":
            ref_avg = ref_s["avg_wall_time_sec"]
            head_avg = head_s["avg_wall_time_sec"]
            print(
                f"  {scenario}: ref {ref_avg:.6e}s/call, head {head_avg:.6e}s/call, "
                f"delta {(head_avg - ref_avg):+.6e}s/call"
            )
            print(
                f"    amp max_abs={comparison['amp']['max_abs_diff']:.6e}, "
                f"amp max_rel={comparison['amp']['max_rel_diff']:.6e}, "
                f"phase max_abs={comparison['phase']['max_abs_diff']:.6e}, "
                f"phase max_rel={comparison['phase']['max_rel_diff']:.6e}"
            )
        else:
            ref_avg = ref_s["avg_wall_time_sec"]
            head_avg = head_s["avg_wall_time_sec"]
            print(
                f"  {scenario}: ref {ref_avg:.6e}s/call, head {head_avg:.6e}s/call, "
                f"delta {(head_avg - ref_avg):+.6e}s/call"
            )
            print(
                f"    max_abs={comparison['full_waveform']['max_abs_diff']:.6e}, "
                f"max_rel={comparison['full_waveform']['max_rel_diff']:.6e}, "
                f"mismatch={comparison['full_waveform']['mismatch']:.6e}, "
                f"time_grid_match={comparison['time_grid_match']}"
            )
    print(
        f"  profile reports: ref={ref_payload['profiling']['function_totals_sec']} "
        f"head={head_payload['profiling']['function_totals_sec']}"
    )


def _worker_main() -> None:
    raise RuntimeError("Worker mode is invoked through the generated child script.")


def _parent_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Compare surrogate outputs and line profiles between two commits."
    )
    parser.add_argument(
        "--ref-commit",
        default="HEAD~1",
        help="Reference git commit to benchmark.",
    )
    parser.add_argument(
        "--head-commit",
        default="HEAD",
        help="Head git commit to benchmark.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=_default_config_path(),
        help="JSON file containing parameter sets.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for artifacts. Defaults to a temp directory.",
    )
    parser.add_argument(
        "--repetitions",
        type=int,
        default=50,
        help="Number of repeated calls to profile for each scenario.",
    )
    parser.add_argument(
        "--keep-worktrees",
        action="store_true",
        help="Keep detached worktrees after the run.",
    )
    args = parser.parse_args(argv)

    repo_root = _repo_root()
    config_path = args.config
    param_sets = _ensure_param_config(config_path)
    output_dir = _ensure_dir(args.output_dir or _default_output_dir(args.ref_commit, args.head_commit))
    worktree_parent = _ensure_dir(output_dir / "worktrees")
    worker_output_root = _ensure_dir(output_dir / "workers")

    ref_worktree = head_worktree = None
    try:
        ref_worktree = _make_worktree(repo_root, "ref", args.ref_commit, worktree_parent)
        head_worktree = _make_worktree(repo_root, "head", args.head_commit, worktree_parent)

        worker_results: dict[tuple[str, int], Path] = {}
        for commit_info in (ref_worktree, head_worktree):
            for index, param_set in enumerate(param_sets):
                worker_results[(commit_info.label, index)] = _run_worker(
                    repo_root=repo_root,
                    worktree=commit_info,
                    param_index=index,
                    param_set=param_set,
                    repetitions=args.repetitions,
                    output_dir=worker_output_root,
                )

        raw_results = {
            "repo_root": str(repo_root),
            "config_path": str(config_path),
            "output_dir": str(output_dir),
            "ref_commit": dataclasses.asdict(ref_worktree),
            "head_commit": dataclasses.asdict(head_worktree),
            "repetitions": args.repetitions,
            "parameter_sets": param_sets,
            "runs": {},
            "comparisons": [],
        }

        rows: list[dict[str, Any]] = []
        for index, param_set in enumerate(param_sets):
            ref_payload = _load_worker_result(worker_results[("ref", index)])
            head_payload = _load_worker_result(worker_results[("head", index)])
            ref_circ_full = ref_payload["scenarios"]["circular_full"]
            head_circ_full = head_payload["scenarios"]["circular_full"]
            ref_circ_amp = ref_payload["scenarios"]["circular_amp_phase"]
            head_circ_amp = head_payload["scenarios"]["circular_amp_phase"]
            ref_ecc_full = ref_payload["scenarios"]["eccentric_full"]
            head_ecc_full = head_payload["scenarios"]["eccentric_full"]
            comparison = {
                "param_index": index,
                "param_set": param_set,
                "circular_full": _comparison_metrics(ref_circ_full, head_circ_full),
                "circular_amp_phase": _comparison_metrics(ref_circ_amp, head_circ_amp),
                "eccentric_full": _comparison_metrics(ref_ecc_full, head_ecc_full),
            }
            raw_results["runs"][f"ref:param_{index:02d}"] = ref_payload
            raw_results["runs"][f"head:param_{index:02d}"] = head_payload
            raw_results["comparisons"].append(comparison)

            def profile_avg(payload: dict[str, Any], name: str) -> float:
                total = payload["profiling"]["function_totals_sec"][name]
                count = payload["profiling"]["call_counts"][name]
                return total / count if count else float("nan")

            _print_summary(
                ref_commit=ref_worktree.resolved_commit,
                head_commit=head_worktree.resolved_commit,
                param_index=index,
                param_set=param_set,
                ref_payload=ref_payload,
                head_payload=head_payload,
                comparison=comparison["circular_full"],
            )

            ref_circ_full_row = {
                "param_index": index,
                "scenario": "circular_full",
                "commit_role": "ref",
                "commit": ref_worktree.resolved_commit,
                "avg_wall_time_sec": ref_circ_full["avg_wall_time_sec"],
                "profile_avg_time_sec": profile_avg(ref_payload, "CircularSurrogate.__call__"),
                "time_grid_match": comparison["circular_full"]["time_grid_match"],
                "max_abs_diff": comparison["circular_full"]["full_waveform"]["max_abs_diff"],
                "max_rel_diff": comparison["circular_full"]["full_waveform"]["max_rel_diff"],
                "mismatch": comparison["circular_full"]["full_waveform"]["mismatch"],
            }
            head_circ_full_row = {
                "param_index": index,
                "scenario": "circular_full",
                "commit_role": "head",
                "commit": head_worktree.resolved_commit,
                "avg_wall_time_sec": head_circ_full["avg_wall_time_sec"],
                "profile_avg_time_sec": profile_avg(head_payload, "CircularSurrogate.__call__"),
                "time_grid_match": comparison["circular_full"]["time_grid_match"],
                "max_abs_diff": comparison["circular_full"]["full_waveform"]["max_abs_diff"],
                "max_rel_diff": comparison["circular_full"]["full_waveform"]["max_rel_diff"],
                "mismatch": comparison["circular_full"]["full_waveform"]["mismatch"],
            }

            ref_amp_phase_row = {
                "param_index": index,
                "scenario": "circular_amp_phase",
                "commit_role": "ref",
                "commit": ref_worktree.resolved_commit,
                "avg_wall_time_sec": ref_circ_amp["avg_wall_time_sec"],
                "profile_avg_time_sec": profile_avg(ref_payload, "CircularSurrogate.__call__"),
                "time_grid_match": comparison["circular_amp_phase"]["time_grid_match"],
                "max_abs_diff_amp": comparison["circular_amp_phase"]["amp"]["max_abs_diff"],
                "max_rel_diff_amp": comparison["circular_amp_phase"]["amp"]["max_rel_diff"],
                "max_abs_diff_phase": comparison["circular_amp_phase"]["phase"]["max_abs_diff"],
                "max_rel_diff_phase": comparison["circular_amp_phase"]["phase"]["max_rel_diff"],
            }
            head_amp_phase_row = {
                "param_index": index,
                "scenario": "circular_amp_phase",
                "commit_role": "head",
                "commit": head_worktree.resolved_commit,
                "avg_wall_time_sec": head_circ_amp["avg_wall_time_sec"],
                "profile_avg_time_sec": profile_avg(head_payload, "CircularSurrogate.__call__"),
                "time_grid_match": comparison["circular_amp_phase"]["time_grid_match"],
                "max_abs_diff_amp": comparison["circular_amp_phase"]["amp"]["max_abs_diff"],
                "max_rel_diff_amp": comparison["circular_amp_phase"]["amp"]["max_rel_diff"],
                "max_abs_diff_phase": comparison["circular_amp_phase"]["phase"]["max_abs_diff"],
                "max_rel_diff_phase": comparison["circular_amp_phase"]["phase"]["max_rel_diff"],
            }

            ref_ecc_full_row = {
                "param_index": index,
                "scenario": "eccentric_full",
                "commit_role": "ref",
                "commit": ref_worktree.resolved_commit,
                "avg_wall_time_sec": ref_ecc_full["avg_wall_time_sec"],
                "profile_avg_time_sec": profile_avg(ref_payload, "EccentricSurrogate.__call__"),
                "time_grid_match": comparison["eccentric_full"]["time_grid_match"],
                "max_abs_diff": comparison["eccentric_full"]["full_waveform"]["max_abs_diff"],
                "max_rel_diff": comparison["eccentric_full"]["full_waveform"]["max_rel_diff"],
                "mismatch": comparison["eccentric_full"]["full_waveform"]["mismatch"],
            }
            head_ecc_full_row = {
                "param_index": index,
                "scenario": "eccentric_full",
                "commit_role": "head",
                "commit": head_worktree.resolved_commit,
                "avg_wall_time_sec": head_ecc_full["avg_wall_time_sec"],
                "profile_avg_time_sec": profile_avg(head_payload, "EccentricSurrogate.__call__"),
                "time_grid_match": comparison["eccentric_full"]["time_grid_match"],
                "max_abs_diff": comparison["eccentric_full"]["full_waveform"]["max_abs_diff"],
                "max_rel_diff": comparison["eccentric_full"]["full_waveform"]["max_rel_diff"],
                "mismatch": comparison["eccentric_full"]["full_waveform"]["mismatch"],
            }

            rows.extend(
                [
                    ref_circ_full_row,
                    head_circ_full_row,
                    ref_amp_phase_row,
                    head_amp_phase_row,
                    ref_ecc_full_row,
                    head_ecc_full_row,
                    {
                        "param_index": index,
                        "scenario": "circular_full",
                        "commit_role": "comparison",
                        "commit": f"{ref_worktree.resolved_commit}..{head_worktree.resolved_commit}",
                        "ref_avg_wall_time_sec": ref_circ_full["avg_wall_time_sec"],
                        "head_avg_wall_time_sec": head_circ_full["avg_wall_time_sec"],
                        "wall_delta_sec": head_circ_full["avg_wall_time_sec"] - ref_circ_full["avg_wall_time_sec"],
                        "wall_delta_pct": (
                            (head_circ_full["avg_wall_time_sec"] - ref_circ_full["avg_wall_time_sec"])
                            / ref_circ_full["avg_wall_time_sec"]
                            if ref_circ_full["avg_wall_time_sec"]
                            else float("nan")
                        ),
                        "ref_profile_avg_time_sec": profile_avg(ref_payload, "CircularSurrogate.__call__"),
                        "head_profile_avg_time_sec": profile_avg(head_payload, "CircularSurrogate.__call__"),
                        "profile_delta_sec": profile_avg(head_payload, "CircularSurrogate.__call__")
                        - profile_avg(ref_payload, "CircularSurrogate.__call__"),
                        "profile_delta_pct": (
                            (profile_avg(head_payload, "CircularSurrogate.__call__")
                             - profile_avg(ref_payload, "CircularSurrogate.__call__"))
                            / profile_avg(ref_payload, "CircularSurrogate.__call__")
                            if profile_avg(ref_payload, "CircularSurrogate.__call__")
                            else float("nan")
                        ),
                        "time_grid_match": comparison["circular_full"]["time_grid_match"],
                        "max_abs_diff": comparison["circular_full"]["full_waveform"]["max_abs_diff"],
                        "max_rel_diff": comparison["circular_full"]["full_waveform"]["max_rel_diff"],
                        "mismatch": comparison["circular_full"]["full_waveform"]["mismatch"],
                    },
                    {
                        "param_index": index,
                        "scenario": "circular_amp_phase",
                        "commit_role": "comparison",
                        "commit": f"{ref_worktree.resolved_commit}..{head_worktree.resolved_commit}",
                        "ref_avg_wall_time_sec": ref_circ_amp["avg_wall_time_sec"],
                        "head_avg_wall_time_sec": head_circ_amp["avg_wall_time_sec"],
                        "wall_delta_sec": head_circ_amp["avg_wall_time_sec"] - ref_circ_amp["avg_wall_time_sec"],
                        "wall_delta_pct": (
                            (head_circ_amp["avg_wall_time_sec"] - ref_circ_amp["avg_wall_time_sec"])
                            / ref_circ_amp["avg_wall_time_sec"]
                            if ref_circ_amp["avg_wall_time_sec"]
                            else float("nan")
                        ),
                        "ref_profile_avg_time_sec": profile_avg(ref_payload, "CircularSurrogate.__call__"),
                        "head_profile_avg_time_sec": profile_avg(head_payload, "CircularSurrogate.__call__"),
                        "profile_delta_sec": profile_avg(head_payload, "CircularSurrogate.__call__")
                        - profile_avg(ref_payload, "CircularSurrogate.__call__"),
                        "profile_delta_pct": (
                            (profile_avg(head_payload, "CircularSurrogate.__call__")
                             - profile_avg(ref_payload, "CircularSurrogate.__call__"))
                            / profile_avg(ref_payload, "CircularSurrogate.__call__")
                            if profile_avg(ref_payload, "CircularSurrogate.__call__")
                            else float("nan")
                        ),
                        "time_grid_match": comparison["circular_amp_phase"]["time_grid_match"],
                        "max_abs_diff_amp": comparison["circular_amp_phase"]["amp"]["max_abs_diff"],
                        "max_rel_diff_amp": comparison["circular_amp_phase"]["amp"]["max_rel_diff"],
                    "max_abs_diff_phase": comparison["circular_amp_phase"]["phase"]["max_abs_diff"],
                    "max_rel_diff_phase": comparison["circular_amp_phase"]["phase"]["max_rel_diff"],
                },




                    {
                        "param_index": index,
                        "scenario": "eccentric_full",
                        "commit_role": "comparison",
                        "commit": f"{ref_worktree.resolved_commit}..{head_worktree.resolved_commit}",
                        "ref_avg_wall_time_sec": ref_ecc_full["avg_wall_time_sec"],
                        "head_avg_wall_time_sec": head_ecc_full["avg_wall_time_sec"],
                        "wall_delta_sec": head_ecc_full["avg_wall_time_sec"] - ref_ecc_full["avg_wall_time_sec"],
                        "wall_delta_pct": (
                            (head_ecc_full["avg_wall_time_sec"] - ref_ecc_full["avg_wall_time_sec"])
                            / ref_ecc_full["avg_wall_time_sec"]
                            if ref_ecc_full["avg_wall_time_sec"]
                            else float("nan")
                        ),
                        "ref_profile_avg_time_sec": profile_avg(ref_payload, "EccentricSurrogate.__call__"),
                        "head_profile_avg_time_sec": profile_avg(head_payload, "EccentricSurrogate.__call__"),
                        "profile_delta_sec": profile_avg(head_payload, "EccentricSurrogate.__call__")
                        - profile_avg(ref_payload, "EccentricSurrogate.__call__"),
                        "profile_delta_pct": (
                            (profile_avg(head_payload, "EccentricSurrogate.__call__")
                             - profile_avg(ref_payload, "EccentricSurrogate.__call__"))
                            / profile_avg(ref_payload, "EccentricSurrogate.__call__")
                            if profile_avg(ref_payload, "EccentricSurrogate.__call__")
                            else float("nan")
                        ),
                        "time_grid_match": comparison["eccentric_full"]["time_grid_match"],
                        "max_abs_diff": comparison["eccentric_full"]["full_waveform"]["max_abs_diff"],
                        "max_rel_diff": comparison["eccentric_full"]["full_waveform"]["max_rel_diff"],
                        "mismatch": comparison["eccentric_full"]["full_waveform"]["mismatch"],
                    },
                ]
            )

        raw_json = output_dir / "raw_results.json"
        with raw_json.open("w", encoding="utf-8") as fh:
            json.dump(raw_results, fh, indent=2, sort_keys=True, default=_json_default)
            fh.write("\n")

        csv_path = output_dir / "comparison_summary.csv"
        fieldnames = sorted({key for row in rows for key in row.keys()})
        with csv_path.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        print(f"Raw JSON: {raw_json}")
        print(f"CSV summary: {csv_path}")
        print(f"Artifacts: {output_dir}")
        return 0
    finally:
        if not args.keep_worktrees:
            for wt in (ref_worktree, head_worktree):
                if wt is not None:
                    _remove_worktree(repo_root, wt.root)


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if "--worker" in argv:
        _worker_main()
    return _parent_main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
