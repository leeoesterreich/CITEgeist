"""Usage: run_bootstrap.py --case C [--inj] --side S [--n N] [--B B] [--procs P] -> synthetic/g6_fix/boot_*.json"""
import argparse
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2] / "scripts"))
import provenance  # noqa: E402
import synthetic_bootstrap as SB  # noqa: E402


PROTOCOL_N_BOOT = 2000
PROTOCOL_SEED = 20261014
ENV_KEYS = ("BOOT_N", "BOOT_B", "BOOT_OUT", "SLURM_JOB_ID", "SLURM_ARRAY_TASK_ID")


def run_settings(n, B, seed0, procs, env=None):
    """Everything that shaped this run, including the sbatch env vars, so the numbers are reproducible from the output alone."""
    env = os.environ if env is None else env
    return {
        "n_datasets": int(n), "n_boot": int(B), "seed0": int(seed0), "procs": int(procs),
        "bootstrap_seed_rule": "each synthetic dataset's own seed (seed0 + dataset index) is used as the block-bootstrap seed",
        "env": {k: env.get(k) for k in ENV_KEYS},
        "protocol_n_boot": PROTOCOL_N_BOOT, "protocol_seed": PROTOCOL_SEED,
        "departs_from_protocol_n_boot": int(B) != PROTOCOL_N_BOOT,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--inj", action="store_true")
    ap.add_argument("--side", type=float, default=6.0)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--B", type=int, default=300)
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--out", type=Path, default=HERE.parent / "synthetic" / "g6_fix")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    res = SB.coverage(a.case, a.side, a.n, a.B, fold=3.0 if a.inj else 1.0, procs=a.procs)
    name = f"boot_{'inj_' if a.inj else ''}{a.case}_s{int(a.side)}"
    provenance.write_stamped_json(a.out / f"{name}.json", {"run": name, "settings": run_settings(a.n, a.B, 1000, a.procs), "result": res}, "Benchmarking/midkine_remediation/src/run_bootstrap.py")
    print(name, res, flush=True)


if __name__ == "__main__":
    main()
