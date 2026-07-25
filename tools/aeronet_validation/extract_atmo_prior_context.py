"""Recover the atmospheric prior consumed by saved validation retrievals.

New retrievals persist this summary directly.  For older experiment records,
the runner log is the authoritative receipt because it records the exact live
MAIAC value returned before the solver ran.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
_MID_RE = re.compile(r"(?:TASK\s+\d+\s+->\s+option=\S+\s+mid=|Array task\s+\S+\s+->\s+)(\S+)")
_LIVE_RE = re.compile(
    r"GEE MAIAC atmo prior .*?aot=(?P<aot>[-+0-9.eE]+)\s+tcwv=(?P<tcwv>[-+0-9.eE]+)"
)
_STAGED_RE = re.compile(r"Prebuilt staged MAIAC backstop: aot=(?P<aot>[-+0-9.eE]+)")


def _from_result(record: dict[str, Any]) -> dict[str, float] | None:
    summary = record.get("atmo_prior")
    if not isinstance(summary, dict):
        return None
    out: dict[str, float] = {}
    for key, value in summary.items():
        if key.endswith("_mean"):
            out[key.removesuffix("_mean")] = float(value)
    return out or None


def _parse_log(path: Path) -> tuple[str, dict[str, float]] | None:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    mid_matches = list(_MID_RE.finditer(text))
    if not mid_matches:
        return None
    matchup_id = mid_matches[-1].group(1)
    live = list(_LIVE_RE.finditer(text))
    if live:
        match = live[-1]
        return matchup_id, {"aot": float(match.group("aot")), "tcwv": float(match.group("tcwv"))}
    staged = list(_STAGED_RE.finditer(text))
    if staged:
        return matchup_id, {"aot": float(staged[-1].group("aot"))}
    return None


def extract(root: Path, result_dir: Path, mids: list[str], output_dir: Path) -> dict[str, int]:
    wanted = set(mids)
    recovered: dict[str, tuple[int, dict[str, float], str]] = {}
    for path in (root / "slurm").glob("lc20s2adapt_*.out"):
        parsed = _parse_log(path)
        if parsed is None or parsed[0] not in wanted:
            continue
        matchup_id, values = parsed
        stamp = path.stat().st_mtime_ns
        previous = recovered.get(matchup_id)
        if previous is None or stamp > previous[0]:
            recovered[matchup_id] = (stamp, values, str(path))

    output_dir.mkdir(parents=True, exist_ok=True)
    counts = {"result": 0, "log": 0, "missing": 0, "skipped_non_ok": 0}
    for matchup_id in mids:
        result_path = result_dir / f"{matchup_id}.json"
        if not result_path.exists():
            counts["missing"] += 1
            continue
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if result.get("status") != "OK":
            counts["skipped_non_ok"] += 1
            continue
        values = _from_result(result)
        source = "saved_result"
        source_log = None
        if values is None:
            receipt = recovered.get(matchup_id)
            if receipt is None:
                counts["missing"] += 1
                continue
            values = receipt[1]
            source = "runner_log"
            source_log = receipt[2]
            counts["log"] += 1
        else:
            counts["result"] += 1
        record: dict[str, Any] = {
            "matchup_id": matchup_id,
            "status": "OK",
            "source": source,
            **values,
        }
        if source_log is not None:
            record["source_log"] = source_log
        (output_dir / f"{matchup_id}.json").write_text(
            json.dumps(record, indent=2) + "\n",
            encoding="utf-8",
        )
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--mids", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    mids = [line.strip() for line in args.mids.read_text(encoding="utf-8").splitlines() if line.strip()]
    counts = extract(args.root, args.result_dir, mids, args.output_dir)
    print(json.dumps(counts, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
