#!/usr/bin/env python3
"""Verify a covering-design certificate."""

from __future__ import annotations

import argparse
import itertools
import math
import re
from pathlib import Path


def read_blocks(path: Path, v: int, k: int) -> list[tuple[int, ...]]:
    rows: list[list[int]] = []
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        values = [int(token) for token in re.split(r"[\s,]+", line) if token]
        if len(values) != k:
            raise ValueError(f"line {line_number}: expected block size {k}, found {len(values)}")
        rows.append(values)

    if not rows:
        raise ValueError("certificate contains no blocks")

    minimum = min(min(row) for row in rows)
    maximum = max(max(row) for row in rows)
    if minimum == 0:
        if maximum >= v:
            raise ValueError(f"zero-based point {maximum} is outside 0..{v - 1}")
        normalized = rows
    else:
        if minimum < 1 or maximum > v:
            raise ValueError(f"one-based points must be in 1..{v}; saw {minimum}..{maximum}")
        normalized = [[point - 1 for point in row] for row in rows]

    blocks: list[tuple[int, ...]] = []
    for index, row in enumerate(normalized, start=1):
        block = tuple(sorted(row))
        if len(set(block)) != k:
            raise ValueError(f"block {index} contains repeated points: {row}")
        blocks.append(block)
    return blocks


def missing_subsets(blocks: list[tuple[int, ...]], v: int, t: int) -> list[tuple[int, ...]]:
    covered: set[tuple[int, ...]] = set()
    for block in blocks:
        covered.update(itertools.combinations(block, t))

    missing: list[tuple[int, ...]] = []
    for subset in itertools.combinations(range(v), t):
        if subset not in covered:
            missing.append(subset)
    return missing


def counting_lower_bound(v: int, k: int, t: int) -> int:
    return math.ceil(math.comb(v, t) / math.comb(k, t))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--v", type=int, required=True, help="number of points")
    parser.add_argument("--k", type=int, required=True, help="block size")
    parser.add_argument("--t", type=int, required=True, help="covered subset size")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not (0 < args.t <= args.k <= args.v):
        raise SystemExit("parameters must satisfy 0 < t <= k <= v")

    blocks = read_blocks(args.certificate, args.v, args.k)
    missing = missing_subsets(blocks, args.v, args.t)
    total = math.comb(args.v, args.t)
    lower = counting_lower_bound(args.v, args.k, args.t)

    if missing:
        print(
            f"FAIL C({args.v},{args.k},{args.t}) blocks={len(blocks)} "
            f"covered={total - len(missing)}/{total} missing={len(missing)} "
            f"counting_lb={lower}"
        )
        for subset in missing[:20]:
            print("missing", " ".join(str(point + 1) for point in subset))
        return 1

    duplicate_count = len(blocks) - len(set(blocks))
    print(
        f"PASS C({args.v},{args.k},{args.t}) <= {len(blocks)} "
        f"covered={total}/{total} counting_lb={lower} duplicates={duplicate_count}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
