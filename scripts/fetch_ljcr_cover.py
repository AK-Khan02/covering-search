#!/usr/bin/env python3
"""Fetch one covering-design certificate from the La Jolla Covering Repository."""

from __future__ import annotations

import argparse
import html
import re
import urllib.request
from pathlib import Path


def parse_pre_blocks(page: str) -> list[str]:
    match = re.search(r"<pre[^>]*>(.*?)</pre>", page, flags=re.IGNORECASE | re.DOTALL)
    if not match:
        raise ValueError("could not find <pre> block in LJCR response")
    body = html.unescape(match.group(1))
    rows: list[str] = []
    for raw_line in body.splitlines():
        line = raw_line.strip()
        if line:
            rows.append(" ".join(line.split()))
    return rows


def extract_header(page: str) -> str:
    match = re.search(r"<h1[^>]*>(.*?)</h1>", page, flags=re.IGNORECASE | re.DOTALL)
    if not match:
        return ""
    text = re.sub(r"<[^>]+>", "", html.unescape(match.group(1)))
    return " ".join(text.split())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--v", type=int, required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--t", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    url = f"https://ljcr.dmgordon.org/cover/show_cover.php?v={args.v}&k={args.k}&t={args.t}"
    with urllib.request.urlopen(url, timeout=30) as response:
        page = response.read().decode("utf-8", errors="replace")

    rows = parse_pre_blocks(page)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        handle.write(f"# source={url}\n")
        header = extract_header(page)
        if header:
            handle.write(f"# ljcr_header={header}\n")
        handle.write(f"# C({args.v},{args.k},{args.t}) blocks={len(rows)}\n")
        for row in rows:
            handle.write(row + "\n")

    print(f"wrote {len(rows)} blocks to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

