#!/usr/bin/env python3
"""
Create mini cytosine-report datasets for comethyl from Bismark CpG reports.

Use case
--------
- Keep a selected set of participants (e.g. 49 samples)
- Optionally restrict to a small genomic subset so the data is lightweight enough
  for GitHub and local comethyl testing
- Copy/update metadata so the mini dataset is reproducible


Supported input
---------------
Input files are expected to look like:
  <sample_id>_...CpG_report.txt.gz

Typical Bismark CpG report columns assumed:
  chrom  position  strand  methylated_count  unmethylated_count  context  trinuc

This script streams gz files line-by-line.

Example of usage
python /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/tests/scripts/make_minidata/00_make_comethyl_minidata.py \
  --input_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/processed/cytosine_reports \
  --output_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/tests/data/processed/ \
  --sample_list /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/processed/sample_ids.txt \
  --keep_chr chr22 \
  --dry_run


"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Make mini CpG_report.txt.gz data for comethyl"
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing original CpG_report.txt.gz files",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory where mini data will be written",
    )
    parser.add_argument(
        "--sample_list",
        required=True,
        help="Text file with one sample ID per line (e.g. 49 participants)",
    )
    parser.add_argument(
        "--metadata",
        default=None,
        help="Optional metadata file (csv/tsv/xlsx not supported here; use csv/tsv)",
    )
    parser.add_argument(
        "--metadata_id_col",
        default="sample_id",
        help="Column in metadata matching sample IDs",
    )
    parser.add_argument(
        "--filename_regex",
        default=r"^([^_]+)",
        help=(
            "Regex used to extract sample ID from filename. "
            "Default grabs text before first underscore."
        ),
    )
    parser.add_argument(
        "--keep_chr",
        nargs="*",
        default=None,
        help="Optional list of chromosomes to keep, e.g. chr22 chrM",
    )
    parser.add_argument(
        "--bed",
        default=None,
        help="Optional BED file of genomic regions to keep (0-based BED assumed)",
    )
    parser.add_argument(
        "--max_lines_per_file",
        type=int,
        default=None,
        help=(
            "Optional hard cap on kept rows per file after genomic filtering. "
            "Useful for making an ultra-tiny test set."
        ),
    )
    parser.add_argument(
        "--copy_full_files",
        action="store_true",
        help=(
            "Copy selected files without genomic filtering. Not recommended for GitHub if files are large."
        ),
    )
    parser.add_argument(
        "--gzip_level",
        type=int,
        default=6,
        help="gzip compression level for output files (default: 6)",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Print what would be done without writing output",
    )
    return parser.parse_args()


def read_sample_list(sample_list_path: Path) -> List[str]:
    samples: List[str] = []
    with open(sample_list_path, "r", encoding="utf-8") as handle:
        for line in handle:
            s = line.strip()
            if s and not s.startswith("#"):
                samples.append(s)
    if not samples:
        raise ValueError("Sample list is empty")
    return samples


def infer_delimiter(path: Path) -> str:
    if path.suffix.lower() == ".tsv":
        return "\t"
    return ","


def filter_metadata(
    metadata_path: Path,
    output_path: Path,
    id_col: str,
    keep_ids: Set[str],
) -> None:
    delim = infer_delimiter(metadata_path)
    with open(metadata_path, "r", encoding="utf-8", newline="") as infile:
        reader = csv.DictReader(infile, delimiter=delim)
        if reader.fieldnames is None:
            raise ValueError("Metadata file has no header")
        if id_col not in reader.fieldnames:
            raise ValueError(f"Column '{id_col}' not found in metadata")
        rows = [row for row in reader if row[id_col] in keep_ids]

    with open(output_path, "w", encoding="utf-8", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames, delimiter=delim)
        writer.writeheader()
        writer.writerows(rows)


BedMap = Dict[str, List[Tuple[int, int]]]


def load_bed(path: Path) -> BedMap:
    regions: BedMap = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            regions.setdefault(chrom, []).append((start, end))

    for chrom in regions:
        regions[chrom].sort()
    return regions


def in_bed(chrom: str, pos_1based: int, bed_map: BedMap) -> bool:
    if chrom not in bed_map:
        return False
    pos_0based = pos_1based - 1
    for start, end in bed_map[chrom]:
        if start <= pos_0based < end:
            return True
        if start > pos_0based:
            return False
    return False


def extract_sample_id(filename: str, pattern: re.Pattern[str]) -> Optional[str]:
    match = pattern.search(filename)
    if not match:
        return None
    return match.group(1)


EXPECTED_SUFFIX = "CpG_report.txt.gz"


def discover_files(
    input_dir: Path,
    selected_samples: Set[str],
    filename_regex: str,
) -> Dict[str, Path]:
    pattern = re.compile(filename_regex)
    matched: Dict[str, Path] = {}

    for path in sorted(input_dir.glob(f"*{EXPECTED_SUFFIX}")):
        sample_id = extract_sample_id(path.name, pattern)
        if sample_id is None:
            continue
        if sample_id in selected_samples:
            matched[sample_id] = path

    return matched


def keep_line(
    chrom: str,
    pos: int,
    keep_chr: Optional[Set[str]],
    bed_map: Optional[BedMap],
) -> bool:
    if keep_chr is not None and chrom not in keep_chr:
        return False
    if bed_map is not None and not in_bed(chrom, pos, bed_map):
        return False
    return True


def filter_cpg_report(
    input_path: Path,
    output_path: Path,
    keep_chr: Optional[Set[str]],
    bed_map: Optional[BedMap],
    max_lines_per_file: Optional[int],
    gzip_level: int,
) -> Tuple[int, int]:
    total = 0
    kept = 0

    with gzip.open(input_path, "rt", encoding="utf-8") as infile, gzip.open(
        output_path, "wt", encoding="utf-8", compresslevel=gzip_level
    ) as outfile:
        for line in infile:
            total += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue

            if keep_line(chrom, pos, keep_chr, bed_map):
                outfile.write(line)
                kept += 1
                if max_lines_per_file is not None and kept >= max_lines_per_file:
                    break

    return total, kept


def copy_file(input_path: Path, output_path: Path) -> int:
    shutil.copy2(input_path, output_path)
    return output_path.stat().st_size


def write_manifest(
    manifest_path: Path,
    rows: List[Dict[str, object]],
) -> None:
    fieldnames = [
        "sample_id",
        "input_file",
        "output_file",
        "mode",
        "input_lines",
        "kept_lines",
        "output_bytes",
    ]
    with open(manifest_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    sample_list_path = Path(args.sample_list).resolve()

    samples = read_sample_list(sample_list_path)
    sample_set = set(samples)

    keep_chr = set(args.keep_chr) if args.keep_chr else None
    bed_map = load_bed(Path(args.bed).resolve()) if args.bed else None

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not sample_list_path.exists():
        raise FileNotFoundError(f"Sample list not found: {sample_list_path}")

    matched = discover_files(
        input_dir=input_dir,
        selected_samples=sample_set,
        filename_regex=args.filename_regex,
    )

    missing = sorted(sample_set - set(matched.keys()))
    if missing:
        print("WARNING: The following sample IDs were not matched to files:")
        for s in missing:
            print(f"  - {s}")

    if args.dry_run:
        print(f"Selected samples requested: {len(samples)}")
        print(f"Matched files: {len(matched)}")
        print(f"Mode: {'copy_full_files' if args.copy_full_files else 'filter'}")
        if keep_chr:
            print(f"Chromosomes kept: {sorted(keep_chr)}")
        if args.bed:
            print(f"BED file: {args.bed}")
        if args.max_lines_per_file is not None:
            print(f"Max kept lines/file: {args.max_lines_per_file}")
        return

    reports_out = output_dir / "cytosine_reports"
    reports_out.mkdir(parents=True, exist_ok=True)

    manifest_rows: List[Dict[str, object]] = []

    for sample_id in samples:
        if sample_id not in matched:
            continue
        input_path = matched[sample_id]
        output_path = reports_out / input_path.name

        if args.copy_full_files:
            output_bytes = copy_file(input_path, output_path)
            row = {
                "sample_id": sample_id,
                "input_file": str(input_path),
                "output_file": str(output_path),
                "mode": "copy_full_files",
                "input_lines": "NA",
                "kept_lines": "NA",
                "output_bytes": output_bytes,
            }
        else:
            total, kept = filter_cpg_report(
                input_path=input_path,
                output_path=output_path,
                keep_chr=keep_chr,
                bed_map=bed_map,
                max_lines_per_file=args.max_lines_per_file,
                gzip_level=args.gzip_level,
            )
            row = {
                "sample_id": sample_id,
                "input_file": str(input_path),
                "output_file": str(output_path),
                "mode": "filtered",
                "input_lines": total,
                "kept_lines": kept,
                "output_bytes": output_path.stat().st_size,
            }

        manifest_rows.append(row)
        print(f"[{sample_id}] wrote {output_path.name}")

    write_manifest(output_dir / "manifest.csv", manifest_rows)

    if args.metadata:
        metadata_path = Path(args.metadata).resolve()
        md_out = output_dir / metadata_path.name
        filter_metadata(
            metadata_path=metadata_path,
            output_path=md_out,
            id_col=args.metadata_id_col,
            keep_ids=set(matched.keys()),
        )
        print(f"Wrote filtered metadata: {md_out}")

    readme_path = output_dir / "README_minidata.md"
    with open(readme_path, "w", encoding="utf-8") as handle:
        handle.write(
            "# Mini cytosine report dataset\n\n"
            "This directory was generated for lightweight local testing of comethyl.\n\n"
            f"- Number of requested samples: {len(samples)}\n"
            f"- Number of matched samples: {len(matched)}\n"
            f"- Mode: {'copy_full_files' if args.copy_full_files else 'filtered'}\n"
            f"- Chromosomes kept: {', '.join(sorted(keep_chr)) if keep_chr else 'all'}\n"
            f"- BED used: {args.bed if args.bed else 'none'}\n"
            f"- Max lines per file: {args.max_lines_per_file if args.max_lines_per_file is not None else 'none'}\n"
        )
    print(f"Wrote README: {readme_path}")
    print(f"Wrote manifest: {output_dir / 'manifest.csv'}")


if __name__ == "__main__":
    main()
