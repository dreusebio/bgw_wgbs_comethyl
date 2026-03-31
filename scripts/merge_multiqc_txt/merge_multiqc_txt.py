#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Safely merge MultiQC .txt files from multiple multiqc_data folders."
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Paths to multiqc_data folders from different runs."
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        required=False,
        help="Optional labels for each input folder. Must match number of inputs."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output directory for merged files."
    )
    parser.add_argument(
        "--sample-column-candidates",
        nargs="+",
        default=["Sample", "sample", "sample_name", "Sample Name", "id", "ID"],
        help="Candidate column names to treat as the sample identifier."
    )
    return parser.parse_args()


def get_label_from_path(path_str):
    path = Path(path_str).resolve()
    if path.name == "multiqc_data":
        return path.parent.name
    return path.name


def detect_sample_column(df, candidates):
    for col in candidates:
        if col in df.columns:
            return col
    return None


def read_table_safely(filepath):
    try:
        df = pd.read_csv(filepath, sep="\t", dtype=str)
        return df, None
    except Exception as e:
        return None, str(e)


def main():
    args = parse_args()

    input_dirs = [Path(x).resolve() for x in args.inputs]
    for d in input_dirs:
        if not d.exists():
            print(f"ERROR: input folder does not exist: {d}", file=sys.stderr)
            sys.exit(1)
        if not d.is_dir():
            print(f"ERROR: input path is not a directory: {d}", file=sys.stderr)
            sys.exit(1)

    if args.labels:
        if len(args.labels) != len(input_dirs):
            print("ERROR: number of --labels must equal number of --inputs", file=sys.stderr)
            sys.exit(1)
        labels = args.labels
    else:
        labels = [get_label_from_path(x) for x in input_dirs]

    if len(set(labels)) != len(labels):
        print("ERROR: labels are not unique. Please provide unique --labels.", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    files_by_name = {}
    for label, folder in zip(labels, input_dirs):
        txt_files = sorted(folder.glob("*.txt"))
        for f in txt_files:
            files_by_name.setdefault(f.name, []).append((label, f))

    summary_rows = []

    for fname, entries in sorted(files_by_name.items()):
        print(f"\nProcessing: {fname}")

        if len(entries) < 2:
            print(f"  Skipping {fname}: present in only one run.")
            summary_rows.append({
                "file": fname,
                "status": "skipped",
                "reason": "present_in_only_one_run",
                "runs_present": ",".join([x[0] for x in entries]),
                "rows_written": 0,
                "duplicates_removed": 0
            })
            continue

        tables = []
        column_reference = None
        compatible = True
        reason = ""
        total_duplicates_removed = 0
        seen_samples = set()

        for label, path in entries:
            df, err = read_table_safely(path)
            if err is not None:
                compatible = False
                reason = f"read_error:{label}:{err}"
                print(f"  WARNING: failed to read {path}: {err}")
                break

            if df.shape[1] < 2:
                compatible = False
                reason = f"not_tabular:{label}"
                print(f"  WARNING: {path} does not appear to be a useful tabular file.")
                break

            if column_reference is None:
                column_reference = list(df.columns)
            else:
                if list(df.columns) != column_reference:
                    compatible = False
                    reason = f"column_mismatch:{label}"
                    print(f"  WARNING: column mismatch in {path}")
                    print(f"    Expected: {column_reference}")
                    print(f"    Found:    {list(df.columns)}")
                    break

            sample_col = detect_sample_column(df, args.sample_column_candidates)

            duplicates_removed_this_file = 0
            if sample_col is not None:
                df[sample_col] = df[sample_col].astype(str)

                before = df.shape[0]
                df = df[~df[sample_col].isin(seen_samples)].copy()
                after = df.shape[0]

                duplicates_removed_this_file = before - after
                total_duplicates_removed += duplicates_removed_this_file

                seen_samples.update(df[sample_col])

                print(
                    f"  {label}: kept {after} rows, removed {duplicates_removed_this_file} duplicate samples"
                )
            else:
                print(
                    f"  NOTE: no sample column detected in {path.name}; "
                    f"cannot deduplicate samples for {label}."
                )

            df["multiqc_run_label"] = label
            tables.append(df)

        if not compatible:
            summary_rows.append({
                "file": fname,
                "status": "skipped",
                "reason": reason,
                "runs_present": ",".join([x[0] for x in entries]),
                "rows_written": 0,
                "duplicates_removed": 0
            })
            continue

        merged = pd.concat(tables, axis=0, ignore_index=True)
        out_file = output_dir / fname
        merged.to_csv(out_file, sep="\t", index=False)

        print(f"  Wrote merged file: {out_file}")
        print(f"  Rows written: {merged.shape[0]}")
        print(f"  Total duplicates removed: {total_duplicates_removed}")

        summary_rows.append({
            "file": fname,
            "status": "merged",
            "reason": "",
            "runs_present": ",".join([x[0] for x in entries]),
            "rows_written": merged.shape[0],
            "duplicates_removed": total_duplicates_removed
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_dir / "merge_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    print("\nDone.")
    print(f"Merged files written to: {output_dir}")
    print(f"Summary written to: {summary_path}")


if __name__ == "__main__":
    main()