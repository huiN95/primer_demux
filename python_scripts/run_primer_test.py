import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Sequence


def run(cmd: Sequence[str]) -> None:
    """Run a command and stream output; raise if non-zero."""
    print("➤", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_capture(cmd: Sequence[str]) -> str:
    """Run command and capture stdout + stderr."""
    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return (p.stdout or "") + "\n" + (p.stderr or "")


def detect_cli_flags(docker_image: str) -> Dict[str, str]:
    """
    Detect which CLI flag style the image supports.
    Supports both --arg-name and --arg_name due to recent changes.
    """
    help_cmd = [
        "docker",
        "run",
        "--rm",
        docker_image,
        "primer_demux",
        "--help",
    ]

    help_text = run_capture(help_cmd)

    # Prefer hyphenated style which is the new standard
    flags = {
        "input_file": "--input-file" if "--input-file" in help_text else "--input_file",
        "output_folder": "--output-folder" if "--output-folder" in help_text else "--output_folder",
        "log_folder": "--log-folder" if "--log-folder" in help_text else "--log_folder",
        "max_distance": "--max-distance" if "--max-distance" in help_text else "--max_distance",
        "pipeline_version": "--pipeline-version" if "--pipeline-version" in help_text else "--pipeline_version",
        "output_format": "--output-format" if "--output-format" in help_text else "--output_format",
        "primer": "--primer" if "--primer" in help_text else "--primer",
    }

    print(f"Detected CLI flags for image {docker_image}:")
    for key, value in flags.items():
        print(f"  {key} = {value}")

    return flags


def docker_primer_demux(
    docker_image: str,
    infile: Path,
    output_dir: Path,
    primer: Path,
    max_distance: int,
    pipeline_version: int,
    output_format: str,
) -> None:
    """
    Run primer_demux in docker.
    """
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    flags = detect_cli_flags(docker_image)

    cmd = [
        "docker",
        "run",
        "--rm",
        "-v", "/data:/data:ro",
        "-v", f"{str(output_dir)}:{str(output_dir)}",
        "--user", f"{os.getuid()}:{os.getgid()}",
        docker_image,
        "primer_demux",
        "-i", str(infile),
        "-o", str(output_dir),
        flags["primer"], str(primer),
        flags["log_folder"], str(output_dir),
        flags["max_distance"], str(max_distance),
        flags["pipeline_version"], str(pipeline_version),
        flags["output_format"], output_format,
    ]

    run(cmd)


def mode_smoke_test(
    test_file_dir: Path,
    output_root: Path,
    docker_image: str,
) -> None:
    """
    Basic smoke test.
    """
    # 1. Discover input files (BAM or FASTQ)
    test_files = sorted(test_file_dir.glob("*.bam"))
    if not test_files:
        # Include .fastq, .fastq.gz, .fq, .fq.gz
        test_files = sorted(list(test_file_dir.glob("*.fastq*")) + list(test_file_dir.glob("*.fq*")))

    # 2. Discover primer file
    primer_file = test_file_dir / "primer40_v2.fasta"
    if not primer_file.exists():
        primer_file = test_file_dir / "primers.fasta"
    
    if not primer_file.exists():
        primer_file = next(test_file_dir.glob("*primer*.fasta"), None)

    if not test_files:
        raise FileNotFoundError(f"No test files (*.bam, *.fastq, *.fq) found under {test_file_dir}")

    if not primer_file or not primer_file.exists():
        raise FileNotFoundError(f"Primer pattern file not found in {test_file_dir}")

    output_dir = output_root / "smoke"

    for file in test_files:
        print(f"\n➡️  Running primer_demux smoke test on {file}")
        sample_output_dir = output_dir / file.stem

        docker_primer_demux(
            docker_image=docker_image,
            infile=file,
            output_dir=sample_output_dir,
            primer=primer_file,
            max_distance=1,
            pipeline_version=1,
            output_format="fasta",
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run primer_demux smoke tests"
    )

    parser.add_argument(
        "--test-file-dir",
        required=True,
        type=Path,
        help="Directory containing test files and primers.fasta",
    )

    parser.add_argument(
        "--output-root",
        required=True,
        type=Path,
        help="Writable root directory for test outputs",
    )

    parser.add_argument(
        "--docker-image",
        required=True,
        help="Current docker image to test",
    )

    args = parser.parse_args()

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    try:
        mode_smoke_test(
            test_file_dir=args.test_file_dir,
            output_root=output_root,
            docker_image=args.docker_image,
        )
        print("\n✅ Smoke test completed successfully!")

    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
