import os
import sys
from typing import Generator, Iterator, Optional, Tuple

import click


def parse_step_input(
    ctx: click.Context, param: click.Parameter, value: str
) -> Optional[int]:
    if not value:
        return None
    try:
        return int(float(value))
    except ValueError:
        raise click.BadParameter(f"The value '{value}' is not a valid number.")


def get_step_from_line(line: str) -> Optional[int]:
    try:
        parts = line.split("=")
        if len(parts) < 2:
            return None
        return int(parts[1].strip())
    except (IndexError, ValueError):
        return None


def iter_steps_with_offset(input_file: str) -> Generator[Tuple[int, int], None, None]:
    """
    Reads the file line by line and yields (step, offset_in_byte).
    Uses readline() + tell() to get the precise position of the start of the line.
    """
    with open(input_file, "r") as f:
        while True:
            # Save the offset BEFORE reading the line
            offset = f.tell()
            line = f.readline()
            if not line:
                break

            if line.startswith("t = "):
                step = get_step_from_line(line)
                if step is not None:
                    yield step, offset


def pick_nearest_offset(
    before: Optional[Tuple[int, int]],
    after: Optional[Tuple[int, int]],
    target_step: int,
) -> Tuple[int, int]:
    """
    Chooses between the 'before' and 'after' candidates based on the distance from the target.
    Returns the winning (step, offset) tuple.
    """
    if before is None and after is None:
        raise ValueError("File appears empty or has no steps.")
    if before is None:
        return after  # type: ignore
    if after is None:
        return before

    step_before, _ = before
    step_after, _ = after

    if (target_step - step_before) <= (step_after - target_step):
        return before
    return after


def scan_nearest(input_file: str, target_step: int) -> Tuple[Optional[int], int, bool]:
    """
    Intelligent scan.
    Returns: (found_step, byte_offset, exact_match_bool)
    Offset is -1 if not found.
    Assumes steps are monotonically increasing in the file.
    """
    before: Optional[Tuple[int, int]] = None

    for step, offset in iter_steps_with_offset(input_file):
        current = (step, offset)

        if step == target_step:
            return step, offset, True

        if step > target_step:
            # We have passed the target. The closest is either 'before' or 'current'.
            best_step, best_offset = pick_nearest_offset(before, current, target_step)
            return best_step, best_offset, False

        before = current

    # End of file without passing the target (target > last step)
    if before is not None:
        return before[0], before[1], False

    return None, -1, False


def scan_steps(
    input_file: str, target_step: int, nearest: bool
) -> Tuple[int, int, bool]:
    """
    Orchestrates the scan.
    Returns (step_to_extract, byte_offset, exact_match).
    For nearest lookups, steps are assumed to be monotonically increasing.
    """
    if nearest:
        step, offset, exact = scan_nearest(input_file, target_step)
        if step is None:
            # Empty or invalid file case
            return target_step, -1, False
        return step, offset, exact

    # Standard scan (optimized early exit)
    for step, offset in iter_steps_with_offset(input_file):
        if step == target_step:
            return step, offset, True
        if step > target_step:
            break

    return target_step, -1, False


def extract_block_fast(input_file: str, output: str, start_offset: int) -> bool:
    """
    Extracts the block by jumping directly to the indicated offset.
    Does not search, just writes until it finds the next 't = '.
    """
    if start_offset < 0:
        return False

    with open(input_file, "r") as infile, open(output, "w") as outfile:
        # HYPERSPACE JUMP: Go straight to the point
        infile.seek(start_offset)

        # Read the first line (which we know is 't = ...')
        header = infile.readline()
        if not header:
            return False
        outfile.write(header)

        # Write everything else until a new block starts
        for line in infile:
            if line.startswith("t = "):
                break
            outfile.write(line)

    return True


@click.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.argument("target_step", callback=parse_step_input)
@click.option(
    "-o", "--output", type=click.Path(), help="Output file. Default: conf_{step}.dat"
)
@click.option(
    "--nearest", is_flag=True, help="Extract closest step if exact match is missing."
)
def extract_trajectory(
    input_file: str, target_step: int, output: Optional[str], nearest: bool
):
    """
    Extracts a configuration from an oxDNA trajectory.
    Optimized with SEEK/TELL for maximum performance on large files.
    """

    # PHASE 1: Scanning (find ONLY the position)
    step_to_extract, offset, found_exact_match = scan_steps(
        input_file, target_step, nearest
    )

    # Error Handling / Nearest
    if nearest:
        if offset == -1:
            click.echo("Error: No valid steps found in file.", err=True)
            sys.exit(1)
        if step_to_extract != target_step:
            click.echo(f"Step {target_step} not found. Closest is {step_to_extract}.")
    elif not found_exact_match:
        click.echo(f"Error: Step {target_step} not found.", err=True)
        click.echo("Tip: use --nearest.", err=True)
        sys.exit(1)

    # PHASE 2: Output setup
    if output is None:
        output = f"conf_{step_to_extract}.dat"

    # PHASE 3: Surgical Extraction
    # Does not reread the entire file, jumps directly to the offset calculated in Phase 1
    success = extract_block_fast(input_file, output, offset)

    if success:
        click.echo(f"Done: step {step_to_extract} -> {output}")
    else:
        click.echo("Error during extraction (file truncated?).", err=True)
        if os.path.exists(output):
            os.remove(output)


if __name__ == "__main__":
    extract_trajectory()

