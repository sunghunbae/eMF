"""
Author: Jessica M. Gonzalez-Delgado
        North Carolina State University

This script creates a file containing relaxation rates that are available
for comparison between two files, e.g. the file will contain spin #1234
only if data for spin #1234 is present in file #1, file #2 and file #3.

This is useful when you have multiple relaxation or eMF digested data
files for the spin system that might not contain data for the same spins
and you want to see the differences in the data files without having to
manually weed out the spins that are not present in both files.

The script supports:
    (1) Input format required for "emf_inputgen.py"
    (2) eMF digested data files

The output file is named 'data4comp.out' and contains:
    - Header with filenames
    - Data rows: spin-number followed by values from both files
"""

import argparse
from typing import Tuple, List


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=f"Extract common spin data for comparison between two NMR relaxation data files."
    )
    parser.add_argument(
        'file_a',
        type=str,
        help='Path to the 1st input file (.dat or .out).'
    )
    parser.add_argument(
        'file_b',
        type=str,
        help='Path to the 2nd input file (.dat or .out).'
    )

    args = parser.parse_args()

    if not args.file_a.endswith(".dat") and not args.file_a.endswith(".out"):
        parser.error("File_a must be *.DAT or *.OUT")
    if not args.file_b.endswith('.dat') and not args.file_b.endswith('.out'):
        parser.error("File_b must be *.DAT or *.OUT")

    return args


def is_emf_format(args: argparse.Namespace) -> bool:
    """Check if files are in EMF digest format.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    ----------
    is_emf : bool
        Indicates if input files are in EMF digest format.
    """
    input_files = [args.file_a, args.file_b]
    is_emf = []

    for input_file in input_files:
        with open(input_file) as file:
            for line in file:
                line = line.strip().split()
                if line[0].startswith("#") or len(line) == 1:
                    continue

                if len(line) == 4:  # NOE, R1, R2 format
                    is_emf.append(False)
                elif len(line) == 16:  # eMF digested format
                    is_emf.append(True)
                else:
                    raise ValueError(f"Error: {input_file} has incorrect format.")

    if len(set(is_emf)) != 1:
        raise ValueError(f"Both input files must have the same format, i.e. 'EMF digest' format or 'NOE R1 R2' format.")

    return is_emf[0]


def get_spin_data(args: argparse.Namespace, is_emf: bool) -> Tuple[List[str], List[str]]:
    """Get relaxation data from spins that are in both input files.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    is_emf : bool
        Indicates if input files are in EMF digest format.

    Returns
    ----------
    data_a : list (str)
        Filtered data from the 1st input file. Filtered data is the
        relaxation data of the common spins between the two input files.
    data_b : list (str)
        Filtered data from the 2nd input file. Filtered data is the
        relaxation data of the common spins between the two input files.
    """
    spins_a, data_a =  set(), []
    spins_b, data_b =  set(), []

    for file, spins, data in [
        (args.file_a, spins_a, data_a), (args.file_b, spins_b, data_b)]:
        with open(file) as fo:
            for line in fo:
                line = line.strip().split()
                if line[0].startswith("#") or len(line) == 1:
                    continue
                entry = [line[0], line[1], line[7], line[9]] if is_emf else line[:2]
                data.append(entry)
                spins.add(int(entry[0]))

    spins_common = spins_a & spins_b
    data_a = [entry for entry in data_a if int(entry[0]) in spins_common]
    data_b = [entry[1:] for entry in data_b if int(entry[0]) in spins_common]

    return data_a, data_b


def write_output_file(args: argparse.Namespace, is_emf: bool, data_a: List[str], data_b: List[str]) -> None:
    """Write the filtered data to an output file.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    is_emf : bool
        Indicates if input files are in EMF digest format.
    data_a : list (str)
        Filtered data from the 1st input file. Filtered data is the
        relaxation data of the common spins between the two input files.
    data_b : list (str)
        Filtered data from the 2nd input file. Filtered data is the
        relaxation data of the common spins between the two input files.

    Returns
    ----------
    None
    """
    output_file = "data_comparison.out"

    with open(output_file, "w") as file:
        file.write(f"# Data files: {args.file_a}\t{args.file_b}\n")
        if is_emf:
            file.write(f"# Res\tS2-1\tte-1\tRex-1\tS2-2\tte-2\tRex-2\n")
            for a, b in zip(data_a, data_b):
                file.write(f"{a[0]}\t{a[1]}\t{a[2]}\t{a[3]}\t{b[0]}\t{b[1]}\t{b[2]}\n")
        else:
            file.write("# Res\tval-1\tval-2\n")
            for a, b in zip(data_a, data_b):
                file.write(f"{a[0]}\t{a[1]}\t{b[0]}\n")

    return None


def main() -> None:
    """Main Program."""
    args = parse_arguments()
    is_emf = is_emf_format(args)
    data_a, data_b = get_spin_data(args, is_emf)
    write_output_file(args, is_emf, data_a, data_b)


if __name__ == "__main__":
    main()
