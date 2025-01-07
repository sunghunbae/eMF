"""
Author: Jessica M. Gonzalez-Delgado
        North Carolina State University

This script creates an input file for the eMF software from individual
relaxation data files.

Run script in the work directory where all the relaxation data files are
present.

The relaxation data files have to have the following format:
    Line 1: Type of relaxation data (values = noe, r1, r2)
    Line 2: Magnetic field (MHz)
    Lines 3+: residue_number value error cluster_number

For the distributed model, the distribution file name must be manually
added after creating the input file.

The relaxation data files must be named as:
    *-r1.dat, *-r2.dat, *-noe.dat
"""

import os
import sys
import textwrap
from typing import List, Tuple


def get_data_files() -> List[str]:
    """Get names of relaxation data files in current directory.

    Returns
    ----------
    data_files : List[str]
        List of relaxation data files.

    Raises
    ----------
    ValueError : If not enough datafiles to run eMF.
    """
    valid_suffixes = ("-r1.dat", "-r2.dat", "-noe.dat")
    data_files = [file for file in os.listdir(
        os.getcwd()) if file.endswith(valid_suffixes)]

    try:
        if len(data_files) < 6:
            raise ValueError("Not enough datafiles to run eMF. Must have NOE, R1, and R2 data in at least 2 fields.")

    except ValueError as e:
        sys.exit(f"Error: {e}")

    return data_files


def get_experiment_field(experiment: str, field: float,
                         line: List[str]) -> Tuple[str, float]:
    """Process line with either experiment name or field.

    Parameters
    ----------
    experiment : str
        Name of experiment (noe, r1, r2).
    field : float
        Experiment field strength (MHz).
    line : str
        Line to process

    Returns
    ----------
    experiment : str
        Name of experiment (noe, r1, r2).
    field : float
        Experiment field strength (MHz).
    """
    if line[0].lower() in {'noe', 'r1', 'r2'}:
        experiment = line[0].lower()

    else:
        try:
            field = float(line[0])
        except ValueError:
            pass

    return experiment, field


def fill_data_arrays(
        noe: List[Tuple[int, int, float, float, float]],
        r1: List[Tuple[int, int, float, float, float]],
        r2: List[Tuple[int, int, float, float, float]],
        experiment: str,
        field: float,
        line: List[str]
) -> None:
    """Fill corresponding experiment data array with line.

    Parameters
    ----------
    noe : List[Tuple[int, int, float, float, float]]
        List with data from NOE experiments.
    r1 : List[Tuple[int, int, float, float, float]]
        List with data from R1 experiments.
    r2 : List[Tuple[int, int, float, float, float]]
        List with data from R2 experiments.
    experiment : str
        Name of experiment (noe, r1, r2).
    field : float
        Experiment field strength (MHz).
    line : list (str)
        Line to process

    Returns
    ----------
    None

    Raises
    ----------
    UserWarning : If NOE value is > 1.0 it gets automatically adjusted
    to 1.0
    """
    spin = int(line.pop(0))
    cluster = int(line.pop())
    value = float(line.pop(0))
    error = float(line[0])

    if experiment == 'noe':
        if value > 1.0:
            print(f"Warning: NOE value > 1.000 at spin {spin}, adjusted to 1.000")
            value = 1.000
        noe.append((spin, cluster, field, value, error))
    elif experiment == 'r1':
        r1.append((spin, cluster, field, value, error))
    elif experiment == 'r2':
        r2.append((spin, cluster, field, value, error))

    return None


def get_data(data_files: List[str]) -> Tuple[
        List[Tuple[int, int, float, float, float]],
        List[Tuple[int, int, float, float, float]],
        List[Tuple[int, int, float, float, float]]
]:
    """Get relaxation data from data files.

    Parameters
    ----------
    data_files : list (str)
        List of relaxation data files.

    Returns
    ----------
    noe : List[Tuple[int, int, float, float, float]]
        List with data from NOE experiments.
    r1 : List[Tuple[int, int, float, float, float]]
        List with data from R1 experiments.
    r2 : List[Tuple[int, int, float, float, float]]
        List with data from R2 experiments.

    Raises
    ----------
    ValueError : If data file has incorrect format.
    """
    noe, r1, r2 = [], [], []

    for data_file in data_files:
        with open(data_file) as file:
            experiment, field = None, None

            for line in file:
                line = line.strip().split()

                try:
                    if line[0].startswith("#"):
                        continue

                    if len(line) not in {1, 4}:
                        raise ValueError(
                            f"{data_file} does not have correct format. Verify that lines are the correct length.")

                    if len(line) == 1:
                        experiment, field = get_experiment_field(experiment, field, line)

                    elif len(line) == 4 and experiment is not None and field is not None:
                        fill_data_arrays(noe, r1, r2, experiment, field, line)

                except ValueError as e:
                    sys.exit(f"Error: {e}")

    return noe, r1, r2


def sort_validate_data(
        noe: List[Tuple[int, int, float, float, float]],
        r1: List[Tuple[int, int, float, float, float]],
        r2: List[Tuple[int, int, float, float, float]]
) -> List[str]:
    """Sort and validate data arrays by veryfing that their length is
    equal and that they contain data for the same spins.

    Parameters
    ----------
    noe : List[Tuple[int, int, float, float, float]]
        List with data from NOE experiments.
    r1 : List[Tuple[int, int, float, float, float]]
        List with data from R1 experiments.
    r2 : List[Tuple[int, int, float, float, float]]
        List with data from R2 experiments.

    Returns
    ----------
    validated_data : List[str]
        List of validated data to be written to input file.

    Raises
    ----------
    ValueError : Raised if data arrays file have incorrect length, they
        do not contain data for the same spins, or the clusters are
        different.
    """
    validated_data = []
    noe.sort(key=lambda x: x[0])
    r1.sort(key=lambda x: x[0])
    r2.sort(key=lambda x: x[0])

    try:
        if not (len(noe) == len(r1) == len(r2)):
            raise ValueError("Amount of spins in NOE, R1, R2 data files do not match.")

        for noe_spin, r1_spin, r2_spin in zip(noe, r1, r2):
            if not (noe_spin[0] == r1_spin[0] == r2_spin[0]):
                raise ValueError("Data files have the same amount of spins but there is a spin mismatch between data files.")

            if not (noe_spin[1] == r1_spin[1] == r2_spin[1]):
                raise ValueError("There is a cluster mismatch between data files.")

            entry = (
                f"{r1_spin[0]:<6d} {r1_spin[1]:<4d} {r1_spin[2]:>6.2f} "
                f"{r1_spin[3]:>10.6f} {r1_spin[4]:>10.6f} "
                f"{r2_spin[3]:>10.6f} {r2_spin[4]:>10.6f} "
                f"{noe_spin[3]:>6.3f} {noe_spin[4]:>6.3f}"
            )

            validated_data.append(entry)

    except ValueError as e:
        sys.exit(f"Error: {e}")

    return validated_data


def write_input_file(validated_data: List[str]) -> None:
    """Write eMF input file.

    Parameters
    ----------
    validated_data : List[str]
        List of validated data to be written to input file.

    Returns
    ----------
    None
    """
    input_file = "input.conf"
    header = textwrap.dedent("""\
        # column 1   : residue number
        # column 2   : cluster number, used for inclusion or exclusion of data
        # column 3   : B0, field strength in MHz
        # column 4-5 : R1 relaxation rate and error (1/s)
        # column 6-7 : R2 relaxation rate and error (1/s)
        # column 8-9 : {1H}-X NOE and error
        # column 10  : (only for distributed) tc distribution file name\n
    """)

    with open(input_file, 'w') as file:
        file.write(header)
        file.write('\n'.join(validated_data))

    return None


def main() -> None:
    """Main program."""
    data_files = get_data_files()
    noe, r1, r2 = get_data(data_files)
    validated_data = sort_validate_data(noe, r1, r2)
    write_input_file(validated_data)

    return None

if __name__ == "__main__":
    main()
