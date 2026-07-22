import argparse
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

# define regexp used to parse the file
_SECTION_RE = re.compile(r"^# Section:\s*(.+)")


@dataclass
class PeriodicLinearizationData:
    time: np.ndarray
    A: np.ndarray
    Phi: np.ndarray
    M: np.ndarray


def parse_periodic_linearization_file(filepath):
    """Parse a MAECOsim periodic_linearization text file into NumPy arrays."""
    filepath = Path(filepath)
    wanted_sections = {
        "Time_Instants",
        "Linearization_A_Matrix",
        "Monodromy_M_matrix",
        "State_transition_Phi(t)_matrix",
    }
    sections = {}
    current_section = None

    with filepath.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue

            section_match = _SECTION_RE.match(line)
            if section_match:
                section = section_match.group(1)
                if section in wanted_sections:
                    current_section = section
                    sections[current_section] = {
                        "data": [],
                    }
                else:
                    current_section = None
                continue

            if current_section is None:
                continue

            if line.startswith("#"):
                continue

            values = [float(token) for token in line.split()]
            sections[current_section]["data"].append(values)

    arrays = {}
    for section, info in sections.items():
        data = np.asarray(info["data"], dtype=float)
        if section == "Time_Instants" and data.ndim == 2 and data.shape[1] == 1:
            data = data.ravel()
        arrays[section] = data

    return arrays


def read_MAECO(filepath) -> PeriodicLinearizationData:
    """Read a MAECOsim periodic_linearization text file and return a data struct."""
    arrays = parse_periodic_linearization_file(filepath)
    return PeriodicLinearizationData(
        time=arrays["Time_Instants"],
        A=arrays["Linearization_A_Matrix"],
        Phi=arrays["State_transition_Phi(t)_matrix"],
        M=arrays["Monodromy_M_matrix"],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse a MAECOsim periodic_linearization output and print time, A(t), Phi(t), and M."
    )
    parser.add_argument("input", nargs="?", default="periodic_linearization6.txt", help="Input text file path")
    args = parser.parse_args()

    data = read_MAECO(args.input)
    print("time:", data.time)
    print("A:")
    print(data.A)
    print("Phi:")
    print(data.Phi)
    print("M:")
    print(data.M)
