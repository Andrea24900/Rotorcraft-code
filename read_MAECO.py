import argparse
from pathlib import Path
import re

import numpy as np

# define regexp used to parse the file
_SECTION_RE = re.compile(r"^# Section:\s*(.+)")
_VARIABLE_RE = re.compile(r'^# Variable:\s*name="([^"]+)"\s*unit="([^"]*)"\s*rank=(\d+)')


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


def save_periodic_linearization_numpy(filepath, output_path=None):
    """Save extracted arrays to a NumPy .npz archive."""
    arrays = parse_periodic_linearization_file(filepath)
    if output_path is None:
        output_path = Path(filepath).with_suffix(".npz")
    else:
        output_path = Path(output_path)

    np.savez(
        output_path,
        time=arrays["Time_Instants"],
        A=arrays["Linearization_A_Matrix"],
        Phi=arrays["State_transition_Phi(t)_matrix"],
        M=arrays["Monodromy_M_matrix"],
    )
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse a MAECOsim periodic_linearization output and save time, A(t), Phi(t), and M as numpy arrays."
    )
    parser.add_argument("input", nargs="?", default="periodic_linearization6.txt", help="Input text file path")
    parser.add_argument("-o", "--output", default=None, help="Output .npz file path")
    args = parser.parse_args()

    output_file = save_periodic_linearization_numpy(args.input, args.output)
    print(f"Saved NumPy arrays to {output_file}")

    data = np.load(output_file)
    print("time:", data["time"])
    print("A:")
    print(data["A"])
    print("Phi:")
    print(data["Phi"])
    print("M:", data["M"])
