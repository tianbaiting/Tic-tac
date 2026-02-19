# Input Profiles

This folder is the canonical location for hand-written solver input files.

## Files

- `input.txt`: baseline profile for routine runs.
- `input_Ay_test.txt`: low-cost 190 MeV/u quick-test profile.
- `input_Ay_fixed.txt`: refined 190 MeV/u profile.
- `input_Ay_nijmegen_legacy.txt`: legacy Nijmegen profile kept for reproducibility.
- `lab_energies.txt`: baseline energy grid.
- `lab_energies_190MeV.txt`: single-point 190 MeV/u dataset used by 190 MeV workflows.

## Usage

Run from repository root:

```bash
./CPP/run CPP/Input/input.txt
```

`energy_input_file`, `output_folder`, and `P123_folder` values in these templates are root-relative.
