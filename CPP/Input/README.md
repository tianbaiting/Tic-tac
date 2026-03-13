# Input Profiles

This folder is the canonical location for hand-written solver input files.

## Files

- `input.txt`: baseline profile for routine runs.
- `input_Ay_test.txt`: low-cost `Tlab = 190 MeV` quick-test profile.
- `input_Ay_fixed.txt`: refined `Tlab = 190 MeV` profile.
- `input_Ay_nijmegen_legacy.txt`: legacy Nijmegen profile kept for reproducibility.
- `lab_energies.txt`: baseline energy grid.
- `tlab_190MeV.txt`: single-point `Tlab = 190 MeV` dataset used by validation workflows.

## Usage

Run from repository root:

```bash
./CPP/run CPP/Input/input.txt
```

`energy_input_file`, `output_folder`, and `P123_folder` values in these templates are root-relative.
