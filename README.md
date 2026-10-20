# Shark and ray ventilation batch analysis

## What it does

This local desktop application recursively finds videos under a Dropbox folder, extracts a motion signal from species specific respiratory regions of interest, estimates ventilation rate using both a power spectrum and NeuroKit2 peak detection, and writes one CSV row per video.

The output CSV includes the root folder, relative folder, immediate subfolder, file name, species, condition, treatment, video quality fields, FFT estimate, peak estimate, final estimate, and error status.

Metadata are inferred from folder names or filenames using these keywords:

- `epaulette_shark` (also recognises `epaulette`, `Hemiscyllium`, and `ocellatum`)
- `shovelnose_ray` (also recognises `shovelnose`, `Aptychotrema`, and `rostrata`)
- `hypoxia` or `normoxia`
- `PBS` or `Poly`

Use a consistent directory layout such as:

```text
Dropbox
└── VentilationVideos
    ├── epaulette_shark
    │   ├── normoxia
    │   │   ├── PBS
    │   │   └── Poly
    │   └── hypoxia
    │       ├── PBS
    │       └── Poly
    └── shovelnose_ray
        └── ...
```

The folders may be in another order. The program searches the whole relative path and the filename.

## Install on your computer

1. Install Python 3.11 or later from https://www.python.org/downloads/.
2. Open Terminal on macOS or Command Prompt on Windows and change into this folder.
3. Create a clean environment rather than using a broken Conda base environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python ventilation_batch_gui.py
```

On Windows, activate the environment with `.venv\\Scripts\\activate`.

NeuroKit2 is used for the contraction peak finding stage. It requires SciPy, so a clean virtual environment is important if an old Conda installation has a damaged SciPy binary.

## Use

1. Click **Browse** beside **Video root folder** and choose the Dropbox folder containing all video subfolders.
2. Click **Choose ROI** for **epaulette shark** and draw a box around a clearly visible buccal or gill movement region. Click **Choose ROI** for **shovelnose ray** and draw a box around a clearly visible spiracle or branchial movement region. Use comparable framing within each species.
3. Set the plausible ventilation range. Start with 20 to 180 breaths per minute, then narrow it after visually checking diagnostic plots and manual counts.
4. Keep **Save one diagnostic plot per video** selected for the first run.
5. Click **Analyse all videos**.
6. Review the CSV and diagnostic plots before using the estimates for analysis.

## Important validation

The tool measures periodic local motion, not ventilation directly. Validate the first batch against blinded manual counts for epaulette sharks and shovelnose rays separately, within each condition, camera angle, and recording setup. Do not pool files with different body orientation or poorly visible respiratory movement without separate validation.
