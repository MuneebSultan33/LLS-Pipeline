# LLS-Pipeline 🔬
### Open-Access 4D Lattice Light-Sheet Microscopy Processing Pipeline

**Built by Muhammad Muneeb Sultan — UVAS Lahore**  
**Developed in collaboration with the Cell Observatory Lab, DKFZ NCT WERA, Würzburg**

---

## What This Tool Does

Lattice Light-Sheet (LLS) microscopes like the Zeiss LLS7 produce massive 4D datasets — stacks of 3D volumes captured over time. Before these datasets can be analysed, three processing steps must happen:

1. **Deskewing** — the microscope stage moves at a 30-degree angle, making raw data appear tilted. This step mathematically straightens it.
2. **Coverslip Rotation** — after deskewing, the data is rotated 30 degrees around the X axis to convert from stage coordinates into true physical coverslip coordinates, so that distances in Z are physically accurate.
3. **Drift Correction** — during long imaging sessions the sample physically drifts by a few pixels. This step detects and corrects that drift between timepoints.

Currently, most LLS labs do these three steps using a chain of three separate programs: **ZEN software → FIJI macros → MATLAB scripts**. This requires expensive software licenses, manual steps between each program, and significant setup time for new labs.

**This pipeline replaces that entire chain with a single Python script. Free. Open-access. Zero paid licenses required.**

An additional **AI Synapse Scoring** layer identifies pixels that are consistently bright across all timepoints — marking high-confidence molecular clustering zones such as immune synapses, GPRC5D or BCMA receptor clusters.

---

## Scientific Background

This pipeline implements the processing workflow described in:

> Ghosh, A. et al. *"Decoding the molecular interplay of CD20 and therapeutic antibodies with fast volumetric nanoscopy."* **Science 387**, eadq4510 (2025). DOI: 10.1126/science.adq4510

The deskewing and coverslip rotation steps follow the methodology described in:

> Iwanski, M.K., Katrukha, E.A., Kapitein, L.C. *"Lattice light-sheet motor-PAINT."* Methods Mol. Biol. **2694**, 151–174 (2024).

---

## Pipeline Steps

```
RAW LLS DATA
     │
     ▼
[Step 2] DESKEWING
         Shear transform: each z-slice shifted in Y by z × tan(30°)
         Corrects the 30-degree stage angle of the Zeiss LLS7
     │
     ▼
[Step 3] COVERSLIP ROTATION  ← Key step
         30-degree rotation around X axis using affine transform
         Converts stage coordinates → true coverslip coordinates
         Z distances become physically accurate after this step
     │
     ▼
[Step 4] DRIFT CORRECTION
         Centre-of-mass tracking across timepoints
         Detects and corrects inter-timepoint sample movement
         Equivalent to the Correlescence FIJI plugin — in Python
     │
     ▼
[Step 5] AI SYNAPSE SCORING
         Temporal persistence metric: mean × (1 - CV)
         Identifies consistently bright pixels = real structures
         Separates true molecular clusters from random noise
     │
     ▼
OUTPUT IMAGE + SUMMARY
```

---

## Requirements

Python 3.8 or higher. No paid software needed.

Install dependencies with one command:

```bash
pip install numpy scipy matplotlib
```

---

## How to Run

**1. Clone or download this repository**

```bash
git clone https://github.com/MuneebSultan33/LLS-Pipeline.git
cd LLS-Pipeline
```

**2. Install requirements**

```bash
pip install numpy scipy matplotlib
```

**3. Run the pipeline**

```bash
python demo.py
```

**4. View the output**

The script saves `lls_pipeline_result.png` in the same folder.  
It also prints a full summary in the terminal.

---

## Using Real Data

To run the pipeline on real LLS localization data:

1. Export your localization CSV from ZEN software (the file containing X, Y, Z coordinates in nanometers)
2. Rename it `real_data.csv`
3. Place it in the same folder as `demo.py`
4. Run `python demo.py`

The pipeline automatically detects the file and loads it. If the file does not contain XYZ localization columns, it falls back to synthetic data with a clear message explaining why.

**Compatible formats:** CSV files with columns containing X, Y, Z coordinates in nanometers — as exported by ZEN Blue software from the Zeiss LLS7.

---

## Output

Running the pipeline produces:

- `lls_pipeline_result.png` — a 5-row visualization grid showing:
  - Row 1: Raw skewed data
  - Row 2: After deskewing
  - Row 3: After coverslip rotation
  - Row 4: After drift correction
  - Row 5: AI synapse score map with high-confidence regions contoured

- Terminal summary showing:
  - Processing parameters used
  - Drift detected at each timepoint
  - Number of high-confidence synapse pixels identified

---

## Current Status

**v2.0 — Proof of Concept**

- [x] Deskewing (30-degree shear transform)
- [x] Coverslip rotation (30-degree affine rotation around X axis)
- [x] Inter-timepoint drift correction (centre-of-mass tracking)
- [x] AI temporal persistence scoring
- [x] Synthetic 4D immune synapse simulation for testing
- [x] Automatic real CSV detection and loading
- [ ] TIFF stack reader for raw Zeiss LLS7 output (in development)
- [ ] Batch processing of multiple datasets (planned)
- [ ] GUI interface (planned)

---

## Why This Exists

New LLS labs face a significant data bottleneck when setting up their processing workflows. The standard pipeline requires MATLAB (expensive license), FIJI (manual macro operation), and ZEN (tied to Zeiss hardware software). A new postdoc joining a lab can spend weeks setting up this chain before processing a single dataset.

This tool aims to be a free, documented, one-command alternative that any lab can install and run immediately — regardless of budget or software availability.

---

## Author

**Muhammad Muneeb Sultan**  
BSc Biotechnology — University of Veterinary and Animal Sciences (UVAS), Lahore, Pakistan  
Email: muneebbioinformatics369@gmail.com  
GitHub: github.com/MuneebSultan33

---

## Acknowledgements

Developed with guidance from **Dr. Arindam Ghosh**, Research Group Leader, Cell Observatory Lab, DKFZ NCT WERA, Würzburg, Germany.

Pipeline methodology based on the LLS-TDI-DNA-PAINT workflow described in Ghosh et al., Science 2025.

---

## License

MIT License — free to use, modify, and distribute with attribution.
