
# clapfcstimulation

## Overview

This repository contains analysis code developed by **H. Atilgan** for studying **mPFC neuronal activity** during **sensory and/or photostimulation**. The analysis focuses on calcium imaging datasets, enabling researchers to evaluate functional responses of cells in relation to stimulation parameters.

---

## 1. System Requirements

### Operating Systems
- macOS (tested on Monterey 12.6.1)
- Windows 10 (should be fine, limited testing)

### Required Software and Packages
- **MATLAB** (R2021a or later)
- **R**
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Image Processing Toolbox
- Parallel Computing Toolbox (for faster processing)

### Hardware Requirements
- Standard desktop or laptop with:
  - At least 16 GB RAM
  - Quad-core processor
- GPU is not required but can speed up certain preprocessing steps if MATLAB GPU support is enabled.

---

## 2. Installation Guide

### Steps
1. Clone the repository:
   ```bash
   git clone https://github.com/huriyeatg/clapfcstimulation.git
   ```
2. Open MATLAB and navigate to the repository directory.
3. Add all subfolders to the MATLAB path:
   ```matlab
   addpath(genpath(pwd));
   savepath;
   ```

### Typical Install Time
- ~2–5 minutes on a standard desktop computer

---

## 3. Demo

Demo is not applicable.

---

## 4. Instructions for Use

### Reproduction Instructions
To reproduce the figures in the upcoming publication:
1. Update the path in main_funcs.py
2. Use the 'clapfcstimulation_master.ipynb' scripts after running the pipeline.
3. Each figure script is labeled to correspond with specific figure panels (e.g., `fig1_activationAffect.ipynb`).

---

## Contact

For questions, bug reports, or contributions, please contact **Huriye Atilgan** via [GitHub](https://github.com/huriyeatg) or huriye.atilgan@gmail.com or submit an issue through this repository.

---

## License

This project is released under the MIT License.
