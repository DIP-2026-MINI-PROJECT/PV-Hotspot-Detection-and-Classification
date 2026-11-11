# PV Hotspot Detection and Analysis

**Project:** Detection and Classification of Hot Spots in Photovoltaic Panels Using Thermal Image Processing Techniques  
**Course:** UE24EC342AC1 – Digital Image Processing  
**Environment:** MATLAB with Image Processing Toolbox

## Overview

This project applies digital image processing techniques to photovoltaic (PV) panel images to identify and analyse potential abnormal intensity regions. The MATLAB implementation covers image fundamentals, transform-domain processing, enhancement and restoration, morphological operations, and colour processing.

The candidate-region analysis uses 2D Discrete Cosine Transform (DCT)-based low-frequency background estimation, residual calculation, adaptive thresholding, panel-region masking, and morphological refinement. Region properties such as area, centroid, and mean residual intensity are extracted for visualization and further inspection.

**Important:** The current implementation identifies candidate intensity deviations. It does not establish that a region is a confirmed PV fault or perform a trained hotspot classification. Pixel intensities are not calibrated temperature measurements, and pseudocolour is an artificial visualization palette.

## Repository contents

- `PV.m` – main MATLAB script covering all four DIP units.
- `pv_dataset_manifest.m` – discovers and organizes dataset images by date and development/evaluation split.
- `pv_validation_summary.m` – summarizes PSNR measurements by experiment groups and collection dates.
- `PV_images_failures.zip` – supplied PV image dataset archive (extract before running).

## Requirements

- MATLAB
- Image Processing Toolbox
- Dataset archive extracted into a folder named `PV_images_failures`

## Project structure

After extracting the dataset, the project should be arranged as follows:

```text
PV-Hotspot-Detection-and-Classification/
├── PV.m
├── pv_dataset_manifest.m
├── pv_validation_summary.m
├── PV_images_failures.zip
└── PV_images_failures/
    ├── 2018_April_27/
    ├── 2018_April_28/
    ├── 2018_May_4/
    ├── 2018_December_20/
    ├── 2018_December_21/
    ├── 2019_January_16/
    └── 2019_January_19/
```

## Running the project

1. Clone or download this repository.
2. Extract `PV_images_failures.zip` into the project root. The extracted folder must be named `PV_images_failures`.
3. Open the project root as the MATLAB Current Folder.
4. Ensure the three `.m` files are together in the project root.
5. Run:

```matlab
PV
```

The script requires the Image Processing Toolbox. It clears the MATLAB workspace, closes open figures, and clears the Command Window at startup.

## Processing workflow

1. Load PV images from the dataset.
2. Convert images to grayscale where required, normalize intensities, and resize to 128 × 128.
3. Apply the 2D DCT and evaluate reconstruction at selected coefficient budgets.
4. Estimate the slowly varying background using low-frequency DCT components.
5. Calculate the residual between the image and estimated background.
6. Apply adaptive thresholding and restrict candidate regions using the PV panel ROI.
7. Refine the binary mask using morphological operations and small-region removal.
8. Extract region measurements and visualize candidate contours.
9. Demonstrate enhancement, noise modelling, restoration, morphology, and colour processing across the four DIP units.

## Outputs

The script creates four unit-wise result figures and prints summary information in the MATLAB Command Window. It also creates result tables in the MATLAB workspace.

The setting `saveResults` is `false` by default. To export the available CSV tables and MAT results, change it to `true` near the top of `PV.m`. The output folder `PV_DCT_results` will then be created by the script.

## Limitations

- Candidate regions are based on image-intensity deviations and are not verified hotspot/fault labels.
- The image intensities are not calibrated temperatures.
- Pseudocolour is used for visualization and does not provide thermal calibration.
- Detection precision, recall, classification accuracy, and hotspot severity cannot be claimed without suitable labelled ground truth and calibrated measurements.
- Results depend on the supplied dataset and have not been established for all PV installations or imaging conditions.

## References

1. Gonzalez, R. C., and Woods, R. E., *Digital Image Processing*, 4th Edition, Pearson, 2018.
2. MathWorks, MATLAB Image Processing Toolbox documentation: https://www.mathworks.com/help/images/
3. MathWorks, `dct2` documentation: https://www.mathworks.com/help/images/ref/dct2.html
4. Altawallbeh, Wejdan, et al. "Detection and Classification of Hot Spots in Photovoltaic Panels Using Thermal Image Processing Techniques." Signals 7.4 (2026): 61.
## Execution status

The source files and project documentation are prepared for repository submission. Successful execution on the target MATLAB installation has **not yet been verified**. Please report any environment-specific dependency or runtime issue before treating the outputs as reproduced.
