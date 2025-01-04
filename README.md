# Single Nuclei Transcription Factor Prediction and Validation

This repository contains the code and data necessary for single nuclei transcription factor prediction and validation using Seurat and Harmony in R.

## Project Structure

- **PTPRF_SingleCell_Validation_Paper_Figure.R**: Main analysis script for performing Seurat-based single-cell analysis.
- **Dockerfile**: Docker container setup for reproducibility.
- **annot_humanAll.csv**: Annotation file used for cell metadata integration.
- **repo/**: Contains essential files and scripts for the project.

## Requirements

- R (4.0 or higher)
- Seurat
- ggplot2
- dplyr
- harmony

## Usage

To run the analysis:

```bash
# Build the Docker container
docker build -t single_nuclei_tf .

# Run the analysis inside the container
docker run --rm -v $(pwd):/workspace single_nuclei_tf
```

## License

This project is licensed under the MIT License.

## Author

Carlos Buss
