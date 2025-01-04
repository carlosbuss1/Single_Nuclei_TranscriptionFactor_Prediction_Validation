FROM rocker/r-ver:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('Seurat', 'ggplot2', 'dplyr', 'harmony', 'BiocManager'))"
RUN R -e "BiocManager::install(c('SeuratData'))"

# Set working directory and copy project files
WORKDIR /home/ulb/Desktop/carlos/Human_Liver_Atlas
COPY ./rawData_human /home/ulb/Desktop/carlos/Human_Liver_Atlas/rawData_human
COPY ./output /home/ulb/Desktop/carlos/Human_Liver_Atlas/output
COPY PTPRF_SingleCell_Validation_Paper_Figure.R /home/ulb/Desktop/carlos/Human_Liver_Atlas/

# Set default command to run the script
CMD ["Rscript", "/home/ulb/Desktop/carlos/Human_Liver_Atlas/PTPRF_SingleCell_Validation_Paper_Figure.R"]

