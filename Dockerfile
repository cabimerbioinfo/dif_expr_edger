# Use Rocker 4.4 as the base image
FROM rocker/r-ver:4.4.0

# Set working directory
WORKDIR /workspace

# Install necessary packages for Rsubread and edgeR
RUN apt-get update && \
    apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    liblzma-dev \
    libbz2-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"

#Install biomaRt dependencies
RUN R -e "install.packages(c('httr', 'xml2', 'jsonlite', 'methods'), repos='http://cran.rstudio.com/')"

# Install Rsubread and edgeR using BiocManager
RUN R -e "BiocManager::install(c('Rsubread', 'edgeR', 'csaw', 'dplyr', 'biomaRt', 'data.table'))"

# Clean up apt-get cache
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy the differential expression analysis scripts using edgeR to the working directory
COPY difexpr.R /scripts/

# Set the default command to execute the R script using ENTRYPOINT
ENTRYPOINT ["Rscript"]
