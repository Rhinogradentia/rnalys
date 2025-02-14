# Use Python base image
FROM python:3.9

# Install R and dependencies
RUN apt-get update && apt-get install -y \
    software-properties-common \
    dirmngr \
    wget

# Add R 4.2 repository
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

# Install R 4.2 and system dependencies
RUN apt-get update && apt-get install -y \
    r-base=4.2.* \
    r-base-dev=4.2.* \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    git

# Clone the repository
WORKDIR /app
RUN git clone https://github.com/Christofferfrisk/rnalys.git .

# Install Python dependencies with exact versions
RUN pip install --no-cache-dir \
    dash==2.6.2 \
    dash-auth==1.4.1 \
    dash-bio==1.0.2 \
    dash-daq==0.5.0 \
    dash-bootstrap-components==1.4.1 \
    dash-core-components==2.0.0 \
    dash-html-components==2.0.0 \
    dash-mantine-components==0.12.1 \
    dash-table==5.0.0 \
    Werkzeug==2.2.2 \
    Flask==2.2.2 \
    Flask-Compress==1.13 \
    Flask-SeaSurf==1.1.1 \
    matplotlib==3.5.3 \
    pandas==1.3.5 \
    plotly==5.10.0 \
    statsmodels==0.13.5

# Create R package installation script
RUN echo '# Set CRAN mirror \n\
    options(repos = c(CRAN = "https://cloud.r-project.org")) \n\
    \n\
    # Install BiocManager \n\
    if (!require("BiocManager", quietly = TRUE)) \n\
        install.packages("BiocManager") \n\
    \n\
    # Install specific Bioconductor version compatible with R 4.2 \n\
    BiocManager::install(version = "3.16", force = TRUE, ask = FALSE) \n\
    \n\
    # Install required packages \n\
    BiocManager::install(c( \n\
        "DESeq2", \n\
        "edgeR", \n\
        "limma", \n\
        "statmod", \n\
        "ggplot2", \n\
        "dplyr", \n\
        "qsmooth" \n\
    ), force = TRUE, ask = FALSE)' > install_r_packages.R

# Run the R package installation script
RUN Rscript install_r_packages.R

# Expose the port the app runs on
EXPOSE 8050

# Command to run the app
CMD ["python", "app.py", "--host", "0.0.0.0"]

#####
# Original source: https://github.com/Christofferfrisk/rnalys
# Dockerfile created with Claude 3.5 Sonnet
