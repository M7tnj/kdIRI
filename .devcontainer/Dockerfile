FROM rocker/rstudio:latest

# Give rstudio user passwordless sudo access
RUN echo "rstudio ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/rstudio

# Set a known password
RUN echo "rstudio:rstudio" | chpasswd

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libglpk-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libsodium-dev \
    libuv1-dev \
    libsqlite3-dev \
    && rm -rf /var/lib/apt/lists/*

# Configure R library path
RUN echo '.libPaths(c("/usr/local/lib/R/site-library", .libPaths()))' >> /home/rstudio/.Rprofile

# Set default working directory
RUN echo 'setwd("/home/rstudio/project")' >> /home/rstudio/.Rprofile

# Set up RStudio config directories
RUN mkdir -p /home/rstudio/.config/rstudio/themes

# Create symlinks for RStudio navigation
RUN ln -sf /workspaces /home/rstudio/workspaces

# Make workspace directories accessible
RUN chmod a+rx /workspaces 2>/dev/null || true

# Ensure rstudio owns their home directory files
RUN chown -R rstudio:rstudio /home/rstudio