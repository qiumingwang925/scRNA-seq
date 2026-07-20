# syntax=docker/dockerfile:1
## ABOUTME: Builds a single multi-arch image that runs all four scRNA-seq Shiny apps.
## ABOUTME: rocker/r-ver base (amd64+arm64); each app is served on its own port by start.sh.

FROM rocker/r-ver:4.4.2

# System libraries required to compile the single-cell + Bioconductor stack
# (hdf5 for Seurat, fftw/gsl for stats deps, glpk/gmp for igraph, geos/udunits
# for spatial deps, cairo/freetype/harfbuzz for graphics devices).
RUN apt-get update && apt-get install -y --no-install-recommends \
      libhdf5-dev \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      libpng-dev \
      libtiff5-dev \
      libjpeg-dev \
      libfreetype6-dev \
      libharfbuzz-dev \
      libfribidi-dev \
      libcairo2-dev \
      libxt-dev \
      libfftw3-dev \
      libgsl-dev \
      libglpk-dev \
      libgmp-dev \
      libgeos-dev \
      libudunits2-dev \
      cmake \
      git \
    && rm -rf /var/lib/apt/lists/*

# Parallelise within-package compilation of the GitHub/Bioc source builds.
ENV MAKEFLAGS="-j4"

# Install packages before copying app code so edits to the apps don't bust the
# (expensive) package layer cache. The library is a BuildKit cache mount: it persists
# across builds (even failed ones), so a re-run resumes and only installs what's still
# missing instead of the whole list. Cache mounts are not part of the image, so copy the
# result into the real site-library once the install succeeds.
COPY docker/install_packages.R /tmp/install_packages.R
RUN --mount=type=cache,target=/rlib,sharing=locked \
    Rscript /tmp/install_packages.R /rlib && \
    cp -an /rlib/. /usr/local/lib/R/site-library/

RUN useradd --create-home --shell /bin/bash appuser

# App code, mirroring the on-disk repo layout: R/ is a sibling of the module dirs
# so each app's source("../R/utils.R") and module-I's shinyFiles "Project" root
# (dirname(getwd()) == /srv/app) resolve unchanged. No data is baked in — it
# arrives through the datasets/ bind mount.
COPY R          /srv/app/R
COPY module-I   /srv/app/module-I
COPY module-II  /srv/app/module-II
COPY module-III /srv/app/module-III
COPY module-IV  /srv/app/module-IV
COPY docker/index.html /srv/app/index.html
COPY docker/start.sh   /usr/local/bin/start.sh

RUN chmod +x /usr/local/bin/start.sh && chown -R appuser:appuser /srv/app

USER appuser

# 3838 = landing page; 3839-3842 = Modules I-IV.
EXPOSE 3838 3839 3840 3841 3842

CMD ["/usr/local/bin/start.sh"]
