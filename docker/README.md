# scNexus Docker Image

A single multi-arch image that runs all four scRNA-seq modules, with every
dependency baked in — no R, no compilers, and no package installs on the user's
machine. Built to sidestep the dependency-install pain (CellChat, LIANA,
biomaRt, Seurat ecosystem).

Base image is `rocker/r-ver` (amd64 + arm64). Each module runs as its own
Shiny process on its own port; a small landing page links to all four. There is
no shiny-server (it is amd64-only), so nothing here is architecture-locked.

## Build

From the **repository root**:

```bash
docker compose build          # or: docker build -f docker/Dockerfile -t scnexus-demo:latest .
```

First build is slow (~30–60 min) and the image is large (~3–5 GB) — normal for
the single-cell + Bioconductor stack. CRAN packages come as Posit Package
Manager binaries; Bioconductor and the GitHub packages compile from source.

## Run

```bash
docker compose up
# or: docker run --rm -p 3838-3842:3838-3842 scnexus-demo:latest
```

| URL | App |
|-----|-----|
| <http://localhost:3838> | Landing page |
| <http://localhost:3839> | Module I — Individual sample |
| <http://localhost:3840> | Module II — Integration |
| <http://localhost:3841> | Module III — Exploration |
| <http://localhost:3842> | Module IV — Cell-cell communication |

If one app crashes it is visible in `docker logs` and the others keep serving.

## Sample data

No data is bundled into the image. The host's `datasets/` folder is bind-mounted
to `/srv/app/datasets` and is the only data source.

- **Module I** uses a server-side file browser restricted to that mount.
- **Modules II–IV** upload from your computer — point them at the same
  `datasets/` folder on the host.

See [data preparation](../docs/README.md#data-preparation) for the demo dataset
download.

## Memory

Single-cell objects and CellChat/LIANA computations are memory-hungry.
`docker-compose.yml` sets a 16 GB limit; raise it (and Docker Desktop's memory
allocation under Settings → Resources) if an app dies mid-computation.

## Publishing a multi-arch image to Docker Hub

```bash
docker buildx build --platform linux/amd64,linux/arm64 \
  -f docker/Dockerfile -t <dockerhub-user>/scnexus:<tag> --push .
```

Building `amd64` on an Apple Silicon Mac runs under QEMU emulation and is very
slow for the source compiles — prefer a native amd64 machine or CI for that arch.

## Reproducibility / pinning

Nothing is version-pinned yet (by design, while components are still in flux).
The build records every resolved package version and GitHub commit SHA to
`/srv/BUILD_VERSIONS.txt` inside the image:

```bash
docker run --rm scnexus-demo:latest cat /srv/BUILD_VERSIONS.txt
```

When the stack stabilizes, pin: the GitHub refs in `docker/install_packages.R`
(`user/repo@<sha>`), a dated P3M snapshot for CRAN, and optionally add an
`renv.lock`.

## What's inside

- `Dockerfile` — `rocker/r-ver` base, apt system libs, package install, app copy
- `install_packages.R` — CRAN → Bioconductor → GitHub; fails the build if any
  package is missing; writes `BUILD_VERSIONS.txt`
- `start.sh` — entrypoint; launches the four apps + landing page
- `index.html` — landing page
