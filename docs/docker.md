# Setup — Running the Platform with Docker

The whole platform ships as a single Docker image with every dependency baked in
(Seurat, CellChat, LIANA, biomaRt, and the rest of the single-cell stack). No R,
no compilers, and no package installs on your machine — just Docker.

Once it's running, follow the [workflow guides](README.md) to use each module.

## Prerequisites

- **Docker Desktop** (or Docker Engine) installed and running.
- **Memory:** single-cell objects and CellChat/LIANA are memory-hungry. Give
  Docker at least **16 GB** (Docker Desktop → Settings → Resources → Memory).

## Pull the prebuilt image (recommended)

Skip the build entirely — the published image runs on both **x86-64
(Intel/AMD)** and **ARM64 (Apple Silicon)**:

```bash
docker pull jixianli/scnexus:latest      # or a pinned release: v0.1.0
```

Docker picks the right architecture for your machine automatically; there is
nothing arch-specific to choose. Single-architecture tags
(`v0.1.0-arm64`, `v0.1.0-amd64`) also exist if you ever need to force one.

Run it with the same ports and data mount the compose file uses, from the
**repository root** (so `datasets/` resolves):

```bash
docker run --rm \
  --name scnexus-demo \
  -p 3838-3842:3838-3842 \
  -v "$(pwd)/datasets:/srv/app/datasets" \
  -e SCNEXUS_DATA_ROOT=/srv/app/datasets \
  --memory 16g \
  jixianli/scnexus:latest
```

Then open the URLs in [Run](#run) below. To use the prebuilt image with
`docker compose` instead, comment out the `build:` block in
`docker-compose.yml` and set `image: jixianli/scnexus:latest`.

## Build

Only needed to run your own changes — the published image already covers both
architectures. From the **repository root**:

```bash
SOURCE_COMMIT=$(git rev-parse HEAD) docker compose build
```

The first build is slow (~30–60 min) and the image is large (~3–5 GB) — normal
for the Bioconductor + single-cell stack. Subsequent builds are cached.

`SOURCE_COMMIT` stamps the image with the commit it came from (see [Recorded
versions](#recorded-versions)). Plain `docker compose build` works too and
records `unknown`.

## Run

```bash
docker compose up
```

Then open the module you want in a browser:

| URL | App |
|-----|-----|
| <http://localhost:3838> | Landing page (links to all four) |
| <http://localhost:3839> | Module I — scNexus-Process |
| <http://localhost:3840> | Module II — scNexus-Integrate |
| <http://localhost:3841> | Module III — scNexus-Explore |
| <http://localhost:3842> | Module IV — scNexus-Interact |

Each module runs as its own process; if one crashes the others keep serving
(check `docker logs scnexus-demo`).

## Sample data

No data is baked into the image. Everything is read from the `datasets/` folder
in the project directory, which is bind-mounted into the container.

Download the demo datasets into `datasets/` before starting — see [data
preparation](README.md#data-preparation) in the user guide.

**Module I**'s server-side file browser is restricted to that mount and cannot
see anything else in the container. **Modules II–IV** upload from your computer,
so point their file pickers at the same `datasets/` folder on the host.

## Troubleshooting

- **Module IV's LIANA analysis runs out of memory** — a known limitation of the
  containerized build. Run that step from RStudio instead; see the
  [RStudio setup guide](shiny.md). Everything else in the platform, Module IV's
  CellChat analysis included, runs in Docker.
- **An app dies mid-computation** — almost always memory. Raise `mem_limit` in
  `docker-compose.yml` and Docker Desktop's memory allocation, then re-run.
- **Check what's running** — `docker logs scnexus-demo` shows each module's
  startup and any crash.

## Recorded versions

**Which source an image came from.** Every image is labelled with the commit it
was built from. This reads the label without downloading or running the image,
so it works for any architecture:

```bash
docker image inspect jixianli/scnexus:latest \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}'
```

`unknown` means the image was built without `SOURCE_COMMIT` set, so its source
cannot be identified — expect that from ad-hoc local builds, not from releases.

**Which package versions an image contains.** Nothing is version-pinned yet, so
two builds made at different times can resolve different package versions. The
build records every resolved version and GitHub SHA inside the image:

```bash
docker run --rm jixianli/scnexus:latest cat /srv/BUILD_VERSIONS.txt
```

The source commit is appended to the end of that file as well.
