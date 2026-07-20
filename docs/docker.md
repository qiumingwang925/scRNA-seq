# Setup — Running the Platform with Docker

The whole platform ships as a single Docker image with every dependency baked in
(Seurat, CellChat, LIANA, biomaRt, and the rest of the single-cell stack). No R,
no compilers, and no package installs on your machine — just Docker.

Once it's running, follow the [workflow guides](README.md) to use each module.

## Prerequisites

- **Docker Desktop** (or Docker Engine) installed and running.
- **Memory:** single-cell objects and CellChat/LIANA are memory-hungry. Give
  Docker at least **16 GB** (Docker Desktop → Settings → Resources → Memory).

## Build

From the **repository root**:

```bash
docker compose build
```

The first build is slow (~30–60 min) and the image is large (~3–5 GB) — normal
for the Bioconductor + single-cell stack. Subsequent builds are cached.

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

Nothing is version-pinned yet. The build records every resolved package version
and GitHub commit SHA inside the image:

```bash
docker run --rm scnexus-demo:latest cat /srv/BUILD_VERSIONS.txt
```
