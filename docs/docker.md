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
| <http://localhost:3839> | Module I — Individual-sample analysis |
| <http://localhost:3840> | Module II — Multi-sample integration |
| <http://localhost:3841> | Module III — Post-annotation exploration |
| <http://localhost:3842> | Module IV — Cell-cell communication |

Each module runs as its own process; if one crashes the others keep serving
(check `docker logs scnexus-demo`).

## Sample data

The bundled sample in `test-data/` is baked into the image.

**Module I** uses a server-side file browser restricted to two roots — it cannot
see anything else in the container:

- **Sample data** (the default root) — the bundled sample. Pick `21401X3` to
  run the demo straight away.
- **Datasets** — your own data. Drop sample folders into the `datasets/` folder
  in the project directory; it is bind-mounted into the container. In the
  browser, switch the root dropdown from *Sample data* to *Datasets* to load
  them.

**Modules II–IV** upload from your computer, so copy the bundled data to the
host first:

```bash
./extract-test-data.sh            # -> ~/Downloads/scNexus-test-data
./extract-test-data.sh /some/dir  # custom destination
```

## Troubleshooting

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
