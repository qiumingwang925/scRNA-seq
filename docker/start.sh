#!/usr/bin/env bash
## ABOUTME: Container entrypoint — launches the four module apps and a landing page.
## ABOUTME: Each app runs as its own R/Shiny process; the landing page is the foreground process.
set -euo pipefail

run_module () {  # $1 = module dir name, $2 = port
  echo "Starting $1 on port $2"
  R --quiet --no-save -e "shiny::runApp('/srv/app/$1', host = '0.0.0.0', port = $2L)" &
}

run_module module-I   3839
run_module module-II  3840
run_module module-III 3841
run_module module-IV  3842

# Landing page (foreground): serves index.html. If this exits the container
# stops; if an individual app crashes the others keep serving (visible in logs).
echo "Starting landing page on port 3838"
exec R --quiet --no-save -e "httpuv::runServer('0.0.0.0', 3838L, list(call = function(req) list(status = 200L, headers = list('Content-Type' = 'text/html; charset=utf-8'), body = paste(readLines('/srv/app/index.html'), collapse = '\n'))))"
