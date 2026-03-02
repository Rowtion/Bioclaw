#!/bin/bash
# Start RStudio Server
/init &

# Start JupyterLab
jupyter lab \
    --ip=0.0.0.0 \
    --port=8888 \
    --no-browser \
    --allow-root \
    --NotebookApp.token="${JUPYTER_TOKEN:-bioclaw}" \
    --notebook-dir=/home/rstudio/data &

wait
