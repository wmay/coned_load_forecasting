#!/bin/bash

# Running conda in a container is not easy! Must change the home directory to
# something writable. Must source conda's bashrc code. `source` is not available
# so use the equivalent `.`. Then finally we can activate
export HOME=/mnt/coe/Will
export CONDARC=/mnt/coe/Will/.condarc
. /mnt/coe/Will/.bashrc && \
    conda activate coned-load-forecasting
# conda info
python3 scripts/get_nwp_data.py && \
    python3 scripts/process_nwp_data.py && \
    Rscript scripts/forecast_tv.R
