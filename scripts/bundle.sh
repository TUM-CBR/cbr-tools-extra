#!/bin/bash

pyinstaller \
    -n cbrtools cbrtools.py \
    --collect-all primer3 \
    --collect-data open3d \
    --hidden-import open3d