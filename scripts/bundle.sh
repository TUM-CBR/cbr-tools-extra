#!/bin/bash

pyinstaller \
    -n cbrtools cbrtools.py \
    --collect-all primer3 \
    --hidden-import open3d \
    --collect-data open3d \
    --exclude-module sklearn \
    --exclude-module scipy