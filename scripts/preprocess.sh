#!/bin/bash

python scanpy_preprocess/preprocess.py \
--dataset E14_hom E13_hom \
--keep_only_highly_variable
# --write_to_file
