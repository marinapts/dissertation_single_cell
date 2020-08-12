#!/bin/bash

python scanpy_preprocess/preprocess.py \
--dataset E14_het E13_het \
--keep_only_highly_variable \
--write_to_file
