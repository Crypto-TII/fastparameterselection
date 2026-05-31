#!/bin/bash

python3 src/estimate.py --param "lambda" --n "8192" --logq "200" --hw "128" --secret "sparse" -v --table
python3 src/estimate.py --param "lambda" --n "8192" --logq "119" --hw "128" --secret "sparse" -v --table
python3 src/estimate.py --param "lambda" --n "8192" --logq "87" --hw "128" --secret "sparse" -v --table
python3 src/estimate.py --param "lambda" --n "8192" --logq "210" --hw "192" --secret "sparse" -v --table
python3 src/estimate.py --param "lambda" --n "8192" --logq "128" --hw "192" --secret "sparse" -v --table
python3 src/estimate.py --param "lambda" --n "8192" --logq "91" --hw "192" --secret "sparse" -v --table

python3 src/estimate.py --param "lambda" --n "32768" --logq "850" --hw "128" --secret "sparse" -v --table
python3 src/estimate.py --param "lambda" --n "32768" --logq "500" --hw "128" --secret "sparse" -v --table
python3 src/estimate.py --param "lambda" --n "32768" --logq "330" --hw "128" --secret "sparse" -v --table
python3 src/estimate.py --param "lambda" --n "32768" --logq "850" --hw "192" --secret "sparse" -v --table
python3 src/estimate.py --param "lambda" --n "32768" --logq "565" --hw "192" --secret "sparse" -v --table
python3 src/estimate.py --param "lambda" --n "32768" --logq "410" --hw "192" --secret "sparse" -v --table