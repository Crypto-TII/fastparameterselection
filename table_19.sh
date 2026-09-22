#!/bin/bash

python3 src/estimate.py --param "logq" --lambda "100" --n "1024" --hw "64" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "100" --n "1024" --hw "128" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "100" --n "1024" --hw "192" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "128" --n "1024" --hw "64" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "128" --n "1024" --hw "128" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "128" --n "1024" --hw "192" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "192" --n "1024" --hw "64" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "192" --n "1024" --hw "128" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "192" --n "1024" --hw "192" --secret "sparse" -v --table

python3 src/estimate.py --param "logq" --lambda "100" --n "32768" --hw "64" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "100" --n "32768" --hw "128" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "100" --n "32768" --hw "192" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "128" --n "32768" --hw "64" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "128" --n "32768" --hw "128" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "128" --n "32768" --hw "192" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "192" --n "32768" --hw "64" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "192" --n "32768" --hw "128" --secret "sparse" -v --table
python3 src/estimate.py --param "logq" --lambda "192" --n "32768" --hw "192" --secret "sparse" -v --table