#!/bin/bash
python3 src/estimate.py --param "lambda" --n "512" --logq "11" --secret "binomial" --error "binomial" --std "1.22" -v
python3 src/estimate.py --param "n" --lambda "128" --logq "11" --secret "binomial" --error "binomial" --std "1.22" -v
python3 src/estimate.py --param "logq" --lambda "128" --n "512" --secret "binomial" --error "binomial" --std "1.22" -v
python3 src/estimate.py --param "std_e" --lambda "128" --n "512" --logq "11" --secret "binomial" -v

python3 src/estimate.py --param "lambda" --n "512" --logq "13" --secret "binomial" --error "binomial" --std "4.05" -v
python3 src/estimate.py --param "n" --lambda "128" --logq "13" --secret "binomial" --error "binomial" --std "4.05" -v
python3 src/estimate.py --param "logq" --lambda "128" --n "512" --secret "binomial" --error "binomial" --std "4.05" -v
python3 src/estimate.py --param "std_e" --lambda "128" --n "512" --logq "13" --secret "binomial" -v