#!/bin/bash

LAMBDAS=("80" "100" "128" "192" "256")
NS=("1024" "2048" "4096" "32768")
SECRETS=("binary" "ternary")
ERROR="gaussian"

for n in "${NS[@]}"; do
  for secret in "${SECRETS[@]}"; do

    # Skip invalid combinations if needed
    if [[ "$n" == "2048" && "$secret" == "ternary" ]]; then
      continue 
    fi
    if [[ "$n" == "32768" && "$secret" == "binary" ]]; then
      continue 
    fi

    for lambda in "${LAMBDAS[@]}"; do
      python3 src/estimate.py --param "logq" --lambda "$lambda" --n "$n" --secret "$secret" --error "$ERROR" --std "3.19" -v --table
    done
  done
done