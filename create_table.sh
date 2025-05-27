#!/bin/bash

# Input file containing the commands
COMMANDS_FILE="tests_commands/std_e.txt"

# Output file to store the compiled results
COMPILED_OUTPUT="compiled_output.csv"

# Debug log file
DEBUG_LOG="debug_log.txt"

# Clear the compiled output and debug log files if they exist
> "$COMPILED_OUTPUT"
> "$DEBUG_LOG"

# Initialize a flag to track if the header has been written
HEADER_WRITTEN=false

# Loop through each line in the commands file
while IFS= read -r command; do
    # Log the command being executed
    echo "Running command: $command" >> "$DEBUG_LOG"

    # Execute the command
    eval "$command"

    # Check if the output.csv file exists
    if [[ -f "output.csv" ]]; then
        # Log the content of output.csv for debugging
        echo "Content of output.csv:" >> "$DEBUG_LOG"
        cat output.csv >> "$DEBUG_LOG"
        echo "----------------------------------------" >> "$DEBUG_LOG"

        # If the header has not been written, write it to the compiled output file
        if [[ "$HEADER_WRITTEN" == false ]]; then
            head -n 1 output.csv >> "$COMPILED_OUTPUT"
            HEADER_WRITTEN=true
        fi

        # Append the content of output.csv (excluding the header) to the compiled output file
        # Format numbers: integers remain as integers, floats limited to 2 decimal places
        tail -n +2 output.csv | awk -F, '{
            for (i=1; i<=NF; i++) {
                if ($i ~ /^[0-9]+$/) {
                    # Integer: print as-is
                    printf "%d", $i
                } else if ($i ~ /^[0-9.-]+$/) {
                    # Float: limit to 2 decimal places
                    printf "%.2f", $i
                } else {
                    # Non-numeric: print as-is
                    printf "%s", $i
                }
                if (i < NF) printf "," # Add a comma between fields
            }
            printf "\n"
        }' >> "$COMPILED_OUTPUT"
    else
        echo "Warning: output.csv not found after running command: $command" >> "$DEBUG_LOG"
    fi
done < "$COMMANDS_FILE"

echo "All commands executed. Results compiled in $COMPILED_OUTPUT."