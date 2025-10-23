#!/bin/bash

# Default values for input and output files
COMMANDS_FILE="tests_commands/lambda_tables.txt"
COMPILED_OUTPUT="compiled_output.csv"
DEBUG_LOG="debug_log.txt"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --commands-file)
            COMMANDS_FILE="$2"
            shift 2
            ;;
        --compiled-output)
            COMPILED_OUTPUT="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Debug log file
> "$DEBUG_LOG"

# Clear the compiled output file if it exists
> "$COMPILED_OUTPUT"

# Initialize a flag to track if the header has been written
HEADER_WRITTEN=false

# Use a while loop with a subshell to ensure the last line is processed
while IFS= read -r command || [[ -n "$command" ]]; do
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
        tail -n +2 output.csv >> "$COMPILED_OUTPUT"
    else
        echo "Warning: output.csv not found after running command: $command" >> "$DEBUG_LOG"
    fi
done < "$COMMANDS_FILE"

echo "All commands executed. Results compiled in $COMPILED_OUTPUT."