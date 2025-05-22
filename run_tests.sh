#!/bin/bash

# Test script for running examples and checking for errors

# Function to run a command and check for errors
run_test() {
    local cmd="$1"
    echo "Running: $cmd"

    case "$output_mode" in
        all)
            eval "$cmd"
            ;;
        failed)
            eval "$cmd" > /dev/null 2>&1
            ;;
        none)
            eval "$cmd" > /dev/null 2>&1
            ;;
        *)
            echo "Invalid output mode: $output_mode"
            exit 1
            ;;
    esac

    # Check the exit status of the command
    if [ $? -ne 0 ]; then
        echo "Test FAILED: $cmd" >> test_results.log
        echo "Error encountered while running: $cmd"
        failed_commands+=("$cmd")  # Add the failed command to the list
        if [ "$output_mode" == "failed" ]; then
            echo "Output of failed test:"
            eval "$cmd"
        fi
    else
        echo "Test PASSED: $cmd" >> test_results.log
    fi
}

# Clear the log file
> test_results.log

# Initialize an array to store failed commands
failed_commands=()

# Default output mode
output_mode="all"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-mode)
            output_mode="$2"
            shift 2
            ;;
        *)
            input="$1"
            shift
            ;;
    esac
done

# Validate output mode
if [[ "$output_mode" != "all" && "$output_mode" != "failed" && "$output_mode" != "none" ]]; then
    echo "Invalid output mode: $output_mode"
    echo "Valid modes are: all, failed, none"
    exit 1
fi

# Check if input is provided
if [ -z "$input" ]; then
    echo "Usage: $0 [--output-mode <all|failed|none>] <commands_file_or_folder>"
    exit 1
fi

# Check if the input is a file or a folder
if [ -f "$input" ]; then
    # Input is a single file
    echo "Running tests from file: $input"
    while IFS= read -r cmd; do
        # Skip empty lines or comments
        [[ -z "$cmd" || "$cmd" =~ ^# ]] && continue
        run_test "$cmd"
    done < "$input"
elif [ -d "$input" ]; then
    # Input is a folder
    echo "Running tests from all .txt files in folder: $input"
    for file in "$input"/*.txt; do
        echo "Processing file: $file"
        while IFS= read -r cmd; do
            # Skip empty lines or comments
            [[ -z "$cmd" || "$cmd" =~ ^# ]] && continue
            run_test "$cmd"
        done < "$file"
    done

    # If there are failed commands, write them to a new file
    if [ ${#failed_commands[@]} -gt 0 ]; then
        failed_tests_file="$input/failed_tests.txt"
        echo "Writing failed commands to $failed_tests_file"
        > "$failed_tests_file"  # Clear or create the file
        for failed_cmd in "${failed_commands[@]}"; do
            echo "$failed_cmd" >> "$failed_tests_file"
        done
        echo "You can rerun the failed tests using: bash run_tests.sh $failed_tests_file"
    fi
else
    echo "Error: '$input' is neither a file nor a folder."
    exit 1
fi

# Print summary
if [ ${#failed_commands[@]} -eq 0 ]; then
    echo "All tests PASSED!"
else
    echo "The following tests FAILED:"
    for failed_cmd in "${failed_commands[@]}"; do
        echo "$failed_cmd"
    done
fi

echo "Test results logged in test_results.log"