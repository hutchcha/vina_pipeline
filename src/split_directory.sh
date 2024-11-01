#!/bin/bash

# Set the source directory (current directory)
src_dir="."

# Set the number of files per directory
files_per_dir=2000

# Set the prefix for the new directories
dir_prefix="part"

# Generate a list of all files in the source directory (excluding directories)
find "$src_dir" -maxdepth 1 -type f > filelist.txt

# Get the total number of files
total_files=$(wc -l < filelist.txt)

# Check if there are any files to process
if [ "$total_files" -eq 0 ]; then
    echo "No files found in the source directory."
    rm filelist.txt
    exit 0
fi

# Initialize counters
file_count=0       # Files moved in the current directory
file_count_total=0 # Total files moved
dir_count=1

# Create the first destination directory
mkdir -p "${dir_prefix}_${dir_count}"

# Start the timer
start_time=$(date +%s)

# Function to display the progress bar
show_progress() {
    local progress=$1
    local total=$2
    local elapsed=$3

    # Calculate percentage
    percent=$(( progress * 100 / total ))

    # Calculate the number of hash marks to display
    bar_length=50
    filled_length=$(( progress * bar_length / total ))

    # Create the bar
    bar=$(printf "%-${filled_length}s" "#" | tr ' ' '#')
    empty=$(printf "%-$((bar_length - filled_length))s" " " | tr ' ' '-')

    # Calculate estimated time remaining
    if [ "$progress" -gt 0 ]; then
        rate=$(( elapsed / progress )) # seconds per file
        remaining=$(( rate * (total - progress) ))
        # Format time as HH:MM:SS
        printf "Progress: |%s%s| %3d%% [%s] ETA: %02d:%02d:%02d\r" \
            "$bar" "$empty" "$percent" "$progress/$total files" \
            $((remaining/3600)) $(( (remaining%3600)/60 )) $((remaining%60))
    else
        printf "Progress: |%s%s| %3d%% [%s] ETA: --:--:--\r" \
            "$bar" "$empty" "$percent" "$progress/$total files"
    fi
}

# Read the list of files and move them
while IFS= read -r file; do
    # Move the file to the current directory
    mv "$file" "${dir_prefix}_${dir_count}/"

    # Increment the file counters
    ((file_count++))
    ((file_count_total++))

    # Show progress every 100 files or on the last file
    if (( file_count % 100 == 0 )) || (( file_count_total == total_files )); then
        current_time=$(date +%s)
        elapsed=$(( current_time - start_time ))
        show_progress "$file_count_total" "$total_files" "$elapsed"
    fi

    # If file_count reaches files_per_dir, reset it and move to next directory
    if (( file_count >= files_per_dir )); then
        # Reset file counter
        file_count=0
        # Increment directory counter
        ((dir_count++))
        # Create the next directory
        mkdir -p "${dir_prefix}_${dir_count}"
    fi

done < filelist.txt

# Ensure the final progress bar shows 100%
current_time=$(date +%s)
elapsed=$(( current_time - start_time ))
show_progress "$file_count_total" "$total_files" "$elapsed"
echo ""

# Remove the temporary file list
rm filelist.txt

# Calculate total elapsed time
end_time=$(date +%s)
total_elapsed=$(( end_time - start_time ))

# Format total elapsed time as HH:MM:SS
hours=$((total_elapsed/3600))
minutes=$(( (total_elapsed%3600)/60 ))
seconds=$((total_elapsed%60))

echo "Files have been split into $dir_count directories in $hours:$minutes:$seconds (HH:MM:SS)."
