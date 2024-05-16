#!/bin/bash

# List of files
files=(

)
# Loop through each file
for file in "${files[@]}"; do
    # Create a new screen session with a unique name based on the file
    screen -d -m  bash -c "sftp beluga << EOF
    get $file
    quit
    EOF
    "
done
