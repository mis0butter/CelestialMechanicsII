#!/bin/bash

# Input file
input_file="starkCoeffs_generic.m"

dos2unix "$input_file"

# Temporary file
temp_file="temp.m" 
> "$temp_file" 

# Initialize previous line variable
prev_line=""

# Read file line by line
while IFS= read -r line
do

    # Check if line starts with hashtag
    if [[ $line == \#* ]]; then

        # Append hashtag line to end of previous line
        prev_line="${prev_line}${line:1}"
        
    else
        # If previous line is not empty, write it to temporary file
        if [[ -n $prev_line ]]; then
            echo "$prev_line" >> "$temp_file"
        fi

        # Update previous line variable
        # prev_line="$line"
        prev_line=$(echo "$line" | sed 's/[ \t]*$//')

    fi

done < "$input_file"

# Write last line to temporary file
echo "$prev_line" >> "$temp_file"

# Use sed to replace "**" with "^"
sed -i 's/\*\*/^/g' "$temp_file"

# % append ; to end of each line 
sed -i 's/$/ ;/g' "$temp_file" 

# remove all # 
# sed -i 's/#//g' "$temp_file"

# Replace input file with temporary file
mv "$temp_file" "$input_file"
