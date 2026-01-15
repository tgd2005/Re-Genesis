#!/bin/bash

# Define sequences.txt path
input_file="sequences/sequences_3to5.txt"

# 1. Calculate totals
total_bases=0
total_strands=0
while read -r seq || [[ -n "$seq" ]]; do
    [[ -z "$seq" || "$seq" == ">"* ]] && continue
    ((total_strands++))
    total_bases=$((total_bases + ${#seq}))
done < "$input_file"

# 2. Write Header
echo "$total_bases $total_strands" > top_3to5.top

# 3. Write Data
current_strand=1
current_offset=0
while read -r seq || [[ -n "$seq" ]]; do
    [[ -z "$seq" || "$seq" == ">"* ]] && continue
    len=${#seq}
    for (( i=0; i<$len; i++ )); do
        base=${seq:$i:1}
        prev=$(( i == 0 ? -1 : current_offset + i - 1 ))
        next=$(( i == len - 1 ? -1 : current_offset + i + 1 ))
        # Use printf for strict spacing: StrandID Base Prev Next
        printf "%d %s %d %d\n" "$current_strand" "$base" "$prev" "$next" >> top_3to5.top
    done
    current_offset=$((current_offset + len))
    ((current_strand++))
done < "$input_file"
