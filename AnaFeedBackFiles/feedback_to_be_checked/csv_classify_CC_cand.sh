#!/bin/bash
# Author: Jialin WU
# Date: 2025-04-08

# delete the columns before "old" "new"
# delete the rows before the first row of vertex information
CSV_FILE="/home/andywu/Desktop/PhysicsAnalysis_FASERnu/EventDisplayChecks/FASERnu_cand_20250408.csv"

for LineNum in {1..32}
do

    EvtType="$(sed -n "${LineNum}p" "$CSV_FILE" | cut -d',' -f2)CC"
    echo "Event Type: $EvtType"

    ModNum="$(sed -n "${LineNum}p" "$CSV_FILE" | cut -d',' -f3)"
    echo "Module Number: $ModNum"

    zone="$(sed -n "${LineNum}p" "$CSV_FILE" | cut -d',' -f4)"
    echo "Zone: $zone"

    # keep one decimal place
    vx="$(sed -n "${LineNum}p" "$CSV_FILE" | cut -d',' -f5 | awk '{printf "%.1f", $0}')"
    vy="$(sed -n "${LineNum}p" "$CSV_FILE" | cut -d',' -f6 | awk '{printf "%.1f", $0}')"
    # keep three figures
    pl="$(sed -n "${LineNum}p" "$CSV_FILE" | cut -d',' -f7 | awk '{printf "%03d", $0}')"

    echo "Vertex Position: ($vx, $vy, $pl)"

    FeedBackFileNamePattern="vtx_F${ModNum}_zone${zone}_reco*_${vx}_${vy}_pl${pl}_v*.feedback"
    echo "Feedback File Name Pattern: $FeedBackFileNamePattern"

    # find the first feedback file matching the pattern
    FeedBackFileName=$(find . -type f -name "$FeedBackFileNamePattern" | head -n 1)
    echo "Feedback File Name: $FeedBackFileName"
    if [[ -n "$FeedBackFileName" && -f "$FeedBackFileName" ]]; then
        echo "Feedback file found: $FeedBackFileName"
        mv "$FeedBackFileName" "./$EvtType/"
    else
        echo "Feedback file not found"
    fi

done

