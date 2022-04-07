#!/bin/bash
FILENAME=$1
FILE_INDEX=1
CONTIG_INDEX=1

while read LINE; do
    if [ `expr $FILE_INDEX % 2` -eq 0 ]; then
        echo "SEQUENCE_TEMPLATE="$LINE
		echo "="
    else
        echo "SEQUENCE_ID=oligo"$CONTIG_INDEX
		CONTIG_INDEX=$(expr $CONTIG_INDEX + 1)
    fi    
	FILE_INDEX=$(expr $FILE_INDEX + 1)
done < $FILENAME

