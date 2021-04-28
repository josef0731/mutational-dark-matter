#!/bin/bash 

while read c1 c2; do
	sbatch -p brc --job-name=$c1 --export=context1=$c1,context2=$c2 zoomvarCDS_32contexts.bash
done < all32contexts
