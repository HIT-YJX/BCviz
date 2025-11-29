#!/bin/bash
datasets=("marvel" "writer" "BookCrossing" "Team" "Actor-movie" "Twitter" "WikiPedia" "DBLP")
algorithms=("BCviz" "BCviz+" "BCviz-")
problems=("MEB" "MVB" "MBB")

echo "/---------SearchTime---(ms)--------/"

for dataname in "${datasets[@]}"; do
	for pro in "${problems[@]}"; do
		for algo in "${algorithms[@]}"; do
			./MBS $dataname $pro ${algo} 4 4
		done
	done
	echo ""
done

