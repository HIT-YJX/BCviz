#!/bin/bash
datasets=("writer" "WikiPedia")
algorithms=("BCviz" "BCviz+" "BCviz-")
problems=("MEB")

echo "/---------VariedSearchTime---(ms)--------/"

for dataname in "${datasets[@]}"; do
	for num in $(seq 3 6); do
		for pro in "${problems[@]}"; do
			for algo in "${algorithms[@]}"; do
				./MBS $dataname $pro ${algo} $num $num
			done
		done
		echo ""
	done
done

