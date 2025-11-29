#!/bin/bash
datasets=("marvel20" "marvel40" "marvel60" "marvel80" "marvel" "Twitter20" "Twitter40" "Twitter60" "Twitter80" "Twitter")
algorithms=("BCviz" "BCviz+" "BCviz-")
problems=("MEB")


echo "/---------IndexConstructionTime-for-Scalability---(s)--------/"

for dataname in "${datasets[@]}"; do
	for pro in "${problems[@]}"; do
		for algo in "${algorithms[@]}"; do
			./BCviz $dataname $pro ${algo}
		done
	done
	echo ""
done

