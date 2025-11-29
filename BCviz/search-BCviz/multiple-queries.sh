#!/bin/bash

problems=(MEB MVB MBB )
datasets=(marvel writer BookCrossing Team Actor-movie Twitter WikiPedia DBLP )
algorithms=(BCviz BCviz+ BCviz- )

echo "/---------MultipleSearchTime---(ms)--------/"

for dataname in ${datasets[*]}; do
    echo "----------$dataname-----------"
    for pro in ${problems[*]}; do
        for s in $(seq 3 5); do
            for t in $(seq 3 5); do
		for algo in ${algorithms[*]}; do
		    ./MBS $dataname $pro $algo $s $t
                done
            done
        done
	echo ""
    done
done