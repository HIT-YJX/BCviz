MBS
==========

This code implements the search algorithms of BCviz called MBS.<br> 
We have compiled the code to 'MBS' with g++ 7.5.0, and you can run it directly.

---

## 1. run.

To run the algorithm, first type ``chmod +x MBS`` or enter ``make`` to complie the code.<br>

Then enter the command in the following format:

```shell
./MBS [dataset] [problem_type] [algo_type] [s] [t] ([BCviz_file])
```

where<br> 
``[dataset]`` is a string denoting the name of a dataset, such as "marvel".<br>
``[problem_type]`` is a string denoting the types of the problem, which can be MEB, MVB, and MBB.<br>
``[algo_type]`` is a string denoting the types of the algorithm, which can be BCviz, BCviz+, and BCviz-.<br>
``[s]`` is an positive integer, which is the threshold on the first vertex set of bipartite graph.<br>
``[t]`` is an positive integer, which is the threshold on the second vertex set of bipartite graph.<br>
``([BCviz_file])`` is an optional string, which is the address of index file. The default value is ``../Index-results/[dataset]_[problem_type]_[algo_type].txt``. <br>


### As an example:

```shell
./MBS writer MEB BCviz- 4 4
./MBS writer MBB BCviz 3 3 ../Index-results/writer_MBB_BCviz_1_1.txt 
```

If you want to save the outputs of the algorithm, please type:

```shell
nohup ./MBS [dataset] [problem_type] [algo_type] [s] [t] ([BCviz_file]) >> ../log/search.txt 2>&1 &
```

Then you can view the output in the file ``../log/search.txt``.

---

## 2. the results.

The algorithm will output the results in the following format:

```
[problem_type] [algo_type] [dataset] [s_min] [t_min] [total_time]
```

where<br> 
``[total_time]`` is the search time (millisecond) of this algorithm.<br>

