BCviz
==========

This code implements the construction algorithms of BCviz.<br> 
We have compiled the code to 'BCviz' with g++ 7.5.0, and you can run it directly.

---

## 1. compile and run
To run the algorithm, first type ``chmod +x BCviz``, or enter ``make`` to complie the code.<br>

Then enter the command in the following format:

```shell
./BCviz [dataset] [problem_type] [algo_type] ([s_min]) ([t_min])
```

where<br> 
``[dataset]`` is a string denoting the name of a dataset, such as `marvel`.
``[problem_type]`` is a string denoting the types of the problem, which can be `MEB`, `MVB`, and `MBB`.<br> 
``[algo_type]`` is a string denoting the types of the algorithm, which can be `BCviz`, `BCviz+`, and `BCviz-`.<br> 
``([s_min])`` is an optional positive integer, which is the minimum threshold on the first vertex set of bipartite graph. The default value is `3`.<br> 
``([t_min])`` is an optional positive integer, which is the minimum threshold on the second vertex set of bipartite graph. The default value is `3`.<br> 


### As an example:

```shell
./BCviz writer MEB BCviz-
./BCviz writer MBB BCviz 1 1
```

If you want to save the output of the algorithm, please type:

```shell
nohup ./BCviz [dataset] [problem_type] [algo_type] ([s_min]) ([t_min]) >> ../log/construct.txt 2>&1 &
```

Then you can view the outputs in the file ``../log/construct.txt``.

---

## 2. the results.
The algorithm will output the results in the following format:

```shell
[problem_type] [algo_type] [dataset] [s_min] [t_min] [total_time]
```

where<br> 
``[total_time]`` is the running time (second) of this algorithm.<br> 

And the index named ``[dataset]_[problem_type]_[algo_type].txt`` or ``[dataset]_[problem_type]_[algo_type]_[s_min]_[t_min].txt`` will be generated in the directory ``~/BCviz/Index-results``<br> 


### As an example:

If you type:
```shell
./BCviz writer MEB BCviz-
./BCviz writer MBB BCviz 1 1
```
Then file "writer_MEB_BCviz-.txt" and file "writer_MBB_BCviz_1_1.txt" will be generated in the directory ``~/BCviz/Index-results``.<br>
