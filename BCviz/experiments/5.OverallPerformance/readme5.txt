5.Overall Performance
==============
Please finish the experiment "2.IndexConstructionTime" first.
And enter the following command in directory "~\BCviz\experiments\5.OverallPerformance" to get index time:
```shell
cp ../2.IndexConstructionTime/index-time.txt ./
```

Here, we can get the following figure of the paper.

Fig.12 
(1) To get the experimental results in a Linux server, please turn to the directory "~/BCviz/search-BCviz".
There are two schemes, and you can select either one to obtain the results:<br>

Scheme 1. (approximate experiment with 2 hours)<br>
First type:
```shell
chmod 775 multiple-queries.sh
nohup ./multiple-queries.sh > multiple-queries.txt 2>&1 &
```

2 hours later, the file "multiple-queries.txt" will be generated in current directory, then type:
```shell
chmod 775 data-collation
nohup ./data-collation > ../experiments/5.OverallPerformance/total-time.txt 2>&1 &
```

Then, the file "total-time.txt" will be generated in directory "~\BCviz\experiments\5.OverallPerformance".<br>

Scheme 2. (total experiment with around 7 days)<br>
First type:
```shell
chmod 775 multiple-queries-total.sh
nohup ./multiple-queries-total.sh > multiple-queries-total.txt 2>&1 &
```

7 days later, the file "multiple-queries-total.txt" will be generated in current directory, then type:
```shell
chmod 775 data-collation-total
nohup ./data-collation-total > ../experiments/5.OverallPerformance/total-time.txt 2>&1 &
```

Then, the file "total-time.txt" will be generated in directory "~\BCviz\experiments\5.OverallPerformance".

(2) To plot the figures, turn to the directory "~\BCviz\experiments\5.OverallPerformance" and run the program "Fig.12.py".
Three figures "MEB-total-time.jpg", "MVB-total-time" and "MBB-total-time.jpg" will be generated in current director.


