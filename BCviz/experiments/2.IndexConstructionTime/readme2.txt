Index Construction Time and Speedup of Construction Time
========

Here, we can get the following figures and tables of the paper.

Fig. 9
(1) To get results of Index Construction Time in Linux server, please turn to directory "~/BCviz/construct-BCviz" and enter:

nohup sh ./index-time.sh > ../experiments/2.IndexConstructionTime/index-time.txt 2>&1 &

The file "index-time.txt" will be generated in directory "~\BCviz\experiments\2.IndexConstructionTime".

(2) To plot the figures, turn to directory "~\BCviz\experiments\2.IndexConstructionTime" and run the program "Fig.9.py".
Three files "MEB-IndexTime.jpg", "MVB-IndexTime.jpg" and "MBB-IndexTime.jpg" will be generated in current directory.


Table 4
After running "Fig.9.py", the speedup ratios will be displayed on the console.