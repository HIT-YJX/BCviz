Search Time
========
Here, we can get the following figures of the paper.

Fig. 10
(1) To get results of Search Time in Linux server, please turn to directory "~/BCviz/search-BCviz" and enter:

chmod 775 search-time.sh
nohup ./search-time.sh > ../experiments/4.SearchTime/search-time.txt 2>&1 &

Then, the file "search-time.txt" will be generated in directory "~\BCviz\experiments\4.SearchTime".

(2) To plot the figures, turn to directory "~\BCviz\experiments\4.SearchTime" and run program "Fig.10.py".
Three figures "MEB-search-time.jpg", "MVB-search-time.jpg" and "MBB-search-time.jpg" will be generated in current directory.


Fig. 11
(1) To get results of Search Time with varied parameters in Linux server, please turn to directory "~/BCviz/search-BCviz" and enter:

chmod 775 varied-search-time.sh
nohup ./varied-search-time.sh > ../experiments/4.SearchTime/varied-search-time.txt 2>&1 &

Then, the file "varied-search-time.txt" will be generated in directory "~\BCviz\experiments\4.SearchTime".

(2) To plot the figures, turn to directory "~\BCviz\experiments\4.SearchTime" and run program "Fig.11.py".
Two figures "Writers-varied-search-time.jpg" and "WikiPedia-varied-search-time.jpg" will be generated in current directory.