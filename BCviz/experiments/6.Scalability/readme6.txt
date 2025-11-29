Scalability
========
Here, we can get the following figure of the paper.

Fig. 13
(1) To get the experimental results in a Linux server, please turn to the directory "~/BCviz/construct-BCviz", and enter the following command:

chmod 775 scalability.sh
nohup ./scalability.sh > ../experiments/6.Scalability/scalability.txt 2>&1 &

Then, the file "scalability.txt" will be generated in directory "~\BCviz\experiments\6.Scalability".

(2) To plot the figures, turn to the directory "~\BCviz\experiments\6.Scalability" and run the program "Fig.13.py".
Two figures "Marvel-scalability.jpg" and "WikiPedia-scalability.jpg" will be generated in current director.


Similar experimental results on scalability are obtained for MVB and MBB, which are omitted to save time.