Datasets in paper "BCviz: A Linear-Space Index for Mining and Visualizing Cohesive Bipartite Subgraphs"
=========

Here are the datasets used in this paper.<br>

To run certain experimental scripts, the datasets **must be strictly named** as follows:<br>  
**marvel**, **writer**, **BookCrossing**, **Team**, **Actor-movie**, **Twitter**, **WikiPedia**, **DBLP**.<br>

Among them, the datasets marvel20, marvel40, marvel60, and marvel80 were obtained by sampling 20%, 40%, 60%, and 80% of the edges from the original marvel dataset. <br>
Similarly, Twitter20, Twitter40, Twitter60, and Twitter80 were derived by sampling from the Twitter dataset.<br>

Please note that the DBLP dataset is compressed as DBLP.tar.xz. To use it, you need to extract it first. 
In a Linux environment, you can extract it with the following command:

```shell
tar -xJvf DBLP.tar.xz
```



All source datasets in the paper can also be found in KONECT (http://konect.cc/).<br>

marvel: http://konect.cc/networks/marvel/ <br>
writer: http://konect.cc/networks/dbpedia-writer/ <br>
BookCrossing: http://konect.cc/networks/bookcrossing_rating/ <br>
Team: http://konect.cc/networks/dbpedia-team/ <br>
Actor-movie: http://konect.cc/networks/actor-movie/ <br>
Twitter: http://konect.cc/networks/munmun_twitterex_ut/ <br>
WikiPedia: https://konect.cc/networks/wiki-en-cat/ <br>
DBLP: http://konect.cc/networks/dblp-author/ <br>
