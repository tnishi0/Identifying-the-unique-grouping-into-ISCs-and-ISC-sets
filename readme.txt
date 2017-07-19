========================================================
 Identifying the unique grouping into ISCs and ISC sets

 Files
-------

The following files are included in the software folder:

Instructions:
readme.txt (this file)

Sample output:
sample_output_Fig1_network.txt
 System requirement
--------------------

Installation of Sagemath software is required. Installation guide for Linux and Mac OS X can be found at



 Installation and usage
------------------------

The quickest way to try the code is to place all the files in the folder where Sagemath is installed, and run it from that folder.

1) Finding all candidate CS patterns for a given network 

From the folder, enter

../sage all_CS_patterns.sage [filename]

where "[filename]" should be replaced by the name of an adjacency matrix file (such as Fig1_network.txt and Fig3_network.txt included in this package; the adjacency matrix is assumed to be binary and symmetric).  This will compute the geometric decomposition and all possible CS patterns for each subgroup component in the decomposition.  


2) Computing the unique grouping of the clusters in a given CS pattern and the associated cluster-based coordinate transformation

From the folder, enter

../local/bin/python grouping_clusters.py [filename]

where "[filename]" is the name of an adjacency matrix file.  It will ask for input on the CS pattern of your choice.  Enter the number of clusters and then the indices of the nodes in each cluster, and it will output a list ISCs and ISC sets, as well as the computed matrices U (transpose) and B.  The matrices will also be written into files U.txt and B.txt.


** To run the code from a folder in a different location, a suitable path for sage and python would be needed.