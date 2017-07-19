Grouping synchronization clusters into ISCs and ISC sets
--------------------------------------------------------

Copyright (C) 2017  Young Sul Cho

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

The full text of the GNU General Public License can be found in the file
“license.txt".


 System requirement
--------------------

Installation of Sagemath software is required. Installation guide for Linux and Mac OS X can be found at

http://doc.sagemath.org/html/en/installation/binary.html#linux-and-os-x

The code has been tested with Sagemath 6.8 under Debian GNU/Linux and with Sagemath 7.5.1 under Mac OS X 10.11.6.


 Installation and usage
------------------------

The quickest way to try the code is to place all the files in the folder where Sagemath is installed, and run it from that folder.

1. *Finding all candidate CS patterns for a given network.* From the folder, enter

   ../sage all_CS_patterns.sage [filename]

   where "[filename]" should be replaced by the name of an adjacency matrix file (such as Fig1_network.txt and Fig3_network.txt included in this package; the adjacency matrix is assumed to be binary and symmetric).  This will compute the geometric decomposition and all possible CS patterns for each subgroup component in the decomposition.  


2. *Computing the unique grouping of the clusters in a given CS pattern and the associated cluster-based coordinate transformation.* From the folder, enter

   ../local/bin/python grouping_clusters.py [filename]

   where "[filename]" is the name of an adjacency matrix file.  It will ask for input on the CS pattern of your choice.  Enter the number of clusters and then the indices of the nodes in each cluster, and it will output a list ISCs and ISC sets, as well as the computed matrices U (transpose) and B.  The matrices will also be written into files U.txt and B.txt.


** To run the code from a folder in a different location, a suitable path for sage and python would be needed.
