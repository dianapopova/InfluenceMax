# InfluenceMax
This repository contains efficient implementations for computing the influence maximization on large graphs using distinct data structures. The details of the implementations are described in the following paper:

Diana Popova, Akshay Khot, Alex Thomo: Data Structures for Efficient Computation of Influence Maximization and Influence Estimation. Accepted by EDBT/ICDT 2018 Joint Conference (Vienna, March 26-29, 2018).

There are three implementations of the algorithm in:

C. Borgs, M. Brautbar, J. Chayes, and B. Lucier. Maximizing social influence
in nearly optimal time. In SODA, pages 946–957, 2014.

using the WebGraph compression in:

P. Boldi and S. Vigna. The webgraph framework I: compression techniques. WWW'04. 

All implementations are coded in Java 8.
# IM_list
Data structure used in this implementation is list of lists (two-dimensional list).
# IM_flat
Data structure used in this implementation is one-dimensional array (flat array).
# IM_flat_compr
Data structures used in this implementation are one-dimensional array (flat array) and custom-compressed flat array.
# IM_flat_parallel
Data structure is the same as in IM_flat, but the hypergraph is built in parallel, independently, by all available cores on the machine, and then combined into one flat array by concatenating the results of each core.
# Borgs
Data structure for the hypergraph in this implementation is Webgraph (compressed adjacency lists). The hypergraph is built in parallel, independently, by all available cores on the machine, and then combined into one file stored on disk.
# Getting Started
a. Download and install Webgraph framework from http://webgraph.di.unimi.it.

b. Get dataset: download a dataset from http://law.di.unimi.it/datasets.php. Our programs use transpose (inverse) graphs, 
which have "-t" after the name. For example, cnr-2000-t.graph. These datasets are already in Webgraph format.

Or, if you would like to use another graph, convert a list of edges (must be sorted) into the Webgraph format. For example:

java -cp "../lib/*":"../bin" it.unimi.dsi.webgraph.BVGraph -1 -g ArcListASCIIGraph dummy cnr2000 < cnr2000-sortedEdges.txt

c. The dataset will be presented by two files, with the extensions .graph and .properties. You need another file, .offset. To get it:

java -cp "../lib/*":"../bin" it.unimi.dsi.webgraph.BVGraph -o -O -L cnr2000

d. The graph should not have self-loops. To eliminate the self-loops, download and run our custom program SelfLoopRemover.

e. Download IM_list.java, IM_flat.java, and/or IM_flat_compr.java from this project, to test different data structures.
# Running the tests
Compile the programs using Webgraph library, for example:

```

javac -cp "lib/*" -d bin src/*.java

```

Run the program with parameters of your choice. For example:

```

java -Xmx16g -cp "../lib/*":"../bin" IM_flat cnr-2000-t 0.1 32 10

```

IM_flat is the program to run; cnr-2000-t is the basename for the graph; 0.1 is p, the probability of edge existence;
32 is the value of beta (coefficient for calculating the weight of hypergraph); and 10 is k, the number of seeds.
# Contacts
If you have a question, please, send e-mail to dpopova@uvic.ca or thomo@uvic.ca
