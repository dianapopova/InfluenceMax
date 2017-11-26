# InfluenceMax
Data structures for an efficient computation of influence maximization and influence estimation.
# Getting Started
Download and install Webgraph framework from http://webgraph.di.unimi.it.

Download a dataset from http://law.di.unimi.it/datasets.php. Our programs use transpose (inverse) graphs, 
which have "-t" after the name. For example, cnr-2000-t.graph. The graph should not have self-loops.

Download IM_list.java, IM_flat.java, and/or IM_flat_compr.java from this project, 
to test different data strucures.
# Running the tests
Compile the programs using Webgraph library, for example:

```

javac -cp "lib/*" -d bin src/*.java

```

Run the program with parameters of your choice. For example:

```

java -Xmx16g -cp "../lib/*":"../bin" IM_flat cnr2000 0.1 32 10

```

IM_flat is the program to run; cnr2000 is the basename for the graph; 0.1 is p, the probability of edge existence;
32 is the value of beta (coefficient for calculating the weight of hypergraph); and 10 is k, the number of seeds.
# References
Diana Popova, Akshay Khot, and Alex Thomo. *Data Structures for an Efficient Computing of Influence Maximization
and Influence Estimation.* Submitted to EDBT/ICDT 2018 Joint Conference.
