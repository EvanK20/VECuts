----------------------------------------
Format of graph files
----------------------------------------

All graph files used by our programs must have the following format:

The first line contains two numbers n, m separated by space, where n = number of nodes, and m = number of edges.
The next m lines contain two numbers x, y separated by space, where (x,y) is a graph edge.
The range of node numbers is [1,...,n]. 

For the programs DFS-VE, DFS-VE(OGDF), and SPQR_VE, no multiple edges are allowed. 
To remove multiple edges, you may use the program "removeMultipleEdges.cpp".

Example of the contents of a graph file:
4 6
1 2
1 3
1 4
2 3
2 4
3 4


----------------------------------------
Graph generation
----------------------------------------

The program "randomGraph.cpp" generates graphs according to the Erdős–Rényi model. 
Input: <number of nodes> (int)
       <number of edges> (int)
       <create Hamiltonian cycle> (int) (0-1)
       <seed> (int)

Example of use: ./randomGraph 10000 15000 1 12345
Output: randGr10000.15000.1.12345


The program "makeBiconnected.cpp" augments an input graph by adding edges (if necessary) so that it is 2-vertex-connected.


----------------------------------------
Modification of OGDF
----------------------------------------

We modified the OGDF source file "OGDF/src/ogdf/graphalp/Triconnectivity.cpp", which contains the implementation of the triconnectivity algorithm 
of Hopcroft and Tarjan, so that we can handle graphs with more than ~45000 vertices.

What we did was to replace the recursive calls in algorithms DFS1, pathSearch, and pathFinder, with stacks.

This modification is contained in file "Triconnectivity.cpp".





