# Quadruple Graph

Algorithm to generate the Quadruple Graph based on 3 different genomes, returning all optimal solutions for the median problem for cycles of lenght up to 4.

## Input
This algorithm has the file "input.txt" as input file must be organized with 4 lines. The first line with the number of genes and the following 3 lines with each genome, where in each genome, each gene must be separated by a space. Note that the genome can have more than one cromossome and that the algorithm does not work properly if the Multiple Breakpoint Graph generated from the input has any double-edge. For example, a viable input file is:

### Example

7

( 1 2 -3 4 5 -6 7 )

( 1 3 2 4 ) ( 5 6 ) ( 7 )

( 1 -2 3 -6 -7 -4 5 )

## Execution
From the 3 genomes used as input, the first step is to create the Multiple Breakpoint Graph (MBG), this will create the corresponding image file but will not show immediately. The next step is to create the Quadruple Graph (QG), to find the solution, the algorithm uses a brute force approach to select the set of n vertex that will have the highest weight, where n is the number of genes. Each optimal solution for the median will output the same graph with nodes colored in different ways.

### Multiple Breakpoint Graph
The Graph will have the number of vertex equal to the number of genes x two, where each gene has a tail and a head. Each adjacency from each genome will be represented by an edge of different color: red, green and blue.

### Quadruple Graph
For each pair of nodes in the Multiple Breapoint Graph, there will be a new vertex in the Quadruple Graph and each edge represent a cycle of lenght 4 if both nodes were selected to be part of the solution. The nodes of color red represents the solution set while the ones in yellow all the others. The edges of color red represents an edge of weight 1 and the color green represents an edge of weight 2, in other words, two different cycles in the MBG.

## Output
There will be created a directory called **solutions** if there is not one. And for each time the algorithm is executed, there will also create a directory called **instance_K**, where k is the number of instances already created inside solutions.
Inside instance_K will have a txt file called **output.txt** with the input, the weight of the solution, the time that took to compute and all viable optimal solutions. There will also have a png file called **MBG.png** to save the Multiple Breakpoint Graph and for each optimal solution, a file called **solution_I.png** to show all the optimal solutions for the problem in the Quadruple Graph, where I is the number of solutions.
