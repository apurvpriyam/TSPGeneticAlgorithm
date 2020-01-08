## Travelling Salesman Problem
The Traveling Salesman Problem (TSP) is a problem arising in op-erations research and computer science. The TSP is widely studieddue to its many applications to practical problems in STEM fields.Despite being a fundamental problem in research since its firstformulation in the 1930’s, the TSP remains to be solved completely- there currently does not exist a method for producing an exactsolution in polynomial time

## Definition
Intuitively, the TSP is the following: given a list of cities and theirlocations, find the shortest tour which visits all cities exactly onceand returns to the starting city.

Formally, this can be stated as, given a graph consisting of edgesand vertices,G=(V,E), if each edgee=(v<sub>i</sub>,v<sub>j</sub>)has corresponding weight w<sub>ij</sub>, find the cycle that minimizes the sum of weights such that each node is visited exactly once. In graph theory, a cycle that visits each node exactly once is referred to as a Hamiltonian Cycle. Thus, the TSP is equivalently to find the minimum Hamiltonian Cycle.

## Genetic Algorithm
A genetic algorithm (GA) is a metaheuristic method for solving op-timization problem inspired by the process of natural selection and 'survival of the fittest'. In a genetic algorithm, a group (called population) of candidate solutions (called individuals or chromosomes) to an optimization problem is evolved in each generation. Weak individuals are killed and stronger ones are given more chance to breed. A typical genetic algorithm has a genetic representation of the solution and a fitness function to score the solution. A standard representation of solution is an array of bits. Traditionally, solutions are represented in binary as strings of 0s and 1s, but encoding of solution can change depending on the problem, just like Travelling Salesman Problem.
