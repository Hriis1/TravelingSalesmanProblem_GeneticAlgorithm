# Traveling Salesman Problem (Genetic Algorithm)

This project solves the Traveling Salesman Problem (TSP) in C++ using a genetic algorithm.

## Main files

- `TravelingSalesmanProblem_GeneticAlgorithm/Source.cpp`  
  Entry point (`main`). It creates test TSP instances, runs the GA solver, compares against nearest neighbor, and prints results.

- `TravelingSalesmanProblem_GeneticAlgorithm/TravelingSalesmanProblem.h`  
  The GA solver itself: population loop, selection, crossover, mutation, elitism, and stopping conditions.

- `TravelingSalesmanProblem_GeneticAlgorithm/Genome.h`  
  Defines one candidate tour (a genome): city order + total tour distance, with helpers to generate and score tours.

- `TravelingSalesmanProblem_GeneticAlgorithm/TSPUtils.h`  
  Utility functions: generate TSP adjacency matrices from points and compute a nearest-neighbor baseline distance.

## How this implementation of GA works

1. Create an initial population of random valid tours.
2. Score each tour by total route distance (including return to start).
3. Keep the best tours (elitism).
4. Build new tours by:
   - selecting parents (tournament selection),
   - combining parents (edge recombination crossover),
   - occasionally mutating a tour (displacement mutation).
5. Repeat for many generations or stop early when there is no improvement.
6. Return the best tour found.
