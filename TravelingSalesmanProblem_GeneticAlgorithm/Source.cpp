#include <iostream>
#include <vector>
#include <iomanip>

#include "TravelingSalesmanProblem.h"
#include "TSPUtils.h"

void outputPath(const std::vector<int>& path)
{
	const int sz = path.size();

	if (sz == 0)
		return;

	//Output path
	for (size_t i = 0; i < sz; i++)
	{
		std::cout << path[i] << " -> ";
	}

	//Return to start city
	if (sz > 1) 
	{
		std::cout << path[0];
	}

	//New line
	std::cout << std::endl;
}

void printMatrix(const std::vector<std::vector<int>>& mat)
{
	for (const auto& row : mat)
	{
		for (const auto& val : row)
		{
			std::cout << std::setw(3) << val << " ";
		}
		std::cout << std::endl;
	}
}

int main() 
{
	//Input adj matrix
	std::vector<std::vector<int>> adjMat;

	//Generate the matrix
	std::cout << "Genereting matrix..." << std::endl;
	adjMat = TSPUtils::generateTspAdjMatrix(20, TSPUtils::TspDatasetType::RandomUniform);
	std::cout << "Matrix generated!" << std::endl;

	//Print the matrix
	std::cout << "Matrix:" << std::endl;
	printMatrix(adjMat);

	//if matrix is not square => invalid input
	if (adjMat.size() == 0 || (adjMat.size() != adjMat[0].size()))
	{
		std::cout << "Invalid input!" << std::endl;
		return -1;
	}

	//Nubmer of cities
	const int nCities = adjMat.size();

	//Set up params;
	int ng = 0; //number of generation
	int npop = 0; //number of children in a generation
	int nnoimpr = 0; //after how many generations with no improvement to stop
	float pc = 0.90f; //probability of crossover
	float pm = 0.05f; //probability of mutation

	//set up params based on num cities
	if (nCities <= 30)
	{
		ng = 500;
		npop = 50;
		nnoimpr = 100;
	}
	else if (nCities >= 30 && nCities <= 100)
	{
		ng = 1000;
		npop = 100;
		nnoimpr = 100;

	}
	else // nCities > 100
	{
		ng = 2000;
		npop = 200;
		nnoimpr = 200;
	}

	//Solve - unseeded
	TravelingSalesmanProblem tsp = TravelingSalesmanProblem(adjMat, ng, npop, nnoimpr, pc, pm);
	tsp.solve();

	//Output path and dist
	std::cout << "Path: " << std::endl;
	outputPath(tsp.getCurrSolutionPath());
	std::cout << "Dist: " << tsp.getCurrSolutionDist() << std::endl;


	std::cin.get();
	return 0;
}