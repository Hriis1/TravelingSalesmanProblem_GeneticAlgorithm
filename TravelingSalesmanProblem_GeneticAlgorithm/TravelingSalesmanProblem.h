#pragma once
#include <climits>

#include "Genome.h"

class TravelingSalesmanProblem
{
public:
	TravelingSalesmanProblem(const std::vector<std::vector<int>>& adjMat, int ng, int npop, float pc, float pm, unsigned int seed = std::random_device{}())
		:_adjMat(&adjMat), _NG(ng), _NPOP(npop), _PC(pc), _PM(pm), _currSolution(adjMat.size()), _gen(seed)
	{
		_currSolution.dist = INT_MAX;
	}

	void reinit(const std::vector<std::vector<int>>& adjMat, int ng, int npop, float pc, float pm)
	{
		_adjMat = &adjMat;
		_NG = ng;
		_NPOP = npop;
		_PC = pc;
		_PM = pm;
		_currSolution = Genome(adjMat.size());
		_currSolution.dist = INT_MAX;
	}

	void reseed(unsigned int seed = std::random_device{}())
	{
		_gen.seed(seed);
	}

	void solve()
	{

	}

	//getters
	int getCurrSolutionDist() const
	{
		return _currSolution.dist;
	}

	const std::vector<int>& getCurrSolutionPath() const
	{
		return _currSolution.path;
	}

private:
	const std::vector<std::vector<int>>* _adjMat = nullptr;
	int _NG = 0;
	int _NPOP = 0;
	float _PC = 0.0f;
	float _PM = 0.0f;
	Genome _currSolution;
	std::mt19937 _gen;
};
