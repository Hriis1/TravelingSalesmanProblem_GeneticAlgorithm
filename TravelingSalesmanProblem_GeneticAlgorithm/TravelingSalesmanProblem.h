#pragma once
#include <climits>

#include "Genome.h"
#include <array>

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

	std::vector<int> edgeRecombinationCrossover(const std::vector<int>& p1, const std::vector<int>& p2) 
	{

		//helper vars for left and right idx
		int leftNbrIdx = 0;
		int rightNbrIdx = 0;

		//Size of parent/child
		int N = p1.size();

		//child vec
		std::vector<int> child(N);

		//Empty parent
		if (N == 0)
			return child;

		//Adjacency matrix and num neighbors for all vertecies
		std::vector<std::array<int, 4>> adjMatrix(N);
		std::vector<uint8_t> deg(N, 0);

		//Build the adj matrix from both parents
		// helper: add neighbor if not already present
		auto addUnique = [&](int city, int neigh)
		{
			if (neigh == city) return;

			for (uint8_t t = 0; t < deg[city]; ++t)
				if (adjMatrix[city][t] == neigh)
					return; // duplicate

			// assumes deg[city] < 4
			adjMatrix[city][deg[city]++] = neigh;
		};

		//p1, dont check for duplicates
		for (size_t i = 0; i < N; i++)
		{
			int city = p1[i];

			//Get the left/right neighbors
			int leftCity = p1[(i - 1 + N) % N];
			int rightCity = p1[(i + 1) % N];

			// Add to adjacency list and increment deg
			adjMatrix[city][deg[city]++] = leftCity;
			adjMatrix[city][deg[city]++] = rightCity;
		}

		// p2, check for duplicates
		for (size_t i = 0; i < N; ++i)
		{
			int city = p2[i];

			//Get the left/right neighbors
			int leftCity = p2[(i - 1 + N) % N];
			int rightCity = p2[(i + 1) % N];

			// Add to adjacency list and increment deg, checking for duplicates
			addUnique(city, leftCity);
			addUnique(city, rightCity);
		}

		// Build the child from adj matrix
		// helper: remove "rem" from neighbor list of "city"
		auto removeNeighbor = [&](int city, int rem)
		{
			for (uint8_t i = 0; i < deg[city]; ++i)
			{
				if (adjMatrix[city][i] == rem)
				{
					// swap with last and shrink
					adjMatrix[city][i] = adjMatrix[city][deg[city] - 1];
					deg[city]--;
					return;
				}
			}
		};

		// helper: after choosing city "chosen", remove it from all its neighbors' lists
		auto removeChosenFromAll = [&](int chosen)
		{
			// chosen currently has up to deg[chosen] neighbors
			uint8_t d = deg[chosen];
			for (uint8_t i = 0; i < d; ++i)
			{
				int u = adjMatrix[chosen][i];
				removeNeighbor(u, chosen);
			}
			deg[chosen] = 0; // clear chosen's list
		};

		//used cities
		std::vector<char> used(N, 0);

		// --- O(1) fallback support: keep a dynamic list of unused cities ---
		std::vector<int> remaining;
		remaining.reserve(N);
		std::vector<int> posInRemaining(N, -1);

		for (int c = 0; c < N; ++c)
		{
			posInRemaining[c] = (int)remaining.size();
			remaining.push_back(c);
		}

		// helper: remove a city from the remaining list in O(1)
		auto removeFromRemaining = [&](int city)
		{
			int idx = posInRemaining[city];
			if (idx < 0) return; // already removed

			int last = remaining.back();
			remaining[idx] = last;
			posInRemaining[last] = idx;

			remaining.pop_back();
			posInRemaining[city] = -1;
		};

		// first city
		child[0] = p1[0];
		used[child[0]] = 1;
		removeFromRemaining(child[0]);
		removeChosenFromAll(child[0]);

		// ERX construction loop
		for (int i = 1; i < N; ++i)
		{
			int lastCity = child[i - 1];

			// collect best candidates among available (not used) neighbors of lastCity
			int bestDeg = INT_MAX;
			std::vector<int> candidates;
			candidates.reserve(4);

			for (uint8_t n = 0; n < deg[lastCity]; ++n)
			{
				int neigh = adjMatrix[lastCity][n];
				if (used[neigh]) continue;

				int d = (int)deg[neigh];
				if (d < bestDeg)
				{
					bestDeg = d;
					candidates.clear();
					candidates.push_back(neigh);
				}
				else if (d == bestDeg)
				{
					candidates.push_back(neigh);
				}
			}

			int nextCity = -1;

			if (!candidates.empty())
			{
				// tie-break randomly among min-degree candidates
				std::uniform_int_distribution<int> pick(0, (int)candidates.size() - 1);
				nextCity = candidates[pick(_gen)];
			}
			else
			{
				// fallback: pick a random unused city (O(1) now)
				std::uniform_int_distribution<int> pick(0, (int)remaining.size() - 1);
				nextCity = remaining[pick(_gen)];
			}

			child[i] = nextCity;
			used[nextCity] = 1;
			removeFromRemaining(nextCity);

			// critical ERX step: remove chosen city from all adjacency lists
			removeChosenFromAll(nextCity);
		}

		return child;
	}

	const std::vector<std::vector<int>>* _adjMat = nullptr;
	int _NG = 0;
	int _NPOP = 0;
	float _PC = 0.0f;
	float _PM = 0.0f;
	Genome _currSolution;
	std::mt19937 _gen;
};
