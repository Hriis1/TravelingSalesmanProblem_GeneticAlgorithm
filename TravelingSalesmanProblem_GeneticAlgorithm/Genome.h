#pragma once

#include <vector>
#include <algorithm>
#include <random>

struct Genome
{
public:

	//vars
	int dist = 0;
	std::vector<int> path;

	//funcs
	Genome(int numCities) 
	{
		path.reserve(numCities);
	}

    void generatePath(const std::vector<std::vector<int>>& adjMat, int startCity)
    {
        int N = adjMat.size();

        path.clear();

        // 1. Fill with all cities
        for (int i = 0; i < N; ++i)
            path.emplace_back(i);

        // 2. Move startCity to front
        auto it = std::find(path.begin(), path.end(), startCity);
        std::iter_swap(path.begin(), it);

        // 3. Shuffle everything except first element
        static std::random_device rd;
        static std::mt19937 gen(rd());

        std::shuffle(path.begin() + 1, path.end(), gen);

        // 4. Compute distance
        dist = 0;

        for (int i = 0; i < N - 1; ++i)
            dist += adjMat[path[i]][path[i + 1]];

        // add return to start
        dist += adjMat[path[N - 1]][path[0]];
    }
};