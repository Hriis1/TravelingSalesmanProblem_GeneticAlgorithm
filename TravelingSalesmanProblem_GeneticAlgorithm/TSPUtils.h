#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <random>
#include <algorithm>

namespace TSPUtils
{

    enum class TspDatasetType {
        RandomUniform,
        Circular,
        ClusteredDeceptive
    };

    int distRoundedInt(double x1, double y1, double x2, double y2)
    {
        const double dx = x1 - x2;
        const double dy = y1 - y2;
        return (int)std::lround(std::sqrt(dx * dx + dy * dy));
    }

    std::vector<std::vector<int>> generateTspAdjMatrix(unsigned int numCities, TspDatasetType type, int planeSize = 1000, unsigned int seed = std::random_device{}())
    {
        if (numCities < 2) return {};

        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> uni(0.0, (double)planeSize);

        std::vector<std::pair<double, double>> pts;
        pts.reserve(numCities);

        if (type == TspDatasetType::RandomUniform)
        {
            for (unsigned int i = 0; i < numCities; i++)
                pts.push_back({ uni(rng), uni(rng) });
        }
        else if (type == TspDatasetType::Circular)
        {
            // Points on a circle + small noise
            const double cx = planeSize * 0.5;
            const double cy = planeSize * 0.5;
            const double R = planeSize * 0.40;

            std::normal_distribution<double> noise(0.0, planeSize * 0.01); // 1% noise

            for (unsigned int i = 0; i < numCities; i++)
            {
                const double angle = (2.0 * M_PI * (double)i) / (double)numCities;
                double x = cx + R * std::cos(angle) + noise(rng);
                double y = cy + R * std::sin(angle) + noise(rng);
                x = std::max(0.0, std::min(x, (double)planeSize)); //clamp
                y = std::max(0.0, std::min(y, (double)planeSize)); //clamp
                pts.push_back({ x, y });
            }

            // Optional: shuffle so the “correct” tour isn't trivially 0..n-1
            std::shuffle(pts.begin(), pts.end(), rng);
        }
        else // ClusteredDeceptive
        {
            // Make K cluster centers far apart, then sample points near centers.
            // This tends to create hard "cluster ordering" decisions.
            const unsigned int K = std::max(2u, (unsigned int)std::sqrt((double)numCities / 2.0));
            std::vector<std::pair<double, double>> centers;
            centers.reserve(K);

            // Place centers in a rough grid to ensure separation
            const unsigned int grid = (unsigned int)std::ceil(std::sqrt((double)K));
            const double step = planeSize / (double)(grid + 1);

            unsigned int placed = 0;
            for (unsigned int gy = 1; gy <= grid && placed < K; gy++)
            {
                for (unsigned int gx = 1; gx <= grid && placed < K; gx++)
                {
                    // jitter center a bit
                    std::normal_distribution<double> jitter(0.0, step * 0.08);
                    double x = gx * step + jitter(rng);
                    double y = gy * step + jitter(rng);
                    x = std::max(0.0, std::min(x, (double)planeSize)); //clamp
                    y = std::max(0.0, std::min(y, (double)planeSize)); //clamp
                    centers.push_back({ x, y });
                    placed++;
                }
            }

            // Tight cluster radius
            const double sigma = planeSize * 0.02; // 2% of plane
            std::normal_distribution<double> clNoise(0.0, sigma);

            for (unsigned int i = 0; i < numCities; i++)
            {
                const auto& c = centers[i % K];
                double x = c.first + clNoise(rng);
                double y = c.second + clNoise(rng);
                x = std::max(0.0, std::min(x, (double)planeSize)); //clamp
                y = std::max(0.0, std::min(y, (double)planeSize)); //clamp
                pts.push_back({ x, y });
            }

            std::shuffle(pts.begin(), pts.end(), rng);
        }

        // Build adjacency matrix
        std::vector<std::vector<int>> adj(numCities, std::vector<int>(numCities, 0));
        for (unsigned int i = 0; i < numCities; i++)
        {
            for (unsigned int j = i + 1; j < numCities; j++)
            {
                int d = distRoundedInt(pts[i].first, pts[i].second, pts[j].first, pts[j].second);
                adj[i][j] = adj[j][i] = d;
            }
        }
        return adj;
    }
}
