#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <vector>

struct Point {
    double x;
    double y;
};

class Solution {
public:
    std::tuple<std::vector<int>, double, double> computeTSP(
        const std::vector<Point>& points, const std::vector<std::vector<int>>& adj,
        int startIndex, int endIndex, double costPerUnit)
    {
        std::vector<int> path = nearestNeighborTour(points, adj, startIndex);
        reorderPathForStartEnd(path, startIndex, endIndex);

        double dist = pathDistance(path, points);
        double totalCost = dist * costPerUnit;

        return {path, dist, totalCost};
    }

private:
    double distance(const Point& a, const Point& b) const
    {
        double dx = a.x - b.x;
        double dy = a.y - b.y;
        return std::sqrt(dx * dx + dy * dy);
    }

    double pathDistance(const std::vector<int>& path,
                        const std::vector<Point>& points) const
    {
        double dist = 0.0;
        for (int i = 0; i + 1 < path.size(); ++i) {
            dist += distance(points[path[i]], points[path[i + 1]]);
        }
        return dist;
    }

    std::vector<int> nearestNeighborTour(const std::vector<Point>& points,
                                         const std::vector<std::vector<int>>& adj,
                                         int startIndex) const
    {
        int n = points.size();
        std::vector<bool> visited(n, false);
        std::vector<int> tour;
        tour.reserve(n);

        int current = startIndex;
        visited[current] = true;
        tour.push_back(current);

        for (int i = 1; i < n; ++i) {
            int next = -1;
            double bestDist = std::numeric_limits<double>::infinity();

            for (int candidate : adj[current]) {
                if (!visited[candidate]) {
                    double d = distance(points[current], points[candidate]);
                    if (d < bestDist) {
                        bestDist = d;
                        next = candidate;
                    }
                }
            }

            if (next == -1) {
                double globalBest = std::numeric_limits<double>::infinity();
                int globalNext = -1;
                for (int u = 0; u < n; ++u) {
                    if (!visited[u]) {
                        double d = distance(points[current], points[u]);
                        if (d < globalBest) {
                            globalBest = d;
                            globalNext = u;
                        }
                    }
                }
                if (globalNext == -1) {
                    break;
                } else {
                    visited[globalNext] = true;
                    tour.push_back(globalNext);
                    current = globalNext;
                }
            } else {
                visited[next] = true;
                tour.push_back(next);
                current = next;
            }
        }

        for (int u = 0; u < n; u++) {
            if (!visited[u]) {
                visited[u] = true;
                tour.push_back(u);
            }
        }

        return tour;
    }

    void reorderPathForStartEnd(std::vector<int>& path, int startIndex,
                                int endIndex) const
    {
        auto itStart = std::find(path.begin(), path.end(), startIndex);
        if (itStart != path.end() && itStart != path.begin()) {
            std::rotate(path.begin(), itStart, path.end());
        }

        auto itEnd = std::find(path.begin(), path.end(), endIndex);
        if (itEnd != path.end() && (itEnd + 1) != path.end()) {
            std::rotate(itEnd + 1, path.end(), path.end());
        }
    }
};

std::tuple<std::vector<int>, double, double> solveTSP(
    const std::vector<Point>& points, const std::vector<std::vector<int>>& adj,
    int startIndex, int endIndex, double costPerUnit)
{
    Solution solver;
    return solver.computeTSP(points, adj, startIndex, endIndex, costPerUnit);
}
