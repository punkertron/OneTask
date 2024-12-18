#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <queue>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include "Solution.hpp"

constexpr int N = 100;
constexpr double RADIUS = 100.0;
constexpr double COST_PER_UNIT = 10.0;

static double distanceBetween(const Point& a, const Point& b)
{
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

std::vector<Point> generateRandomPointsInsideCircle(int n, double radius)
{
    std::vector<Point> points(n);
    static std::mt19937_64 rng(std::random_device{}());
    std::uniform_real_distribution<double> angleDist(0.0, 2.0 * std::numbers::pi);
    std::uniform_real_distribution<double> radiusDist(0.0, 1.0);

    for (int i = 0; i < n; ++i) {
        double t = angleDist(rng);
        double r = radius * std::sqrt(radiusDist(rng));
        points[i].x = r * std::cos(t);
        points[i].y = r * std::sin(t);
    }
    return points;
}

std::vector<int> findComponents(const std::vector<std::vector<int>>& adj)
{
    int n = adj.size();
    std::vector<int> comp(n, -1);
    int c = 0;
    for (int i = 0; i < n; ++i) {
        if (comp[i] == -1) {
            std::queue<int> q;
            q.push(i);
            comp[i] = c;
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                for (auto w : adj[u]) {
                    if (comp[w] == -1) {
                        comp[w] = c;
                        q.push(w);
                    }
                }
            }
            c++;
        }
    }
    return comp;
}

bool connectAllComponents(std::vector<std::vector<int>>& adj,
                          const std::vector<Point>& points, int minDeg = 2,
                          int maxDeg = 6)
{
    int n = adj.size();
    auto degree = std::vector<int>(n, 0);
    for (int i = 0; i < n; i++) {
        degree[i] = adj[i].size();
    }

    auto comp = findComponents(adj);
    int cCount = 0;
    for (auto x : comp) {
        if (x > cCount)
            cCount = x;
    }
    cCount++;

    if (cCount == 1)
        return true;

    std::vector<std::vector<int>> comps(cCount);
    for (int i = 0; i < n; i++) {
        comps[comp[i]].push_back(i);
    }

    for (int cid = 1; cid < cCount; cid++) {
        std::vector<int> leftSide, rightSide;
        for (int i = 0; i < cid; i++) {
            for (auto v : comps[i]) {
                if (degree[v] < maxDeg) {
                    leftSide.push_back(v);
                }
            }
        }
        for (auto v : comps[cid]) {
            if (degree[v] < maxDeg) {
                rightSide.push_back(v);
            }
        }

        if (leftSide.empty() || rightSide.empty()) {
            return false;
        }

        double bestDist = std::numeric_limits<double>::infinity();
        int uBest = -1, vBest = -1;
        for (auto u : leftSide) {
            for (auto w : rightSide) {
                double d = distanceBetween(points[u], points[w]);
                if (d < bestDist) {
                    bestDist = d;
                    uBest = u;
                    vBest = w;
                }
            }
        }
        adj[uBest].push_back(vBest);
        adj[vBest].push_back(uBest);
        degree[uBest]++;
        degree[vBest]++;
        comp = findComponents(adj);
    }

    comp = findComponents(adj);
    int finalC = *std::max_element(comp.begin(), comp.end());
    return finalC == 0;
}

std::vector<std::vector<int>> buildRestrictedDegreeGraph(const std::vector<Point>& points)
{
    int n = points.size();
    std::vector<std::vector<std::pair<double, int>>> neighbors(n);
    for (int i = 0; i < n; ++i) {
        neighbors[i].reserve(n - 1);
        for (int j = 0; j < n; ++j) {
            if (i == j)
                continue;
            double d = distanceBetween(points[i], points[j]);
            neighbors[i].push_back({d, j});
        }
        std::sort(neighbors[i].begin(), neighbors[i].end(), [](auto& a, auto& b) {
            return a.first < b.first;
        });
    }

    std::vector<std::vector<int>> adj(n);
    std::vector<int> degree(n, 0);

    for (int i = 0; i < n; ++i) {
        int added = 0;
        for (auto& pr : neighbors[i]) {
            int nb = pr.second;
            if (degree[i] < 6 && degree[nb] < 6) {
                if (std::find(adj[i].begin(), adj[i].end(), nb) == adj[i].end()) {
                    adj[i].push_back(nb);
                    adj[nb].push_back(i);
                    degree[i]++;
                    degree[nb]++;
                    added++;
                }
            }
            if (added == 2)
                break;
        }
    }

    for (int i = 0; i < n; ++i) {
        if (degree[i] < 2) {
            int idx = 0;
            while (degree[i] < 2 && idx < neighbors[i].size()) {
                int nb = neighbors[i][idx].second;
                if (nb != i && degree[i] < 6 && degree[nb] < 6 &&
                    std::find(adj[i].begin(), adj[i].end(), nb) == adj[i].end()) {
                    adj[i].push_back(nb);
                    adj[nb].push_back(i);
                    degree[i]++;
                    degree[nb]++;
                }
                idx++;
            }
        }
    }

    if (!connectAllComponents(adj, points, 2, 6)) {
        std::cerr << "Error in graph construction\n";
    }

    return adj;
}

int main()
{
    auto points = generateRandomPointsInsideCircle(N, RADIUS);

    auto adj = buildRestrictedDegreeGraph(points);

    int endIndex = 0;
    std::cout << "Enter endIndex:\n";
    std::cin >> endIndex;

    auto [path, dist, cost] = solveTSP(points, adj, 0, endIndex, COST_PER_UNIT);

    if (path.size() == N) {
        std::cout << "All nodes visited successfully.\n";
    } else {
        std::cout << "Not all nodes visited. Path size: " << path.size() << " out of "
                  << N << "\n";
    }

    std::cout << "Final path:\n";
    for (int i = 0; i < path.size(); ++i) {
        std::cout << path[i] << (i + 1 < path.size() ? " -> " : "\n");
    }
    std::cout << "Distance: " << dist << "\n";
    std::cout << "Total cost: " << cost << " USD\n";

    return 0;
}
