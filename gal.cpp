#include <vector>
#include <queue>
#include <numeric>
#include <functional>
#include <chrono>

#include <omp.h>

#include <boost/range.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/graph/graph_utility.hpp>

using namespace boost;

struct Vertex {
    std::string name;
};

struct Edge {
    double weight;
};

struct Path {
    std::vector<int> vertices;
    std::vector<double> costs;

    double totalCost;

    Path() : totalCost{0.0} {}
};

// Define comparison operators for Path
bool operator>(const Path &p1, const Path &p2) {
    return (p1.totalCost > p2.totalCost);
}

bool operator<(const Path &p1, const Path &p2) {
    return (p1.totalCost < p2.totalCost);
}

typedef adjacency_list<vecS, vecS, directedS, Vertex, Edge> graph_t;

graph_t getGraph(const std::string &filename);

bool getStartAndEndNodeIndexes(const graph_t &graph, const std::string &startNode, const std::string &endNode,
                               int &startNodeIndex, int &endNodeIndex);

void findShortestPaths(const graph_t &graph, int startVertexIndex, int endVertexIndex, int pathCount,
                       std::vector<Path> &shortestPaths);

void findShortestPathsLoopless(const graph_t &graph, int startVertexIndex, int endVertexIndex, int pathCount,
                               std::vector<Path> &shortestPaths);

void printShortestPaths(const std::vector<Path> &shortestPaths, const graph_t &graph);

void
printPathQueue(std::priority_queue<Path, std::vector<Path>, std::greater<>> q, const graph_t &graph, int pathCount);

void findShortestPathsParallel(const graph_t &graph, int startVertexIndex, int endVertexIndex, int pathCount,
                               std::vector<Path> &shortestPaths);

void findShortestPathsHighIQ(const graph_t &graph, int startVertexIndex, int endVertexIndex, int pathCount,
                             std::priority_queue<Path, std::vector<Path>, std::greater<>> &accumulatedShortestPaths);

template<class T>
void fastInsert(std::deque<T> &q, T &&item);

int main(int argc, char *argv[]) {
    std::string filename, startNode, endNode;
    int startNodeIndex{};
    int endNodeIndex{};
    int pathCount;

    if (argc != 6) {
        std::cerr << "./gal19 graph start_node end_node path_count parallel_threads" << std::endl
                  << "\tgraph - path to graph file in .dot format" << std::endl
                  << "\tstart_node - identifier of start node" << std::endl
                  << "\tend_node - identifier of end node" << std::endl
                  << "\tpath_count - maximum number of shortest paths to find" << std::endl
                  << "\tparallel_threads - number of parallel threads to run (0 will run the reference sequential version)"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    graph_t graph;

    filename = argv[1];
    startNode = argv[2];
    endNode = argv[3];
    pathCount = std::stoi(argv[4]);
    auto threadNum = std::stoi(argv[5]);

    graph = getGraph(filename);
    if (!getStartAndEndNodeIndexes(graph, startNode, endNode, startNodeIndex, endNodeIndex)) {
        std::cerr << "Start or end node is not in graph" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (threadNum > 0) {
        omp_set_num_threads(threadNum);

        std::priority_queue<Path, std::vector<Path>, std::greater<>> accumulatedShortestPaths;
        auto start = std::chrono::high_resolution_clock::now();
        findShortestPathsHighIQ(graph, startNodeIndex, endNodeIndex, pathCount, accumulatedShortestPaths);
        auto stop = std::chrono::high_resolution_clock::now();
        printPathQueue(accumulatedShortestPaths, graph, pathCount);
        std::cout << std::endl << "Runtime: "
                  << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << " usec"
                  << std::endl;
    } else {
        std::vector<Path> shortestPaths;
        auto start = std::chrono::high_resolution_clock::now();
        findShortestPaths(graph, startNodeIndex, endNodeIndex, pathCount, shortestPaths);
        auto stop = std::chrono::high_resolution_clock::now();
        printShortestPaths(shortestPaths, graph);
        std::cout << std::endl << "Runtime: "
                  << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << " usec"
                  << std::endl;
    }
}

graph_t getGraph(const std::string &filename) {
    // Construct an empty graph and prepare the dynamic_property_maps.
    graph_t graph(0);

    dynamic_properties dp;
    dp.property("node_id", get(&Vertex::name, graph));
    dp.property("weight", get(&Edge::weight, graph));

    std::ifstream dot(filename);

    if (!read_graphviz(dot, graph, dp)) {
        std::cerr << "Can't read input graph" << std::endl;
        exit(EXIT_FAILURE);
    }

    write_graphviz_dp(std::cout, graph, dp);

    return graph;
}

bool getStartAndEndNodeIndexes(const graph_t &graph, const std::string &startNode, const std::string &endNode,
                               int &startNodeIndex, int &endNodeIndex) {
    int index = 0;
    bool startNodeFound, endNodeFound;
    startNodeFound = endNodeFound = false;
    for (const auto &vertex : graph.m_vertices) {
        if (vertex.m_property.name == startNode) {
            startNodeIndex = index;
            startNodeFound = true;
        } else if (vertex.m_property.name == endNode) {
            endNodeIndex = index;
            endNodeFound = true;
        }
        index++;
        if (startNodeFound && endNodeFound) {
            break;
        }
    }

    return (startNodeFound && endNodeFound);
}

void findShortestPaths(const graph_t &graph, int startVertexIndex, int endVertexIndex, int pathCount,
                       std::vector<Path> &shortestPaths) {
    // Allocate shortest path priority queue
    std::priority_queue<Path, std::vector<Path>, std::greater<>> pathQueue;
    //std::deque<Path> pathQueue;
    // Allocate per-node path count vector
    std::vector<int> pathCounts(graph.m_vertices.size(), 0);

    shortestPaths.clear();

    // Insert initial path
    auto initialPath = Path{};
    initialPath.vertices.push_back(startVertexIndex);
    initialPath.costs.push_back(0.0);
    pathQueue.push(initialPath);
    //pathQueue.push_back(initialPath);

    while (!pathQueue.empty() && pathCounts[endVertexIndex] < pathCount) {
        Path currentShortestPath{pathQueue.top()};
        pathQueue.pop();
        //Path currentShortestPath{pathQueue[0]};
        //pathQueue.erase(pathQueue.begin());

        auto pathEndVertex = currentShortestPath.vertices.back();
        pathCounts[pathEndVertex]++;

        if (pathEndVertex == endVertexIndex) {
            shortestPaths.push_back(currentShortestPath);
        }

        if (pathCounts[pathEndVertex] <= pathCount) {
            for (const auto &edge : graph.out_edge_list(pathEndVertex)) {
                auto targetIndex = edge.m_target;
                auto newPath = Path(currentShortestPath);
                newPath.vertices.push_back(targetIndex);

                auto edgeCost = edge.get_property().weight;
                newPath.costs.push_back(edgeCost);
                newPath.totalCost += edgeCost;

                pathQueue.push(newPath);
                //pathQueue.push_back(newPath);
                //fastInsert(pathQueue, std::move(newPath));
            }
            //std::sort(pathQueue.begin(), pathQueue.end());
        }
    }
}

void findShortestPathsHighIQ(const graph_t &graph, int startVertexIndex, int endVertexIndex, int pathCount,
                             std::priority_queue<Path, std::vector<Path>, std::greater<>> &accumulatedShortestPaths) {
    auto threadCount = omp_get_max_threads();

    // Allocate shortest path priority queue
    std::priority_queue<Path, std::vector<Path>, std::greater<>> pathQueue;
    // Allocate per-node path count vector
    std::vector<int> pathCounts(graph.m_vertices.size(), 0);

    // Insert initial path
    auto initialPath = Path{};
    initialPath.vertices.push_back(startVertexIndex);
    initialPath.costs.push_back(0.0);
    pathQueue.push(initialPath);

    while (!pathQueue.empty() && pathCounts[endVertexIndex] < pathCount && pathQueue.size() < threadCount) {
        Path currentShortestPath{pathQueue.top()};
        pathQueue.pop();

        auto pathEndVertex = currentShortestPath.vertices.back();
        pathCounts[pathEndVertex]++;

        if (pathEndVertex == endVertexIndex) {
            accumulatedShortestPaths.push(currentShortestPath);
        }

        if (pathCounts[pathEndVertex] <= pathCount) {
            for (const auto &edge : graph.out_edge_list(pathEndVertex)) {
                auto targetIndex = edge.m_target;
                auto newPath = Path(currentShortestPath);
                newPath.vertices.push_back(targetIndex);

                auto edgeCost = edge.get_property().weight;
                newPath.costs.push_back(edgeCost);
                newPath.totalCost += edgeCost;

                pathQueue.push(newPath);
            }
        }
    }

    std::vector<double> threadLastPathCosts(threadCount);
    double sharedMaxCost = 0.0;
    bool stop = false;
#pragma omp parallel default(none) firstprivate(pathQueue, pathCounts) shared(stop, threadCount, sharedMaxCost, threadLastPathCosts, graph, startVertexIndex, endVertexIndex, pathCount, accumulatedShortestPaths)
    {
        auto threadId = omp_get_thread_num();
        auto itemCount = pathQueue.size() / threadCount;
        auto itemModulo = pathQueue.size() % threadCount;

        auto removeCount = itemCount * threadId + std::min<int>(threadId, itemModulo);

        for (int i = 0; i < removeCount; i++)
            pathQueue.pop();

        std::priority_queue<Path, std::vector<Path>, std::greater<>> localQueue;
        itemCount = threadId < itemModulo ? itemCount + 1 : itemCount;
        for (int i = 0; i < itemCount; i++) {
            localQueue.push(pathQueue.top());
            pathQueue.pop();
        }

        while (!localQueue.empty() && pathCounts[endVertexIndex] < pathCount && !stop) {
            Path currentShortestPath{localQueue.top()};
            localQueue.pop();

            auto pathEndVertex = currentShortestPath.vertices.back();
            pathCounts[pathEndVertex]++;

            if (pathEndVertex == endVertexIndex) {
#pragma omp critical
                {
                    threadLastPathCosts[threadId] = currentShortestPath.totalCost;
                    if (currentShortestPath.totalCost > sharedMaxCost)
                        sharedMaxCost = currentShortestPath.totalCost;
                    accumulatedShortestPaths.push(currentShortestPath);
                };
            }

            if (pathCounts[pathEndVertex] <= pathCount) {
                for (const auto &edge : graph.out_edge_list(pathEndVertex)) {
                    auto targetIndex = edge.m_target;
                    auto newPath = Path(currentShortestPath);
                    newPath.vertices.push_back(targetIndex);

                    auto edgeCost = edge.get_property().weight;
                    newPath.costs.push_back(edgeCost);
                    newPath.totalCost += edgeCost;

                    localQueue.push(newPath);
                }
            }

#pragma omp master
            {
                bool allCostsGreater = true;
                for (const auto &cost : threadLastPathCosts) {
                    if (cost < sharedMaxCost) {
                        allCostsGreater = false;
                        break;
                    }
                }

                if (allCostsGreater && accumulatedShortestPaths.size() >= pathCount)
                    stop = true;
            };
        }
    };
}

void printShortestPaths(const std::vector<Path> &shortestPaths, const graph_t &graph) {
    for (const auto &path : shortestPaths) {
        std::cout << "Cost: " << path.totalCost << "\tPath: ";
        int loopIndex = 0; // To check for last vertex in path (cannot check by value in case of loopy paths)
        for (const auto &vertexIndex : path.vertices) {
            std::cout << graph[vertexIndex].name;
            if (loopIndex != path.vertices.size() - 1)
                std::cout << " -> ";

            loopIndex++;
        }
        std::cout << std::endl;
    }
}

void
printPathQueue(std::priority_queue<Path, std::vector<Path>, std::greater<>> q, const graph_t &graph, int pathCount) {
    auto end = std::min<int>(q.size(), pathCount);
    for (int i = 0; i < end; i++) {
        auto path = q.top();
        q.pop();
        std::cout << "Cost: " << path.totalCost << "\tPath: ";
        int loopIndex = 0; // To check for last vertex in path (cannot check by value in case of loopy paths)
        for (const auto &vertexIndex : path.vertices) {
            std::cout << graph[vertexIndex].name;
            if (loopIndex != path.vertices.size() - 1)
                std::cout << " -> ";

            loopIndex++;
        }
        std::cout << std::endl;
    }
}

template<class T>
void fastInsert(std::deque<T> &q, T &&item) {
    for (auto it = q.rbegin(); it != q.rend(); it++) {
        if (item > *it) {
            q.insert(it.base(), item);
            return;
        }
    }

    q.push_back(item);
}