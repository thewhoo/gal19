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

void findShortestPathsParallel(const graph_t &graph, int startVertexIndex, int endVertexIndex, int pathCount,
                               std::vector<Path> &shortestPaths) {
    auto threadCount = omp_get_max_threads();
    // Allocate shortest path priority queue
    //std::vector<Path> currentPaths;
    //currentPaths.reserve(10000);
    //std::priority_queue<Path, std::vector<Path>, std::greater<>> currentPaths;
    std::deque<Path> currentPaths;
    std::vector<std::vector<Path>> prospectivePaths(threadCount);
    std::vector<double> prospectivePathMinCosts(threadCount);
    // Allocate per-node path count vector
    std::vector<int> pathCounts(graph.m_vertices.size(), 0);

    shortestPaths.clear();

    // Insert initial path
    auto initialPath = Path{};
    initialPath.vertices.push_back(startVertexIndex);
    initialPath.costs.push_back(0.0);
    //currentPaths.push_back(initialPath);
    currentPaths.push_front(initialPath);

    std::cout << "thr count: " << threadCount << std::endl;

#pragma omp parallel default(none) shared(graph, currentPaths, pathCounts, endVertexIndex, pathCount, prospectivePaths, prospectivePathMinCosts, shortestPaths, threadCount, std::cout)
    while (!currentPaths.empty() && pathCounts[endVertexIndex] < pathCount) {
        auto threadId = omp_get_thread_num();

        if (threadId < currentPaths.size()) {
            auto precalcPath = currentPaths[threadId];
            auto pathEndVertex = precalcPath.vertices.back();

            std::vector<Path> generatedPaths;
            double minCost = std::numeric_limits<double>::max();
            if (pathCounts[pathEndVertex] <= pathCount) {
                for (const auto &edge : graph.out_edge_list(pathEndVertex)) {
                    auto targetIndex = edge.m_target;
                    auto &&newPath = Path(precalcPath);

                    newPath.vertices.push_back(targetIndex);

                    auto edgeCost = edge.get_property().weight;
                    newPath.costs.push_back(edgeCost);
                    newPath.totalCost += edgeCost;

                    if (newPath.totalCost < minCost)
                        minCost = newPath.totalCost;

                    generatedPaths.push_back(newPath);
                }
            }

            //prospectivePaths.emplace(prospectivePaths.begin() + threadId, std::move(generatedPaths));
            prospectivePaths[threadId] = std::move(generatedPaths);
            prospectivePathMinCosts[threadId] = minCost;
        }
#pragma omp barrier

#pragma omp single
        {
            auto precalcResultCount = std::min<int>(threadCount, currentPaths.size());
            for (auto i = 0; i < precalcResultCount; i++) {
                if (i == 0 || currentPaths[0].totalCost < prospectivePathMinCosts[i]) {
                    //if(true) {
                    Path currentShortestPath{currentPaths[0]};
                    currentPaths.erase(currentPaths.begin());
                    //currentPaths.pop();

                    auto pathEndVertex = currentShortestPath.vertices.back();
                    pathCounts[pathEndVertex]++;

                    if (pathEndVertex == endVertexIndex) {
                        shortestPaths.push_back(currentShortestPath);
                    }

                    //currentPaths.insert(currentPaths.end(), prospectivePaths[i].begin(), prospectivePaths[i].end());
                    //std::sort(currentPaths.begin(), currentPaths.end());
                    for (auto &p : prospectivePaths[i]) {
                        //currentPaths.push_back(std::move(p));
                        fastInsert(currentPaths, std::move(p));
                    }
                    //std::sort(currentPaths.begin(), currentPaths.end());
                } else
                    break;
            }
            //std::sort(currentPaths.begin(), currentPaths.end());
        };
    }
}

void findShortestPathsLoopless(const graph_t &graph, int startVertexIndex, int endVertexIndex, int pathCount,
                               std::vector<Path> &shortestPaths) {
    shortestPaths.clear();

    // Priority queue of candidate best paths
    std::priority_queue<Path, std::vector<Path>, std::greater<>> candidateQueue;

    // Get initial shortest path from start to end
    findShortestPaths(graph, startVertexIndex, endVertexIndex, 1, shortestPaths);
    if (shortestPaths.empty())
        return;

    for (int i = 1; i < pathCount; i++) {
        //const auto &previousPath = shortestPaths[i - 1];
        const auto &previousPath = shortestPaths.back();
        // Spur node from first node to next-to-last node in previous path
        for (int j = 0; j < previousPath.vertices.size() - 2; j++) {

            auto rootPath = previousPath;
            // Get subpath of previous shortest path from root node to current spur node
            auto spurNode = rootPath.vertices[j];
            rootPath.vertices.erase(rootPath.vertices.begin() + j + 1, rootPath.vertices.end());
            rootPath.costs.erase(rootPath.costs.begin() + j + 1, rootPath.costs.end());
            // Recalculate subpath cost
            rootPath.totalCost = std::accumulate(rootPath.costs.begin(), rootPath.costs.end(), 0.0);

            std::cout << "i = " << i << " j = " << j << " spur: " << spurNode << " nsp size: "
                      << rootPath.vertices.size() << " new sp cost: " << rootPath.totalCost << std::endl;


            // Copy of graph should be O(n)
            auto graphCopy = graph;

            // Disconnect all previous shortest paths that share the same root path
            for (const auto &path : shortestPaths) {
                auto currentSubpath = path.vertices;
                currentSubpath.erase(currentSubpath.begin() + j + 1, currentSubpath.end());
                if (currentSubpath == rootPath.vertices) {
                    remove_edge(path.vertices[j], path.vertices[j + 1], graphCopy);
                }
            }

            // Remove disconnected nodes
            for (const auto &vertex : rootPath.vertices) {
                if (vertex != spurNode) {
                    //remove_edge(graphCopy.out_edge_list(vertex), graphCopy);
                    for (int v = 0; i < num_vertices(graphCopy); i++) {
                        remove_edge(vertex, v, graphCopy);
                        remove_edge(v, vertex, graphCopy);
                    }
                }
            }

            //int indexOffset = 0;
            //for(const auto &vertex : rootPath.vertices) {
            //    if(vertex != spurNode) {
            //        remove_vertex(vertex, graphCopy);
            //        indexOffset++;
            //    }
            //}

            std::vector<Path> spurPathVec;
            //std::cout << "offset: " << indexOffset << " fsp spur: " << spurNode - indexOffset << " fsp end: " <<endVertexIndex - indexOffset << std::endl;
            findShortestPaths(graphCopy, spurNode, endVertexIndex, 1, spurPathVec);
            auto spurPath = spurPathVec.front();

            auto itVertices = spurPath.vertices.begin();
            auto itCosts = spurPath.costs.begin();
            while (itVertices != spurPath.vertices.end() && itCosts != spurPath.costs.end()) {
                rootPath.vertices.push_back((*itVertices));
                rootPath.costs.push_back(*itCosts);
                itVertices++;
                itCosts++;
            }
            rootPath.totalCost += spurPath.totalCost;

            candidateQueue.push(rootPath);
        }

        if (candidateQueue.empty())
            break;

        shortestPaths.push_back(candidateQueue.top());
        candidateQueue.pop();
    }
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