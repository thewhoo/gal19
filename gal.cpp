#include <vector>
#include <queue>

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
    double cost;

    Path() : cost{0.0} {}
};

// Define comparison operators for Path
bool operator >(const Path &p1, const Path &p2)
{
    return (p1.cost > p2.cost);
}

typedef adjacency_list<vecS, vecS, directedS, Vertex, Edge> graph_t;

graph_t getGraph(const std::string &filename);
bool getStartAndEndNodeIndexes(const graph_t &graph, const std::string &startNode, const std::string &endNode, int &startNodeIndex, int &endNodeIndex);
void findShortestPaths(const graph_t &graph, int startVertexIndex, int endVertexIndex, int pathCount, std::vector<Path> &shortestPaths);
void printShortestPaths(const std::vector<Path> &shortestPaths, const graph_t &graph);

int main(int argc, char *argv[]) {
    std::string filename, startNode, endNode;
    int startNodeIndex{};
    int endNodeIndex{};
    int pathCount;

    if(argc != 5)
    {
        std::cerr << "Wrong number of arguments" << std::endl;
        exit(EXIT_FAILURE);
    }

    graph_t graph;

    filename = argv[1];
    startNode = argv[2];
    endNode = argv[3];
    pathCount = std::stoi(argv[4]);


    graph = getGraph(filename);
    if(!getStartAndEndNodeIndexes(graph, startNode, endNode, startNodeIndex, endNodeIndex))
    {
        std::cerr << "Start or end node is not in graph" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<Path> shortestPaths;
    findShortestPaths(graph, startNodeIndex, endNodeIndex, pathCount, shortestPaths);
    printShortestPaths(shortestPaths, graph);
}

graph_t getGraph(const std::string &filename)
{
    // Construct an empty graph and prepare the dynamic_property_maps.
    graph_t graph(0);

    dynamic_properties dp;
    dp.property("node_id", get(&Vertex::name, graph));
    dp.property("weight", get(&Edge::weight, graph));

    std::ifstream dot(filename);

    if(!read_graphviz(dot, graph, dp)) 
    {
        std::cerr << "Can't read input graph" << std::endl;
        exit(EXIT_FAILURE);
    }

    write_graphviz_dp(std::cout, graph, dp);

    return graph;
}

bool getStartAndEndNodeIndexes(const graph_t &graph, const std::string &startNode, const std::string &endNode, int &startNodeIndex, int &endNodeIndex)
{
    int index = 0;
    bool startNodeFound, endNodeFound;
    startNodeFound = endNodeFound = false;
    for(const auto &vertex : graph.m_vertices)
    {
        if(vertex.m_property.name == startNode)
        {
            startNodeIndex = index;
            startNodeFound = true;
        }
        else if(vertex.m_property.name == endNode)
        {
            endNodeIndex = index;
            endNodeFound = true;
        }
        index++;
        if(startNodeFound && endNodeFound)
        {
            break;
        }
    }

    return (startNodeFound && endNodeFound);
}

void findShortestPaths(const graph_t &graph, int startVertexIndex, int endVertexIndex, int pathCount, std::vector<Path> &shortestPaths)
{
    // Allocate shortest path priority queue
    std::priority_queue<Path, std::vector<Path>, std::greater<>> pathQueue;
    // Allocate per-node path count vector
    std::vector<int> pathCounts(graph.m_vertices.size(), 0);

    shortestPaths.clear();

    // Insert initial path
    auto initialPath = Path{};
    initialPath.vertices.push_back(startVertexIndex);
    pathQueue.push(initialPath);

    while(!pathQueue.empty() && pathCounts[endVertexIndex] < pathCount) {
        Path currentShortestPath{pathQueue.top()};
        pathQueue.pop();

        auto pathEndVertex = currentShortestPath.vertices.back();
        pathCounts[pathEndVertex]++;

        if(pathEndVertex == endVertexIndex) {
            shortestPaths.push_back(currentShortestPath);
        }

        if(pathCounts[pathEndVertex] <= pathCount) {
            for(const auto &edge : graph.out_edge_list(pathEndVertex)) {
                auto targetIndex = edge.m_target;
                auto &&newPath = Path(currentShortestPath);
                newPath.vertices.push_back(targetIndex);
                newPath.cost += edge.get_property().weight;
                pathQueue.push(newPath);
            }
        }
    }
}

void printShortestPaths(const std::vector<Path> &shortestPaths, const graph_t &graph)
{
    for(const auto &path : shortestPaths) {
        std::cout << "Cost: " << path.cost << "\tPath: ";
        int loopIndex = 0; // To check for last vertex in path (cannot check by value in case of loopy paths)
        for(const auto &vertexIndex : path.vertices) {
            std::cout << graph[vertexIndex].name;
            if(loopIndex != path.vertices.size() - 1)
                std::cout << " -> ";

            loopIndex++;
        }
        std::cout << std::endl;
    }
}