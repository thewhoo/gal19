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

typedef adjacency_list<vecS, vecS, directedS, Vertex, Edge> graph_t;

graph_t getGraph(const std::string &filename);
bool getStartAndEndNodeIndexes(const graph_t &graph, const std::string &startNode, const std::string &endNode, int &startNodeIndex, int &endNodeIndex);

int main(int argc, char *argv[]) {
    std::string filename, startNode, endNode;
    int startNodeIndex{};
    int endNodeIndex{};

    if(argc != 4)
    {
        std::cerr << "Wrong number of arguments" << std::endl;
        exit(EXIT_FAILURE);
    }

    graph_t graph;

    filename = argv[1];
    startNode = argv[2];
    endNode = argv[3];

    graph = getGraph(filename);
    if(!getStartAndEndNodeIndexes(graph, startNode, endNode, startNodeIndex, endNodeIndex))
    {
        std::cerr << "Start or end node is not graph" << std::endl;
    }
    
    std::cout << "start node " << startNodeIndex << " end node " << endNodeIndex << std::endl;

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