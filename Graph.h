#ifndef DM_EUROPE_GRAPH_H
#define DM_EUROPE_GRAPH_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stack>
#include <set>
#include <algorithm>

using namespace std;

const int INF = INT_MAX;


struct NodeInfo {
    int saturation;
    int degree;
    int vertex_num;
};

struct MaxSaturation {
    bool operator()(const NodeInfo& lhs, const NodeInfo& rhs) const {
        return std::tie(lhs.saturation, lhs.degree, lhs.vertex_num) > std::tie(rhs.saturation, rhs.degree, rhs.vertex_num);
    }
};

class Graph {
private:
    int edge_count;
    int vertex_count;
    vector<vector<int>> adjacency_list;
    vector<vector<int>> adjacency_matrix;
    vector<pair<int, int>> edges_list;

    int components_count;

    int diam;
    int rad;
    vector<int> center;

    void ReadGraph(const string& filename);
    void CreateAdjacencyMatrix();
    void CreateEdgesList();
    int DFSComp();

    void DFS(int u, vector<int>& visited);

    void FloydWarshall(vector<vector<int>>& distance_matrix);

    void BronKerbosch(set<int> clique, set<int> candidates, set<int> excluded, const vector<set<int>>& adjacency_list);
    vector<set<int>> BuildAdjacencyList();

    int DFSVcomp(int v, int predok);


    int DFSEcomp(int v, int predok);

public:
    void Reset();
    Graph(const string& filename);

    void GetEdgeCount();
    void GetVertexCount();
    void DegreesCount();
    void CountCyclomaticNum();

    void EccentricitiesCount();

    void CountChromaticNum();

    void FindMaximalCliques();

    void Vcomp();
    void Ecomp();

    void KruskalMST();
};


#endif //DM_EUROPE_GRAPH_H
