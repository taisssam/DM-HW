#ifndef DM_EUROPE_GRAPH_H
#define DM_EUROPE_GRAPH_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <set>
#include <algorithm>
#include <map>

using namespace std;

const int INF = INT_MAX;

struct NodeInfo {
    int saturation;
    int degree;
    int vertex_num;
};

struct MaxSaturation {
    bool operator()(const NodeInfo& lhs, const NodeInfo& rhs) const {
        return tie(lhs.saturation, lhs.degree, lhs.vertex_num) > tie(rhs.saturation, rhs.degree, rhs.vertex_num);
    }
};

struct WeightedEdge {
    int from;
    int to;
    double weight;

    WeightedEdge(int from, int to, double weight) : from(from), to(to), weight(weight) {}

    bool operator<(const WeightedEdge& other) const {
        return weight < other.weight;
    }
};

class DisjointSetUnion {
private:
    vector<int> parent, rank;

public:
    DisjointSetUnion(int n) : parent(n+1), rank(n+1, 0) {
        for (int i = 0; i <= n; ++i) {
            parent[i] = i;
        }
    }

    int FindSet(int u) {
        if (u != parent[u]) parent[u] = FindSet(parent[u]);
        return parent[u];
    }

    void UnionSets(int u, int v) {
        u = FindSet(u);
        v = FindSet(v);
        if (u != v) {
            if (rank[u] < rank[v]) swap(u, v);
            parent[v] = u;
            if (rank[u] == rank[v]) ++rank[u];
        }
    }
};

class Graph {
private:
    int edge_count;
    int vertex_count;
    vector<vector<int>> adjacency_list;
    vector<vector<int>> adjacency_matrix;
    vector<pair<int, int>> edges_list;
    vector<pair<int, int>> mstEdges;
    map<int, vector<int>> adjacencyList;


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

    void BronKerbosch(set<int> R, set<int> P, set<int> X, vector<set<int>>& cliques);
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

    void AddWeightedEdge(int from, int to, double weight);
    void AddWeightedEdges();

    void KruskalMST();

    vector<int> generatePruferCode();
    vector<int> generateBinaryCode();
    void dfsForBinCode(int vertex, map<int, bool>& visited, vector<int>& ans);
};


#endif //DM_EUROPE_GRAPH_H
