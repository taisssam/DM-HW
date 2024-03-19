#include "Graph.h"


vector<int> visited;
vector<int> time_enter;
int timer;

vector<vector<pair<int, int>>> incident_edges;
vector<int> art_point;
vector<pair<int, int>> bridges;

Graph::Graph(const string& filename){
    edge_count = 0;
    vertex_count = 0;
    ReadGraph(filename);
    CreateAdjacencyMatrix();
    CreateEdgesList();

    diam = 0;
    rad = 0;

}

void Graph::ReadGraph(const string& filename) {
    ifstream file(filename);
    string line;
    int country_id = 0;

    while (getline(file, line)) {
        vertex_count++;
        vector<int> neighbors;
        int pos = line.find(':');

        if (pos != string::npos) {
            string neighbors_str = line.substr(pos + 1);
            int start = 0;
            int end = neighbors_str.find(' ');

            while (end != string::npos) {
                string neighbor_str = neighbors_str.substr(start, end - start);
                if (!neighbor_str.empty()) {
                    neighbors.push_back(stol(neighbor_str) - 1);
                    ++edge_count;
                }

                start = end + 1;
                end = neighbors_str.find(' ', start);
            }
            string last_neighbor_str = neighbors_str.substr(start);

            if (!last_neighbor_str.empty()) {
                neighbors.push_back(stol(last_neighbor_str) - 1);
                ++edge_count;
            }
        }
        adjacency_list.push_back(neighbors);
        ++country_id;
    }
}


void Graph::GetEdgeCount() {
    cout << "Edges count: " << edge_count / 2;
}

void Graph::GetVertexCount() {
    cout << "Vertexes count: " << vertex_count;
}


void Graph::Reset() {
    visited.clear();
    time_enter.clear();
    timer = 0;
}

void Graph::CreateAdjacencyMatrix() {
    adjacency_matrix.resize(adjacency_list.size(), vector<int>(adjacency_list.size(), 0));

    for (int i = 0; i < adjacency_list.size(); ++i) {
        const vector<int>& neighbors = adjacency_list[i];
        for (int j = 0; j < neighbors.size(); ++j) {
            int to = neighbors[j];
            adjacency_matrix[i][to] = 1;
        }
    }
}

void Graph::DegreesCount() {
    vector<int> degrees(vertex_count, 0);

    for (int i = 0; i < vertex_count; ++i) {
        for (int j = 0; j < vertex_count; ++j) {
            degrees[i] += adjacency_matrix[i][j];
        }
    }

    int min_degree = INF;
    int max_degree = 0;

    for (int i = 0; i < degrees.size(); ++i) {
        if (degrees[i] > max_degree) {
            max_degree = degrees[i];
        }
        if (degrees[i] < min_degree) {
            if (degrees[i] != 0) {
                min_degree = degrees[i];
            }

        }
    }

    cout << "Max degree: " << max_degree << '\n';
    cout << "Min degree: " << min_degree;
}

void Graph::FloydWarshall(vector<vector<int>>& distance_matrix) {
    int n = adjacency_matrix.size();

    distance_matrix.resize(n, vector<int>(n, INF));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (adjacency_matrix[i][j] == 1) {
                distance_matrix[i][j] = 1;
                distance_matrix[j][i] = 1;
            }
        }
    }

    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (distance_matrix[i][k] != INF && distance_matrix[k][j] != INF &&
                    distance_matrix[i][k] + distance_matrix[k][j] < distance_matrix[i][j]) {
                    distance_matrix[i][j] = distance_matrix[i][k] + distance_matrix[k][j];
                }
            }
        }
    }

    cout << "Distance matrix: " << '\n';

    for (int i = 0; i < n; ++i) {
        distance_matrix[i][i] = 0;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {

            if (distance_matrix[i][j] == INF) {
                cout << '-' << ' ';
            } else {
                cout << distance_matrix[i][j] << " ";
            }

        }
        cout << '\n';
    }

    for (int i = 0; i < n; ++i) {
        distance_matrix[i][i] = 0;
        for (int j = 0; j < n; ++j) {
            if (distance_matrix[i][j] == INF) {
                distance_matrix[i][j] = 0;
                distance_matrix[j][i] = 0;
            }

        }
    }

}

void Graph::EccentricitiesCount() {

    int n = vertex_count;
    vector<vector<int>> distance_matrix(n, vector<int>(n, INF));
    FloydWarshall(distance_matrix);


    vector<int> eccentricities(vertex_count, 0);
    diam = 0;
    rad = INF;

    for (int i = 0; i < n; ++i) {
        int max_distance = 0;
        for (int j = 0; j < n; ++j) {
            if (distance_matrix[i][j] > max_distance) {
                max_distance = distance_matrix[i][j];
            }
        }
        eccentricities[i] = max_distance;
    }

    for (int i = 0; i < n; ++i) {
        if (eccentricities[i] > diam) {
            diam = eccentricities[i];
        }
        if (eccentricities[i] < rad) {
            if (eccentricities [i] > 1) {
                rad = eccentricities[i];
            }
        }

    }
    cout << "Diameter: " << diam << ", " << "Radius: " << rad << '\n';

    for (int i = 0; i < n; ++i) {
        if (eccentricities[i] == rad) {
            center.push_back(i + 1);
        }
    }

    cout << "Center: ";
    for (int i = 0; i < center.size(); ++i) {
        cout << center[i] << ' ';
    }
}


void Graph::DFS(int v, vector<int>& visited) {
    visited[v] = 1;

    for (int u : adjacency_list[v]) {
        if (!visited[u]) {
            DFS(u, visited);
        }
    }
}

int Graph::DFSComp() {
    visited.resize(vertex_count, 0);

    for (int i = 0; i < vertex_count; ++i) {
        if (!visited[i]) {
            DFS(i, visited);
            components_count++;
        }
    }
    return components_count;
}

void Graph::CountCyclomaticNum() {
    components_count = DFSComp();
    int vertices_count = vertex_count;
    int cyclomatic_number = edge_count / 2 - vertices_count + components_count;
    cout << "Cyclomatic number: " << cyclomatic_number << endl;
}

void Graph::BronKerbosch(set<int> clique, set<int> candidates, set<int> excluded, const vector<set<int>>& adjacency_list) {
    if (candidates.empty() && excluded.empty()) {
        bool maximal = true;
        for (int vertex : clique) {
            for (int neighbor : clique) {
                if (vertex != neighbor && adjacency_list[vertex].find(neighbor) == adjacency_list[vertex].end()) {
                    maximal = false;
                    break;
                }
            }
            if (!maximal) break;
        }
        if (maximal) {
            for (int vertex : clique) {
                cout << vertex << ' ';
            }
            cout << '\n';
        }
        return;
    }

    set<int> candidates_copy = candidates;
    for (int vertex : candidates_copy) {
        set<int> new_clique = clique;
        new_clique.insert(vertex);

        set<int> new_candidates;
        for (int neighbor : adjacency_list[vertex]) {
            if (candidates.find(neighbor) != candidates.end()) {
                new_candidates.insert(neighbor);
            }
        }

        set<int> new_excluded;
        for (int neighbor : adjacency_list[vertex]) {
            if (excluded.find(neighbor) != excluded.end()) {
                new_excluded.insert(neighbor);
            }
        }

        BronKerbosch(new_clique, new_candidates, new_excluded, adjacency_list);
        candidates.erase(vertex);
        excluded.insert(vertex);
    }
}

vector<set<int>> Graph::BuildAdjacencyList() {
    vector<set<int>> adjacency_list(adjacency_matrix.size());
    for (int i = 0; i < adjacency_matrix.size(); ++i) {
        for (int j = 0; j < adjacency_matrix[i].size(); ++j) {
            if (adjacency_matrix[i][j] == 1) {
                adjacency_list[i].insert(j);
                adjacency_list[j].insert(i);
            }
        }
    }
    return adjacency_list;
}



void Graph::FindMaximalCliques() {
    vector<set<int>> adjacency_list = BuildAdjacencyList();
    set<int> clique, candidates, excluded;
    for (int i = 0; i < adjacency_matrix.size(); ++i) {
        candidates.insert(i);
    }
    BronKerbosch(clique, candidates, excluded, adjacency_list);
}

void Graph::Vcomp() {
    incident_edges.resize(vertex_count);

    for (const auto& edge : edges_list) {
        int u = edge.first;
        int v = edge.second;
        if (u != v) {
            incident_edges[u].push_back(edge);
            incident_edges[v].push_back(edge);
        }
    }

    visited.resize(vertex_count, 0);
    time_enter.resize(vertex_count, 0);
    timer = 0;

    DFSVcomp(0, -1);

    cout << "Articulation points: ";
    for (int i = 0; i < art_point.size(); ++i) {
        if (i != art_point.size() - 1) {
            cout << art_point[i] + 1 << ", ";
        } else {
            cout << art_point[i] + 1;
        }

    }
}

int Graph::DFSVcomp(int v, int predok) {
    time_enter[v] = timer++;
    int min_time = time_enter[v];
    int child_count = 0;

    for (auto p : incident_edges[v]) {
        int u = p.first;
        if (u != predok) {
            if (time_enter[u] == 0) {
                int t = DFSVcomp(u, v);
                if (t >= time_enter[v]) {
                    ++child_count;
                }
                min_time = min(min_time, t);
            } else {
                min_time = min(min_time, time_enter[u]);
            }
        }
    }

    if ((predok != -1 && child_count > 0) || (predok == -1 && child_count > 1)) {
        if (min_time < time_enter[v]) {
            art_point.push_back(v);
        }
    }

    return min_time;
}

void Graph::CreateEdgesList() {
    for (int u = 0; u < adjacency_list.size(); ++u) {
        for (int v : adjacency_list[u]) {
            edges_list.push_back(make_pair(u, v));
        }
    }
}

int Graph::DFSEcomp(int v, int predok) {
    time_enter[v] = timer++;
    int min_time = time_enter[v];

    for (auto u : adjacency_list[v]) {
        if (u != predok) {
            int t;
            if (time_enter[u] == 0) {
                t = DFSEcomp(u, v);
                min_time = min(min_time, t);
                if (t > time_enter[v]) {
                    bridges.push_back(make_pair(min(v, u) + 1, max(v, u) + 1));
                }
            } else {
                min_time = min(min_time, time_enter[u]);
            }
        }
    }
    return min_time;
}

void Graph::Ecomp() {
    visited.resize(vertex_count, 0);
    time_enter.resize(vertex_count, 0);
    timer = 0;

    for (int i = 0; i < vertex_count; ++i) {
        if (!visited[i]) {
            DFSEcomp(i, -1);
        }

    }

    cout << "Bridges: ";
    for (auto bridge: bridges) {
        cout << '(' << bridge.first << ", " << bridge.second  << ')' << ", ";
    }
}

void Graph::CountChromaticNum() {
    vector<int> colors(vertex_count + 1, -1);
    vector<int> d(vertex_count + 1);
    vector<set<int>> adj_colors(vertex_count + 1, set<int>());
    vector<bool> visited(vertex_count + 1, false);
    set<NodeInfo, MaxSaturation> Q;
    set<NodeInfo, MaxSaturation>::iterator max_ptr;

    int u = 0, i = 0;

    for (u = 1; u <= vertex_count; ++u) {
        d[u] = adjacency_list[u].size();
        Q.emplace(NodeInfo{0, d[u], u});
    }
    while (!Q.empty()) {
        max_ptr = Q.begin();
        u = (*max_ptr).vertex_num;
        Q.erase(max_ptr);

        for (auto v : adjacency_list[u]) {
            if (colors[v] != -1) {
                visited[colors[v]] = true;
            }
        }

        for (i = 1; i <= visited.size(); ++i) {
            if (!visited[i]) {
                break;
            }
        }

        for (auto v : adjacency_list[u]) {
            if (colors[v] != -1) {
                visited[colors[v]] = false;
            }
        }

        colors[u] = i;

        for (auto v : adjacency_list[u]) {
            if (colors[v] == -1) {
                Q.erase({static_cast<int>(adj_colors[v].size()), d[v], v});
                adj_colors[v].insert(i);
                --d[v];
                Q.emplace(NodeInfo {static_cast<int>(adj_colors[v].size()), d[v], v});
            }
        }
    }
    set<int> uniq_colors;
    for (int j = 1; j <= vertex_count; ++j) {
        uniq_colors.insert(colors[j]);
    }
    cout << "Chromatic number: " << uniq_colors.size();
}


struct Edge {
    int u, v;
    double weight;

    Edge(int u, int v, double weight) : u(u), v(v), weight(weight) {}

    // Оператор для сравнения ребер по весу
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

class DisjointSetUnion {
private:
    vector<int> parents;
    vector<int> rank;

public:
    DisjointSetUnion(int n) {
        parents.resize(n + 1);
        rank.resize(n + 1);

        for (int i = 1; i <= n; ++i) {
            parents[i] = i;
            rank[i] = 0;
        }
    }

    int FindSet(int u) {
        if (u != parents[u]) {
            parents[u] = FindSet(parents[u]);
        }
        return parents[u];
    }

    void UnionSets(int u, int v) {
        u = FindSet(u);
        v = FindSet(v);

        if (u != v) {
            if (rank[u] < rank[v]) {
                swap(u, v);
            }
            parents[v] = u;
            if (rank[u] == rank[v]) {
                ++rank[u];
            }
        }
    }
};

//void Graph::KruskalMST() {
//    vector<Edge> sorted_edges;
//
//
//    for (int u = 0; u < adjacency_list.size(); ++u) {
//        for (Edge edge : adjacency_list[u]) {
//            sorted_edges.push_back(edge);
//        }
//    }
//
//
//    sort(sorted_edges.begin(), sorted_edges.end());
//
//    DisjointSetUnion dsu(vertex_count);
//
//    vector<Edge> result_edges;
//
//    for (const auto& edge : sorted_edges) {
//        int u = edge.u;
//        int v = edge.v;
//        double weight = edge.weight;
//
//        if (dsu.FindSet(u) != dsu.FindSet(v)) {
//            result_edges.push_back(edge);
//            dsu.UnionSets(u, v);
//        }
//    }
//
//    for (const auto& edge : result_edges) {
//        cout << edge.u << " - " << edge.v << ": " << edge.weight << endl;
//    }
//}