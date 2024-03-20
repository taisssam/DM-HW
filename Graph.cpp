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

    Graph::AddWeightedEdges();
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

void Graph::BronKerbosch(set<int> R, set<int> P, set<int> X, vector<set<int>>& cliques) {
    if (P.empty() && X.empty()) {
        cliques.push_back(R);
        return;
    }

    set<int> P_copy = P;
    for (int vertex : P_copy) {
        set<int> newR = R;
        newR.insert(vertex);

        set<int> newP;
        set<int> newX;

        for (int neighbor : adjacency_list[vertex]) {
            if (P.find(neighbor) != P.end()) {
                newP.insert(neighbor);
            }
        }

        for (int neighbor : adjacency_list[vertex]) {
            if (X.find(neighbor) != X.end()) {
                newX.insert(neighbor);
            }
        }

        BronKerbosch(newR, newP, newX, cliques);

        P.erase(vertex);
        X.insert(vertex);
    }
}

void Graph::FindMaximalCliques() {
    set<int> R, P, X;
    vector<set<int>> cliques;

    for (int i = 0; i < adjacency_list.size(); ++i) {
        P.insert(i);
    }

    BronKerbosch(R, P, X, cliques);

    size_t max_size = 0;
    for (const auto& clique : cliques) {
        if (clique.size() > max_size) {
            max_size = clique.size();
        }
    }

    for (const auto& clique : cliques) {
        if (clique.size() == max_size) {
            for (int vertex : clique) {
                cout << vertex + 1 << ' ';
            }
            cout << '\n';
        }
    }
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

vector<WeightedEdge> weighted_edges_list;
void Graph::AddWeightedEdge(int from, int to, double weight) {
    weighted_edges_list.push_back(WeightedEdge(from - 1, to - 1, weight));
}

void Graph::AddWeightedEdges() {
    AddWeightedEdge(1, 19, 501.4);
    AddWeightedEdge(1, 33, 153.5);
    AddWeightedEdge(1, 31, 132);
    AddWeightedEdge(1, 40, 391.2);
    AddWeightedEdge(2, 43, 494.3);
    AddWeightedEdge(2, 16, 708.3);
    AddWeightedEdge(3, 5, 452.9);
    AddWeightedEdge(3, 18, 170);
    AddWeightedEdge(3, 46, 994.5);
    AddWeightedEdge(4, 12, 252.8);
    AddWeightedEdge(4, 17, 524);
    AddWeightedEdge(4, 45, 685.1);
    AddWeightedEdge(4, 25, 527.7);
    AddWeightedEdge(4, 23, 765.5);
    AddWeightedEdge(4, 42, 276.6);
    AddWeightedEdge(4, 20, 214.7);
    AddWeightedEdge(4, 41, 57.3);
    AddWeightedEdge(5, 18, 447.6);
    AddWeightedEdge(5, 38, 1930.3);
    AddWeightedEdge(5, 46, 1445.4);
    AddWeightedEdge(6, 38, 677.5);
    AddWeightedEdge(6, 24, 403.6);
    AddWeightedEdge(6, 26, 172.1);
    AddWeightedEdge(6, 35, 476.9);
    AddWeightedEdge(6, 47, 434.2);
    AddWeightedEdge(7, 32, 173);
    AddWeightedEdge(7, 16, 265.1);
    AddWeightedEdge(7, 17, 651.2);
    AddWeightedEdge(7, 27, 187.7);
    AddWeightedEdge(8, 10, 289.4);
    AddWeightedEdge(8, 31, 172.1);
    AddWeightedEdge(8, 40, 196.9);
    AddWeightedEdge(9, 37, 296.2);
    AddWeightedEdge(9, 40, 329.7);
    AddWeightedEdge(9, 33, 174.1);
    AddWeightedEdge(9, 19, 525.5);
    AddWeightedEdge(9, 46, 855.1);
    AddWeightedEdge(10, 42, 117);
    AddWeightedEdge(10, 20, 300);
    AddWeightedEdge(10, 40, 368.2);
    AddWeightedEdge(10, 31, 457.6);
    AddWeightedEdge(12, 17, 279.7);
    AddWeightedEdge(12, 41, 292.1);
    AddWeightedEdge(12, 35, 518.5);
    AddWeightedEdge(13, 17, 356.8);
    AddWeightedEdge(24, 14, 279.6);
    AddWeightedEdge(38, 14, 870.4);
    AddWeightedEdge(15, 38, 895);
    AddWeightedEdge(15, 44, 397.4);
    AddWeightedEdge(15, 34, 790.6);
    AddWeightedEdge(16, 17, 879);
    AddWeightedEdge(16, 27, 288);
    AddWeightedEdge(16, 45, 436.3);
    AddWeightedEdge(16, 23, 1106.6);
    AddWeightedEdge(16, 43, 1052.5);
    AddWeightedEdge(17, 32, 577.6);
    AddWeightedEdge(17, 27, 602.4);
    AddWeightedEdge(17, 45, 752.3);
    AddWeightedEdge(17, 35, 519.5);
    AddWeightedEdge(18, 38, 1648);
    AddWeightedEdge(18, 46, 1026);
    AddWeightedEdge(19, 46, 819);
    AddWeightedEdge(19, 33, 487.8);
    AddWeightedEdge(20, 41, 160.1);
    AddWeightedEdge(20, 42, 381.7);
    AddWeightedEdge(20, 40, 317.3);
    AddWeightedEdge(20, 37, 644.1);
    AddWeightedEdge(23, 49, 2.7);
    AddWeightedEdge(23, 39, 227);
    AddWeightedEdge(23, 45, 689.6);
    AddWeightedEdge(23, 42, 489.5);
    AddWeightedEdge(24, 38, 844.6);
    AddWeightedEdge(24, 26, 262.4);
    AddWeightedEdge(25, 45, 158.8);
    AddWeightedEdge(26, 38, 792.8);
    AddWeightedEdge(26, 35, 394);
    AddWeightedEdge(29, 47, 400.7);
    AddWeightedEdge(29, 37, 357.6);
    AddWeightedEdge(30, 16, 690.2);
    AddWeightedEdge(31, 40, 281);
    AddWeightedEdge(33, 40, 323.2);
    AddWeightedEdge(34, 44, 418.8);
    AddWeightedEdge(34, 38, 1649.8);
    AddWeightedEdge(35, 38, 1154.3);
    AddWeightedEdge(35, 47, 691.6);
    AddWeightedEdge(35, 41, 530.4);
    AddWeightedEdge(36, 43, 503.9);
    AddWeightedEdge(37, 47, 746.8);
    AddWeightedEdge(37, 40, 450);
    AddWeightedEdge(38, 47, 756.6);
    AddWeightedEdge(41, 47, 1004.6);
    AddWeightedEdge(47, 20, 901.4);
}

void Graph::KruskalMST() {
    sort(weighted_edges_list.begin(), weighted_edges_list.end());

    DisjointSetUnion dsu(vertex_count + 1);

    double total_weight = 0;
    mstEdges.clear();

    for (const auto& edge : weighted_edges_list) {
        if (dsu.FindSet(edge.from) != dsu.FindSet(edge.to)) {
            dsu.UnionSets(edge.from, edge.to);
            total_weight += edge.weight;
            mstEdges.push_back({edge.from, edge.to});
            cout << "Edge (" << edge.from + 1 << ", " << edge.to + 1 << ") with weight: " << edge.weight << endl;
        }
    }

    cout << "Total weight of MST: " << total_weight << endl;

    adjacencyList.clear();
    for (const auto& edge : mstEdges) {
        adjacencyList[edge.first].push_back(edge.second);
        adjacencyList[edge.second].push_back(edge.first);
    }
}

vector<int> Graph::generatePruferCode() {
    map<int, vector<int>> adjacencyListForPrufer;
    for (const auto& edge : mstEdges) {
        adjacencyListForPrufer[edge.first].push_back(edge.second);
        adjacencyListForPrufer[edge.second].push_back(edge.first);
    }

    map<int, int> vertexDegree;
    for (const auto& item : adjacencyListForPrufer) {
        vertexDegree[item.first] = item.second.size();
    }

    vector<int> leaves;
    for (const auto& item : vertexDegree) {
        if (item.second == 1) {
            leaves.push_back(item.first);
        }
    }
    sort(leaves.begin(), leaves.end());

    vector<int> pruferCode;
    while (adjacencyListForPrufer.size() > 2) {
        int leaf = leaves.front();
        leaves.erase(leaves.begin());

        int neighbor = adjacencyListForPrufer[leaf][0];

        pruferCode.push_back(neighbor);

        --vertexDegree[neighbor];

        adjacencyListForPrufer[neighbor].erase(remove(adjacencyListForPrufer[neighbor].begin(), adjacencyListForPrufer[neighbor].end(), leaf), adjacencyListForPrufer[neighbor].end());

        if (vertexDegree[neighbor] == 1) {
            leaves.push_back(neighbor);
            sort(leaves.begin(), leaves.end());
        }

        adjacencyListForPrufer.erase(leaf);
        vertexDegree.erase(leaf);
    }

    return pruferCode;
}

void dfsForBinaryCode(int vertex, map<int, bool>& visited, vector<int>& ans, const map<int, vector<int>>& graph) {
    visited[vertex] = true;
    for (int neighbor : graph.at(vertex)) {
        if (!visited[neighbor]) {
            ans.push_back(1);
            dfsForBinaryCode(neighbor, visited, ans, graph);
        }
    }
    ans.push_back(0);
}

vector<int> Graph::generateBinaryCode() {
    map<int, bool> visited;
    for (const auto& item : adjacencyList) {
        visited[item.first] = false;
    }

    vector<int> binaryCode;
    int startVertex = adjacencyList.begin()->first;
    dfsForBinaryCode(startVertex, visited, binaryCode, adjacencyList);

    return binaryCode;
}

struct Edge {
    int u, v;
    double weight;

    Edge(int u, int v, double weight) : u(u), v(v), weight(weight) {}
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};
