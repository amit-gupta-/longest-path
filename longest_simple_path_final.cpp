#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <iostream>
#include <iterator>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class Graph {
  public:
    using Vertex = size_t;
    using VertexSet = std::unordered_set<Vertex>;
    using AdjacencyList = std::unordered_map<Vertex, VertexSet>;

    void AddVertex(Vertex v) {
        adjacency_list_[v];
    }

    void AddEdge(Vertex u, Vertex v) {
        adjacency_list_[u].insert(v);
        adjacency_list_[v].insert(u);
    }

    const VertexSet& AdjacentVertices(Vertex v) const {
        const auto it = adjacency_list_.find(v);
        if (it != adjacency_list_.end()) {
            return it->second;
        } else {
            return empty_set_;
        }
    }

    VertexSet AllVertices() const {
        VertexSet vs;
        vs.reserve(adjacency_list_.size());
        for (const auto& pair : adjacency_list_) {
            const auto& vertex = pair.first;
            vs.insert(vertex);
        }
        return vs;
    }

    const AdjacencyList& AsAdjacencyList() const {
        return adjacency_list_;
    }

  private:
    AdjacencyList adjacency_list_;
    static const VertexSet empty_set_;
};

const Graph::VertexSet Graph::empty_set_;

class LongestPath {
  public:
    explicit LongestPath(const Graph& graph)
        : graph_(graph)  { visited_.resize(graph.AllVertices().size(), false); }

    void GetStartVertex() {
        // Получаем случайную вершину из графа
        const auto vertices = graph_.AllVertices();
        auto iter = vertices.begin();
        std::advance(iter, rand() % vertices.size());
        // Помечаем ее посещенной и засовываем в путь
        visited_[*iter] = true;
        deque_.push_front(*iter);
    }

    Graph::VertexSet CandidatesToRemove() const {
        Graph::VertexSet candidates;
        if (deque_.size() < 1)
            return candidates;
        candidates.insert(deque_.front());
        if (deque_.size() > 1)
            candidates.insert(deque_.back());
        return candidates;
    }

    const Graph::VertexSet CandidatesToAdd() const {
        Graph::VertexSet candidates;
        for (auto vertex : graph_.AdjacentVertices(deque_.front())) {
            if (visited_[vertex] == false)
                candidates.insert(vertex);
        }
        for (auto vertex : graph_.AdjacentVertices(deque_.back())) {
            if (visited_[vertex] == false)
                candidates.insert(vertex);
        }
        return candidates;
    }

    void Add(Graph::Vertex v) {
        auto front_adj = graph_.AdjacentVertices(deque_.front());
        if (std::find(front_adj.begin(), front_adj.end(), v) != front_adj.end()) {
            deque_.push_front(v);
            visited_[v] = true;
        } else {
            deque_.push_back(v);
            visited_[v] = true;
        }
    }

    void Remove(Graph::Vertex v) {
        if (deque_.front() == v) {
            deque_.pop_front();
            visited_[v] = false;
        } else if (deque_.back() == v) {
            deque_.pop_back();
            visited_[v] = false;
        }
    }

    size_t Cost() const {
        return deque_.size();
    }

    const Graph& GetGraph() const {
        return graph_;
    }

    const std::deque<Graph::Vertex> GetPath() const {
        return deque_;
    }

  private:
    const Graph& graph_;
    std::vector<bool> visited_;
    std::deque<Graph::Vertex> deque_;
};

void GraphEdges(std::ostream& out, const Graph::AdjacencyList& adjacency_list, const LongestPath& longest_path) {
    // Конструируем из очереди список смежности
    std::deque<size_t> path = longest_path.GetPath();
    // Выводим первую вершину, потом оставшийся путь, затем последнюю
    auto iter = path.begin();
    auto last_elem_iter = path.end();
    out << "\t" << *iter << " -- ";
    --last_elem_iter;
    ++iter;
    for (; iter != last_elem_iter && iter != path.end(); ++iter) {
        out << *iter << " [color=\"red\" penwidth=3]\n" << "\t" << *iter << " -- ";
    }
    out << *last_elem_iter << " [color=\"red\" penwidth=3]\n";

    for (const auto& pair : adjacency_list) {
        const auto& vertex = pair.first;
        const auto& neighbours = pair.second;
        for (const auto adj_vertex : neighbours) {
            out << "\t" << vertex << " -- " << adj_vertex << "\n";
        }
    }
}

void GraphViz(std::ostream& out, const LongestPath& longest_path) {
    out << "THE SIZE OF DEQUE IS " << longest_path.Cost() << "\n";
    out << "strict graph {\n";
    for (const auto& pair : longest_path.GetGraph().AsAdjacencyList()) {
        const auto& vertex = pair.first;
        out << "\t" << vertex << "\n";
    }
    GraphEdges(out, longest_path.GetGraph().AsAdjacencyList(), longest_path);
    out << "}\n";
}

class LongestPathSolver {
  public:
    virtual LongestPath Solve(const Graph& graph) const = 0;
    virtual ~LongestPathSolver() = default;
};

class GradientDescent final: public LongestPathSolver {
    LongestPath Solve(const Graph& graph) const {
        LongestPath longest_path(graph);
        longest_path.GetStartVertex();
        while (1) {
            const auto candidates = longest_path.CandidatesToAdd();
            if (candidates.empty()) {
                return longest_path;
            }
            auto random_candidate = candidates.begin();
            std::advance(random_candidate, rand() % candidates.size());
            longest_path.Add(*random_candidate);
        }
    }
};

class Metropolis final: public LongestPathSolver {
  public:
    Metropolis(double k, double t, bool annealing=true, size_t iterations=150)
        : k_(k), t_(t), annealing_(annealing), iterations_(iterations) {
    }

    LongestPath Solve(const Graph& graph) const {
        double t = t_;
        LongestPath longest_path(graph);
        longest_path.GetStartVertex();
        for (size_t i = 0; i < iterations_; ++i) {
            const auto remove_candidates = longest_path.CandidatesToRemove();
            const auto add_candidates = longest_path.CandidatesToAdd();
            if (add_candidates.size() > 0) {
                auto random_candidate = add_candidates.begin();
                std::advance(random_candidate, rand() % add_candidates.size());
                longest_path.Add(*random_candidate);
            } else if (static_cast<double>(rand()) / RAND_MAX >= exp(-1. / k_ / t) && remove_candidates.size() > 1) {
                auto random_candidate = remove_candidates.begin();
                std::advance(random_candidate, rand() % remove_candidates.size());
                longest_path.Remove(*random_candidate);
            }
            if (annealing_) {
                t /= 2;
            }
        }
        return longest_path;
    }

  private:
    double k_;
    double t_;
    bool annealing_;
    size_t iterations_;
};

Graph RandomGraph(size_t size, double edge_probability) {
    Graph graph;
    for (Graph::Vertex v = 1; v <= size; ++v) {
        graph.AddVertex(v);
    }
    for (Graph::Vertex v = 1; v <= size; ++v) {
        for (Graph::Vertex u = v + 1; u <= size; ++u) {
            if (double(rand()) / RAND_MAX <= edge_probability) {
                graph.AddEdge(v, u);
            }
        }
    }
    return graph;
}

Graph StarGraph(size_t size) {
    Graph graph;
    for (Graph::Vertex v = 2; v <= size; ++v) {
        graph.AddEdge(1, v);
    }
    return graph;
}

int InitRandSeed(int argc, const char* argv[]) {
    int rand_seed;
    if (argc >= 2) {
        rand_seed = atoi(argv[1]);
    } else {
        rand_seed = time(nullptr);
    }
    srand(rand_seed);
    return rand_seed;
}

void TrySolver(const LongestPathSolver& solver, const Graph& graph, std::vector<size_t>& bests) {
    auto best_cost = 0;
    for (int attempt = 1; attempt < 100; ++attempt) {
        const auto longest_path = solver.Solve(graph);
        auto cost = longest_path.Cost();
        if (cost > best_cost) {
            best_cost = cost;
            // GraphViz(std::cout, longest_path);
        }
    }
    bests.push_back(best_cost);
}

int main(int argc, const char* argv[]) {
    int number_of_vertices = 30;
    std::cout << "Choose number of vertices: ";
    std::cin >> number_of_vertices;
    double edge_probability = 0.5;
    std::cout << "Choose probability of edge between two vertices: ";
    std::cin >> edge_probability;

    std::vector<size_t> gd_bests;
    std::vector<size_t> m_bests;
    std::vector<size_t> mwa_bests;

    for (int i = 0; i != 30; ++i) {
        std::cout << "\n";
        // InitRandSeed(argc, argv);
        // std::cout << "Using rand seed: " << InitRandSeed(argc, argv) << "\n";
        const auto graph = RandomGraph(number_of_vertices, edge_probability);
        GradientDescent gradient_descent;
        TrySolver(gradient_descent, graph, gd_bests);
        std::cout << "GD " << gd_bests.back() << " ";
        Metropolis metropolis(1, 10000, false, 100);
        TrySolver(metropolis, graph, m_bests);
        std::cout << "M " << m_bests.back() << " ";
        Metropolis metropolis2(1, 10000, true, 100);
        TrySolver(metropolis2, graph, mwa_bests);
        std::cout << "MwA " << mwa_bests.back();
    }
    size_t count = 0;
    for (size_t j = 0; j != gd_bests.size(); ++j) {
        int longest = std::max(gd_bests[j], m_bests[j]);
        if (longest <= mwa_bests[j]) {
            ++count;
        }
    }
    std::cout << "\nAnnealing was at least not worse than the other in " << static_cast<double>(count) * 100 / gd_bests.size() << "% of cases";
    return 0;
}
