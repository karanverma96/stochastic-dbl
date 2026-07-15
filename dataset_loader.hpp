#pragma once
#include "graph.hpp"
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>

// Loads a SNAP-format directed edge list:
//   lines starting with '#' are comments
//   remaining lines are "FromNodeId<whitespace>ToNodeId"
// Original node ids are remapped to a dense [0, n) range since SNAP ids
// are often sparse/non-contiguous.
struct LoadedEdgeList {
    int numVertices;
    std::vector<std::pair<int,int>> edges;
    std::unordered_map<long long, int> idMap; // original id -> dense id
};

inline LoadedEdgeList loadSnapEdgeList(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("could not open dataset file: " + path);

    LoadedEdgeList result;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        long long a, b;
        if (!(iss >> a >> b)) continue;

        auto remap = [&](long long id) -> int {
            auto it = result.idMap.find(id);
            if (it != result.idMap.end()) return it->second;
            int dense = static_cast<int>(result.idMap.size());
            result.idMap[id] = dense;
            return dense;
        };
        int u = remap(a), v = remap(b);
        result.edges.emplace_back(u, v);
    }
    result.numVertices = static_cast<int>(result.idMap.size());
    return result;
}
