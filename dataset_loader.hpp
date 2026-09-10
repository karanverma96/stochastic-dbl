#pragma once
#include "graph.hpp"
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <stdexcept>

// Loads a SNAP-format directed edge list: '#' starts a comment, every other
// line is "FromNodeId<whitespace>ToNodeId". SNAP ids are sparse, so they get
// remapped to a dense [0, n) range.
//
// Repeated edges are dropped, because degree here is just adjacency-list
// length -- a duplicate would be counted twice, and those counts pick the
// landmarks and set the WeightedCascade probabilities.
struct LoadedEdgeList {
    int numVertices;
    std::vector<std::pair<int,int>> edges;
    std::unordered_map<long long, int> idMap; // original id -> dense id
    long long duplicatesDropped = 0;
    long long malformedLines = 0;
};

inline LoadedEdgeList loadSnapEdgeList(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("could not open dataset file: " + path);

    LoadedEdgeList result;
    std::unordered_set<long long> seen;

    auto remap = [&](long long id) -> int {
        auto it = result.idMap.find(id);
        if (it != result.idMap.end()) return it->second;
        int dense = static_cast<int>(result.idMap.size());
        result.idMap[id] = dense;
        return dense;
    };

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        long long a, b;
        if (!(iss >> a >> b)) { ++result.malformedLines; continue; }

        // Two declarators in one declaration are sequenced left to right, so
        // the ids come out the same on every compiler. As f(remap(a),
        // remap(b)) this would have been unspecified before C++17.
        int u = remap(a), v = remap(b);

        long long key = (static_cast<long long>(u) << 32)
                      | static_cast<unsigned int>(v);
        if (!seen.insert(key).second) { ++result.duplicatesDropped; continue; }

        result.edges.emplace_back(u, v);
    }
    result.numVertices = static_cast<int>(result.idMap.size());
    return result;
}
