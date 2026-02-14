#pragma once
#include "Topology_enhanced.h"
#include "TopoLineCompact_enhanced.hpp"  // Use separate file
#include <string>
#include <vector>
#include <functional>

class TopologyDB_enhanced {
public:
    struct Record {
        std::string name;
        Topology_enhanced topo;
    };

    explicit TopologyDB_enhanced(std::string path);

    // Basic operations
    bool append(const Topology_enhanced& T) const;
    std::vector<Record> loadAll() const;
    bool loadByName(const std::string& name, Topology_enhanced& out) const;

    // Deduplication
    int dedupeByContentHash(bool keep_last = false) const;
    int dedupeByName(bool keep_last = false) const;

    // ✨ NEW: External-specific queries
    std::vector<Record> loadWithExternals() const;  // Load only topologies with external curves
    std::vector<Record> loadWithoutExternals() const;  // Load only topologies without external curves
    std::vector<Record> loadWithExternalsOnPort(int parent_type, int port) const;
    
    // ✨ NEW: Statistics
    struct ExternalStats {
        int total_topologies;
        int topologies_with_externals;
        int total_externals;
        int externals_on_blocks;
        int externals_on_sidelinks;
        int externals_on_instantons;
        int max_externals_per_topology;
        int max_port_index_used;
    };
    ExternalStats getExternalStatistics() const;

    // Serialization
    static std::string serializeCanonical(const Topology_enhanced& T);
    static bool deserializeCanonical(std::istream& in, int nLines, Topology_enhanced& out);
    
    // Migration from basic Topology
    static Topology_enhanced upgradeFromBasic(const struct Topology& basic);
    bool importFromBasicDB(const std::string& basic_db_path);

private:
    std::string path_;
    
    // Internal utilities
    static int countLines(const std::string& s);
    static std::string cheapHashHex(const std::string& s);
    static bool writeFileAtomic(const std::string& path, const std::string& content);
};

// TopoLineCompact_enhanced is now in its own file
