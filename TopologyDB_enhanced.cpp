#include "TopologyDB_enhanced.hpp"
#include "TopoLineCompact_enhanced.hpp"  // Include separate implementation
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <filesystem>
#include <algorithm>
#include <cctype>

// ===== Helper functions for LKind conversion =====
static inline int kindToInt(LKind k) {
    switch (k) {
        case LKind::g: return 0;
        case LKind::L: return 1;
        case LKind::S: return 2;
        case LKind::I: return 3;
        case LKind::E: return 4;  // Enhanced: External
    }
    return -1;
}

static inline LKind intToKind(int k) {
    switch (k) {
        case 0: return LKind::g;
        case 1: return LKind::L;
        case 2: return LKind::S;
        case 3: return LKind::I;
        case 4: return LKind::E;  // Enhanced: External
    }
    return LKind::g;
}

// ===== TopologyDB_enhanced implementation =====
TopologyDB_enhanced::TopologyDB_enhanced(std::string path) : path_(std::move(path)) {}

int TopologyDB_enhanced::countLines(const std::string& s) {
    int n = 0;
    for (char c : s) if (c == '\n') ++n;
    return n;
}

std::string TopologyDB_enhanced::cheapHashHex(const std::string& s) {
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s) {
        h ^= c;
        h *= 1099511628211ULL;
    }
    std::ostringstream os;
    os << std::hex << h;
    return os.str();
}

bool TopologyDB_enhanced::writeFileAtomic(const std::string& path, const std::string& content) {
    const auto p = std::filesystem::path(path);
    std::filesystem::create_directories(p.parent_path());
    const auto tmp = (p.parent_path() / (p.filename().string() + ".tmp")).string();
    
    {
        std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
        if (!out) return false;
        out.write(content.data(), (std::streamsize)content.size());
        if (!out.good()) return false;
    }
    
    std::error_code ec;
    std::filesystem::rename(tmp, path, ec);
    if (ec) {
        std::filesystem::remove(tmp);
        return false;
    }
    return true;
}

// ===== Serialization =====
std::string TopologyDB_enhanced::serializeCanonical(const Topology_enhanced& T) {
    std::ostringstream os;
    os << "name:" << T.name << "\n";
    os << "blocks:" << T.block.size() << "\n";
    for (const auto& b : T.block) {
        os << "  " << kindToInt(b.kind) << "," << b.param << "\n";
    }
    
    os << "side_links:" << T.side_links.size() << "\n";
    for (const auto& s : T.side_links) {
        os << "  " << s.param << "\n";
    }
    
    os << "instantons:" << T.instantons.size() << "\n";
    for (const auto& i : T.instantons) {
        os << "  " << i.param << "\n";
    }
    
    // ✨ NEW: Externals
    os << "externals:" << T.externals.size() << "\n";
    for (const auto& e : T.externals) {
        os << "  " << e.param << "\n";
    }
    
    os << "l_conn:" << T.l_connection.size() << "\n";
    for (const auto& e : T.l_connection) {
        os << "  " << e.u << "," << e.v << "\n";
    }
    
    os << "s_conn:" << T.s_connection.size() << "\n";
    for (const auto& e : T.s_connection) {
        os << "  " << e.u << "," << e.v << "\n";
    }
    
    os << "i_conn:" << T.i_connection.size() << "\n";
    for (const auto& e : T.i_connection) {
        os << "  " << e.u << "," << e.v << "\n";
    }
    
    // ✨ NEW: External connections
    os << "e_conn:" << T.e_connection.size() << "\n";
    for (const auto& e : T.e_connection) {
        os << "  " << e.parent_id << "," << e.parent_type << "," 
           << e.port_idx << "," << e.external_id << "\n";
    }
    
    return os.str();
}

static void trim_end(std::string& s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n')) s.pop_back();
}

static int parseCount(const std::string& s, const char* key) {
    if (s.rfind(key, 0) != 0) return 0;
    return std::stoi(s.substr(std::string(key).size()));
}

bool TopologyDB_enhanced::deserializeCanonical(std::istream& in, int /*nLines*/, Topology_enhanced& T) {
    T.Initialize();
    std::string line;

    if (!std::getline(in, line) || line.rfind("name:", 0) != 0) return false;
    T.name = line.substr(5);

    // Blocks
    if (!std::getline(in, line)) return false;
    trim_end(line);
    {
        int n = parseCount(line, "blocks:");
        for (int i = 0; i < n; i++) {
            if (!std::getline(in, line)) return false;
            trim_end(line);
            std::istringstream ls(line);
            int k, p;
            char comma;
            if (!(ls >> k >> comma >> p)) return false;
            T.addBlock(intToKind(k), p);
        }
    }

    // Side links
    if (!std::getline(in, line)) return false;
    trim_end(line);
    {
        int n = parseCount(line, "side_links:");
        for (int i = 0; i < n; i++) {
            if (!std::getline(in, line)) return false;
            trim_end(line);
            int param = std::stoi(line);
            T.addSideLink(param);
        }
    }

    // Instantons
    if (!std::getline(in, line)) return false;
    trim_end(line);
    {
        int n = parseCount(line, "instantons:");
        for (int i = 0; i < n; i++) {
            if (!std::getline(in, line)) return false;
            trim_end(line);
            int param = std::stoi(line);
            T.addInstanton(param);
        }
    }

    // ✨ NEW: Externals
    if (!std::getline(in, line)) return false;
    trim_end(line);
    {
        int n = parseCount(line, "externals:");
        for (int i = 0; i < n; i++) {
            if (!std::getline(in, line)) return false;
            trim_end(line);
            int param = std::stoi(line);
            T.addExternal(param);
        }
    }

    // L connections
    if (!std::getline(in, line)) return false;
    trim_end(line);
    {
        int n = parseCount(line, "l_conn:");
        std::vector<InteriorStructure> L;
        for (int i = 0; i < n; i++) {
            if (!std::getline(in, line)) return false;
            trim_end(line);
            std::istringstream ls(line);
            int u, v;
            char comma;
            if (ls >> u >> comma >> v) L.push_back({u, v});
        }
        T.LinkingBlocks(std::move(L));
    }

    // S connections
    if (!std::getline(in, line)) return false;
    trim_end(line);
    {
        int n = parseCount(line, "s_conn:");
        std::vector<SideLinkStructure> S;
        for (int i = 0; i < n; i++) {
            if (!std::getline(in, line)) return false;
            trim_end(line);
            std::istringstream ls(line);
            int u, v;
            char comma;
            if (ls >> u >> comma >> v) S.push_back({u, v});
        }
        T.LinkingSideLinks(std::move(S));
    }

    // I connections
    if (!std::getline(in, line)) return false;
    trim_end(line);
    {
        int n = parseCount(line, "i_conn:");
        std::vector<InstantonStructure> I;
        for (int i = 0; i < n; i++) {
            if (!std::getline(in, line)) return false;
            trim_end(line);
            std::istringstream ls(line);
            int u, v;
            char comma;
            if (ls >> u >> comma >> v) I.push_back({u, v});
        }
        T.LinkingInstantonStructure(std::move(I));
    }

    // ✨ NEW: E connections
    if (!std::getline(in, line)) return false;
    trim_end(line);
    {
        int n = parseCount(line, "e_conn:");
        std::vector<ExternalStructure> E;
        for (int i = 0; i < n; i++) {
            if (!std::getline(in, line)) return false;
            trim_end(line);
            std::istringstream ls(line);
            int parent_id, parent_type, port_idx, external_id;
            char c1, c2, c3;
            if (ls >> parent_id >> c1 >> parent_type >> c2 >> port_idx >> c3 >> external_id) {
                E.push_back({parent_id, parent_type, port_idx, external_id});
            }
        }
        T.LinkingExternals(std::move(E));
    }

    return true;
}

// ===== Basic I/O =====
bool TopologyDB_enhanced::append(const Topology_enhanced& T) const {
    const std::string payload = serializeCanonical(T);
    const int lines = countLines(payload);
    std::ofstream out(path_, std::ios::app);
    if (!out) return false;
    out << T.name << "\t" << lines << "\n" << payload;
    return true;
}

std::vector<TopologyDB_enhanced::Record> TopologyDB_enhanced::loadAll() const {
    std::vector<Record> out;
    std::ifstream in(path_);
    if (!in) return out;

    std::string header;
    while (std::getline(in, header)) {
        if (header.empty()) continue;
        std::istringstream hs(header);
        std::string nm, cnt;
        if (!std::getline(hs, nm, '\t')) continue;
        if (!std::getline(hs, cnt)) continue;
        int n = std::stoi(cnt);

        std::ostringstream payload;
        for (int i = 0; i < n; i++) {
            std::string l;
            if (!std::getline(in, l)) break;
            payload << l << "\n";
        }

        std::istringstream pin(payload.str());
        Topology_enhanced T;
        if (!deserializeCanonical(pin, n, T)) continue;
        if (T.name.empty()) T.name = nm;

        out.push_back(Record{nm, std::move(T)});
    }
    return out;
}

bool TopologyDB_enhanced::loadByName(const std::string& name, Topology_enhanced& out) const {
    auto all = loadAll();
    for (auto& r : all) {
        const std::string nm = r.topo.name.empty() ? r.name : r.topo.name;
        if (nm == name) {
            out = r.topo;
            return true;
        }
    }
    return false;
}

// ✨ NEW: External-specific queries
std::vector<TopologyDB_enhanced::Record> TopologyDB_enhanced::loadWithExternals() const {
    std::vector<Record> result;
    auto all = loadAll();
    for (auto& r : all) {
        if (r.topo.hasExternalCurves()) {
            result.push_back(std::move(r));
        }
    }
    return result;
}

std::vector<TopologyDB_enhanced::Record> TopologyDB_enhanced::loadWithoutExternals() const {
    std::vector<Record> result;
    auto all = loadAll();
    for (auto& r : all) {
        if (!r.topo.hasExternalCurves()) {
            result.push_back(std::move(r));
        }
    }
    return result;
}

// ✨ NEW: Statistics
TopologyDB_enhanced::ExternalStats TopologyDB_enhanced::getExternalStatistics() const {
    ExternalStats stats = {0, 0, 0, 0, 0, 0, 0, 0};
    auto all = loadAll();
    
    stats.total_topologies = static_cast<int>(all.size());
    
    for (const auto& r : all) {
        const auto& T = r.topo;
        if (T.hasExternalCurves()) {
            stats.topologies_with_externals++;
            stats.total_externals += static_cast<int>(T.externals.size());
            stats.max_externals_per_topology = std::max(
                stats.max_externals_per_topology, 
                static_cast<int>(T.externals.size())
            );
            
            for (const auto& ec : T.e_connection) {
                stats.max_port_index_used = std::max(stats.max_port_index_used, ec.port_idx);
                switch (ec.parent_type) {
                    case 0: stats.externals_on_blocks++; break;
                    case 1: stats.externals_on_sidelinks++; break;
                    case 2: stats.externals_on_instantons++; break;
                }
            }
        }
    }
    
    return stats;
}

// ===== Deduplication =====
int TopologyDB_enhanced::dedupeByContentHash(bool keep_last) const {
    auto all = loadAll();
    if (all.empty()) return 0;

    std::unordered_map<std::string, int> seen;
    std::vector<Record> kept;
    kept.reserve(all.size());

    if (!keep_last) {
        for (auto& r : all) {
            const std::string payload = serializeCanonical(r.topo);
            const std::string h = cheapHashHex(payload);
            if (seen.emplace(h, 1).second) kept.push_back(std::move(r));
        }
    } else {
        std::vector<Record> rev(all.rbegin(), all.rend());
        for (auto& r : rev) {
            const std::string payload = serializeCanonical(r.topo);
            const std::string h = cheapHashHex(payload);
            if (seen.emplace(h, 1).second) kept.push_back(std::move(r));
        }
        std::reverse(kept.begin(), kept.end());
    }

    std::ostringstream buf;
    for (auto& r : kept) {
        const std::string payload = serializeCanonical(r.topo);
        buf << r.topo.name << "\t" << countLines(payload) << "\n" << payload;
    }
    if (!writeFileAtomic(path_, buf.str())) return 0;
    return (int)all.size() - (int)kept.size();
}

int TopologyDB_enhanced::dedupeByName(bool keep_last) const {
    auto all = loadAll();
    if (all.empty()) return 0;

    std::unordered_map<std::string, int> seen;
    std::vector<Record> kept;
    kept.reserve(all.size());

    if (!keep_last) {
        for (auto& r : all) {
            const std::string nm = r.topo.name.empty() ? r.name : r.topo.name;
            if (seen.emplace(nm, 1).second) kept.push_back(std::move(r));
        }
    } else {
        std::vector<Record> rev(all.rbegin(), all.rend());
        for (auto& r : rev) {
            const std::string nm = r.topo.name.empty() ? r.name : r.topo.name;
            if (seen.emplace(nm, 1).second) kept.push_back(std::move(r));
        }
        std::reverse(kept.begin(), kept.end());
    }

    std::ostringstream buf;
    for (auto& r : kept) {
        const std::string payload = serializeCanonical(r.topo);
        buf << r.topo.name << "\t" << countLines(payload) << "\n" << payload;
    }
    if (!writeFileAtomic(path_, buf.str())) return 0;
    return (int)all.size() - (int)kept.size();
}
