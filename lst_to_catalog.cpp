// lst_to_catalog.cpp
// Converts existing LST_dedup topocompactline files to unified catalog format
// Usage: ./lst_to_catalog <input_file_or_dir> <output_catalog>
//
// Unified format: id | type | T | topo
//   type = Nk (k = node count from topocompactline)
//   topo = original topocompactline (verbatim)
//   T = sig_neg of IF

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <filesystem>
#include <cmath>
#include <functional>
#include <Eigen/Dense>

#include "Tensor.h"
#include "Theory_enhanced.h"
#include "Topology_enhanced.h"
#include "TopoLineCompact_enhanced.hpp"
#include "anomaly_tables.h"

namespace fs = std::filesystem;

// ============================================================================
// Build IF from Topology_enhanced (simplified from sugra_generator_v2)
// ============================================================================

Spec make_spec(LKind kind, int param) {
    switch (kind) {
        case LKind::g: return n(param);
        case LKind::L: return i(param);
        case LKind::S: return s(param);
        case LKind::I: return s(param);
        case LKind::E: return e(param);
        default: return n(param);
    }
}

Eigen::MatrixXi build_IF_from_topology(const Topology_enhanced& T) {
    TheoryGraph G;
    std::vector<NodeRef> blockNodes;
    std::vector<NodeRef> sideNodes;
    std::vector<NodeRef> instNodes;
    std::vector<NodeRef> extNodes;
    
    // 1. Blocks (chain backbone)
    for (size_t i = 0; i < T.block.size(); ++i) {
        const auto& b = T.block[i];
        blockNodes.push_back(G.add(make_spec(b.kind, b.param)));
    }
    
    // 2. SideLinks
    for (size_t i = 0; i < T.side_links.size(); ++i) {
        sideNodes.push_back(G.add(s(T.side_links[i].param)));
    }
    
    // 3. Instantons
    for (size_t i = 0; i < T.instantons.size(); ++i) {
        instNodes.push_back(G.add(s(T.instantons[i].param)));
    }
    
    // 4. Externals
    for (size_t i = 0; i < T.externals.size(); ++i) {
        extNodes.push_back(G.add(e(T.externals[i].param)));
    }
    
    // 5. Connections
    // Block-block (chain links: implicit from deserialize via addBlockRight)
    for (const auto& conn : T.l_connection) {
        if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
            conn.v >= 0 && conn.v < (int)blockNodes.size()) {
            G.connect(blockNodes[conn.u], blockNodes[conn.v]);
        }
    }
    
    // SideLink-Block
    for (const auto& conn : T.s_connection) {
        if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
            conn.v >= 0 && conn.v < (int)sideNodes.size()) {
            G.connect(sideNodes[conn.v], blockNodes[conn.u]);
        }
    }
    
    // Instanton-Block
    for (const auto& conn : T.i_connection) {
        if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
            conn.v >= 0 && conn.v < (int)instNodes.size()) {
            G.connect(instNodes[conn.v], blockNodes[conn.u]);
        }
    }
    
    // External connections
    for (const auto& conn : T.e_connection) {
        if (conn.external_id < 0 || conn.external_id >= (int)extNodes.size()) 
            continue;
        
        NodeRef* parentRef = nullptr;
        switch (conn.parent_type) {
            case 0:
                if (conn.parent_id >= 0 && conn.parent_id < (int)blockNodes.size())
                    parentRef = &blockNodes[conn.parent_id];
                break;
            case 1:
                if (conn.parent_id >= 0 && conn.parent_id < (int)sideNodes.size())
                    parentRef = &sideNodes[conn.parent_id];
                break;
            case 2:
                if (conn.parent_id >= 0 && conn.parent_id < (int)instNodes.size())
                    parentRef = &instNodes[conn.parent_id];
                break;
        }
        
        if (parentRef) {
            AttachmentPoint ap(conn.port_idx);
            G.connect(extNodes[conn.external_id], AttachmentPoint(-1), *parentRef, ap);
        }
    }
    
    return G.ComposeIF_Gluing();
}

// ============================================================================
// Compute signature (T = sig_neg)
// ============================================================================

struct SigResult {
    int sig_pos, sig_neg, sig_zero;
};

SigResult compute_signature(const Eigen::MatrixXi& IF) {
    SigResult r = {0, 0, 0};
    if (IF.rows() == 0) return r;
    
    Eigen::MatrixXd IFd = IF.cast<double>();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(IFd);
    auto eigs = solver.eigenvalues();
    
    const double eps = 1e-6;
    for (int i = 0; i < eigs.size(); ++i) {
        if (eigs(i) > eps) r.sig_pos++;
        else if (eigs(i) < -eps) r.sig_neg++;
        else r.sig_zero++;
    }
    return r;
}

// ============================================================================
// Compute anomaly from Topology_enhanced
// ============================================================================

struct AnomalyResult {
    int T, H, V;
    int residual; // H - V + 29T - 273
};

AnomalyResult compute_anomaly(const Topology_enhanced& T, int sig_neg) {
    AnomalyResult result;
    result.H = 0;
    result.V = 0;
    result.T = sig_neg;
    
    for (const auto& b : T.block) {
        AnomalyContrib c;
        switch (b.kind) {
            case LKind::g: c = get_node_contrib(b.param); break;
            case LKind::L: c = get_interior_contrib(b.param); break;
            case LKind::S: c = get_sidelink_contrib(b.param); break;
            case LKind::I: c = get_instanton_contrib(b.param); break;
            case LKind::E: c = get_external_contrib(b.param); break;
            default: c = {0, 0}; break;
        }
        result.H += c.H;
        result.V += c.V;
    }
    
    for (const auto& sl : T.side_links) {
        auto c = get_sidelink_contrib(sl.param);
        result.H += c.H;
        result.V += c.V;
    }
    
    for (const auto& inst : T.instantons) {
        auto c = get_instanton_contrib(inst.param);
        result.H += c.H;
        result.V += c.V;
    }
    
    for (const auto& ext : T.externals) {
        auto c = get_external_contrib(ext.param);
        result.H += c.H;
        result.V += c.V;
    }
    
    result.residual = result.H - result.V + 29 * result.T - 273;
    return result;
}

// ============================================================================
// Count g-nodes in topology (for Nk type)
// ============================================================================

int count_nodes(const Topology_enhanced& T) {
    int count = 0;
    for (const auto& b : T.block) {
        if (b.kind == LKind::g) count++;
    }
    return count;
}

// ============================================================================
// Process one topocompactline
// ============================================================================

struct CatalogEntry {
    int node_count;
    int T;
    std::string topo_line;  // original topocompactline (verbatim)
    bool valid;
};

CatalogEntry process_line(const std::string& line) {
    CatalogEntry entry;
    entry.valid = false;
    entry.topo_line = line;
    
    // Deserialize
    Topology_enhanced topo;
    if (!TopoLineCompact_enhanced::deserialize(line, topo)) {
        return entry;
    }
    
    // Count nodes
    entry.node_count = count_nodes(topo);
    
    // Build IF
    Eigen::MatrixXi IF = build_IF_from_topology(topo);
    if (IF.rows() == 0) {
        return entry;
    }
    
    // Signature → T
    SigResult sig = compute_signature(IF);
    entry.T = sig.sig_neg;
    entry.valid = true;
    
    return entry;
}

// ============================================================================
// Find all dedup files in a directory (recursively)
// ============================================================================

std::vector<std::string> find_dedup_files(const std::string& dir) {
    std::vector<std::string> files;
    
    for (const auto& entry : fs::recursive_directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        std::string name = entry.path().filename().string();
        if (name.find("_dedup") != std::string::npos && name.find("._") == std::string::npos) {
            files.push_back(entry.path().string());
        }
    }
    
    std::sort(files.begin(), files.end());
    return files;
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file_or_dir> <output_catalog> [--verbose]\n";
        std::cerr << "\n  input: single topocompactline file, or directory with *_dedup files\n";
        std::cerr << "  output: unified catalog file\n";
        return 1;
    }
    
    std::string input = argv[1];
    std::string output = argv[2];
    bool verbose = (argc > 3 && std::string(argv[3]) == "--verbose");
    
    // Collect input files
    std::vector<std::string> files;
    if (fs::is_directory(input)) {
        files = find_dedup_files(input);
        std::cout << "Found " << files.size() << " dedup files in " << input << "\n";
    } else {
        files.push_back(input);
    }
    
    // Process all files
    std::ofstream out(output);
    if (!out) {
        std::cerr << "Error: Cannot open " << output << " for writing\n";
        return 1;
    }
    
    out << "# Unified LST Catalog\n";
    out << "# id | type | T | topo\n";
    out << "# NOTE: For Nk entries, topo field contains '|' (topocompactline format)\n";
    out << "#   First 3 fields are pipe-delimited, topo is everything after 3rd '|'\n";
    out << "#   Parser: split on '|', fields[0..2] = id,type,T, topo = join(fields[3:], '|')\n";
    
    int total_lines = 0, success = 0, fail = 0, skipped_dup = 0;
    std::map<int, int> node_counts;  // Nk → count
    std::set<std::string> seen_lines;  // dedup across files
    
    for (const auto& file : files) {
        std::ifstream fin(file);
        if (!fin) {
            std::cerr << "Warning: Cannot open " << file << "\n";
            continue;
        }
        
        std::string line;
        int file_ok = 0, file_fail = 0;
        
        while (std::getline(fin, line)) {
            // Skip empty lines and comments
            if (line.empty() || line[0] == '#') continue;
            
            // Dedup across files
            if (seen_lines.count(line)) {
                skipped_dup++;
                continue;
            }
            seen_lines.insert(line);
            
            total_lines++;
            CatalogEntry entry = process_line(line);
            
            if (entry.valid) {
                out << success << " | N" << entry.node_count << " | "
                    << entry.T << " | "
                    << entry.topo_line << "\n";
                success++;
                file_ok++;
                node_counts[entry.node_count]++;
            } else {
                fail++;
                file_fail++;
                if (verbose) {
                    std::cerr << "  FAIL: " << line.substr(0, 80) << "...\n";
                }
            }
        }
        
        if (verbose) {
            std::cout << "  " << fs::path(file).filename().string() 
                      << ": " << file_ok << " ok, " << file_fail << " fail\n";
        }
    }
    
    // Summary
    out << "# Summary: " << success << " Nk entries\n";
    for (const auto& [k, cnt] : node_counts) {
        out << "# N" << k << ": " << cnt << "\n";
    }
    
    out.close();
    
    std::cout << "\n=== Conversion Complete ===\n";
    std::cout << "Unique lines: " << total_lines << "\n";
    std::cout << "Skipped duplicates: " << skipped_dup << "\n";
    std::cout << "Success: " << success << "\n";
    std::cout << "Failed: " << fail << "\n";
    std::cout << "\nBreakdown by node count:\n";
    for (const auto& [k, cnt] : node_counts) {
        std::cout << "  N" << k << ": " << cnt << "\n";
    }
    std::cout << "\nOutput: " << output << "\n";
    
    return 0;
}
