// nhc_filter.cpp
// Read topocompactline files → build IF → run NHC check → output pass/fail
//
// Usage:
//   ./nhc_filter <input_dir_or_file> [output.cat] [--verbose] [--pass-only] [--fail-only]
//
// Input: directory of *_dedup files OR single file with topocompactline entries
// Output: catalog with NHC pass/fail annotation
//
// Compile:
//   g++ -std=c++17 -O2 -o nhc_filter nhc_filter.cpp Tensor.C \
//       TopoLineCompact_enhanced.cpp Topology_enhanced.cpp \
//       -I$(brew --prefix eigen)/include/eigen3 -I.

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
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Tensor.h"
#include "Theory_enhanced.h"
#include "TopoLineCompact_enhanced.hpp"
#include "Topology_enhanced.h"
#include "anomaly_tables.h"

namespace fs = std::filesystem;

// ============================================================================
// Eigenvalue Spectrum Computation (from no_node_theory_v2.cpp)
// ============================================================================

std::vector<double> compute_eigenvalue_spectrum(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    if (n == 0) return {};
    Eigen::MatrixXd IF_double = IF.cast<double>();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(IF_double);
    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    std::vector<double> spectrum(n);
    for (int i = 0; i < n; ++i) spectrum[i] = eigenvalues(i);
    std::sort(spectrum.begin(), spectrum.end());
    return spectrum;
}

bool are_spectra_equal(const std::vector<double>& s1, const std::vector<double>& s2, double tol = 1e-9) {
    if (s1.size() != s2.size()) return false;
    for (size_t i = 0; i < s1.size(); ++i)
        if (std::abs(s1[i] - s2[i]) > tol) return false;
    return true;
}

// ============================================================================
// Gauge info for port gauge lookup
// ============================================================================

struct PortGaugeInfo {
    const char* gauge_name;
    int dim;        // dimension of gauge algebra
    int h_dual;     // dual Coxeter number
};

// Central charge: c(k) = k * dim / (k + h∨)
// k = intersection number between -1 curve and the neighbor curve
double central_charge(const PortGaugeInfo& g, int k) {
    if (g.dim == 0) return 0.0;
    return (double)(k * g.dim) / (double)(k + g.h_dual);
}

const PortGaugeInfo GAUGE_NONE  = {"none",    0,  0};
const PortGaugeInfo GAUGE_SP1   = {"sp1",     3,  2};
const PortGaugeInfo GAUGE_SU2   = {"su2",     3,  2};
const PortGaugeInfo GAUGE_SU3   = {"su3",     8,  3};
const PortGaugeInfo GAUGE_G2    = {"g2",     14,  4};
const PortGaugeInfo GAUGE_SO7   = {"so7",    21,  5};
const PortGaugeInfo GAUGE_SO8   = {"so8",    28,  6};
const PortGaugeInfo GAUGE_SO9   = {"so9",    36,  7};
const PortGaugeInfo GAUGE_SO10  = {"so10",   45,  8};
const PortGaugeInfo GAUGE_SO11  = {"so11",   55,  9};
const PortGaugeInfo GAUGE_SO12  = {"so12",   66, 10};
const PortGaugeInfo GAUGE_F4    = {"f4",     52,  9};
const PortGaugeInfo GAUGE_E6    = {"e6",     78, 12};
const PortGaugeInfo GAUGE_E7    = {"e7",    133, 18};
const PortGaugeInfo GAUGE_E8    = {"e8",    248, 30};

// ============================================================================
// NHC Definition with Port Gauges
// ============================================================================

struct NHCDefinition {
    std::string name;
    bool is_linear;
    std::vector<int> self_ints;
    std::vector<PortGaugeInfo> gauges;
    std::vector<double> spectrum;
    Eigen::MatrixXi IF_template;
};

std::vector<NHCDefinition> build_nhc_table() {
    std::vector<NHCDefinition> table;

    auto add_linear = [&](const std::string& name,
                          std::vector<int> self_ints,
                          std::vector<PortGaugeInfo> gauges) {
        NHCDefinition def;
        def.name = name;
        def.is_linear = true;
        def.self_ints = self_ints;
        def.gauges = gauges;

        int n = self_ints.size();
        Eigen::MatrixXi IF = Eigen::MatrixXi::Zero(n, n);
        for (int i = 0; i < n; ++i) {
            IF(i, i) = self_ints[i];
            if (i + 1 < n) { IF(i, i+1) = 1; IF(i+1, i) = 1; }
        }
        def.IF_template = IF;
        def.spectrum = compute_eigenvalue_spectrum(IF);
        table.push_back(def);
    };

    // Single curves
    add_linear("-3",   {-3},   {GAUGE_SU3});
    add_linear("-4",   {-4},   {GAUGE_SO8});
    add_linear("-5",   {-5},   {GAUGE_F4});
    add_linear("-6",   {-6},   {GAUGE_E6});
    add_linear("-7",   {-7},   {GAUGE_E7});
    add_linear("-8",   {-8},   {GAUGE_E7});
    add_linear("-12",  {-12},  {GAUGE_E8});

    // Multi-curve NHCs
    add_linear("-2-3",   {-2, -3},       {GAUGE_SU2, GAUGE_G2});
    add_linear("-2-3-2", {-2, -3, -2},   {GAUGE_SU2, GAUGE_SO7, GAUGE_SU2});
    add_linear("-2-2-3", {-2, -2, -3},   {GAUGE_NONE, GAUGE_SP1, GAUGE_G2});

    return table;
}

const std::vector<NHCDefinition>& get_nhc_table() {
    static std::vector<NHCDefinition> cache = build_nhc_table();
    return cache;
}

bool is_pure_minus2_cluster(const std::vector<int>& diag) {
    for (int d : diag) if (d != -2) return false;
    return true;
}

const NHCDefinition* identify_nhc(const std::vector<double>& spectrum) {
    for (const auto& def : get_nhc_table())
        if (are_spectra_equal(spectrum, def.spectrum)) return &def;
    return nullptr;
}

// ============================================================================
// Chain walk and port gauge lookup
// ============================================================================

std::vector<int> walk_linear_chain(const Eigen::MatrixXi& IF,
                                   const std::vector<int>& cluster_curves) {
    int m = cluster_curves.size();
    if (m <= 1) return cluster_curves;

    std::map<int, std::vector<int>> adj;
    for (int a : cluster_curves)
        for (int b : cluster_curves)
            if (a != b && IF(a, b) != 0) adj[a].push_back(b);

    int start = cluster_curves[0];
    for (int c : cluster_curves)
        if ((int)adj[c].size() <= 1) { start = c; break; }

    std::vector<int> ordered;
    std::set<int> visited_set;
    int cur = start;
    while (true) {
        ordered.push_back(cur);
        visited_set.insert(cur);
        bool found_next = false;
        for (int nb : adj[cur]) {
            if (visited_set.find(nb) == visited_set.end()) {
                cur = nb;
                found_next = true;
                break;
            }
        }
        if (!found_next) break;
    }
    return ordered;
}

PortGaugeInfo lookup_port_gauge(const Eigen::MatrixXi& IF,
                                const std::vector<int>& cluster_curves,
                                int target_global_idx,
                                const NHCDefinition* nhc_def) {
    if (!nhc_def) return GAUGE_NONE;

    if (nhc_def->is_linear) {
        auto ordered = walk_linear_chain(IF, cluster_curves);
        std::vector<int> chain_si;
        for (int c : ordered) chain_si.push_back(IF(c, c));

        // Forward match
        if (chain_si == nhc_def->self_ints) {
            for (size_t i = 0; i < ordered.size(); ++i)
                if (ordered[i] == target_global_idx) return nhc_def->gauges[i];
        }

        // Reverse match
        std::vector<int> rev_si(chain_si.rbegin(), chain_si.rend());
        if (rev_si == nhc_def->self_ints) {
            int m2 = ordered.size();
            for (size_t i = 0; i < ordered.size(); ++i)
                if (ordered[i] == target_global_idx) return nhc_def->gauges[m2 - 1 - i];
        }
    }

    int si = IF(target_global_idx, target_global_idx);
    for (size_t i = 0; i < nhc_def->self_ints.size(); ++i)
        if (nhc_def->self_ints[i] == si) return nhc_def->gauges[i];

    return GAUGE_NONE;
}

// ============================================================================
// NHC Check Result
// ============================================================================

struct NHCCheckResult {
    bool passes;
    std::string failure_reason;
    bool layer1_pass = true;
    bool layer2_pass = true;
    bool layer3_pass = true;
    int n_clusters = 0;
    std::vector<std::vector<int>> clusters;
};

// ============================================================================
// Layer 1: Cluster Recognition
// ============================================================================

bool check_layer1_clusters(const Eigen::MatrixXi& IF,
                           const std::vector<std::vector<int>>& clusters,
                           std::map<int, const NHCDefinition*>& curve_to_nhc,
                           std::string& fail_reason) {
    for (const auto& cluster : clusters) {
        int m = cluster.size();
        Eigen::MatrixXi sub(m, m);
        for (int a = 0; a < m; ++a)
            for (int b = 0; b < m; ++b)
                sub(a, b) = IF(cluster[a], cluster[b]);

        auto spec = compute_eigenvalue_spectrum(sub);

        std::vector<int> diag;
        for (int c : cluster) diag.push_back(IF(c, c));

        if (is_pure_minus2_cluster(diag)) {
            for (int c : cluster) curve_to_nhc[c] = nullptr;
            continue;
        }

        const NHCDefinition* nhc = identify_nhc(spec);
        if (!nhc) {
            std::ostringstream ss;
            ss << "unrecognized cluster {";
            for (int a = 0; a < m; ++a) {
                if (a > 0) ss << ",";
                ss << IF(cluster[a], cluster[a]);
            }
            ss << "} spectrum [";
            for (size_t a = 0; a < spec.size(); ++a) {
                if (a > 0) ss << ",";
                ss << std::fixed << std::setprecision(4) << spec[a];
            }
            ss << "]";
            fail_reason = ss.str();
            return false;
        }

        for (int c : cluster) curve_to_nhc[c] = nhc;
    }
    return true;
}

// ============================================================================
// E8 Subalgebra Pair Table
// ============================================================================

bool is_e8_subalgebra_pair(const std::string& g1, const std::string& g2) {
    if (g1 == "none" || g2 == "none") return true;

    static const std::set<std::pair<std::string,std::string>> allowed = {
        {"sp1", "sp1"}, {"sp1", "su2"}, {"sp1", "su3"},
        {"sp1", "g2"}, {"sp1", "so7"}, {"sp1", "so8"},
        {"sp1", "so9"}, {"sp1", "so10"}, {"sp1", "so11"},
        {"sp1", "so12"}, {"sp1", "f4"}, {"sp1", "e6"}, {"sp1", "e7"},

        {"su2", "su2"}, {"su2", "su3"}, {"su2", "g2"},
        {"su2", "so7"}, {"su2", "so8"}, {"su2", "so9"},
        {"su2", "so10"}, {"su2", "so11"}, {"su2", "so12"},
        {"su2", "f4"}, {"su2", "e6"}, {"su2", "e7"},

        {"su3", "su3"}, {"su3", "g2"}, {"su3", "so7"},
        {"su3", "so8"}, {"su3", "so9"}, {"su3", "so10"},
        {"su3", "f4"}, {"su3", "e6"},

        {"g2", "g2"}, {"g2", "so7"}, {"g2", "so8"},
        {"g2", "so9"}, {"g2", "f4"},

        {"so7", "su3"},
        {"so8", "su3"},
        {"so10", "su3"},
        {"f4", "su3"}, {"f4", "g2"},
        {"e6", "su3"},
    };

    return allowed.count({g1, g2}) > 0 || allowed.count({g2, g1}) > 0;
}

// ============================================================================
// Layer 2: E8 Gauge Constraint
// ============================================================================

bool check_layer2_e8_constraint(const Eigen::MatrixXi& IF,
                                const std::vector<int>& minus1_curves,
                                const std::vector<std::vector<int>>& clusters,
                                const std::map<int, const NHCDefinition*>& curve_to_nhc,
                                int T_tensor,
                                std::string& fail_reason) {
    int n = IF.rows();
    // E8 level-1 central charge = 248/31 = 8
    const double E8_BUDGET = 8.0;

    std::map<int, int> curve_to_cluster_idx;
    for (size_t ci = 0; ci < clusters.size(); ++ci)
        for (int c : clusters[ci])
            curve_to_cluster_idx[c] = ci;

    for (int m1 : minus1_curves) {
        std::vector<std::pair<PortGaugeInfo, int>> neighbor_info;  // (gauge, k)
        double total_c = 0.0;

        for (int j = 0; j < n; ++j) {
            if (j == m1 || IF(m1, j) == 0) continue;
            if (IF(j, j) == -1) continue;

            // k = intersection number between -1 curve and neighbor
            int k = std::abs(IF(m1, j));

            auto it = curve_to_cluster_idx.find(j);
            if (it == curve_to_cluster_idx.end()) {
                neighbor_info.push_back({GAUGE_NONE, k});
                continue;
            }

            int ci = it->second;
            const NHCDefinition* nhc = nullptr;
            auto nhc_it = curve_to_nhc.find(j);
            if (nhc_it != curve_to_nhc.end()) nhc = nhc_it->second;

            auto pg = lookup_port_gauge(IF, clusters[ci], j, nhc);
            neighbor_info.push_back({pg, k});
            total_c += central_charge(pg, k);
        }

        // E8 subalgebra pair check
        for (size_t i = 0; i < neighbor_info.size(); ++i) {
            for (size_t j = i + 1; j < neighbor_info.size(); ++j) {
                if (!is_e8_subalgebra_pair(neighbor_info[i].first.gauge_name,
                                           neighbor_info[j].first.gauge_name)) {
                    std::ostringstream ss;
                    ss << "E8 subalgebra violation at -1 curve " << m1
                       << ": " << neighbor_info[i].first.gauge_name
                       << " x " << neighbor_info[j].first.gauge_name
                       << " not in E8";
                    fail_reason = ss.str();
                    return false;
                }
            }
        }

        // Central charge budget: sum c_g(k) <= 8
        if (total_c > E8_BUDGET + 1e-9) {
            std::ostringstream ss;
            ss << "central charge budget at -1 curve " << m1
               << ": sum = " << std::fixed << std::setprecision(2) << total_c 
               << " > " << E8_BUDGET << " (gauges:";
            for (const auto& [pg, k] : neighbor_info) 
                ss << " " << pg.gauge_name << "(k=" << k 
                   << ",c=" << std::setprecision(2) << central_charge(pg, k) << ")";
            ss << ")";
            fail_reason = ss.str();
            return false;
        }
    }
    return true;
}

// ============================================================================
// Layer 3: Flavor Compatibility (placeholder)
// ============================================================================

bool check_layer3_flavor_compat(const Eigen::MatrixXi& IF,
                                const std::vector<int>& minus1_curves,
                                const std::vector<std::vector<int>>& clusters,
                                const std::map<int, const NHCDefinition*>& curve_to_nhc,
                                std::string& fail_reason) {
    (void)IF; (void)minus1_curves; (void)clusters; (void)curve_to_nhc; (void)fail_reason;
    return true;
}

// ============================================================================
// Main NHC Check (from no_node_theory_v2.cpp)
// ============================================================================

NHCCheckResult check_nhc_clusters(const Eigen::MatrixXi& IF,
                                  bool run_l2 = true, bool run_l3 = true) {
    NHCCheckResult result;
    result.passes = true;

    int n = IF.rows();
    if (n == 0) return result;

    // Find -1 curves and non-(-1) curves
    std::vector<int> minus1_curves;
    std::vector<int> non_minus1;
    for (int i = 0; i < n; ++i) {
        if (IF(i, i) == -1) minus1_curves.push_back(i);
        else non_minus1.push_back(i);
    }

    if (non_minus1.empty()) return result;

    // BFS: connected components of non-(-1) curves
    std::vector<bool> visited(n, false);
    for (int start_idx : non_minus1) {
        if (visited[start_idx]) continue;
        std::vector<int> component;
        std::vector<int> queue = {start_idx};
        visited[start_idx] = true;
        while (!queue.empty()) {
            int cur = queue.back();
            queue.pop_back();
            component.push_back(cur);
            for (int j : non_minus1) {
                if (!visited[j] && IF(cur, j) != 0) {
                    visited[j] = true;
                    queue.push_back(j);
                }
            }
        }
        std::sort(component.begin(), component.end());
        result.clusters.push_back(component);
    }
    result.n_clusters = result.clusters.size();

    // Layer 1
    std::map<int, const NHCDefinition*> curve_to_nhc;
    std::string fail1;
    result.layer1_pass = check_layer1_clusters(IF, result.clusters, curve_to_nhc, fail1);
    if (!result.layer1_pass) {
        result.passes = false;
        result.failure_reason = "L1: " + fail1;
        return result;
    }

    // Layer 2
    if (run_l2) {
        int T_tensor = n;
        std::string fail2;
        result.layer2_pass = check_layer2_e8_constraint(IF, minus1_curves, result.clusters,
                                                        curve_to_nhc, T_tensor, fail2);
        if (!result.layer2_pass) {
            result.passes = false;
            result.failure_reason = "L2: " + fail2;
        }
    }

    // Layer 3
    if (run_l3) {
        std::string fail3;
        result.layer3_pass = check_layer3_flavor_compat(IF, minus1_curves, result.clusters,
                                                         curve_to_nhc, fail3);
        if (!result.layer3_pass) {
            result.passes = false;
            if (!result.failure_reason.empty()) result.failure_reason += "; ";
            result.failure_reason += "L3: " + fail3;
        }
    }

    return result;
}

// ============================================================================
// Build IF from Topology (from lst_to_catalog.cpp)
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
    std::vector<NodeRef> blockNodes, sideNodes, instNodes, extNodes;

    for (size_t i = 0; i < T.block.size(); ++i)
        blockNodes.push_back(G.add(make_spec(T.block[i].kind, T.block[i].param)));
    for (size_t i = 0; i < T.side_links.size(); ++i)
        sideNodes.push_back(G.add(s(T.side_links[i].param)));
    for (size_t i = 0; i < T.instantons.size(); ++i)
        instNodes.push_back(G.add(s(T.instantons[i].param)));
    for (size_t i = 0; i < T.externals.size(); ++i)
        extNodes.push_back(G.add(e(T.externals[i].param)));

    for (const auto& conn : T.l_connection)
        if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
            conn.v >= 0 && conn.v < (int)blockNodes.size())
            G.connect(blockNodes[conn.u], blockNodes[conn.v]);
    for (const auto& conn : T.s_connection)
        if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
            conn.v >= 0 && conn.v < (int)sideNodes.size())
            G.connect(sideNodes[conn.v], blockNodes[conn.u]);
    for (const auto& conn : T.i_connection)
        if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
            conn.v >= 0 && conn.v < (int)instNodes.size())
            G.connect(instNodes[conn.v], blockNodes[conn.u]);
    for (const auto& conn : T.e_connection) {
        if (conn.external_id < 0 || conn.external_id >= (int)extNodes.size()) continue;
        NodeRef* p = nullptr;
        if (conn.parent_type == 0 && conn.parent_id >= 0 && conn.parent_id < (int)blockNodes.size())
            p = &blockNodes[conn.parent_id];
        else if (conn.parent_type == 1 && conn.parent_id >= 0 && conn.parent_id < (int)sideNodes.size())
            p = &sideNodes[conn.parent_id];
        else if (conn.parent_type == 2 && conn.parent_id >= 0 && conn.parent_id < (int)instNodes.size())
            p = &instNodes[conn.parent_id];
        if (p) G.connect(extNodes[conn.external_id], AttachmentPoint(-1), *p, AttachmentPoint(conn.port_idx));
    }
    return G.ComposeIF_Gluing();
}

// ============================================================================
// Signature computation
// ============================================================================

int compute_sig_neg(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    if (n == 0) return 0;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
    int neg = 0;
    for (int i = 0; i < n; ++i)
        if (es.eigenvalues()(i) < -1e-8) neg++;
    return neg;
}

// ============================================================================
// File collection (recursive directory scan)
// ============================================================================

std::vector<std::string> collect_dedup_files(const std::string& path) {
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto& entry : fs::recursive_directory_iterator(path)) {
            if (entry.is_regular_file()) {
                std::string name = entry.path().filename().string();
                if (name.find("_dedup") != std::string::npos && name.find(".") == std::string::npos)
                    files.push_back(entry.path().string());
            }
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

// ============================================================================
// Count nodes
// ============================================================================

int count_nodes(const Topology_enhanced& T) {
    int count = 0;
    for (const auto& b : T.block)
        if (b.kind == LKind::g) count++;
    return count;
}

// ============================================================================
// Main
// ============================================================================

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <input_dir_or_file> [output.cat] [options]\n"
              << "Options:\n"
              << "  --verbose       Show per-file progress\n"
              << "  --pass-only     Output only NHC-passing entries\n"
              << "  --fail-only     Output only NHC-failing entries\n"
              << "  --summary       Print summary statistics only\n"
              << "  --l12           Layer 1+2 only (skip L3)\n"
              << "  --l13           Layer 1+3 only (skip L2)\n"
              << "  --no-l2         Disable Layer 2\n"
              << "  --no-l3         Disable Layer 3\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string input_path = argv[1];
    std::string output_file = "nhc_filtered.cat";
    bool verbose = false, pass_only = false, fail_only = false, summary_only = false;
    bool use_l1 = true, use_l2 = true, use_l3 = true;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--verbose") verbose = true;
        else if (arg == "--pass-only") pass_only = true;
        else if (arg == "--fail-only") fail_only = true;
        else if (arg == "--summary") summary_only = true;
        else if (arg == "--no-l2") use_l2 = false;
        else if (arg == "--no-l3") use_l3 = false;
        else if (arg == "--l12") { use_l1 = true; use_l2 = true; use_l3 = false; }
        else if (arg == "--l13") { use_l1 = true; use_l2 = false; use_l3 = true; }
        else output_file = arg;
    }

    auto files = collect_dedup_files(input_path);
    std::cout << "Found " << files.size() << " dedup file(s)\n";

    std::ofstream out;
    if (!summary_only) {
        out.open(output_file);
        out << "# NHC-Filtered topocompactline entries\n";
        out << "# Pass-only: " << (pass_only ? "yes" : "no") 
            << ", Fail-only: " << (fail_only ? "yes" : "no") << "\n";
    }

    int total = 0, nhc_pass = 0, nhc_fail = 0;
    int l1_fail = 0, l2_fail = 0, l3_fail = 0;
    int id = 0;
    std::set<std::string> seen;

    // Failure reason statistics
    std::map<std::string, int> fail_reasons;

    for (const auto& filepath : files) {
        std::ifstream fin(filepath);
        if (!fin) continue;

        std::string line;
        int file_ok = 0, file_fail = 0;

        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;

            // Dedup
            if (seen.count(line)) continue;
            seen.insert(line);

            total++;

            Topology_enhanced topo;
            if (!TopoLineCompact_enhanced::deserialize(line, topo)) continue;

            // Build IF
            Eigen::MatrixXi IF = build_IF_from_topology(topo);
            int T = compute_sig_neg(IF);
            int nk = count_nodes(topo);
            std::string type = "N" + std::to_string(nk);

            // Compute H, V from anomaly tables
            int H = 0, V = 0;
            for (const auto& b : topo.block) {
                if (b.kind == LKind::g) {
                    auto c = get_node_contrib(b.param);
                    H += c.H; V += c.V;
                }
            }
            for (const auto& sl : topo.side_links) {
                auto c = get_sidelink_contrib(sl.param);
                H += c.H; V += c.V;
            }
            for (const auto& inst : topo.instantons) {
                auto c = get_sidelink_contrib(inst.param);
                H += c.H; V += c.V;
            }
            for (const auto& ext : topo.externals) {
                auto c = get_external_contrib(ext.param);
                H += c.H; V += c.V;
            }

            // NHC check
            auto nhc_result = check_nhc_clusters(IF, use_l2, use_l3);

            if (nhc_result.passes) {
                nhc_pass++;
                file_ok++;
            } else {
                nhc_fail++;
                file_fail++;
                if (!nhc_result.layer1_pass) l1_fail++;
                if (!nhc_result.layer2_pass) l2_fail++;
                if (!nhc_result.layer3_pass) l3_fail++;
                fail_reasons[nhc_result.failure_reason]++;
            }

            // Output
            if (!summary_only) {
                bool should_write = true;
                if (pass_only && !nhc_result.passes) should_write = false;
                if (fail_only && nhc_result.passes) should_write = false;

                if (should_write) {
                    out << line << "\n";
                    id++;
                }
            }
        }

        if (verbose) {
            std::string fname = fs::path(filepath).filename().string();
            std::cout << "  " << fname << ": " << file_ok << " pass, " << file_fail << " fail\n";
        }
    }

    // Summary
    std::cout << "\n=== NHC Filter Summary ===\n";
    std::cout << "Total unique entries: " << total << "\n";
    std::cout << "NHC pass: " << nhc_pass << " (" << std::fixed << std::setprecision(1)
              << (total > 0 ? 100.0 * nhc_pass / total : 0) << "%)\n";
    std::cout << "NHC fail: " << nhc_fail << " (" << std::fixed << std::setprecision(1)
              << (total > 0 ? 100.0 * nhc_fail / total : 0) << "%)\n";
    std::cout << "  Layer 1 (cluster recognition) fail: " << l1_fail << "\n";
    std::cout << "  Layer 2 (E8 constraint) fail: " << l2_fail << "\n";
    std::cout << "  Layer 3 (flavor compat) fail: " << l3_fail << "\n";

    if (!fail_reasons.empty()) {
        std::cout << "\nTop failure reasons:\n";
        std::vector<std::pair<int, std::string>> sorted_reasons;
        for (auto& [reason, count] : fail_reasons)
            sorted_reasons.push_back({count, reason});
        std::sort(sorted_reasons.rbegin(), sorted_reasons.rend());
        for (int i = 0; i < std::min(20, (int)sorted_reasons.size()); i++)
            std::cout << "  " << sorted_reasons[i].first << "x: " << sorted_reasons[i].second << "\n";
    }

    if (!summary_only) {
        out.close();
        std::cout << "\nOutput: " << output_file << " (" << id << " entries)\n";
    }

    return 0;
}
