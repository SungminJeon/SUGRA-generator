// sugra_generator.h
// SUGRA base construction from LST catalog
//
// Pipeline:
//   1. Load LST catalog → reconstruct IF (via theory_sugra.h)
//   2. Attach external curves (hardcoded rules)
//   3. Optionally glue dummy LSTs (A-type: 1222...221)
//   4. NHC check: all non-(-1) clusters must be recognized NHCs
//   5. H, V from NHC table → gravitational anomaly: H - V + 29T - 273 = 0

#pragma once

#include "theory_sugra.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <functional>

// ############################################################################
// PART 1: GAUGE DATA + CENTRAL CHARGE
// ############################################################################

struct GaugeInfo {
    std::string name;
    int dim;
    int h_dual;
    int rank;
};

inline const GaugeInfo GAUGE_NONE = {"none", 0, 0, 0};
inline const GaugeInfo GAUGE_SU2 = {"su2",  3, 2, 1};
inline const GaugeInfo GAUGE_SP1 = {"sp1",  3, 2, 1};
inline const GaugeInfo GAUGE_SU3 = {"su3",  8, 3, 2};
inline const GaugeInfo GAUGE_G2  = {"g2",  14, 4, 2};
inline const GaugeInfo GAUGE_SO7 = {"so7", 21, 5, 3};
inline const GaugeInfo GAUGE_SO8 = {"so8", 28, 6, 4};
inline const GaugeInfo GAUGE_F4  = {"f4",  52, 9, 4};
inline const GaugeInfo GAUGE_E6  = {"e6",  78, 12, 6};
inline const GaugeInfo GAUGE_E7  = {"e7", 133, 18, 7};
inline const GaugeInfo GAUGE_E8  = {"e8", 248, 30, 8};
inline const GaugeInfo GAUGE_SO16 = {"so16", 120, 14, 8};
inline const GaugeInfo GAUGE_SU8 = {"su8", 63, 8, 7};
inline const GaugeInfo GAUGE_SU16 = {"su16", 255, 16, 15};

inline double central_charge(const GaugeInfo& g, int k) {
    if (g.dim == 0) return 0.0;
    return (double)(k * g.dim) / (double)(k + g.h_dual);
}

inline GaugeInfo gauge_from_si(int si) {
    switch (si) {
        case -3:  return GAUGE_SU3;
        case -4:  return GAUGE_SO8;
        case -5:  return GAUGE_F4;
        case -6:  return GAUGE_E6;
        case -7: case -8: return GAUGE_E7;
        case -12: return GAUGE_E8;
        default:  return GAUGE_NONE;
    }
}

// ############################################################################
// PART 2: NHC TABLE
// ############################################################################
//
// Each NHC has:
//   - self_ints: the chain of self-intersections
//   - gauges: gauge algebra per curve
//   - H, V: total hypermultiplet / vector multiplet contribution
//   - spectrum: eigenvalues for matching (auto-computed)
//
// After removing -1 curves, every connected component of non-(-1)
// curves must match one of these NHCs.
// Total anomaly: sum(H_i) - sum(V_i) + 29T - 273 = 0

struct NHCDef {
    std::string name;
    std::vector<int> self_ints;
    std::vector<GaugeInfo> gauges;
    int H;                              // hypermultiplet contribution
    int V;                              // vector multiplet contribution
    std::vector<double> spectrum;       // auto-computed
};

inline std::vector<double> eigenvalue_spectrum(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    if (n == 0) return {};
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
    std::vector<double> s(n);
    for (int i = 0; i < n; i++) s[i] = es.eigenvalues()(i);
    std::sort(s.begin(), s.end());
    return s;
}

inline bool spectra_equal(const std::vector<double>& a, const std::vector<double>& b, double tol = 1e-9) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++)
        if (std::abs(a[i] - b[i]) > tol) return false;
    return true;
}

inline const std::vector<NHCDef>& get_nhc_table() {
    static std::vector<NHCDef> table = []() {
        std::vector<NHCDef> t;
        
        auto add = [&](const std::string& name,
                       std::vector<int> si,
                       std::vector<GaugeInfo> g,
                       int H, int V,
                       std::vector<int> ints = {}) {
            NHCDef def;
            def.name = name;
            def.self_ints = si;
            def.gauges = g;
            def.H = H;
            def.V = V;
            int n = si.size();
            Eigen::MatrixXi IF = Eigen::MatrixXi::Zero(n, n);
            for (int i = 0; i < n; i++) {
                IF(i, i) = si[i];
                if (i + 1 < n) {
                    int k = (i < (int)ints.size()) ? ints[i] : 1;
                    IF(i, i+1) = k; IF(i+1, i) = k;
                }
            }
            def.spectrum = eigenvalue_spectrum(IF);
            t.push_back(def);
        };
        
        // ┌─────────────────────────────────────────────────────────┐
        // │ NHC Table: name, {self_ints}, {gauges}, H, V           │
        // │                                                         │
        // │ USER: fill in correct H, V for each NHC                │
        // │ H=0, V=0 are placeholders                              │
        // └─────────────────────────────────────────────────────────┘
        
        // Single curves
        add("-3",     {-3},         {GAUGE_SU3},                    0, 8);
        add("-4",     {-4},         {GAUGE_SO8},                    0, 28);
        add("-5",     {-5},         {GAUGE_F4},                     0, 52);
        add("-6",     {-6},         {GAUGE_E6},                     0, 78);
        add("-7",     {-7},         {GAUGE_E7},                     28, 133);
        add("-8",     {-8},         {GAUGE_E7},                     0, 133);
        add("-12",    {-12},        {GAUGE_E8},                     0, 248);
        
        // Multi-curve chains
        add("-2-3",   {-2, -3},     {GAUGE_SU2, GAUGE_G2},         8, 17);
        add("-2-3-2", {-2, -3, -2}, {GAUGE_SU2, GAUGE_SO7, GAUGE_SU2}, 16, 27);
        add("-2-2-3", {-2, -2, -3}, {GAUGE_NONE, GAUGE_SP1, GAUGE_G2}, 8, 17);
        
        // USER: add more NHC chains as needed
        // add("name", {self_ints...}, {gauges...}, H, V);
	add("-2-4", {-2, -4}, {GAUGE_SU8, GAUGE_SO16},             256,183, {2});

        
        return t;
    }();
    return table;
}

// ############################################################################
// PART 3: NHC CHECK → CLUSTER ASSIGNMENT
// ############################################################################

struct ClusterInfo {
    int cluster_id;
    std::vector<int> curve_indices;
    const NHCDef* nhc;             // nullptr for pure -2 clusters
    int H, V;
};

struct NHCResult {
    bool passes;
    std::string fail_reason;
    std::vector<ClusterInfo> clusters;
    
    // Per-curve info
    std::vector<GaugeInfo> curve_gauges;    // indexed by curve in IF
    std::vector<int> curve_cluster;         // cluster_id per curve (-1 for -1 curves)
    
    int total_H() const {
        int h = 0;
        for (auto& c : clusters) h += c.H;
        return h;
    }
    int total_V() const {
        int v = 0;
        for (auto& c : clusters) v += c.V;
        return v;
    }
};

inline std::vector<int> walk_chain(const Eigen::MatrixXi& IF, const std::vector<int>& cluster) {
    int m = cluster.size();
    if (m <= 1) return cluster;
    std::map<int, std::vector<int>> adj;
    for (int a : cluster)
        for (int b : cluster)
            if (a != b && IF(a, b) != 0) adj[a].push_back(b);
    int start = cluster[0];
    for (int c : cluster) if ((int)adj[c].size() <= 1) { start = c; break; }
    std::vector<int> ordered;
    std::set<int> vis;
    int cur = start;
    while (true) {
        ordered.push_back(cur);
        vis.insert(cur);
        bool found = false;
        for (int nb : adj[cur])
            if (vis.find(nb) == vis.end()) { cur = nb; found = true; break; }
        if (!found) break;
    }
    return ordered;
}

inline NHCResult check_nhc(const Eigen::MatrixXi& IF, double cc_budget = 8.0) {
    NHCResult result;
    result.passes = true;
    int n = IF.rows();
    if (n == 0) return result;
    
    result.curve_gauges.assign(n, GAUGE_NONE);
    result.curve_cluster.assign(n, -1);
    
    // Separate -1 and non-(-1) curves
    std::vector<int> minus1, non_minus1;
    for (int i = 0; i < n; i++) {
        if (IF(i, i) == -1) minus1.push_back(i);
        else non_minus1.push_back(i);
    }
    if (non_minus1.empty()) return result;
    
    // BFS: connected components of non-(-1) curves
    std::vector<std::vector<int>> comps;
    std::vector<bool> visited(n, false);
    for (int start : non_minus1) {
        if (visited[start]) continue;
        std::vector<int> comp;
        std::vector<int> queue = {start};
        visited[start] = true;
        while (!queue.empty()) {
            int cur = queue.back(); queue.pop_back();
            comp.push_back(cur);
            for (int j : non_minus1)
                if (!visited[j] && IF(cur, j) != 0) { visited[j] = true; queue.push_back(j); }
        }
        std::sort(comp.begin(), comp.end());
        comps.push_back(comp);
    }
    
    // Match each component to NHC table
    for (int ci = 0; ci < (int)comps.size(); ci++) {
        auto& comp = comps[ci];
        int m = comp.size();
        
        // Extract submatrix
        Eigen::MatrixXi sub(m, m);
        for (int a = 0; a < m; a++)
            for (int b = 0; b < m; b++)
                sub(a, b) = IF(comp[a], comp[b]);
        
        auto spec = eigenvalue_spectrum(sub);
        
        // Check pure -2
        bool pure2 = true;
        for (int c : comp) if (IF(c, c) != -2) { pure2 = false; break; }
        
        if (pure2) {
            ClusterInfo cl;
            cl.cluster_id = ci;
            cl.curve_indices = comp;
            cl.nhc = nullptr;
            cl.H = 0;
            cl.V = 0;
            result.clusters.push_back(cl);
            for (int c : comp) result.curve_cluster[c] = ci;
            continue;
        }
        
        // NHC lookup
        const NHCDef* nhc = nullptr;
        for (const auto& def : get_nhc_table())
            if (spectra_equal(spec, def.spectrum)) { nhc = &def; break; }
        
        if (!nhc) {
            result.passes = false;
            std::ostringstream ss;
            ss << "L1: unknown cluster {";
            for (int a = 0; a < m; a++) { if (a) ss << ","; ss << IF(comp[a], comp[a]); }
            ss << "}";
            result.fail_reason = ss.str();
            return result;
        }
        
        // Assign gauges per curve
        auto ordered = walk_chain(IF, comp);
        std::vector<int> chain_si;
        for (int c : ordered) chain_si.push_back(IF(c, c));
        
        bool forward = (chain_si == nhc->self_ints);
        std::vector<int> rev(chain_si.rbegin(), chain_si.rend());
        bool reverse = (!forward && rev == nhc->self_ints);
        
        if (forward) {
            for (size_t i = 0; i < ordered.size(); i++) {
                result.curve_gauges[ordered[i]] = nhc->gauges[i];
                result.curve_cluster[ordered[i]] = ci;
            }
        } else if (reverse) {
            int sz = ordered.size();
            for (size_t i = 0; i < ordered.size(); i++) {
                result.curve_gauges[ordered[i]] = nhc->gauges[sz - 1 - i];
                result.curve_cluster[ordered[i]] = ci;
            }
        } else {
            // Fallback: assign by self-int
            for (int c : comp) {
                result.curve_gauges[c] = gauge_from_si(IF(c, c));
                result.curve_cluster[c] = ci;
            }
        }
        
        ClusterInfo cl;
        cl.cluster_id = ci;
        cl.curve_indices = comp;
        cl.nhc = nhc;
        cl.H = nhc->H;
        cl.V = nhc->V;
        result.clusters.push_back(cl);
    }
    
    // Central charge check on each -1 curve
    for (int m1 : minus1) {
        double total_c = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == m1 || IF(m1, j) == 0 || IF(j, j) == -1) continue;
            int k = std::abs(IF(m1, j));
            total_c += central_charge(result.curve_gauges[j], k);
        }
        if (total_c > cc_budget + 1e-9) {
            result.passes = false;
            std::ostringstream ss;
            ss << "L2: c.c.=" << std::fixed << std::setprecision(2)
               << total_c << " > " << cc_budget << " at -1 curve " << m1;
            result.fail_reason = ss.str();
            return result;
        }
    }
    
    return result;
}

// ── External -2 → -1 gauge enhancement ──
// When an external -2 curve is attached to a -1 curve, it carries su(2)
// gauge algebra: V=3, H=4 (4 half-hypers in fundamental).
// This must be applied AFTER check_nhc and BEFORE compute_anomaly.
inline void enhance_external_m2_gauge(NHCResult& nhc,
                                       const Eigen::MatrixXi& IF,
                                       int ext_idx, int target_si,
                                       double cc_budget = 8.0) {
    if (IF(ext_idx, ext_idx) != -2 || target_si != -1) return;
    
    // Set gauge on the external -2
    nhc.curve_gauges[ext_idx] = GAUGE_SU2;
    
    // Update cluster H, V
    for (auto& cl : nhc.clusters) {
        for (int c : cl.curve_indices) {
            if (c == ext_idx) {
                cl.V += 3;
                cl.H += 4;
                goto done_cluster;
            }
        }
    }
    done_cluster:
    
    // Re-check central charge on neighboring -1 curves
    int n = IF.rows();
    for (int i = 0; i < n; i++) {
        if (IF(i, i) != -1 || IF(i, ext_idx) == 0) continue;
        double total_c = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i || IF(i, j) == 0 || IF(j, j) == -1) continue;
            int k = std::abs(IF(i, j));
            total_c += central_charge(nhc.curve_gauges[j], k);
        }
        if (total_c > cc_budget + 1e-9) {
            nhc.passes = false;
            nhc.fail_reason = "L2(ext): c.c. exceeded after su(2) enhancement";
            return;
        }
    }
}

// ── hat{1} gauge enhancement ──
// hat{1} (C²=-1, b₀·C=-1) carries gauge algebra depending on target:
//   → -1 target: su(8),  V=63,  H=36  (1 symmetric hyper)
//   → -2 target: su(16), V=255, H=264 (1 symmetric + 8 fundamental hypers)
// This must be applied AFTER check_nhc and BEFORE compute_anomaly.
inline void enhance_hat1_gauge(NHCResult& nhc,
                                const Eigen::MatrixXi& IF,
                                int hat1_idx, int target_si,
                                double cc_budget = 8.0) {
    if (IF(hat1_idx, hat1_idx) != -1) return;
    
    GaugeInfo gauge;
    int dV, dH;
    if (target_si == -1) {
        gauge = GAUGE_SU8; dV = 63; dH = 36;
    } else if (target_si == -2) {
        gauge = GAUGE_SU16; dV = 255; dH = 264;
    } else {
        return;
    }
    
    nhc.curve_gauges[hat1_idx] = gauge;
    
    // Update cluster H, V
    for (auto& cl : nhc.clusters) {
        for (int c : cl.curve_indices) {
            if (c == hat1_idx) {
                cl.V += dV;
                cl.H += dH;
                goto done_hat1;
            }
        }
    }
    // hat1 is a -1 curve, not in any cluster — add a new one
    {
        ClusterInfo cl;
        cl.cluster_id = (int)nhc.clusters.size();
        cl.curve_indices = {hat1_idx};
        cl.nhc = nullptr;
        cl.H = dH;
        cl.V = dV;
        nhc.clusters.push_back(cl);
        nhc.curve_cluster[hat1_idx] = cl.cluster_id;
    }
    done_hat1:
    
    // Re-check central charge on neighboring -1 curves
    int n = IF.rows();
    for (int i = 0; i < n; i++) {
        if (i == hat1_idx || IF(i, i) != -1 || IF(i, hat1_idx) == 0) continue;
        double total_c = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i || IF(i, j) == 0) continue;
            if (IF(j, j) == -1 && j != hat1_idx) continue;
            int k = std::abs(IF(i, j));
            total_c += central_charge(nhc.curve_gauges[j], k);
        }
        if (total_c > cc_budget + 1e-9) {
            nhc.passes = false;
            nhc.fail_reason = "L2(hat1): c.c. exceeded after gauge enhancement";
            return;
        }
    }
}

// ############################################################################
// PART 4: SIGNATURE + DETERMINANT
// ############################################################################

struct SigInfo {
    int sig_pos, sig_neg, sig_zero, det;
};

inline SigInfo compute_sig(const Eigen::MatrixXi& IF) {
    SigInfo s = {0, 0, 0, 0};
    int n = IF.rows();
    if (n == 0) return s;
    Eigen::MatrixXd Fd = IF.cast<double>();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Fd);
    for (int i = 0; i < n; i++) {
        double ev = es.eigenvalues()(i);
        if (std::abs(ev) < 1e-8) s.sig_zero++;
        else if (ev > 0) s.sig_pos++;
        else s.sig_neg++;
    }
    s.det = (int)std::round(Fd.determinant());
    return s;
}

// ############################################################################
// PART 5: ANOMALY
// ############################################################################

struct AnomalyResult {
    int T, H_charged, V;
    int H_neutral;      // 273 + V - H_charged - 29T >= 0
};

inline AnomalyResult compute_anomaly(const Eigen::MatrixXi& IF, const NHCResult& nhc) {
    AnomalyResult a;
    
    // T = sig_neg
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
    a.T = 0;
    for (int i = 0; i < IF.rows(); i++)
        if (es.eigenvalues()(i) < -1e-8) a.T++;
    
    a.H_charged = nhc.total_H();
    a.V = nhc.total_V();
    a.H_neutral = 273 + a.V - a.H_charged - 29 * a.T;
    return a;
}

// ############################################################################
// PART 6: ATTACHMENT RULES
// ############################################################################

struct AttachmentRule {
    int ext_self_int;
    int target_self_int;
    int intersection_num;
    std::vector<int> allowed_sidelinks;         // empty = all
    std::vector<std::string> allowed_nobles;    // empty = all
    std::vector<int> allowed_block_params;      // empty = all (for Nk)
};

// ┌───────────────────────────────────────────────────────────────┐
// │ USER: Fill allowed_sidelinks / allowed_nobles per rule        │
// │ Example:                                                      │
// │   {-2, -3, 1, {33, 331}, {"313"}, {}},                      │
// └───────────────────────────────────────────────────────────────┘

inline std::vector<AttachmentRule> build_attachment_rules() {
    std::vector<AttachmentRule> rules;
    
    rules.push_back({-2, -1, 1, {}, {}, {}});
    rules.push_back({-2, -3, 1, {}, {}, {}});
    rules.push_back({-2, -4, 2, {}, {}, {}});
    
    rules.push_back({-3, -1, 1, {}, {}, {}});
    rules.push_back({-3, -1, 2, {}, {}, {}});
    rules.push_back({-3, -1, 3, {}, {}, {}});
    rules.push_back({-3, -1, 4, {}, {}, {}});
    rules.push_back({-3, -1, 5, {}, {}, {}});
    rules.push_back({-3, -1, 6, {}, {}, {}});
    rules.push_back({-3, -1, 7, {}, {}, {}});
    rules.push_back({-3, -1, 8, {}, {}, {}});
    rules.push_back({-3, -1, 9, {}, {}, {}});









    rules.push_back({-3, -2, 1, {}, {}, {}});
    
    rules.push_back({-4, -1, 1, {}, {}, {}});
    rules.push_back({-4, -2, 2, {}, {}, {}});
    
    rules.push_back({-5,  -1, 1, {}, {}, {}});
    rules.push_back({-6,  -1, 1, {}, {}, {}});
    rules.push_back({-7,  -1, 1, {}, {}, {}});
    rules.push_back({-8,  -1, 1, {}, {}, {}});
    rules.push_back({-12, -1, 1, {}, {}, {}});
    
    return rules;
}

inline bool rule_allows_base(const AttachmentRule& rule, const CatalogEntry& entry) {
    if (entry.is_nn()) {
        auto info = nn_detail::parse_nn_topo(entry.topo);
        if (!info.valid) return false;
        if (info.is_sidelink && !rule.allowed_sidelinks.empty()) {
            bool ok = false;
            for (int p : rule.allowed_sidelinks) if (p == info.sidelink_param) { ok = true; break; }
            if (!ok) return false;
        }
        if (!info.is_sidelink && !rule.allowed_nobles.empty()) {
            bool ok = false;
            for (const auto& nm : rule.allowed_nobles) if (nm == info.noble_name) { ok = true; break; }
            if (!ok) return false;
        }
    }
    return true;
}

// ############################################################################
// PART 7: DUMMY LST + IF OPERATIONS
// ############################################################################

struct DummyLST {
    int length;
    Eigen::MatrixXi IF;
    int T;
};

inline DummyLST make_dummy_lst(int length) {
    DummyLST d;
    d.length = length;
    d.T = length - 1;
    d.IF = Eigen::MatrixXi::Zero(length, length);
    d.IF(0, 0) = -1;
    d.IF(length - 1, length - 1) = -1;
    for (int i = 1; i < length - 1; i++) d.IF(i, i) = -2;
    for (int i = 0; i < length - 1; i++) {
        d.IF(i, i + 1) = 1; d.IF(i + 1, i) = 1;
    }
    return d;
}

inline std::vector<int> find_curves(const Eigen::MatrixXi& IF, int si) {
    std::vector<int> r;
    for (int i = 0; i < IF.rows(); i++) if (IF(i, i) == si) r.push_back(i);
    return r;
}

inline Eigen::MatrixXi attach_curve(const Eigen::MatrixXi& IF,
                                     int target, int ext_si, int int_num) {
    int n = IF.rows();
    Eigen::MatrixXi nIF = Eigen::MatrixXi::Zero(n + 1, n + 1);
    nIF.block(0, 0, n, n) = IF;
    nIF(n, n) = ext_si;
    nIF(n, target) = int_num;
    nIF(target, n) = int_num;
    return nIF;
}

inline Eigen::MatrixXi glue_dummy(const Eigen::MatrixXi& base,
                                   int target, const DummyLST& dummy, int ep = 0) {
    int nb = base.rows(), nd = dummy.length;
    Eigen::MatrixXi nIF = Eigen::MatrixXi::Zero(nb + nd - 1, nb + nd - 1);
    nIF.block(0, 0, nb, nb) = base;
    auto map_d = [&](int j) -> int {
        if (j == ep) return target;
        return nb + (j < ep ? j : j - 1);
    };
    for (int j = 0; j < nd; j++) {
        if (j == ep) continue;
        nIF(map_d(j), map_d(j)) = dummy.IF(j, j);
    }
    for (int j = 0; j < nd; j++)
        for (int k = j + 1; k < nd; k++)
            if (dummy.IF(j, k) != 0) {
                nIF(map_d(j), map_d(k)) += dummy.IF(j, k);
                nIF(map_d(k), map_d(j)) += dummy.IF(k, j);
            }
    return nIF;
}

// ############################################################################
// PART 8: SUGRA RESULT + CONFIG
// ############################################################################

struct SUGRAResult {
    int catalog_id;
    std::string catalog_type;
    Eigen::MatrixXi base_IF;
    Eigen::MatrixXi final_IF;
    
    struct Attachment { int ext_si, target_idx, int_num; };
    std::vector<Attachment> externals;
    
    struct DummyGlue { int length, target_idx, endpoint; };
    std::vector<DummyGlue> dummies;
    
    NHCResult nhc;
    AnomalyResult anomaly;
    SigInfo sig;
    
    bool valid;
    std::string fail_reason;
};

struct SUGRAConfig {
    int T_max = 193;
    int max_externals = 10;
    int max_dummy = 0;
    int dummy_min_length = 2;
    int dummy_max_length = 10;
    
    double cc_budget = 8.0;
    
    bool check_nhc = true;
    bool check_anomaly = false;
    bool check_determinant = false;
    int required_det = 0;
    
    bool include_nn = true;
    bool include_nk = true;
    bool include_dummy = true;
    int catalog_T_min = 0;
    int catalog_T_max = 9999;
    
    bool verbose = false;
};

// ############################################################################
// PART 9: CANDIDATE FINDING
// ############################################################################

inline double used_cc_on_minus1(const Eigen::MatrixXi& IF, int m1) {
    double used = 0;
    for (int j = 0; j < IF.rows(); j++) {
        if (j == m1 || IF(m1, j) == 0 || IF(j, j) == -1) continue;
        used += central_charge(gauge_from_si(IF(j, j)), std::abs(IF(m1, j)));
    }
    return used;
}

struct AttachCandidate {
    int rule_idx, target_idx, ext_si, int_num;
};

inline std::vector<AttachCandidate> find_candidates(
    const Eigen::MatrixXi& IF,
    const std::vector<AttachmentRule>& rules,
    double budget = 8.0)
{
    std::vector<AttachCandidate> cands;
    for (int ri = 0; ri < (int)rules.size(); ri++) {
        const auto& rule = rules[ri];
        for (int tidx : find_curves(IF, rule.target_self_int)) {
            if (rule.target_self_int == -1) {
                auto g = gauge_from_si(rule.ext_self_int);
                if (g.dim > 0) {
                    double rem = budget - used_cc_on_minus1(IF, tidx);
                    if (central_charge(g, rule.intersection_num) > rem + 1e-9) continue;
                }
            }
            cands.push_back({ri, tidx, rule.ext_self_int, rule.intersection_num});
        }
    }
    return cands;
}

// ############################################################################
// PART 10: GENERATOR
// ############################################################################

// Dummy LST IF builders
// A(n): 1-2-2-...-2-1, n curves total (n >= 2)
inline Eigen::MatrixXi build_dummy_A(int n) {
    Eigen::MatrixXi IF = Eigen::MatrixXi::Zero(n, n);
    IF(0, 0) = -1;
    IF(n-1, n-1) = -1;
    for (int i = 1; i < n-1; i++) IF(i, i) = -2;
    for (int i = 0; i < n-1; i++) { IF(i, i+1) = 1; IF(i+1, i) = 1; }
    return IF;
}

// D(n): D-type dummy LST
// n=3: linear chain [-2, -1, -2] (degenerate, no branch)
// n>=4: 2{2,2}2...21
//   Main chain: [0]=-2, [1]=-2, ..., [n-4]=-2, [n-3]=-1  (n-2 curves)
//   Two branches: [n-2]=-2, [n-1]=-2, both connected to [0]
//   Node 0 has degree: 2 (branch) + 1 (chain) = 3  (for n=4: chain is just -1)

inline Eigen::MatrixXi build_dummy_D(int n) {
    if (n == 3) {
        // Linear chain: -2, -1, -2
        Eigen::MatrixXi IF = Eigen::MatrixXi::Zero(3, 3);
        IF(0, 0) = -2; IF(1, 1) = -1; IF(2, 2) = -2;
        IF(0, 1) = 1; IF(1, 0) = 1;
        IF(1, 2) = 1; IF(2, 1) = 1;
        return IF;
    }
    // n >= 4: total n curves
    // Chain: indices 0..n-3, lengths n-2
    //   [0]=-2 (branch point), [1]=-2, ..., [n-4]=-2, [n-3]=-1
    // Two branches: [n-2]=-2, [n-1]=-2, connected to [0]
    Eigen::MatrixXi IF = Eigen::MatrixXi::Zero(n, n);
    for (int i = 0; i <= n - 4; i++) IF(i, i) = -2;
    IF(n - 3, n - 3) = -1;
    IF(n - 2, n - 2) = -2;
    IF(n - 1, n - 1) = -2;
    // Chain adjacency
    for (int i = 0; i < n - 3; i++) { IF(i, i + 1) = 1; IF(i + 1, i) = 1; }
    // Two branches at node 0
    IF(0, n - 2) = 1; IF(n - 2, 0) = 1;
    IF(0, n - 1) = 1; IF(n - 1, 0) = 1;
    return IF;
}

// Generate from raw IF (for dummy LSTs)
inline std::vector<SUGRAResult> generate_from_IF(
    const Eigen::MatrixXi& base_IF,
    const std::string& label,
    int label_id,
    const std::vector<AttachmentRule>& rules,
    const SUGRAConfig& config)
{
    std::vector<SUGRAResult> results;
    
    auto base_sig = compute_sig(base_IF);
    if (base_sig.sig_neg > config.T_max) return results;
    
    auto cands = find_candidates(base_IF, rules, config.cc_budget);
    
    for (const auto& c : cands) {
        // No rule_allows_base check for dummy LSTs (all rules allowed)
        
        Eigen::MatrixXi new_IF = attach_curve(base_IF, c.target_idx, c.ext_si, c.int_num);
        auto new_sig = compute_sig(new_IF);
        
        if (new_sig.sig_neg > config.T_max) continue;
        if (config.check_determinant && std::abs(new_sig.det) != 1) continue;
        
        NHCResult nhc = {true, "", {}, {}, {}};
        if (config.check_nhc) {
            nhc = check_nhc(new_IF, config.cc_budget);
            if (!nhc.passes) continue;
        }
        
        // External -2 → -1: su(2) gauge enhancement
        enhance_external_m2_gauge(nhc, new_IF, new_IF.rows() - 1,
                                  base_IF(c.target_idx, c.target_idx), config.cc_budget);
        if (!nhc.passes) continue;
        
        AnomalyResult anom = {0, 0, 0, 0};
        if (nhc.passes) {
            anom = compute_anomaly(new_IF, nhc);
            if (config.check_anomaly && anom.H_neutral < 0) continue;
        }
        
        SUGRAResult r;
        r.catalog_id = label_id;
        r.catalog_type = label;
        r.base_IF = base_IF;
        r.final_IF = new_IF;
        r.externals.push_back({c.ext_si, c.target_idx, c.int_num});
        r.nhc = nhc;
        r.anomaly = anom;
        r.sig = new_sig;
        r.valid = true;
        results.push_back(r);
    }
    
    return results;
}

inline std::vector<SUGRAResult> generate_from_entry(
    const CatalogEntry& entry,
    const std::vector<AttachmentRule>& rules,
    const SUGRAConfig& config)
{
    std::vector<SUGRAResult> results;
    
    Eigen::MatrixXi base_IF = reconstruct_IF(entry);
    if (base_IF.rows() == 0) return results;
    
    auto base_sig = compute_sig(base_IF);
    if (base_sig.sig_neg > config.T_max) return results;
    
    auto cands = find_candidates(base_IF, rules, config.cc_budget);
    
    for (const auto& c : cands) {
        if (!rule_allows_base(rules[c.rule_idx], entry)) continue;
        
        Eigen::MatrixXi new_IF = attach_curve(base_IF, c.target_idx, c.ext_si, c.int_num);
        auto new_sig = compute_sig(new_IF);
        
        if (new_sig.sig_neg > config.T_max) continue;
        if (config.check_determinant && std::abs(new_sig.det) != 1) continue;
        
        NHCResult nhc = {true, "", {}, {}, {}};
        if (config.check_nhc) {
            nhc = check_nhc(new_IF, config.cc_budget);
            if (!nhc.passes) continue;
        }
        
        // External -2 → -1: su(2) gauge enhancement
        enhance_external_m2_gauge(nhc, new_IF, new_IF.rows() - 1,
                                  base_IF(c.target_idx, c.target_idx), config.cc_budget);
        if (!nhc.passes) continue;
        
        AnomalyResult anom = {0, 0, 0, 0};
        if (nhc.passes) {
            anom = compute_anomaly(new_IF, nhc);
            if (config.check_anomaly && anom.H_neutral < 0) continue;
        }
        
        SUGRAResult r;
        r.catalog_id = entry.id;
        r.catalog_type = entry.type;
        r.base_IF = base_IF;
        r.final_IF = new_IF;
        r.externals.push_back({c.ext_si, c.target_idx, c.int_num});
        r.nhc = nhc;
        r.anomaly = anom;
        r.sig = new_sig;
        r.valid = true;
        results.push_back(r);
    }
    
    return results;
}

inline std::vector<SUGRAResult> generate_sugra(
    const std::vector<CatalogEntry>& catalog,
    const SUGRAConfig& config)
{
    auto rules = build_attachment_rules();
    std::vector<SUGRAResult> all;
    
    // ── Catalog entries ──
    int processed = 0;
    for (const auto& entry : catalog) {
        if (entry.is_nn() && !config.include_nn) continue;
        if (!entry.is_nn() && !config.include_nk) continue;
        if (entry.T < config.catalog_T_min || entry.T > config.catalog_T_max) continue;
        
        auto res = generate_from_entry(entry, rules, config);
        for (auto& r : res) all.push_back(std::move(r));
        
        processed++;
        if (config.verbose && processed % 1000 == 0)
            std::cout << processed << " → " << all.size() << " bases\n";
    }
    if (config.verbose)
        std::cout << "Catalog: " << processed << " → " << all.size() << " bases\n";
    
    // ── Dummy LSTs ──
    if (config.include_dummy) {
        int dummy_count = 0;
        
        // A-type: 1-2-...-2-1, n curves → T_base = n-1
        // Cap: catalog_T_min ≤ T_base ≤ catalog_T_max
        int a_min = std::max(2, config.catalog_T_min + 1);  // n = T_base + 1
        int a_max = config.catalog_T_max + 1;
        for (int n = a_min; n <= a_max; n++) {
            auto IF = build_dummy_A(n);
            auto res = generate_from_IF(IF, "DM:A(" + std::to_string(n) + ")", -n, rules, config);
            for (auto& r : res) all.push_back(std::move(r));
            dummy_count += res.size();
        }
        
        // D-type: D(3)=212, D(4)=2{2,2}1, D(5)=2{2,2}21, ...
        int d_min = std::max(3, config.catalog_T_min + 1);
        int d_max = config.catalog_T_max + 1;
        for (int n = d_min; n <= d_max; n++) {
            auto IF = build_dummy_D(n);
            auto res = generate_from_IF(IF, "DM:D(" + std::to_string(n) + ")", -(1000+n), rules, config);
            for (auto& r : res) all.push_back(std::move(r));
            dummy_count += res.size();
        }
        
        if (config.verbose)
            std::cout << "Dummy LSTs: " << dummy_count << " bases\n";
    }
    
    if (config.verbose)
        std::cout << "Total: " << all.size() << " bases\n";
    return all;
}

// ############################################################################
// PART 11: OUTPUT
// ############################################################################

inline void print_sugra_summary(const std::vector<SUGRAResult>& results) {
    std::map<int, int> T_dist, ext_dist;
    int anom_pass = 0;
    for (const auto& r : results) {
        T_dist[r.anomaly.T]++;
        for (const auto& e : r.externals) ext_dist[e.ext_si]++;
        if (r.anomaly.H_neutral >= 0) anom_pass++;
    }
    std::cout << "SUGRA: " << results.size() << " bases | H_neutral>=0: " << anom_pass << "\n";
    std::cout << "T distribution:\n";
    for (auto& [T, cnt] : T_dist) std::cout << "  T=" << T << ": " << cnt << "\n";
    std::cout << "Externals:\n";
    for (auto& [si, cnt] : ext_dist) std::cout << "  " << si << ": " << cnt << "\n";
}

inline void print_gauge_assignment(const SUGRAResult& r) {
    std::cout << "ID=" << r.catalog_id << " " << r.catalog_type
              << " T=" << r.anomaly.T << " Hc=" << r.anomaly.H_charged
              << " V=" << r.anomaly.V << " Hn=" << r.anomaly.H_neutral << "\n";
    std::cout << "  IF " << r.final_IF.rows() << "x" << r.final_IF.cols()
              << " det=" << r.sig.det << "\n";
    for (int ci = 0; ci < (int)r.nhc.clusters.size(); ci++) {
        auto& cl = r.nhc.clusters[ci];
        std::cout << "  cluster " << ci << ": ";
        if (cl.nhc) std::cout << cl.nhc->name;
        else std::cout << "pure-2";
        std::cout << " H=" << cl.H << " V=" << cl.V << " curves={";
        for (int i = 0; i < (int)cl.curve_indices.size(); i++) {
            if (i) std::cout << ",";
            int c = cl.curve_indices[i];
            std::cout << c << "(" << r.final_IF(c, c) << ")";
        }
        std::cout << "}\n";
    }
}

// end of sugra_generator.h
