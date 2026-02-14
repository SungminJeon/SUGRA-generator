// no_node_theory.cpp
// Generates all no-node LSTs (side links + noble molecules with extra curve)
//
// Original output: LaTeX file with quiver-style notation
// NEW: --catalog mode outputs structured catalog for unified_catalog integration
//
// Algorithm: Intersection matrix → Quiver LaTeX directly
// - Follow a_{i,i+1} for linear chain
// - If a_{i,i+1}=0, find j where a_{i,j}!=0, stack curves (i+1 ~ j-1) on j
// - Extra curve (last row) attached with red overset
// 
// Deduplication: Compare eigenvalue spectra of intersection forms
//
// NEW features:
//   --catalog <file>     Output catalog-format lines (pipe-delimited)
//   --no-nhc-filter      Include NHC-failing entries (tagged NHC_FAIL)
//   --nhc-only           Only output NHC-passing entries (default)
//   --verbose / -v        Print progress

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>
#include <set>
#include <map>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "Tensor.h"
#include "Theory_enhanced.h"
#include "anomaly_tables.h"

// ============================================================================
// All Side Link Params (from Theory_enhanced.h)
// ============================================================================

const std::vector<int> ALL_SIDELINK_PARAMS = {
    // instantons
//    1, 882, 883, 884, 885, 886, 887, 8881, 889, 8810, 8811,
       
    // interiors
    22, 33, 44, 55, 331,
    23, 24, 34, 35, 45,
    
    // alkali 2 links (no -5)
    991, 9920, 993,
    
    // alkali 1 links (no -5)
    91, 92, 93, 94, 95, 96, 97, 98, 99, 910, 911, 912, 913, 914, 915, 916, 917,
    
    // alkali 3 links (one -5)
    99910, 99920, 99930,
    
    // alkali 2 links (one -5)
    994, 995, 996, 997, 998, 999, 9910, 9911, 9912, 9913, 9914,
    
    // alkali 1 links (one -5)
    918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930,
    931, 932, 933, 934, 935, 936, 937, 938, 939, 940, 941, 942, 943, 944, 945,
    958,
    
    // alkali 2 links (two -5)
    9915, 9916, 9917,
    
    // alkali 1 links (two -5)
    946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 956, 957
};

// Instanton sidelinks (excluded from LST count with --no-instanton)
const std::set<int> INSTANTON_PARAMS = {
    1, 882, 883, 884, 885, 886, 887, 8881, 889, 8810, 8811
};

// ============================================================================
// Gluing Rules
// ============================================================================

// Allowed intersection numbers for an extra curve
// For LST/SCFT: always 1. Higher intersection numbers only in SUGRA.
std::vector<int> get_allowed_intersections(int self_int) {
    (void)self_int;
    return {1};
}

const std::vector<int> EXTRA_CURVES = {-1, -2, -3, -5};

std::vector<int> get_allowed_extra_curves(int base_self_int) {
    switch (base_self_int) {
        case -1: return {-1, -2, -3, -5};
        case -2: return {-1, -3};
        case -3: return {-1, -2};
        case -5: return {-1};
        default: return {-1};
    }
}

// ============================================================================
// LST Check
// ============================================================================

bool is_LST(const Tensor& t) {
    return t.NullDirection() == 1;
}

// ============================================================================
// P-type Endpoint Check (UNCHANGED from original)
// ============================================================================

bool has_P_type_endpoints(Tensor t) {
    if (!is_LST(t)) return false;
    
    t.Setb0Q();
    t.ForcedBlowdown();
    
    Eigen::MatrixXi IF_final = t.GetIntersectionForm();
    std::vector<int> b0_final = t.Getb0Q();
    
    if (IF_final.rows() == 0) return false;
    
    const int n = IF_final.rows();
    
    if (n == 1) {
        return (IF_final(0, 0) == 0) && (b0_final.size() > 0 && b0_final[0] == 2);
    }
    
    std::vector<int> endpoints;
    for (int i = 0; i < n; ++i) {
        int degree = 0;
        for (int j = 0; j < n; ++j) {
            if (i != j && IF_final(i, j) != 0) degree++;
        }
        if (degree == 1) {
            endpoints.push_back(i);
        }
    }
    
    if (endpoints.empty()) return false;
    
    for (int ep : endpoints) {
        if (ep >= n || ep >= (int)b0_final.size()) return false;
        if (IF_final(ep, ep) != 0) return false;
        if (b0_final[ep] != 2) return false;
    }
    
    return true;
}

bool is_valid_LST(const Tensor& t) {
    if (!is_LST(t)) return false;
    
    Tensor t_copy;
    t_copy.SetIF(t.GetIntersectionForm());
    
    return has_P_type_endpoints(t_copy);
}

// ============================================================================
// Eigenvalue Spectrum Computation and Comparison
// (여기서 먼저 정의 — NHC check에서 사용)
// ============================================================================

std::vector<double> compute_eigenvalue_spectrum(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    if (n == 0) return {};
    
    Eigen::MatrixXd IF_double = IF.cast<double>();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(IF_double);
    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    
    std::vector<double> spectrum(n);
    for (int i = 0; i < n; ++i) {
        spectrum[i] = eigenvalues(i);
    }
    std::sort(spectrum.begin(), spectrum.end());
    
    return spectrum;
}

bool are_spectra_equal(const std::vector<double>& s1, const std::vector<double>& s2, double tol = 1e-9) {
    if (s1.size() != s2.size()) return false;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (std::abs(s1[i] - s2[i]) > tol) return false;
    }
    return true;
}

// ============================================================================
// Signature Computation
// ============================================================================

struct SignatureInfo {
    int sig_pos, sig_neg, sig_zero;
    int det;
};

SignatureInfo compute_signature(const Eigen::MatrixXi& IF) {
    SignatureInfo info = {0, 0, 0, 0};
    int n = IF.rows();
    if (n == 0) return info;
    
    Eigen::MatrixXd IF_d = IF.cast<double>();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(IF_d);
    auto evals = solver.eigenvalues();
    
    const double tol = 1e-8;
    for (int i = 0; i < evals.size(); ++i) {
        if (std::abs(evals(i)) < tol) info.sig_zero++;
        else if (evals(i) > tol) info.sig_pos++;
        else info.sig_neg++;
    }
    
    info.det = static_cast<int>(std::round(IF_d.determinant()));
    return info;
}


// ============================================================================
// NHC Check: Three-Layer Validation
// ============================================================================
//
// 용어:
//   "extra curve" = LST를 만들기 위해 붙이는 추가 curve (빨간색)
//   "external curve" = SUGRA를 만들 때 나중에 붙이는 curve (다른 단계)
//
// Three layers:
//   Layer 1: Cluster recognition (eigenvalue spectrum comparison)
//   Layer 2: E8 central charge bound (sum on each -1 ≤ 8)
//   Layer 3: -1 flavor gauge compatibility
//
// ┌──────────────────────────────────────────────────────────────────┐
// │ 예시: base = s(23) = [-1,-2,-3,-1] 에 extra -2를 port 4에 부착 │
// │                                                                  │
// │ Full IF:                                                         │
// │   [-1, 1, 0, 0, 0]   curve 0: -1                               │
// │   [ 1,-2, 1, 0, 0]   curve 1: -2  ┐                            │
// │   [ 0, 1,-3, 1, 0]   curve 2: -3  ├─ cluster A: [-2,-3]       │
// │   [ 0, 0, 1,-1, 1]   curve 3: -1                               │
// │   [ 0, 0, 0, 1,-2]   curve 4: -2  ── cluster B: [-2]          │
// │                                                                  │
// │ Layer 1: cluster A = [-2,-3] → spectrum match → recognized NHC  │
// │          cluster B = [-2] alone → pure -2, always OK            │
// │                                                                  │
// │ Layer 2: curve 0 (-1): neighbor = curve 1 (-2, cluster A)      │
// │            port gauge of -2 in [-2,-3] = su2, c = 1.0          │
// │            total c = 1.0 ≤ 8 ✓                                  │
// │          curve 3 (-1): neighbors = curve 2 (-3, cluster A)     │
// │                                  + curve 4 (-2, cluster B)     │
// │            port gauge of -3 in [-2,-3] = su3, c = 2.0          │
// │            port gauge of -2 alone = none, c = 0                 │
// │            total c = 2.0 ≤ 8 ✓                                  │
// │                                                                  │
// │ Layer 3: curve 0: {su2} → compatible ✓                          │
// │          curve 3: {su3, none} → compatible ✓                    │
// └──────────────────────────────────────────────────────────────────┘

// ── Gauge info for port gauge lookup ──

struct PortGaugeInfo {
    const char* gauge_name;
    double central_charge;
};

const PortGaugeInfo GAUGE_NONE  = {"none",  0.0};
const PortGaugeInfo GAUGE_SP1   = {"sp1",   2.0};    // h∨(sp1) = 2
const PortGaugeInfo GAUGE_SU2   = {"su2",   2.0};    // h∨(su2) = 2
const PortGaugeInfo GAUGE_SU3   = {"su3",   3.0};    // h∨(su3) = 3
const PortGaugeInfo GAUGE_G2    = {"g2",    4.0};    // h∨(g2) = 4
const PortGaugeInfo GAUGE_SO7   = {"so7",   5.0};    // h∨(so7) = 5
const PortGaugeInfo GAUGE_SO8   = {"so8",   6.0};    // h∨(so8) = 6
const PortGaugeInfo GAUGE_SO9   = {"so9",   7.0};    // h∨(so9) = 7
const PortGaugeInfo GAUGE_SO10  = {"so10",  8.0};    // h∨(so10) = 8
const PortGaugeInfo GAUGE_SO11  = {"so11",  9.0};    // h∨(so11) = 9
const PortGaugeInfo GAUGE_SO12  = {"so12", 10.0};    // h∨(so12) = 10
const PortGaugeInfo GAUGE_F4    = {"f4",    9.0};    // h∨(f4) = 9
const PortGaugeInfo GAUGE_E6    = {"e6",   12.0};    // h∨(e6) = 12
const PortGaugeInfo GAUGE_E7    = {"e7",   18.0};    // h∨(e7) = 18
const PortGaugeInfo GAUGE_E8    = {"e8",   30.0};    // h∨(e8) = 30

// NOTE: central charge here = (1/6) * dim(adj) 기준이 아님.
// F-theory에서 사용하는 convention에 맞춰야 함.
// 
// USER: 아래 값들을 논문의 convention에 맞게 확인/수정할 것.
// 현재 값은 c_g = (rank of gauge) 또는 dual Coxeter number 기반 추정.
// 정확한 공식: "A에 의한 E8 central charge 소모" = ?
//
// E8 budget 계산에서 사용하는 central charge:
//   E8 total = 8 (혹은 다른 convention?)
//   su2 → 1, su3 → 2, g2 → 2, so7 → 3.5, so8 → 4, f4 → 4.5, ...

// ── NHC Definition with Port Gauges ──

struct NHCDefinition {
    std::string name;
    bool is_linear;
    
    // For linear NHC: ordered self-intersections of the chain
    // 예: {-2, -3, -2} means -2 — -3 — -2
    std::vector<int> self_ints;
    
    // Gauge per curve (same ordering as self_ints)
    // gauge[i] = gauge algebra on curve i within this NHC
    std::vector<PortGaugeInfo> gauges;
    
    // Eigenvalue spectrum (auto-computed at init)
    std::vector<double> spectrum;
    
    // For branching NHC: full IF template
    Eigen::MatrixXi IF_template;
};

// ── Build the NHC table ──
//
// 각 NHC의 self-int chain과 curve별 gauge를 hardcode.
// spectrum은 자동 계산.
//
// USER: 여기에 허용할 NHC와 port gauge를 추가/수정.
//
// ┌────────────────────────────────────────────────────────────┐
// │ 추가하는 법:                                                │
// │                                                            │
// │ 선형:                                                      │
// │ add_linear("name",                                         │
// │            {-2, -3, -2},                // self_ints       │
// │            {GAUGE_SU2, GAUGE_SO7, GAUGE_SU2});  // gauges │
// └────────────────────────────────────────────────────────────┘

std::vector<NHCDefinition> build_nhc_table() {
    std::vector<NHCDefinition> table;
    
    // Helper: linear NHC (NHC는 linear chain만 존재)
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
    
    // ══════════════════════════════════════════════════════════════
    // Standard NHCs (Morrison-Taylor classification)
    //
    // Format: (name, {self_ints}, {gauge per curve})
    //
    // 각 curve의 gauge는 해당 curve가 -1과 만날 때 
    // -1이 "보는" gauge algebra.
    // 내부 curve는 -1과 직접 안 만나므로 실질적 의미 없지만,
    // 완전성을 위해 기록.
    // ══════════════════════════════════════════════════════════════
    
    // ── Single curves ──
    add_linear("-3",   {-3},   {GAUGE_SU3});
    add_linear("-4",   {-4},   {GAUGE_SO8});
    add_linear("-5",   {-5},   {GAUGE_F4});
    add_linear("-6",   {-6},   {GAUGE_E6});
    add_linear("-7",   {-7},   {GAUGE_E7});
    add_linear("-8",   {-8},   {GAUGE_E7});
    add_linear("-12",  {-12},  {GAUGE_E8});
    
    // ── Two-curve chains ──
    //
    //   -2 — -3  : -2 carries su2, -3 carries g2
    //   (주의: standalone -3은 su3이지만, -2와 붙으면 -3의 gauge가 g2로 enhance)
    //
    // NHC는 한 방향만 등록. lookup_port_gauge()에서 forward/reverse walk으로
    // 어느 방향이든 자동으로 gauge 위치를 찾아냄.
    
    // ── Single curves (mandatory gauge) ──
    add_linear("-3",   {-3},   {GAUGE_SU3});
    add_linear("-4",   {-4},   {GAUGE_SO8});
    add_linear("-5",   {-5},   {GAUGE_F4});
    add_linear("-6",   {-6},   {GAUGE_E6});
    add_linear("-7",   {-7},   {GAUGE_E7});
    add_linear("-8",   {-8},   {GAUGE_E7});
    add_linear("-12",  {-12},  {GAUGE_E8});
    
    // ── Multi-curve NHCs (이게 전부) ──
    add_linear("-2-3",   {-2, -3},       {GAUGE_SU2, GAUGE_G2});
    add_linear("-2-3-2", {-2, -3, -2},   {GAUGE_SU2, GAUGE_SO7, GAUGE_SU2});
    add_linear("-2-2-3", {-2, -2, -3},   {GAUGE_NONE, GAUGE_SP1, GAUGE_G2});
    
    // USER: 더 필요한 linear chain 추가
    // add_linear("-2-3-2-5", {-2,-3,-2,-5}, {?, ?, ?, ?});
    
    // ══════════════════════════════════════════════════════════════
    // NHC는 linear chain만 존재. Branching topology는 NHC가 아님.
    // NHC table에 없는 cluster → L1 fail → reject.
    //
    // USER: 여기에 추가 linear NHC를 넣을 것.
    //   add_linear("name", {self_ints...}, {gauges...});
    //
    // 프로그램 실행 시 spectrum이 자동 계산됨.
    // fail 메시지에 cluster의 self-int과 spectrum이 출력되므로,
    // 그걸 보고 빠진 NHC를 추가하면 됨.
    // ══════════════════════════════════════════════════════════════
    
    return table;
}

// Global NHC table (built once)
const std::vector<NHCDefinition>& get_nhc_table() {
    static std::vector<NHCDefinition> cache = build_nhc_table();
    return cache;
}

// ── NHC Identification ──
//
// Given a cluster's eigenvalue spectrum, find the matching NHC definition.
// Returns nullptr if no match (and the cluster is not a pure-(-2) ADE cluster).

bool is_pure_minus2_cluster(const std::vector<int>& diag) {
    for (int d : diag) if (d != -2) return false;
    return true;
}

const NHCDefinition* identify_nhc(const std::vector<double>& spectrum) {
    for (const auto& def : get_nhc_table()) {
        if (are_spectra_equal(spectrum, def.spectrum)) return &def;
    }
    return nullptr;
}

// ── Cluster → Port Gauge Lookup ──
//
// 주어진 cluster 내에서 global_idx 위치의 curve가 무슨 gauge를 갖는지 반환.
//
// Linear NHC의 경우:
//   cluster의 curve들을 chain walk으로 정렬 → NHC definition의 self_ints와 비교
//   → 매칭되면 gauges[위치] 반환
//   → reverse도 시도 (chain 방향 모호할 수 있으므로)
//
// Branching NHC의 경우:
//   IF submatrix를 NHC definition의 IF_template과 비교 (permutation 탐색)

// Helper: walk a linear chain from one endpoint, return ordered curve indices
std::vector<int> walk_linear_chain(const Eigen::MatrixXi& IF,
                                   const std::vector<int>& cluster_curves) {
    int m = cluster_curves.size();
    if (m <= 1) return cluster_curves;
    
    // adjacency within cluster
    std::map<int, std::vector<int>> adj;
    for (int a : cluster_curves) {
        for (int b : cluster_curves) {
            if (a != b && IF(a, b) != 0) adj[a].push_back(b);
        }
    }
    
    // find endpoint (degree 1 within cluster)
    int start = cluster_curves[0];
    for (int c : cluster_curves) {
        if ((int)adj[c].size() <= 1) { start = c; break; }
    }
    
    // walk
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

// Look up the port gauge of a specific curve in a cluster
PortGaugeInfo lookup_port_gauge(const Eigen::MatrixXi& IF,
                                const std::vector<int>& cluster_curves,
                                int target_global_idx,
                                const NHCDefinition* nhc_def) {
    if (!nhc_def) return GAUGE_NONE;  // pure -2 cluster → no gauge
    
    if (nhc_def->is_linear) {
        // Walk the chain
        auto ordered = walk_linear_chain(IF, cluster_curves);
        
        // Read off self-ints
        std::vector<int> chain_si;
        for (int c : ordered) chain_si.push_back(IF(c, c));
        
        // Try forward match
        if (chain_si == nhc_def->self_ints) {
            for (size_t i = 0; i < ordered.size(); ++i) {
                if (ordered[i] == target_global_idx)
                    return nhc_def->gauges[i];
            }
        }
        
        // Try reverse match
        std::vector<int> rev_si(chain_si.rbegin(), chain_si.rend());
        if (rev_si == nhc_def->self_ints) {
            int m = ordered.size();
            for (size_t i = 0; i < ordered.size(); ++i) {
                if (ordered[i] == target_global_idx)
                    return nhc_def->gauges[m - 1 - i];
            }
        }
        
        // self-int sequence mismatch (shouldn't happen if spectrum matched)
        // fallback: return gauge based on self-int alone
    }
    
    // Branching or fallback: try to match by self-int position
    // (crude but usable until full permutation matching is implemented)
    int si = IF(target_global_idx, target_global_idx);
    for (size_t i = 0; i < nhc_def->self_ints.size(); ++i) {
        if (nhc_def->self_ints[i] == si) return nhc_def->gauges[i];
    }
    
    return GAUGE_NONE;
}

// ── NHC Check Result ──

struct NHCCheckResult {
    bool passes;
    std::string failure_reason;
    
    // Per-layer results
    bool layer1_pass = true;  // cluster recognition
    bool layer2_pass = true;  // E8 central charge bound
    bool layer3_pass = true;  // flavor compatibility
    
    // Debug
    int n_clusters = 0;
    std::vector<std::vector<int>> clusters;
};

// ── Layer 1: Cluster Recognition ──

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
        
        // pure -2 cluster → always OK, no gauge constraint
        std::vector<int> diag;
        for (int c : cluster) diag.push_back(IF(c, c));
        
        if (is_pure_minus2_cluster(diag)) {
            for (int c : cluster) curve_to_nhc[c] = nullptr;
            continue;
        }
        
        // NHC lookup
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

// ── E8 Subalgebra Pair Table ──
//
// -1 curve 양쪽의 gauge가 G₁ × G₂ ⊂ E₈ 인지 체크.
// 
// E₈의 maximal subalgebra 분해:
//   A₁ × E₇,  A₂ × E₆,  G₂ × F₄,  D₅ × A₃,  A₄ × A₄
//
// 각 factor의 sub-algebra도 허용 (G₁' ⊂ G₁ → G₁' × G₂ ⊂ E₈).
//
// 테이블: (g1, g2) 쌍이 허용되는지 (순서 무관)

bool is_e8_subalgebra_pair(const std::string& g1, const std::string& g2) {
    if (g1 == "none" || g2 == "none") return true;
    
    static const std::set<std::pair<std::string,std::string>> allowed = {
        // sp1 pairs (⊂ A₁, partner ⊂ E₇)
        {"sp1", "sp1"}, {"sp1", "su2"}, {"sp1", "su3"},
        {"sp1", "g2"}, {"sp1", "so7"}, {"sp1", "so8"},
        {"sp1", "so9"}, {"sp1", "so10"}, {"sp1", "so11"},
        {"sp1", "so12"}, {"sp1", "f4"}, {"sp1", "e6"}, {"sp1", "e7"},
        
        // su2 pairs (= A₁, partner ⊂ E₇)
        {"su2", "su2"}, {"su2", "su3"}, {"su2", "g2"},
        {"su2", "so7"}, {"su2", "so8"}, {"su2", "so9"},
        {"su2", "so10"}, {"su2", "so11"}, {"su2", "so12"},
        {"su2", "f4"}, {"su2", "e6"}, {"su2", "e7"},
        
        // su3 pairs (⊂ A₂, partner ⊂ E₆)
        {"su3", "su3"}, {"su3", "g2"}, {"su3", "so7"},
        {"su3", "so8"}, {"su3", "so9"}, {"su3", "so10"},
        {"su3", "f4"}, {"su3", "e6"},
        
        // g2 pairs (G₂, partner ⊂ F₄)
        {"g2", "g2"}, {"g2", "so7"}, {"g2", "so8"},
        {"g2", "so9"}, {"g2", "f4"},
        
        // so7 pairs (⊂ F₄, commutant = A₁ × G₂)
        {"so7", "su3"},
        
        // so8 pairs (D₄ ⊂ D₅, commutant ⊂ A₁ × A₃)
        {"so8", "su3"},
        
        // so10 pairs (= D₅, partner ⊂ A₃)
        {"so10", "su3"},
        
        // f4 pairs (F₄, partner ⊂ G₂)
        {"f4", "su3"}, {"f4", "g2"},
        
        // e6 pairs (E₆, partner ⊂ A₂)
        {"e6", "su3"},
    };
    
    // 양방향 체크
    return allowed.count({g1, g2}) > 0 || allowed.count({g2, g1}) > 0;
}

// ── Layer 2: E8 Gauge Constraint (subalgebra + h∨ budget) ──
//
// 각 -1 curve에 대해:
//   (a) 이웃 non-(-1) gauge들의 직접곱이 E₈의 subalgebra인지
//   (b) sum(h∨) ≤ bound (20 for T=1, 16 for T>1)
//
// 2개 이웃: pair check
// 3+ 이웃: 모든 쌍이 호환 + h∨ 합 체크 (더 보수적)

bool check_layer2_e8_constraint(const Eigen::MatrixXi& IF,
                                const std::vector<int>& minus1_curves,
                                const std::vector<std::vector<int>>& clusters,
                                const std::map<int, const NHCDefinition*>& curve_to_nhc,
                                int T_tensor,
                                std::string& fail_reason) {
    int n = IF.rows();
    double hv_bound = (T_tensor == 1) ? 20.0 : 16.0;
    
    // Build: curve → which cluster it belongs to
    std::map<int, int> curve_to_cluster_idx;
    for (size_t ci = 0; ci < clusters.size(); ++ci)
        for (int c : clusters[ci])
            curve_to_cluster_idx[c] = ci;
    
    for (int m1 : minus1_curves) {
        // Collect all non-(-1) neighbors and their gauges
        std::vector<PortGaugeInfo> neighbor_gauges;
        double total_hv = 0.0;
        
        for (int j = 0; j < n; ++j) {
            if (j == m1 || IF(m1, j) == 0) continue;
            if (IF(j, j) == -1) continue;
            
            auto it = curve_to_cluster_idx.find(j);
            if (it == curve_to_cluster_idx.end()) {
                neighbor_gauges.push_back(GAUGE_NONE);
                continue;
            }
            
            int ci = it->second;
            const NHCDefinition* nhc = nullptr;
            auto nhc_it = curve_to_nhc.find(j);
            if (nhc_it != curve_to_nhc.end()) nhc = nhc_it->second;
            
            auto pg = lookup_port_gauge(IF, clusters[ci], j, nhc);
            neighbor_gauges.push_back(pg);
            total_hv += pg.central_charge;
        }
        
        // Check (a): E8 subalgebra pair
        // For 2 neighbors: direct pair check
        // For 3+ neighbors: check all pairs (necessary but not sufficient)
        for (size_t i = 0; i < neighbor_gauges.size(); ++i) {
            for (size_t j = i + 1; j < neighbor_gauges.size(); ++j) {
                if (!is_e8_subalgebra_pair(neighbor_gauges[i].gauge_name, 
                                           neighbor_gauges[j].gauge_name)) {
                    std::ostringstream ss;
                    ss << "E8 subalgebra violation at -1 curve " << m1 
                       << ": " << neighbor_gauges[i].gauge_name 
                       << " × " << neighbor_gauges[j].gauge_name 
                       << " ⊄ E8";
                    fail_reason = ss.str();
                    return false;
                }
            }
        }
        
        // Check (b): h∨ budget
        if (total_hv > hv_bound + 1e-9) {
            std::ostringstream ss;
            ss << "h∨ budget violation at -1 curve " << m1 
               << ": sum = " << total_hv << " > " << hv_bound 
               << " (T=" << T_tensor << ", gauges:";
            for (const auto& pg : neighbor_gauges) ss << " " << pg.gauge_name;
            ss << ")";
            fail_reason = ss.str();
            return false;
        }
    }
    return true;
}

// ── Layer 3: Flavor Compatibility ──
//
// -1에 연결된 gauge의 matter representation이 호환되는지 체크.
// 
// 구체적으로: -1 위의 flavor symmetry가 양쪽 gauge에서 오는
// matter content를 수용할 수 있는지.
//
// 현재 구현: E8 subalgebra + h∨ budget이 대부분 커버.
// 추가 규칙이 필요한 경우 아래에 추가.

bool check_layer3_flavor_compat(const Eigen::MatrixXi& IF,
                                const std::vector<int>& minus1_curves,
                                const std::vector<std::vector<int>>& clusters,
                                const std::map<int, const NHCDefinition*>& curve_to_nhc,
                                std::string& fail_reason) {
    // USER: 추가 flavor rule 구현 시 여기에 작성
    // 
    // 예: 특정 gauge 조합에서 matter content가 inconsistent한 경우
    //     half-hyper 조건 위반 등
    (void)IF; (void)minus1_curves; (void)clusters; (void)curve_to_nhc; (void)fail_reason;
    return true;
}

// ── Main NHC Check: Three Layers Combined ──

NHCCheckResult check_nhc_clusters(const Tensor& t_in) {
    NHCCheckResult result;
    result.passes = true;
    
    Eigen::MatrixXi IF = t_in.GetIntersectionForm();
    int n = IF.rows();
    if (n == 0) return result;
    
    // ── Find -1 curves and non-(-1) curves ──
    std::vector<int> minus1_curves;
    std::vector<int> non_minus1;
    for (int i = 0; i < n; ++i) {
        if (IF(i, i) == -1) minus1_curves.push_back(i);
        else non_minus1.push_back(i);
    }
    
    if (non_minus1.empty()) return result;  // 전부 -1 → trivially pass
    
    // ── BFS to find connected components of non-(-1) curves ──
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
    
    // ── Layer 1: Cluster Recognition ──
    std::map<int, const NHCDefinition*> curve_to_nhc;
    std::string fail1;
    result.layer1_pass = check_layer1_clusters(IF, result.clusters, curve_to_nhc, fail1);
    if (!result.layer1_pass) {
        result.passes = false;
        result.failure_reason = "L1: " + fail1;
        return result;  // can't do L2/L3 without cluster identification
    }
    
    // ── Layer 2: E8 Central Charge Bound ──
    // T = number of curves in the base
    int T_tensor = n;
    std::string fail2;
    result.layer2_pass = check_layer2_e8_constraint(IF, minus1_curves, result.clusters, 
                                                curve_to_nhc, T_tensor, fail2);
    if (!result.layer2_pass) {
        result.passes = false;
        result.failure_reason = "L2: " + fail2;
    }
    
    // ── Layer 3: Flavor Compatibility ──
    std::string fail3;
    result.layer3_pass = check_layer3_flavor_compat(IF, minus1_curves, result.clusters,
                                                     curve_to_nhc, fail3);
    if (!result.layer3_pass) {
        result.passes = false;
        if (!result.failure_reason.empty()) result.failure_reason += "; ";
        result.failure_reason += "L3: " + fail3;
    }
    
    return result;
}

// Simple boolean wrapper
bool passes_nhc_check(const Tensor& t) {
    return check_nhc_clusters(t).passes;
}

// ============================================================================
// NEW: Noble Molecule Anomaly Table
// ============================================================================
// 
// Each noble molecule has intrinsic (H, V) contributions from its internal 
// gauge structure. These must be computed separately because noble molecules
// are NOT in the sidelink param system — they're built directly via AT/ATS.
//
// The table below has placeholder values {0, 0}. 
// USER MUST FILL THESE from Mathematica/paper calculations.
// 
// After filling, the anomaly for a noble-molecule-based LST is:
//   H_total = noble_H + extra_H
//   V_total = noble_V + extra_V
//   grav_anom = H_total - V_total + 29*T - 273

struct NobleMoleculeAnomaly {
    const char* name;
    int H;
    int V;
};

// clang-format off
const std::map<std::string, AnomalyContrib> NOBLE_MOLECULE_ANOMALY = {
    // ── noble 2-Molecules with no -5 curves ──
    {"1{1,3}22",              {0, 0}},
    
    // ── noble 0-Molecules with no -5 curves ──
    {"2{2,3}13",              {0, 0}},
    {"2{2,3}132",             {0, 0}},
    {"2{2,3}1322",            {0, 0}},
    {"23132",                 {0, 0}},
    {"223132",                {0, 0}},
    {"3132",                  {0, 0}},
    {"2132",                  {0, 0}},
    {"3123",                  {0, 0}},
    {"2231322",               {0, 0}},
    {"31322",                 {0, 0}},
    {"21322",                 {0, 0}},
    {"313",                   {0, 0}},
    {"23213",                 {0, 0}},
    {"213",                   {0, 0}},
    
    // ── noble 4-Molecules with one -5 curve ──
    {"1{1,{1,5}}1",           {0, 0}},
    {"1{1{1,{1,5}}}13",       {0, 0}},
    {"1{1{1,{1,5}}}132",      {0, 0}},
    {"1{1,5}1",               {0, 0}},
    
    // ── noble 2-Molecules with one -5 curve ──
    {"151{1,3}2",             {0, 0}},
    {"1{1,5}12",              {0, 0}},
    {"31{1,5}1",              {0, 0}},
    {"231{1,5}1",             {0, 0}},
    {"2231{1,5}1",            {0, 0}},
    {"2231{1,5}131",          {0, 0}},
    
    // ── noble 1-Molecules with one -5 curve ──
    {"51{1,3}2",              {0, 0}},
    {"3151{1,3}2",            {0, 0}},
    {"23151{1,3}2",           {0, 0}},
    {"223151{1,3}2",          {0, 0}},
    {"51{1,3}22",             {0, 0}},
    {"512{1,3}2",             {0, 0}},
    {"31{1,5}12",             {0, 0}},
    {"231{1,5}12",            {0, 0}},
    {"2231{1,5}12",           {0, 0}},
    {"31{1,5}13",             {0, 0}},
    {"231{1,5}13",            {0, 0}},
    {"2231{1,5}13",           {0, 0}},
    {"231{1,5}132",           {0, 0}},
    {"2231{1,5}132",          {0, 0}},
    {"2231{1,5}1322",         {0, 0}},
    {"13215132",              {0, 0}},
    {"223151231",             {0, 0}},
    {"3151231",               {0, 0}},
    
    // ── noble 0-Molecules with one -5 curve ──
    {"31{3{1,5}13}",          {0, 0}},
    {"231{3{1,5}13}",         {0, 0}},
    {"2231{3{1,5}13}",        {0, 0}},
    {"231{3{1,5}132}",        {0, 0}},
    {"2231{3{1,5}132}",       {0, 0}},
    {"3215",                  {0, 0}},
    {"2315",                  {0, 0}},
    {"32215",                 {0, 0}},
    {"22315",                 {0, 0}},
    {"315",                   {0, 0}},
    {"23215",                 {0, 0}},
    {"215",                   {0, 0}},
    {"2215",                  {0, 0}},
    {"22215",                 {0, 0}},
    {"313215",                {0, 0}},
    {"231315",                {0, 0}},
    {"31315",                 {0, 0}},
    
    // ── noble molecules with one -5 (longer chains) ──
    {"3215132",               {0, 0}},
    {"2315132",               {0, 0}},
    {"22315132",              {0, 0}},
    {"315132",                {0, 0}},
    {"23215132",              {0, 0}},
    {"215132",                {0, 0}},
    {"2215132",               {0, 0}},
    {"22315123",              {0, 0}},
    {"315123",                {0, 0}},
    {"215123",                {0, 0}},
    {"223151322",             {0, 0}},
    {"3151322",               {0, 0}},
    {"232151322",             {0, 0}},
    {"2151322",               {0, 0}},
    {"22151322",              {0, 0}},
    {"31513",                 {0, 0}},
    {"2321513",               {0, 0}},
    {"21513",                 {0, 0}},
    {"221513",                {0, 0}},
    {"2151232",               {0, 0}},
    {"21512",                 {0, 0}},
    {"31315132",              {0, 0}},
    {"313151322",             {0, 0}},
    {"3131513",               {0, 0}},
    
    // ── noble 1-Molecules with two -5 curves ──
    {"2231513151321",         {0, 0}},
    {"13151315132",           {0, 0}},
    {"131513151322",          {0, 0}},
    
    // ── noble 0-Molecules with two -5 curves ──
    {"513215",                {0, 0}},
    {"5132215",               {0, 0}},
    {"51315",                 {0, 0}},
    {"5123215",               {0, 0}},
    {"231513215",             {0, 0}},
    {"2231513215",            {0, 0}},
    {"31513215",              {0, 0}},
    {"21513215",              {0, 0}},
    {"31512315",              {0, 0}},
    {"2315132215",            {0, 0}},
    {"22315132215",           {0, 0}},
    {"315132215",             {0, 0}},
    {"32151315",              {0, 0}},
    {"23151315",              {0, 0}},
    {"223151315",             {0, 0}},
    {"3151315",               {0, 0}},
    {"2151315",               {0, 0}},
    {"23151315132",           {0, 0}},
    {"223151315132",          {0, 0}},
    {"3151315132",            {0, 0}},
    {"2151315132",            {0, 0}},
    {"22315131513222",        {0, 0}},
    {"31513151322",           {0, 0}},
    {"21513151322",           {0, 0}},
    {"315131513",             {0, 0}},
    {"215131513",             {0, 0}},
    
    // ── noble 0-Molecules with three -5 curves ──
    {"5131513215",            {0, 0}},
    {"513151315",             {0, 0}},
    
    // ── additional noble molecules ──
    {"1{1,5}13151",           {0, 0}},
    {"15132151",              {0, 0}},
    {"21513151",              {0, 0}},
    {"1513151315",            {0, 0}},
    {"51232151",              {0, 0}},
};
// clang-format on

AnomalyContrib get_noble_molecule_contrib(const std::string& name) {
    auto it = NOBLE_MOLECULE_ANOMALY.find(name);
    if (it != NOBLE_MOLECULE_ANOMALY.end()) {
        return it->second;
    }
    std::cerr << "WARNING: no anomaly entry for noble molecule '" << name << "'\n";
    return {0, 0};
}

// ============================================================================
// Result Structures
// ============================================================================

struct LSTResult {
    Eigen::MatrixXi IF;
    std::string latex;
    std::string description;
    std::vector<double> spectrum;
    
    void compute_spectrum() {
        spectrum = compute_eigenvalue_spectrum(IF);
    }
};

// NEW: Extended result for catalog mode
struct CatalogResult {
    // From LSTResult
    Eigen::MatrixXi IF;
    std::string latex;
    std::string description;
    std::vector<double> spectrum;
    
    // Source info
    enum class Source { SideLink, NobleMolecule };
    Source source;
    int sidelink_param;        // if source == SideLink
    std::string noble_name;    // if source == NobleMolecule
    int port;
    int extra_curve;
    int int_num;
    
    // Physics
    int T;           // tensor multiplets (= sig_neg)
    int H, V;        // anomaly contributions
    int grav_anom;   // H - V + 29T - 273
    
    // Geometry
    SignatureInfo sig;
    int n_curves;
    
    // NHC
    bool nhc_ok;
    std::string nhc_reason;
    
    void compute_spectrum() {
        spectrum = compute_eigenvalue_spectrum(IF);
    }
    
    // Structure descriptor for catalog
    std::string get_structure_desc() const {
        std::ostringstream ss;
        if (source == Source::SideLink) {
            ss << "SL:" << sidelink_param << ":" << port 
               << ":" << extra_curve << ":" << int_num;
        } else {
            ss << "NM:" << noble_name << ":" << port 
               << ":" << extra_curve << ":" << int_num;
        }
        return ss.str();
    }
};

// ============================================================================
// Deduplication (works on both types)
// ============================================================================

std::vector<LSTResult> remove_duplicates_by_spectrum(std::vector<LSTResult>& results) {
    for (auto& r : results) r.compute_spectrum();
    
    std::vector<LSTResult> unique_results;
    for (const auto& r : results) {
        bool is_duplicate = false;
        for (const auto& u : unique_results) {
            if (are_spectra_equal(r.spectrum, u.spectrum)) {
                is_duplicate = true;
                break;
            }
        }
        if (!is_duplicate) unique_results.push_back(r);
    }
    return unique_results;
}

std::vector<CatalogResult> remove_duplicates_catalog(std::vector<CatalogResult>& results) {
    for (auto& r : results) r.compute_spectrum();
    
    std::vector<CatalogResult> unique;
    for (const auto& r : results) {
        bool is_dup = false;
        for (const auto& u : unique) {
            if (are_spectra_equal(r.spectrum, u.spectrum)) {
                is_dup = true;
                break;
            }
        }
        if (!is_dup) unique.push_back(r);
    }
    return unique;
}

// ============================================================================
// Matrix to Quiver LaTeX (UNCHANGED)
// ============================================================================

std::string make_extra_marker(int extra_self_int, int extra_int_num) {
    std::string extra_str = std::to_string(extra_self_int);
    if (extra_int_num == 1) {
        return "\\textcolor{red}{" + extra_str + "}";
    } else {
        return "\\textcolor{red}{" + extra_str + "^{" + std::to_string(extra_int_num) + "}}";
    }
}

std::string matrix_to_quiver_latex(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    if (n == 0) return "\\emptyset";
    if (n == 1) return "{" + std::to_string(-IF(0,0)) + "}";
    
    int extra_attached_to = -1;
    int extra_int_num = 0;
    int extra_self_int = -IF(n-1, n-1);
    int base_n = n - 1;
    
    for (int j = 0; j < base_n; ++j) {
        if (IF(n-1, j) != 0) {
            extra_attached_to = j;
            extra_int_num = IF(n-1, j);
            break;
        }
    }
    
    std::vector<std::string> result_parts;
    std::vector<bool> used(base_n, false);
    
    int i = 0;
    while (i < base_n) {
        if (used[i]) { i++; continue; }
        
        used[i] = true;
        int diag_i = -IF(i, i);
        
        if (i + 1 < base_n && IF(i, i+1) != 0 && !used[i+1]) {
            std::string part = std::to_string(diag_i);
            if (extra_attached_to == i) {
                part = "\\overset{" + make_extra_marker(extra_self_int, extra_int_num) + "}{" + part + "}";
            }
            result_parts.push_back(part);
            i = i + 1;
        }
        else {
            int j = -1;
            for (int k = i + 2; k < base_n; ++k) {
                if (IF(i, k) != 0 && !used[k]) { j = k; break; }
            }
            
            if (j == -1) {
                std::string part = std::to_string(diag_i);
                if (extra_attached_to == i) {
                    part = "\\overset{" + make_extra_marker(extra_self_int, extra_int_num) + "}{" + part + "}";
                }
                result_parts.push_back(part);
                i++;
                continue;
            }
            
            std::string part_i = std::to_string(diag_i);
            if (extra_attached_to == i) {
                part_i = "\\overset{" + make_extra_marker(extra_self_int, extra_int_num) + "}{" + part_i + "}";
            }
            result_parts.push_back(part_i);
            
            std::string stacked = "";
            for (int m = i + 1; m < j; ++m) {
                used[m] = true;
                std::string mpart = std::to_string(-IF(m, m));
                if (extra_attached_to == m) {
                    mpart = "\\overset{" + make_extra_marker(extra_self_int, extra_int_num) + "}{" + mpart + "}";
                }
                stacked += mpart;
            }
            
            int diag_j = -IF(j, j);
            std::string part_j = std::to_string(diag_j);
            
            if (!stacked.empty()) {
                part_j = "\\overset{" + stacked + "}{" + part_j + "}";
            }
            if (extra_attached_to == j) {
                if (!stacked.empty()) {
                    part_j = "\\underset{" + make_extra_marker(extra_self_int, extra_int_num) + "}{" + part_j + "}";
                } else {
                    part_j = "\\overset{" + make_extra_marker(extra_self_int, extra_int_num) + "}{" + part_j + "}";
                }
            }
            
            result_parts.push_back(part_j);
            used[j] = true;
            i = j + 1;
        }
    }
    
    std::ostringstream ss;
    ss << "{";
    for (size_t k = 0; k < result_parts.size(); ++k) {
        ss << result_parts[k];
    }
    ss << "}";
    return ss.str();
}

// ============================================================================
// Noble Molecules (UNCHANGED)
// ============================================================================

struct NobleMolecule {
    std::string name;
    std::function<Tensor()> builder;
};

std::vector<NobleMolecule> get_noble_molecules() {
    return {
        // === Noble molecules from Appendix D ===
	
	// noble 2-Molecules with no -5 curves 
	{"1{1,3}22", []{ Tensor t; t.AT(-1); t.ATS(-1,-3); t.AT(-2); t.AT(-2); return t; }},
        
	// noble 0-Molecules with no -5 curves
        {"2{2,3}13", []{ Tensor t; t.AT(-2); t.ATS(-2,-3); t.AT(-1); t.AT(-3); return t; }},
        {"2{2,3}132", []{ Tensor t; t.AT(-2); t.ATS(-2,-3); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"2{2,3}1322", []{ Tensor t; t.AT(-2); t.ATS(-2,-3); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"23132", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"223132", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"3132", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"2132", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"3123", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-2); t.AT(-3); return t; }},
	{"2231322", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"31322", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"21322", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"313", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); return t; }},
	{"23213", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-3); return t; }},
	{"213", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-3); return t; }},

	// noble 4-Molecules with one -5 curve
	{"1{1,{1,5}}1", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-1); t.not_intersect(4,5); t.intersect(5,3); return t; }},
	{"1{1{1,{1,5}}}13", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-1); t.not_intersect(4,5); t.intersect(5,3); t.AT(-3); t.not_intersect(6,5); t.intersect(6,3); return t; }},
	{"1{1{1,{1,5}}}132", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-1); t.not_intersect(4,5); t.intersect(5,3); t.AT(-3); t.not_intersect(6,5); t.intersect(6,3); t.AT(-2); return t; }},
        {"1{1,5}1", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); return t; }},

	// noble 2-Molecules with one -5 curve
	{"151{1,3}2", []{ Tensor t; t.AT(-1); t.AT(-5); t.AT(-1); t.ATS(-1,-3); t.AT(-2); return t; }},
	{"1{1,5}12", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-2); return t; }},
	{"31{1,5}1", []{ Tensor t; t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); return t; }},
	{"231{1,5}1", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); return t; }},
	{"2231{1,5}1", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); return t; }},
	{"2231{1,5}131", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-1); return t; }},

	// noble 1-Molecules with one -5 curve 
	{"51{1,3}2", []{ Tensor t; t.AT(-5); t.AT(-1); t.ATS(-1,-3); t.AT(-2); return t; }},
	{"3151{1,3}2", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.ATS(-1,-3); t.AT(-2); return t; }},
	{"23151{1,3}2", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.ATS(-1,-3); t.AT(-2); return t; }},
	{"223151{1,3}2", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.ATS(-1,-3); t.AT(-2); return t; }},
	{"51{1,3}22", []{ Tensor t; t.AT(-5); t.AT(-1); t.ATS(-1,-3); t.AT(-2); t.AT(-2); return t; }},
	{"512{1,3}2", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-2); t.ATS(-1,-3); t.AT(-2); return t; }},
	{"31{1,5}12", []{ Tensor t; t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-2); return t; }},
	{"231{1,5}12", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-2); return t; }},
	{"2231{1,5}12", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-2); return t; }},
	{"31{1,5}13", []{ Tensor t; t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); return t; }},
	{"231{1,5}13", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); return t; }},
	{"2231{1,5}13", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); return t; }},
	{"231{1,5}132", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"2231{1,5}132", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"2231{1,5}1322", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
	{"13215132", []{ Tensor t; t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"223151231", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-1); return t; }},
	{"3151231", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-1); return t; }},

	// noble 0-Molecules with one -5 curve
	{"31{3{1,5}13}", []{ Tensor t; t.AT(-3); t.AT(-1); t.ATS2(-3,-1,-5); t.AT(-1); t.AT(-3); return t; }},
	{"231{3{1,5}13}", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.ATS2(-3,-1,-5); t.AT(-1); t.AT(-3); return t; }},
	{"2231{3{1,5}13}", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS2(-3,-1,-5); t.AT(-1); t.AT(-3); return t; }},
	{"231{3{1,5}132}", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.ATS2(-3,-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"2231{3{1,5}132}", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS2(-3,-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"3215", []{ Tensor t; t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"2315", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"32215", []{ Tensor t; t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"22315", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"315", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"23215", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"215", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"2215", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"22215", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"313215", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"231315", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"31315", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},

	// noble molecules with -5 curves (longer chains)
	{"3215132", []{ Tensor t; t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"2315132", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"22315132", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"315132", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"23215132", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"215132", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"2215132", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"22315123", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); return t; }},
	{"315123", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); return t; }},
	{"215123", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); return t; }},
	{"223151322", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
	{"3151322", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
	{"232151322", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
	{"2151322", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
	{"22151322", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
	{"31513", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},
	{"2321513", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},
	{"21513", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},
	{"221513", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},
	{"2151232", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-2); return t; }},
	{"21512", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); return t; }},
	{"31315132", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"313151322", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
	{"3131513", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},

	// noble 1-Molecules with two -5 curves
	{"2231513151321", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); return t; }},
	{"13151315132", []{ Tensor t; t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"131513151322", []{ Tensor t; t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},

	// noble 0-Molecules with two -5 curves
	{"513215", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"5132215", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"51315", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"5123215", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"231513215", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"2231513215", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"31513215", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"21513215", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"31512315", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"2315132215", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"22315132215", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"315132215", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"32151315", []{ Tensor t; t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"23151315", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"223151315", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"3151315", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"2151315", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"23151315132", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"223151315132", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"3151315132", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"2151315132", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
	{"22315131513222", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
	{"31513151322", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
	{"21513151322", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
	{"315131513", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},
	{"215131513", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},

	// noble 0-Molecules with three -5 curves
	{"5131513215", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
	{"513151315", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},

	// additional noble molecules
	{"1{1,5}13151", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); return t; }},
	{"15132151", []{ Tensor t; t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); return t; }},
	{"21513151", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); return t; }},
	{"1513151315", []{ Tensor t; t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"51232151", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); return t; }},
    };
}

// ============================================================================
// Collect functions: original (for backward compatibility)
// ============================================================================

std::vector<LSTResult> collect_sidelink_results(bool skip_instanton = false) {
    std::vector<LSTResult> results;
    
    for (int param : ALL_SIDELINK_PARAMS) {
        if (skip_instanton && INSTANTON_PARAMS.count(param)) continue;
        Tensor base = build_tensor(s(param));
        if (base.GetT() == 0) continue;
        
        auto IF = base.GetIntersectionForm();
        int n = IF.rows();
        
        for (int port = 1; port <= n; ++port) {
            int base_self_int = IF(port - 1, port - 1);
            auto allowed_extras = get_allowed_extra_curves(base_self_int);
            
            for (int extra_curve : allowed_extras) {
                auto allowed = get_allowed_intersections(extra_curve);
                
                for (int int_num : allowed) {
                    Tensor t_ext;
                    t_ext.SetIF(IF);
                    t_ext.ATE(extra_curve, port, int_num);
                    
                    if (is_valid_LST(t_ext) && passes_nhc_check(t_ext)) {
                        LSTResult result;
                        result.IF = t_ext.GetIntersectionForm();
                        result.latex = matrix_to_quiver_latex(result.IF);
                        
                        std::ostringstream desc;
                        desc << "sidelink=" << param 
                             << ", port=" << port 
                             << ", extra=" << extra_curve
                             << ", int=" << int_num;
                        result.description = desc.str();
                        results.push_back(result);
                    }
                }
            }
        }
    }
    return results;
}

std::vector<LSTResult> collect_noble_molecule_results() {
    std::vector<LSTResult> results;
    
    for (const auto& nm : get_noble_molecules()) {
        Tensor base = nm.builder();
        if (base.GetT() == 0) continue;
        
        auto IF = base.GetIntersectionForm();
        int n = IF.rows();
        
        for (int port = 1; port <= n; ++port) {
            int base_self_int = IF(port - 1, port - 1);
            auto allowed_extras = get_allowed_extra_curves(base_self_int);
            
            for (int extra_curve : allowed_extras) {
                auto allowed = get_allowed_intersections(extra_curve);
                
                for (int int_num : allowed) {
                    Tensor t_ext;
                    t_ext.SetIF(IF);
                    t_ext.ATE(extra_curve, port, int_num);
                    
                    if (is_valid_LST(t_ext) && passes_nhc_check(t_ext)) {
                        LSTResult result;
                        result.IF = t_ext.GetIntersectionForm();
                        result.latex = matrix_to_quiver_latex(result.IF);
                        
                        std::ostringstream desc;
                        desc << "noble=" << nm.name
                             << ", port=" << port
                             << ", extra=" << extra_curve
                             << ", int=" << int_num;
                        result.description = desc.str();
                        results.push_back(result);
                    }
                }
            }
        }
    }
    return results;
}

// ============================================================================
// NEW: Catalog-mode collect functions
// ============================================================================
//
// Same iteration as originals, but compute:
//   - Anomaly (H, V) from anomaly_tables.h / noble_molecule table
//   - Signature (sig+, sig-, sig0, det) from eigenvalues  
//   - NHC check via Tensor blowdown
//   - LaTeX from matrix_to_quiver_latex

std::vector<CatalogResult> collect_sidelink_results_catalog(bool verbose, bool skip_instanton = false) {
    std::vector<CatalogResult> results;
    int tried = 0;
    
    for (int param : ALL_SIDELINK_PARAMS) {
        if (skip_instanton && INSTANTON_PARAMS.count(param)) continue;
        Tensor base = build_tensor(s(param));
        if (base.GetT() == 0) continue;
        
        auto IF = base.GetIntersectionForm();
        int n = IF.rows();
        
        // Get base anomaly contribution
        AnomalyContrib base_anom = get_sidelink_contrib(param);
        
        for (int port = 1; port <= n; ++port) {
            int base_self_int = IF(port - 1, port - 1);
            auto allowed_extras = get_allowed_extra_curves(base_self_int);
            
            for (int extra_curve : allowed_extras) {
                auto allowed = get_allowed_intersections(extra_curve);
                
                for (int int_num : allowed) {
                    tried++;
                    
                    Tensor t_ext;
                    t_ext.SetIF(IF);
                    t_ext.ATE(extra_curve, port, int_num);
                    
                    if (!is_valid_LST(t_ext)) continue;
                    
                    CatalogResult cr;
                    cr.IF = t_ext.GetIntersectionForm();
                    cr.n_curves = cr.IF.rows();
                    cr.latex = matrix_to_quiver_latex(cr.IF);
                    
                    // Source info
                    cr.source = CatalogResult::Source::SideLink;
                    cr.sidelink_param = param;
                    cr.port = port;
                    cr.extra_curve = extra_curve;
                    cr.int_num = int_num;
                    
                    // Description (for backward compat)
                    std::ostringstream desc;
                    desc << "sidelink=" << param << ", port=" << port 
                         << ", extra=" << extra_curve << ", int=" << int_num;
                    cr.description = desc.str();
                    
                    // Signature
                    cr.sig = compute_signature(cr.IF);
                    cr.T = cr.sig.sig_neg;
                    
                    // Anomaly: base + extra
                    AnomalyContrib extra_anom = get_external_contrib(std::abs(extra_curve));
                    cr.H = base_anom.H + extra_anom.H;
                    cr.V = base_anom.V + extra_anom.V;
                    cr.grav_anom = cr.H - cr.V + 29 * cr.T - 273;
                    
                    // NHC check
                    auto nhc = check_nhc_clusters(t_ext);
                    cr.nhc_ok = nhc.passes;
                    cr.nhc_reason = nhc.failure_reason;
                    
                    results.push_back(cr);
                }
            }
        }
    }
    
    if (verbose) {
        std::cout << "  Sidelink: tried " << tried 
                  << ", valid LSTs " << results.size() << "\n";
    }
    return results;
}

std::vector<CatalogResult> collect_noble_molecule_results_catalog(bool verbose) {
    std::vector<CatalogResult> results;
    int tried = 0;
    
    for (const auto& nm : get_noble_molecules()) {
        Tensor base = nm.builder();
        if (base.GetT() == 0) continue;
        
        auto IF = base.GetIntersectionForm();
        int n = IF.rows();
        
        // Get base anomaly contribution from noble molecule table
        AnomalyContrib base_anom = get_noble_molecule_contrib(nm.name);
        
        for (int port = 1; port <= n; ++port) {
            int base_self_int = IF(port - 1, port - 1);
            auto allowed_extras = get_allowed_extra_curves(base_self_int);
            
            for (int extra_curve : allowed_extras) {
                auto allowed = get_allowed_intersections(extra_curve);
                
                for (int int_num : allowed) {
                    tried++;
                    
                    Tensor t_ext;
                    t_ext.SetIF(IF);
                    t_ext.ATE(extra_curve, port, int_num);
                    
                    if (!is_valid_LST(t_ext)) continue;
                    
                    CatalogResult cr;
                    cr.IF = t_ext.GetIntersectionForm();
                    cr.n_curves = cr.IF.rows();
                    cr.latex = matrix_to_quiver_latex(cr.IF);
                    
                    // Source info
                    cr.source = CatalogResult::Source::NobleMolecule;
                    cr.noble_name = nm.name;
                    cr.port = port;
                    cr.extra_curve = extra_curve;
                    cr.int_num = int_num;
                    
                    // Description
                    std::ostringstream desc;
                    desc << "noble=" << nm.name << ", port=" << port 
                         << ", extra=" << extra_curve << ", int=" << int_num;
                    cr.description = desc.str();
                    
                    // Signature
                    cr.sig = compute_signature(cr.IF);
                    cr.T = cr.sig.sig_neg;
                    
                    // Anomaly: base + extra
                    AnomalyContrib extra_anom = get_external_contrib(std::abs(extra_curve));
                    cr.H = base_anom.H + extra_anom.H;
                    cr.V = base_anom.V + extra_anom.V;
                    cr.grav_anom = cr.H - cr.V + 29 * cr.T - 273;
                    
                    // NHC check
                    auto nhc = check_nhc_clusters(t_ext);
                    cr.nhc_ok = nhc.passes;
                    cr.nhc_reason = nhc.failure_reason;
                    
                    results.push_back(cr);
                }
            }
        }
    }
    
    if (verbose) {
        std::cout << "  Noble: tried " << tried 
                  << ", valid LSTs " << results.size() << "\n";
    }
    return results;
}

// ============================================================================
// NEW: Catalog Output
// ============================================================================
//
// Format:
//   id | NN | structure | T | H | V | grav_anom | det | (s+,s-,s0) | n | IF_flat | latex | NHC_OK/NHC_FAIL
//
// IF_flat: row-major comma-separated
// structure: "SL:param:port:ext:int" or "NM:name:port:ext:int"

void write_catalog(const std::string& filename, 
                   const std::vector<CatalogResult>& results,
                   bool nhc_filter) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Cannot open " << filename << " for writing\n";
        return;
    }
    
    // Unified catalog format (compatible with lst_to_catalog Nk output)
    out << "# Unified LST Catalog (No-Node entries)\n";
    out << "# id | type | T | H | V | topo\n";
    
    int nhc_pass = 0, nhc_fail = 0;
    int written = 0;
    
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        
        if (r.nhc_ok) nhc_pass++;
        else nhc_fail++;
        
        if (nhc_filter && !r.nhc_ok) continue;
        
        out << written << " | NN | " 
            << r.T << " | " << r.H << " | " << r.V << " | " 
            << r.get_structure_desc() << "\n";
        
        written++;
    }
    
    out << "# Summary: " << written << " entries written\n";
    out << "# NHC pass: " << nhc_pass << ", NHC fail: " << nhc_fail << "\n";
    
    out.close();
    
    std::cout << "Catalog written to: " << filename << "\n";
    std::cout << "  Entries: " << written << "\n";
    std::cout << "  NHC pass: " << nhc_pass << ", NHC fail: " << nhc_fail << "\n";
}

// ============================================================================
// NEW: NHC Statistics
// ============================================================================

void print_nhc_statistics(const std::vector<CatalogResult>& results) {
    int pass = 0, fail = 0;
    std::map<std::string, int> failure_reasons;
    
    for (const auto& r : results) {
        if (r.nhc_ok) {
            pass++;
        } else {
            fail++;
            if (!r.nhc_reason.empty()) {
                failure_reasons[r.nhc_reason]++;
            }
        }
    }
    
    std::cout << "\n=== NHC Check Statistics ===\n";
    std::cout << "Pass: " << pass << " / " << results.size() << "\n";
    std::cout << "Fail: " << fail << " / " << results.size() << "\n";
    
    if (!failure_reasons.empty()) {
        std::cout << "\nFailure reasons:\n";
        for (const auto& [reason, count] : failure_reasons) {
            std::cout << "  " << reason << ": " << count << "\n";
        }
    }
}

// ============================================================================
// Output Functions (UNCHANGED from original)
// ============================================================================

void write_intersection_form(std::ostream& out, const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j > 0) out << " ";
            out << IF(i, j);
        }
        out << "\n";
    }
}

void write_intersection_forms_file(const std::string& filename, 
                                   const std::vector<LSTResult>& results) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Cannot open " << filename << " for writing\n";
        return;
    }
    
    for (size_t i = 0; i < results.size(); ++i) {
        write_intersection_form(out, results[i].IF);
        if (i + 1 < results.size()) out << "\n";
    }
    
    out.close();
    std::cout << "Intersection forms written to: " << filename << "\n";
}

void write_results_to_files(std::ofstream& tex_out, std::ofstream& data_out,
                            const std::vector<LSTResult>& results,
                            const std::string& section_title,
                            int& count) {
    tex_out << "% === " << section_title << " ===\n\n";
    
    std::vector<std::string> latex_results;
    
    for (const auto& r : results) {
        ++count;
        latex_results.push_back(r.latex);
        
        data_out << count << ": " << r.description << "\n";
        data_out << "  LaTeX: " << r.latex << "\n";
        data_out << "  Eigenvalues: [";
        for (size_t i = 0; i < r.spectrum.size(); ++i) {
            if (i > 0) data_out << ", ";
            data_out << r.spectrum[i];
        }
        data_out << "]\n";
    }
    
    if (!latex_results.empty()) {
        tex_out << "\\begin{align*}\n";
        for (size_t i = 0; i < latex_results.size(); ++i) {
            if (i % 3 == 0) {
                if (i > 0) tex_out << " \\\\\n";
                tex_out << "&" << latex_results[i];
            } else {
                tex_out << ", &&" << latex_results[i];
            }
        }
        tex_out << "\n\\end{align*}\n\n";
    }
}

// ============================================================================
// Main
// ============================================================================

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n\n"
              << "Original mode (default):\n"
              << "  " << prog << " [tex_file] [data_file] [if_file]\n\n"
              << "Catalog mode:\n"
              << "  " << prog << " --catalog <output_file> [options]\n\n"
              << "Options:\n"
              << "  --catalog FILE    Output catalog-format file\n"
              << "  --no-nhc-filter   Include NHC-failing entries (tagged NHC_FAIL)\n"
              << "  --no-instanton    Exclude instanton sidelinks\n"
              << "  --nhc-stats       Print NHC check statistics\n"
              << "  -v, --verbose     Verbose output\n"
              << "  -h, --help        Show this help\n";
}

int main(int argc, char** argv) {
    // ===== Parse arguments =====
    std::string tex_filename = "no_node_LSTs.tex";
    std::string data_filename = "no_node_LSTs.txt";
    std::string if_filename = "no_node_LSTs_IF.txt";
    std::string catalog_filename;
    bool catalog_mode = false;
    bool nhc_filter = true;    // default: filter out NHC failures
    bool nhc_stats = false;
    bool verbose = false;
    bool skip_instanton = false;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--catalog" && i + 1 < argc) {
            catalog_filename = argv[++i];
            catalog_mode = true;
        } else if (arg == "--no-nhc-filter") {
            nhc_filter = false;
        } else if (arg == "--nhc-stats") {
            nhc_stats = true;
        } else if (arg == "--no-instanton") {
            skip_instanton = true;
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else {
            // Legacy positional arguments
            if (i == 1) tex_filename = arg;
            else if (i == 2) data_filename = arg;
            else if (i == 3) if_filename = arg;
        }
    }
    
    // ===== CATALOG MODE =====
    if (catalog_mode) {
        std::cout << "=== No-Node LST Catalog Builder ===\n";
        std::cout << "Output: " << catalog_filename << "\n";
        std::cout << "NHC filter: " << (nhc_filter ? "ON" : "OFF") << "\n";
        std::cout << "Instanton: " << (skip_instanton ? "EXCLUDED" : "INCLUDED") << "\n\n";
        
        // Collect with catalog data
        std::cout << "Collecting sidelink results...\n";
        auto sl_results = collect_sidelink_results_catalog(verbose, skip_instanton);
        std::cout << "Found " << sl_results.size() << " from side links (before dedup)\n";
        
        std::cout << "Collecting noble molecule results...\n";
        auto nm_results = collect_noble_molecule_results_catalog(verbose);
        std::cout << "Found " << nm_results.size() << " from noble molecules (before dedup)\n";
        
        // Combine
        std::vector<CatalogResult> all_results;
        all_results.reserve(sl_results.size() + nm_results.size());
        all_results.insert(all_results.end(), sl_results.begin(), sl_results.end());
        all_results.insert(all_results.end(), nm_results.begin(), nm_results.end());
        
        std::cout << "\nTotal before deduplication: " << all_results.size() << "\n";
        
        // Dedup
        std::cout << "Removing duplicates by eigenvalue spectrum...\n";
        auto unique_results = remove_duplicates_catalog(all_results);
        std::cout << "After deduplication: " << unique_results.size() << "\n";
        std::cout << "Removed " << (all_results.size() - unique_results.size()) << " duplicates\n";
        
        // NHC stats
        if (nhc_stats || verbose) {
            print_nhc_statistics(unique_results);
        }
        
        // Write catalog
        write_catalog(catalog_filename, unique_results, nhc_filter);
        
        std::cout << "\n=== DONE ===\n";
        return 0;
    }
    
    // ===== ORIGINAL MODE (UNCHANGED) =====
    std::ofstream tex_out(tex_filename);
    std::ofstream data_out(data_filename);
    
    if (!tex_out || !data_out) {
        std::cerr << "Error: Cannot open output files\n";
        return 1;
    }
    
    std::cout << "Collecting side link results...\n";
    std::vector<LSTResult> sidelink_results = collect_sidelink_results(skip_instanton);
    std::cout << "Found " << sidelink_results.size() << " from side links (before dedup)\n";
    
    std::cout << "Collecting noble molecule results...\n";
    std::vector<LSTResult> noble_results = collect_noble_molecule_results();
    std::cout << "Found " << noble_results.size() << " from noble molecules (before dedup)\n";
    
    std::vector<LSTResult> all_results;
    all_results.reserve(sidelink_results.size() + noble_results.size());
    all_results.insert(all_results.end(), sidelink_results.begin(), sidelink_results.end());
    all_results.insert(all_results.end(), noble_results.begin(), noble_results.end());
    
    std::cout << "\nTotal before deduplication: " << all_results.size() << "\n";
    
    std::cout << "Removing duplicates by eigenvalue spectrum...\n";
    std::vector<LSTResult> unique_results = remove_duplicates_by_spectrum(all_results);
    std::cout << "After deduplication: " << unique_results.size() << "\n";
    std::cout << "Removed " << (all_results.size() - unique_results.size()) << " duplicates\n\n";
    
    write_intersection_forms_file(if_filename, unique_results);
    
    tex_out << "% No-Node LST Classification (P-type endpoints only)\n";
    tex_out << "% Generated automatically\n";
    tex_out << "% Extra curves shown in RED\n";
    tex_out << "% Duplicates removed by eigenvalue spectrum comparison\n";
    tex_out << "% Usage: \\input{no_node_LSTs.tex} in your main document\n";
    tex_out << "% Requires: \\usepackage{amsmath, xcolor}\n\n";
    
    tex_out << "\\subsubsection{Unique P-type LSTs (Side Links + Noble Molecules with extra curve)}\n\n";
    
    int count = 0;
    write_results_to_files(tex_out, data_out, unique_results, 
                          "UNIQUE P-TYPE LSTs (DEDUPLICATED BY EIGENVALUE SPECTRUM)", count);
    
    tex_out << "\n% Summary:\n";
    tex_out << "% Total before deduplication: " << all_results.size() << "\n";
    tex_out << "% Total after deduplication: " << unique_results.size() << "\n";
    tex_out << "% Duplicates removed: " << (all_results.size() - unique_results.size()) << "\n";
    
    tex_out.close();
    data_out.close();
    
    std::cout << "=== DONE ===\n";
    std::cout << "Total unique P-type LSTs: " << unique_results.size() << "\n";
    std::cout << "Output files:\n";
    std::cout << "  LaTeX: " << tex_filename << "\n";
    std::cout << "  Data:  " << data_filename << "\n";
    std::cout << "  IF:    " << if_filename << "\n";
    
    return 0;
}
