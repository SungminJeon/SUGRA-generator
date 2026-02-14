// gen_sugra.cpp
// SUGRA base generator: catalog → attach external → filter → LaTeX output
//
// Algorithm:
//   1. Load unified.cat (18,883 LST entries: 119 NN + 18,764 Nk)
//   2. For each catalog entry with T in [T_min, T_max]:
//      a. Reconstruct IF from topo string
//         - NN: parse "SL:param:port:ext:int" → build_tensor(s(param)) → ATE
//         - Nk: deserialize topocompactline → TheoryGraph → ComposeIF_Gluing
//      b. Find attachment candidates:
//         - For each rule in attachment table, find target curves in IF
//         - Check central charge budget on -1 curves: sum c_i(k) ≤ 8
//      c. For each candidate, attach external curve (append row/col to IF)
//      d. Filter:
//         - NHC check: all non-(-1) connected components must be known NHCs
//         - Central charge: WZW c(k) = k·dim(g)/(k+h∨) ≤ 8 per -1 curve
//         - H_neutral = 273 + V - H_charged - 29T ≥ 0
//         - |det(IF)| = 1 (unimodular)
//         - Signature: automatically (1, T) from construction
//   3. Output LaTeX with quiver topology + physics data
//
// Compile:
//   g++ -std=c++17 -O2 -o gen_sugra gen_sugra.cpp Tensor.C \
//       TopoLineCompact_enhanced.cpp Topology_enhanced.cpp \
//       -I/usr/include/eigen3 -I.

#include "sugra_generator.h"
#include <fstream>
#include <cstdlib>
#include <set>
#include <iomanip>
#include <numeric>

// ============================================================================
// Blowdown algorithm
// ============================================================================
//
// Repeatedly blow down (-1)-curves until no (-1)-curves remain.
// Returns the minimal surface type.

struct BlowdownResult {
    Eigen::MatrixXi minimal_IF;    // curve-only IF after all blowdowns
    Eigen::MatrixXi extended_IF;   // full (curves+b₀) matrix after blowdowns
    int num_blowdowns;
    std::string surface;           // "P^2", "F_n", or "UNKNOWN"
    int F_n;                       // n for F_n, -1 for P^2
    bool is_P2 = false;
    int b0_sq = 0;                 // b₀² from extended matrix after blowdown
};

inline Eigen::MatrixXi blowdown_curve(const Eigen::MatrixXi& IF, int idx) {
    int n = IF.rows();
    // Modify: for all j,k != idx, IF(j,k) += IF(j,idx)*IF(k,idx)
    Eigen::MatrixXi M = IF;
    for (int j = 0; j < n; j++) {
        if (j == idx) continue;
        for (int k = 0; k < n; k++) {
            if (k == idx) continue;
            M(j, k) += IF(j, idx) * IF(k, idx);
        }
    }
    // Remove row/col idx
    Eigen::MatrixXi result(n - 1, n - 1);
    int ri = 0;
    for (int i = 0; i < n; i++) {
        if (i == idx) continue;
        int ci = 0;
        for (int j = 0; j < n; j++) {
            if (j == idx) continue;
            result(ri, ci) = M(i, j);
            ci++;
        }
        ri++;
    }
    return result;
}

// Format matrix as LaTeX pmatrix (small)
std::string matrix_to_latex(const Eigen::MatrixXi& M) {
    int n = M.rows();
    std::string s = "\\scriptstyle\\begin{pmatrix}";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            s += std::to_string(M(i, j));
            if (j < n - 1) s += "&";
        }
        if (i < n - 1) s += "\\\\";
    }
    s += "\\end{pmatrix}";
    return s;
}

// Build extended (T+1)×(T+1) matrix with b₀Q row/col
// Following Tensor::Setb0Q() and GetIFb0Q()
inline Eigen::MatrixXi build_b0Q_matrix(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    // b₀² = 9 - T where T = number of negative eigenvalues (not matrix dimension)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
    int sig_neg = 0;
    for (int i = 0; i < es.eigenvalues().size(); i++)
        if (es.eigenvalues()(i) < -1e-8) sig_neg++;
    
    Eigen::MatrixXi m = Eigen::MatrixXi::Zero(n + 1, n + 1);
    // Copy IF
    m.block(0, 0, n, n) = IF;
    // b₀·C_i = IF(i,i) + 2
    for (int i = 0; i < n; i++) {
        m(i, n) = IF(i, i) + 2;
        m(n, i) = IF(i, i) + 2;
    }
    // b₀² = 9 - T_SUGRA (T_SUGRA = sig_neg)
    m(n, n) = 9 - sig_neg;
    return m;
}

// Compute T_max for non-unimodular bases:
// Build b₀Q extended matrix (no blowdown), sweep b₀² = 9-T,
// find T where det = 0 → T_max = T_det0 - 1
struct TmaxResult {
    int T_sugra;          // original T of the SUGRA base
    int T_max;            // maximum T before det vanishes
    double T_crit;        // exact T where det = 0
    bool exact;           // true if T_crit is integer
    Eigen::MatrixXi ext;  // extended matrix at original T
};

inline TmaxResult compute_Tmax(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    TmaxResult res;
    
    // Compute sig_neg = T_SUGRA
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es0(IF.cast<double>());
    res.T_sugra = 0;
    for (int i = 0; i < es0.eigenvalues().size(); i++)
        if (es0.eigenvalues()(i) < -1e-8) res.T_sugra++;
    
    // Build extended matrix with b₀·C_i, b₀² as variable
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n + 1, n + 1);
    M.block(0, 0, n, n) = IF.cast<double>();
    for (int i = 0; i < n; i++) {
        M(i, n) = IF(i, i) + 2;
        M(n, i) = IF(i, i) + 2;
    }
    
    // det(extended) is linear in b₀²: det = A * b₀² + B
    M(n, n) = 0.0;
    double B = M.determinant();
    M(n, n) = 1.0;
    double A = M.determinant() - B;
    
    // det = 0 at b₀² = -B/A, i.e. T_crit = 9 - (-B/A) = 9 + B/A
    if (std::abs(A) < 1e-10) {
        res.T_max = -1; res.T_crit = 1e9; res.exact = false;
    } else {
        double b0sq_crit = -B / A;
        double T_crit = 9.0 - b0sq_crit;
        res.T_crit = T_crit;
        int T_crit_int = (int)std::round(T_crit);
        res.exact = (std::abs(T_crit - T_crit_int) < 1e-6);
        if (T_crit <= res.T_sugra) {
            res.T_max = res.T_sugra;
        } else if (res.exact) {
            res.T_max = T_crit_int - 1;
        } else {
            res.T_max = (int)std::floor(T_crit);
        }
    }
    
    // Store extended matrix at original T
    Eigen::MatrixXi ext = Eigen::MatrixXi::Zero(n + 1, n + 1);
    ext.block(0, 0, n, n) = IF;
    for (int i = 0; i < n; i++) {
        ext(i, n) = IF(i, i) + 2;
        ext(n, i) = IF(i, i) + 2;
    }
    ext(n, n) = 9 - res.T_sugra;
    res.ext = ext;
    
    return res;
}

BlowdownResult blowdown(const Eigen::MatrixXi& IF) {
    // Build extended matrix with b₀ as last row/col
    Eigen::MatrixXi cur = build_b0Q_matrix(IF);
    int bd = 0;
    
    while (true) {
        int n = cur.rows();
        // Find a (-1)-curve among curve indices (all except last = b₀)
        int m1 = -1;
        for (int i = 0; i < n - 1; i++) {  // skip last index (b₀)
            if (cur(i, i) == -1) { m1 = i; break; }
        }
        if (m1 < 0) break;
        cur = blowdown_curve(cur, m1);
        bd++;
    }
    
    BlowdownResult res;
    int n = cur.rows();  // includes b₀ row
    int n_curves = n - 1;  // number of remaining curve generators
    res.num_blowdowns = bd;
    res.F_n = -1;
    res.b0_sq = cur(n - 1, n - 1);  // last diagonal = final b₀²
    res.extended_IF = cur;  // full matrix with b₀ row/col
    
    // Extract curve-only IF (without b₀ row/col) for surface identification
    res.minimal_IF = cur.block(0, 0, n_curves, n_curves);
    
    if (n_curves == 1 && cur(0, 0) > 0) {
        res.surface = "\\mathbb{P}^2";
        res.F_n = -1;
        res.is_P2 = true;
    } else if (n_curves == 2) {
        int a = cur(0, 0), b = cur(1, 1);
        if (b == 0) {
            res.F_n = std::abs(a);
            res.surface = "\\mathbb{F}_{" + std::to_string(res.F_n) + "}";
        } else if (a == 0) {
            res.F_n = std::abs(b);
            res.surface = "\\mathbb{F}_{" + std::to_string(res.F_n) + "}";
        } else {
            res.surface = "\\text{UNKNOWN}_{2\\times 2}";
        }
    } else if (n_curves == 0) {
        res.surface = "\\text{pt}";
    } else {
        res.surface = "\\text{UNKNOWN}_{" + std::to_string(n_curves) + "}";
    }
    
    return res;
}

// ============================================================================
// IF → Quiver LaTeX: proper graph-walk based
// ============================================================================
//
// Walk adjacency graph from an endpoint, produce chain notation.
// Branches shown with \overset{branch}{node}.
// The SUGRA external (last curve index) is shown in red.

std::string if_to_latex(const Eigen::MatrixXi& IF, int ext_idx = -1) {
    int n = IF.rows();
    if (n == 0) return "\\emptyset";
    if (n == 1) return std::to_string(-IF(0,0));
    
    if (ext_idx == -1) ext_idx = n - 1;  // default: last curve is external
    // ext_idx == -2 means no external (for minimal IF display)
    
    // Build adjacency list
    std::vector<std::vector<std::pair<int,int>>> adj(n); // adj[i] = {(j, int_num), ...}
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            if (IF(i,j) != 0) {
                adj[i].push_back({j, IF(i,j)});
                adj[j].push_back({i, IF(j,i)});
            }
    
    // Find endpoint: node with degree 1, prefer not the external
    int start = -1;
    for (int i = 0; i < n; i++)
        if ((int)adj[i].size() == 1 && i != ext_idx) { start = i; break; }
    if (start < 0) {
        for (int i = 0; i < n; i++)
            if ((int)adj[i].size() == 1) { start = i; break; }
    }
    if (start < 0) start = 0;
    
    // DFS walk: main chain + side branches
    std::vector<bool> visited(n, false);
    
    std::function<std::string(int, int)> walk = [&](int cur, int from) -> std::string {
        visited[cur] = true;
        int si = -IF(cur, cur);
        std::string si_str = std::to_string(si);
        // Wrap multi-digit self-intersections in parentheses for readability
        if (si >= 10) si_str = "(" + si_str + ")";
        std::string label = si_str;
        
        // Color external red
        if (cur == ext_idx)
            label = "\\textcolor{red}{" + si_str + "}";
        
        // Find unvisited neighbors
        std::vector<int> children;
        for (auto& [nb, intnum] : adj[cur])
            if (!visited[nb]) children.push_back(nb);
        
        if (children.empty()) {
            return label;
        } else if (children.size() == 1) {
            // Simple chain continuation
            return label + walk(children[0], cur);
        } else {
            // Branch point: main chain = longest path, others = overset branches
            // Pick the non-external child as main if possible
            int main_child = -1;
            std::vector<int> branch_children;
            for (int c : children) {
                if (main_child < 0 && c != ext_idx) main_child = c;
                else branch_children.push_back(c);
            }
            if (main_child < 0) { main_child = children[0]; branch_children.clear(); for (size_t i=1;i<children.size();i++) branch_children.push_back(children[i]); }
            
            // Build branch strings: collect chain, reverse (tip→branch point)
            std::string branches;
            for (int i = (int)branch_children.size() - 1; i >= 0; i--) {
                int bc = branch_children[i];
                // Collect linear chain from branch point outward
                std::vector<std::pair<int,bool>> chain; // {self_int, is_ext}
                int node = bc;
                while (node >= 0) {
                    visited[node] = true;
                    chain.push_back({-IF(node, node), node == ext_idx});
                    int next = -1;
                    for (auto& [nb, intnum] : adj[node])
                        if (!visited[nb]) { next = nb; break; }
                    node = next;
                }
                // Reverse: tip first
                std::reverse(chain.begin(), chain.end());
                std::string br;
                for (auto& [s, is_ext] : chain) {
                    std::string s_str = std::to_string(s);
                    if (s >= 10) s_str = "(" + s_str + ")";
                    if (is_ext)
                        br += "\\textcolor{red}{" + s_str + "}";
                    else
                        br += s_str;
                }
                if (!branches.empty()) branches += ",";
                branches += br;
            }
            
            std::string main_part = walk(main_child, cur);
            
            if (!branches.empty())
                return "\\overset{" + branches + "}{" + label + "}" + main_part;
            else
                return label + main_part;
        }
    };
    
    return walk(start, -1);
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv) {
    // Defaults
    std::string catalog_file = "unified.cat";
    int T_min = 0, T_max = 20;
    bool det_check = true;
    
    if (argc > 1) catalog_file = argv[1];
    if (argc > 2) T_max = std::atoi(argv[2]);
    if (argc > 3) T_min = std::atoi(argv[3]);
    
    std::string suffix = "_T" + std::to_string(T_min) + "_" + std::to_string(T_max);
    
    auto catalog = load_catalog(catalog_file);
    std::cout << "Catalog: " << catalog.size() << " entries\n";
    
    SUGRAConfig config;
    config.T_max = 193;
    config.check_nhc = true;
    config.check_anomaly = true;       // H_neutral >= 0
    config.check_determinant = false;
    config.required_det = 0;
    config.catalog_T_min = T_min;
    config.catalog_T_max = T_max;
    config.verbose = true;
    
    auto results = generate_sugra(catalog, config);
    
    // No det filter — keep all
    std::vector<SUGRAResult>& all_results = results;
    
    std::cout << "Total bases (before dedup): " << all_results.size() << "\n";
    
    // Deduplicate by structural invariants (not just eigenvalues)
    auto dedup_key = [](const Eigen::MatrixXi& IF) -> std::string {
        int n = IF.rows();
        // 1. Sorted eigenvalues
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
        auto ev = es.eigenvalues();
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(8);
        ss << "E:";
        for (int i = 0; i < ev.size(); i++) {
            double v = ev(i);
            if (std::abs(v) < 1e-10) v = 0.0;
            if (i) ss << ",";
            ss << v;
        }
        // 2. Sorted vertex descriptors: (self_int, degree, sorted_neighbor_self_ints)
        struct VDesc {
            int si, deg;
            std::vector<int> nbr_si;
            bool operator<(const VDesc& o) const {
                if (si != o.si) return si < o.si;
                if (deg != o.deg) return deg < o.deg;
                return nbr_si < o.nbr_si;
            }
        };
        std::vector<VDesc> vds(n);
        for (int i = 0; i < n; i++) {
            vds[i].si = IF(i,i);
            vds[i].deg = 0;
            for (int j = 0; j < n; j++) {
                if (j != i && IF(i,j) != 0) {
                    vds[i].deg++;
                    vds[i].nbr_si.push_back(IF(j,j));
                }
            }
            std::sort(vds[i].nbr_si.begin(), vds[i].nbr_si.end());
        }
        std::sort(vds.begin(), vds.end());
        ss << "|V:";
        for (auto& vd : vds) {
            ss << "(" << vd.si << "," << vd.deg << ",{";
            for (int k = 0; k < (int)vd.nbr_si.size(); k++) {
                if (k) ss << ",";
                ss << vd.nbr_si[k];
            }
            ss << "})";
        }
        // 3. Sorted edge multiset: (si_i, si_j, int_num) with si_i <= si_j
        std::vector<std::tuple<int,int,int>> edges;
        for (int i = 0; i < n; i++)
            for (int j = i+1; j < n; j++)
                if (IF(i,j) != 0) {
                    int a = IF(i,i), b = IF(j,j);
                    if (a > b) std::swap(a, b);
                    edges.push_back({a, b, IF(i,j)});
                }
        std::sort(edges.begin(), edges.end());
        ss << "|E:";
        for (auto& [a,b,w] : edges) ss << "(" << a << "," << b << "," << w << ")";
        return ss.str();
    };
    
    std::set<std::string> seen;
    std::vector<SUGRAResult> filtered;
    for (auto& r : all_results) {
        std::string key = dedup_key(r.final_IF);
        if (seen.insert(key).second)
            filtered.push_back(std::move(r));
    }
    
    std::cout << "After dedup: " << filtered.size() << "\n";
    
    // Sort by T, then H_neutral
    std::sort(filtered.begin(), filtered.end(), [](const SUGRAResult& a, const SUGRAResult& b){
        if (a.anomaly.T != b.anomaly.T) return a.anomaly.T < b.anomaly.T;
        if (a.anomaly.H_neutral != b.anomaly.H_neutral) return a.anomaly.H_neutral < b.anomaly.H_neutral;
        return a.catalog_id < b.catalog_id;
    });
    
    // Split into unimodular and non-unimodular
    std::vector<SUGRAResult> uni, nonuni;
    for (auto& r : filtered) {
        if (std::abs(r.sig.det) == 1)
            uni.push_back(r);
        else
            nonuni.push_back(r);
    }
    
    std::cout << "Unimodular (|det|=1): " << uni.size() << "\n";
    std::cout << "Non-unimodular: " << nonuni.size() << "\n";
    
    // ─── Helper: write LaTeX preamble ───
    auto write_preamble = [&](std::ofstream& tex) {
        tex << "\\documentclass[10pt]{article}\n";
        tex << "\\usepackage[a4paper,margin=0.6in]{geometry}\n";
        tex << "\\usepackage{longtable}\n";
        tex << "\\usepackage{booktabs}\n";
        tex << "\\usepackage{amsmath,amssymb}\n";
        tex << "\\usepackage{xcolor}\n";
        tex << "\\usepackage{array}\n";
        tex << "\\setlength{\\parindent}{0pt}\n";
        tex << "\\setlength{\\parskip}{2pt}\n";
        tex << "\\begin{document}\n\n";
    };
    
    // ─── Helper: spectrum key for grouping ───
    auto base_spec_key = [](const Eigen::MatrixXi& IF) -> std::string {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
        auto ev = es.eigenvalues();
        std::vector<double> vals(ev.data(), ev.data() + ev.size());
        std::sort(vals.begin(), vals.end());
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(6);
        for (size_t i = 0; i < vals.size(); i++) {
            if (i) ss << ",";
            ss << vals[i];
        }
        return ss.str();
    };
    
    // ─── Sort entries by T, then by base LST spectrum ───
    auto sort_entries = [&](std::vector<SUGRAResult>& v) {
        // Precompute base spectrum keys
        std::map<int, std::string> keys;
        for (int i = 0; i < (int)v.size(); i++)
            keys[i] = base_spec_key(v[i].base_IF);
        
        std::vector<int> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            if (v[a].anomaly.T != v[b].anomaly.T) return v[a].anomaly.T < v[b].anomaly.T;
            return keys[a] < keys[b];
        });
        std::vector<SUGRAResult> sorted;
        sorted.reserve(v.size());
        for (int i : idx) sorted.push_back(std::move(v[i]));
        v = std::move(sorted);
    };
    sort_entries(uni);
    sort_entries(nonuni);
    
    // ─── Unimodular: with blowdown ───
    {
        std::ofstream tex(("sugra_unimodular" + suffix + ".tex").c_str());
        write_preamble(tex);
        
        // === Section 1: Statistics ===
        tex << "\\section*{6D SUGRA Bases --- Unimodular ($|\\det|=1$)}\n\n";
        
        // Count distinct LST bases used
        std::set<std::string> uni_lst_bases;
        for (auto& r : uni) uni_lst_bases.insert(base_spec_key(r.base_IF));
        
        tex << "Total SUGRA bases: " << uni.size() 
            << " (from " << uni_lst_bases.size() << " distinct LST bases, "
            << "$T_{\\text{LST}} \\in [" << T_min << ", " << T_max << "]$).\n\n";
        tex << "Filters: NHC, $c(k) \\leq 8$, $H_n \\geq 0$. Deduplicated by graph invariants.\n\n";
        tex << "\\textcolor{red}{Red} = attached external curve. $\\scriptstyle [\\#=n]$ = intersection number $n$.\n";
        tex << "$\\Delta = H_{\\text{charged}} - V + 29T - 273 \\leq 0$.\n\n";
        
        // T distribution
        std::map<int,int> T_count;
        for (auto& r : uni) T_count[r.anomaly.T]++;
        
        // Surface distribution
        std::map<int, int> fn_count;
        int p2_count = 0, unk_count = 0;
        for (auto& r : uni) {
            auto bd = blowdown(r.final_IF);
            if (bd.F_n >= 0) fn_count[bd.F_n]++;
            else if (bd.is_P2) p2_count++;
            else unk_count++;
        }
        
        // Two tables centered together
        tex << "\\begin{center}\n";
        tex << "\\begin{tabular}{c|c}\n\\toprule\n$T$ & Count \\\\\n\\midrule\n";
        for (auto& [T, cnt] : T_count) tex << T << " & " << cnt << " \\\\\n";
        tex << "\\bottomrule\n\\end{tabular}\n\\qquad\\qquad\n";
        tex << "\\begin{tabular}{c|c}\n\\toprule\nMinimal Surface & Count \\\\\n\\midrule\n";
        if (p2_count) tex << "$\\mathbb{P}^2$ & " << p2_count << " \\\\\n";
        for (auto& [n, cnt] : fn_count) tex << "$\\mathbb{F}_{" << n << "}$ & " << cnt << " \\\\\n";
        if (unk_count) tex << "UNKNOWN & " << unk_count << " \\\\\n";
        tex << "\\bottomrule\n\\end{tabular}\n";
        tex << "\\end{center}\n\n";
        
        // === Section 2: Entries by T, grouped by base LST ===
        tex << "\\newpage\n";
        int cur_T = -1;
        std::string cur_base_key = "";
        for (auto& r : uni) {
            if (r.anomaly.T != cur_T) {
                cur_T = r.anomaly.T;
                cur_base_key = "";
                tex << "\n\\subsection*{$T = " << cur_T << "$ \\quad ("
                    << T_count[cur_T] << " bases)}\n";
            }
            std::string bk = base_spec_key(r.base_IF);
            if (bk != cur_base_key) {
                cur_base_key = bk;
                std::string base_latex = if_to_latex(r.base_IF, -2);  // -2 = no external highlight
                tex << "\\smallskip{\\small\\textbf{Base:} $" << base_latex << "$}\\smallskip\n\n";
            }
            std::string topo_latex = if_to_latex(r.final_IF);
            int delta = r.anomaly.H_charged - r.anomaly.V + 29 * r.anomaly.T - 273;
            int int_num = r.externals.empty() ? 1 : r.externals[0].int_num;
            auto bd = blowdown(r.final_IF);
            
            tex << "$" << topo_latex << "$";
            if (int_num > 1)
                tex << " $\\scriptstyle [\\#=" << int_num << "]$";
            tex << " \\quad $\\Delta="
                << delta << ",\\; n_+="
                << r.sig.sig_pos << ",\\; n_-="
                << r.sig.sig_neg << ",\\; n_0="
                << r.sig.sig_zero
                << "$ \\quad $\\to " << bd.surface << "$\n\n";
        }
        
        tex << "\\end{document}\n";
        tex.close();
        std::cout << "Written " << uni.size() << " entries to sugra_unimodular.tex\n";
    }
    
    // ─── Non-unimodular: split by ext type, show T_max ───
    {
        struct NonuniEntry {
            const SUGRAResult* r;
            TmaxResult tm;
        };
        
        // Group by ext type
        std::vector<NonuniEntry> ext3, ext24, ext_rest;
        for (auto& r : nonuni) {
            auto tm = compute_Tmax(r.final_IF);
            NonuniEntry e{&r, tm};
            int ext_si = r.externals.empty() ? 0 : r.externals[0].ext_si;
            if (ext_si == -3) ext3.push_back(e);
            else if (ext_si == -2 || ext_si == -4) ext24.push_back(e);
            else ext_rest.push_back(e);
        }
        
        std::cout << "Non-unimodular split: ext=-3: " << ext3.size()
                  << ", ext=-2/-4: " << ext24.size()
                  << ", ext=-5/-6/-7/-8/-12: " << ext_rest.size() << "\n";
        
        // Helper lambda to write a nonunimodular file
        auto write_nonuni_file = [&](const std::string& filename,
                                     const std::string& title,
                                     std::vector<NonuniEntry>& entries) {
            std::ofstream tex(filename);
            write_preamble(tex);
            tex << "\\section*{" << title << "}\n\n";
            
            std::set<std::string> lst_bases;
            for (auto& e : entries) lst_bases.insert(base_spec_key(e.r->base_IF));
            
            // Count exact vs frac
            int n_exact = 0, n_frac = 0;
            for (auto& e : entries) { if (e.tm.exact) n_exact++; else n_frac++; }
            
            tex << "Total: " << entries.size()
                << " bases (from " << lst_bases.size() << " distinct LST bases, "
                << "$T_{\\text{LST}} \\in [" << T_min << ", " << T_max << "]$).\n\n"
                << "Exact $\\det=0$: " << n_exact
                << ", fractional $T_{\\text{crit}}$: " << n_frac << ".\n\n";
            tex << "\\textcolor{red}{Red} = attached external curve. "
                << "$\\scriptstyle [\\#=n]$ = intersection number.\n\n";
            
            // Statistics tables
            std::map<int,int> T_count, Tmax_count;
            std::map<int,int> intnum_exact, intnum_frac, det_exact, det_frac;
            for (auto& e : entries) {
                T_count[e.r->anomaly.T]++;
                Tmax_count[e.tm.T_max]++;
                int num = e.r->externals.empty() ? 1 : e.r->externals[0].int_num;
                int da = std::abs(e.r->sig.det);
                if (e.tm.exact) { intnum_exact[num]++; det_exact[da]++; }
                else { intnum_frac[num]++; det_frac[da]++; }
            }
            
            // T and T_max tables
            tex << "\\begin{center}\n";
            tex << "\\begin{tabular}{c|c}\n\\toprule\n$T$ & Count \\\\\n\\midrule\n";
            for (auto& [T, cnt] : T_count) tex << T << " & " << cnt << " \\\\\n";
            tex << "\\bottomrule\n\\end{tabular}\n\\qquad\\qquad\n";
            tex << "\\begin{tabular}{c|c}\n\\toprule\n$T_{\\max}$ & Count \\\\\n\\midrule\n";
            for (auto& [T, cnt] : Tmax_count) tex << T << " & " << cnt << " \\\\\n";
            tex << "\\bottomrule\n\\end{tabular}\n";
            tex << "\\end{center}\n\n";
            
            // Intersection number table
            std::set<int> all_nums;
            for (auto& [n,_] : intnum_exact) all_nums.insert(n);
            for (auto& [n,_] : intnum_frac) all_nums.insert(n);
            tex << "\\begin{center}\n";
            tex << "\\begin{tabular}{c|r r r}\n\\toprule\n$\\#$ & Exact & Frac & Total \\\\\n\\midrule\n";
            for (int n : all_nums) {
                int ex = intnum_exact.count(n) ? intnum_exact[n] : 0;
                int fr = intnum_frac.count(n) ? intnum_frac[n] : 0;
                tex << n << " & " << ex << " & " << fr << " & " << ex + fr << " \\\\\n";
            }
            tex << "\\bottomrule\n\\end{tabular}\n\\qquad\\qquad\n";
            
            // |det| table
            std::set<int> all_dets;
            for (auto& [d,_] : det_exact) all_dets.insert(d);
            for (auto& [d,_] : det_frac) all_dets.insert(d);
            tex << "\\begin{tabular}{c|r r r}\n\\toprule\n$|\\det|$ & Exact & Frac & Total \\\\\n\\midrule\n";
            for (int d : all_dets) {
                int ex = det_exact.count(d) ? det_exact[d] : 0;
                int fr = det_frac.count(d) ? det_frac[d] : 0;
                tex << d << " & " << ex << " & " << fr << " & " << ex + fr << " \\\\\n";
            }
            tex << "\\bottomrule\n\\end{tabular}\n";
            tex << "\\end{center}\n\n";
            
            // Entries
            tex << "\\newpage\n";
            int cur_T = -1;
            std::string cur_base_key = "";
            for (auto& e : entries) {
                auto& r = *e.r;
                if (r.anomaly.T != cur_T) {
                    cur_T = r.anomaly.T;
                    cur_base_key = "";
                    tex << "\n\\subsection*{$T = " << cur_T << "$ \\quad ("
                        << T_count[cur_T] << " bases)}\n";
                }
                std::string bk = base_spec_key(r.base_IF);
                if (bk != cur_base_key) {
                    cur_base_key = bk;
                    tex << "\\smallskip{\\small\\textbf{Base:} $" << if_to_latex(r.base_IF, -2) << "$}\\smallskip\n\n";
                }
                int delta = r.anomaly.H_charged - r.anomaly.V + 29 * r.anomaly.T - 273;
                int int_num = r.externals.empty() ? 1 : r.externals[0].int_num;
                tex << "$" << if_to_latex(r.final_IF) << "$";
                if (int_num > 1) tex << " $\\scriptstyle [\\#=" << int_num << "]$";
                tex << " \\quad $\\Delta=" << delta
                    << ",\\; |\\det|=" << std::abs(r.sig.det)
                    << ",\\; (n_+,n_-,n_0)=(" << r.sig.sig_pos << "," << r.sig.sig_neg << "," << r.sig.sig_zero << ")";
                if (e.tm.exact) {
                    tex << ",\\; T_{\\max}=" << e.tm.T_max
                        << ",\\; T_{\\text{crit}}=" << e.tm.T_max + 1
                        << "$ \\textcolor{blue}{(exact)}";
                } else {
                    tex << ",\\; T_{\\text{crit}}=" << std::fixed << std::setprecision(2) << e.tm.T_crit
                        << ",\\; T_{\\max}=" << e.tm.T_max << "$";
                }
                tex << "\n\n";
            }
            tex << "\\end{document}\n";
            tex.close();
            std::cout << "Written " << entries.size() << " entries to " << filename << "\n";
        };
        
        write_nonuni_file(("sugra_nonuni_ext3" + suffix + ".tex").c_str(),
            "Non-unimodular SUGRA Bases: ext $= -3$", ext3);
        write_nonuni_file(("sugra_nonuni_ext24" + suffix + ".tex").c_str(),
            "Non-unimodular SUGRA Bases: ext $= -2, -4$", ext24);
        write_nonuni_file(("sugra_nonuni_ext_rest" + suffix + ".tex").c_str(),
            "Non-unimodular SUGRA Bases: ext $= -5, -6, -7, -8, -12$", ext_rest);
    }
    
    // Print summary
    print_sugra_summary(filtered);
    
    return 0;
}
