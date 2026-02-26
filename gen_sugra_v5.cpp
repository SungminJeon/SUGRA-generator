// gen_sugra_v5.cpp
// SUGRA base generator: \hat{1} external curve attachment
//
// \hat{1} external properties:
//   C^2 = -1  (self-intersection)
//   b_0 · C = -1  (NOT the standard C^2 + 2 = 1)
//   Attaches to: -1 curves (int_num=1), -2 curves (int_num=1)
//
// Key differences from v1-v3:
//   - \hat{1} has ext_si = -1 (same as internal -1 curves in IF)
//   - Must track which -1 curves are \hat{1} externals for b_0Q construction
//   - Modified b_0Q matrix: b_0 · C_{\hat{1}} = -1 instead of C^2+2 = 1
//   - LaTeX: \hat{1} shown as $\hat{1}$ in quiver notation
//
// Compile:
//   g++ -std=c++17 -O2 -o gen_sugra_v5 gen_sugra_v5.cpp Tensor.C \
//       TopoLineCompact_enhanced.cpp Topology_enhanced.cpp \
//       -I/usr/include/eigen3 -I.
//
// Usage:
//   ./gen_sugra_v5 unified.cat [T_max=193] [T_min=0]

#include "sugra_generator.h"
#include <fstream>
#include <cstdlib>
#include <set>
#include <iomanip>
#include <numeric>
#include <map>
#include <algorithm>

// ============================================================================
// Block-level topology info (reused from v2/v3)
// ============================================================================

enum class BlockType { Node, InteriorLink, SideLink, Instanton, External, Unknown };

struct BlockInfo {
    BlockType type;
    int param;
    int start, end;
};

struct TopoBlockMap {
    std::vector<BlockInfo> blocks;
    std::vector<int> curve_to_block;

    bool is_node(int c) const {
        if (c < 0 || c >= (int)curve_to_block.size()) return false;
        return blocks[curve_to_block[c]].type == BlockType::Node;
    }
    bool is_interior_link(int c) const {
        if (c < 0 || c >= (int)curve_to_block.size()) return false;
        return blocks[curve_to_block[c]].type == BlockType::InteriorLink;
    }
    int block_id(int c) const {
        if (c < 0 || c >= (int)curve_to_block.size()) return -1;
        return curve_to_block[c];
    }
    const BlockInfo& block_of(int c) const { return blocks[curve_to_block[c]]; }
};

TopoBlockMap get_block_map(const CatalogEntry& entry) {
    TopoBlockMap result;
    if (entry.is_nn()) return result;
    Topology_enhanced topo;
    if (!TopoLineCompact_enhanced::deserialize(entry.topo, topo)) return result;
    int offset = 0;
    auto add = [&](BlockType t, int p, Kind k, int sp) {
        int nc = getCurveCount(sp, k);
        result.blocks.push_back({t, p, offset, offset + nc});
        offset += nc;
    };
    for (auto& b : topo.block) {
        Spec sp = nk_detail::make_spec(b.kind, b.param);
        add(b.kind == LKind::g ? BlockType::Node : BlockType::InteriorLink,
            b.param, sp.kind, sp.param);
    }
    for (auto& s : topo.side_links)
        add(BlockType::SideLink, s.param, Kind::SideLink, s.param);
    for (auto& i : topo.instantons)
        add(BlockType::Instanton, i.param, Kind::SideLink, i.param);
    for (auto& e : topo.externals)
        add(BlockType::External, e.param, Kind::External, e.param);
    result.curve_to_block.resize(offset, -1);
    for (int b = 0; b < (int)result.blocks.size(); b++)
        for (int c = result.blocks[b].start; c < result.blocks[b].end; c++)
            result.curve_to_block[c] = b;
    return result;
}

// ============================================================================
// \hat{1} SUGRA result: extends SUGRAResult with hat1 tracking
// ============================================================================

struct Hat1Result {
    int catalog_id;
    std::string catalog_type;
    Eigen::MatrixXi base_IF;
    Eigen::MatrixXi final_IF;

    // \hat{1} attachment info
    int hat1_idx;                    // index of \hat{1} curve in final_IF
    int target_idx;                  // primary target curve in base_IF (first target)
    int target_si;                   // self-intersection of target(s): -1 or -2
    int int_num;                     // intersection number per target (typically 1)
    std::vector<int> all_targets;    // all target indices in base_IF (size≥1)

    NHCResult nhc;
    AnomalyResult anomaly;
    SigInfo sig;
    bool valid;
    
    int num_targets() const { return (int)all_targets.size(); }
};

// ============================================================================
// Modified b_0Q matrix: accounts for \hat{1} having b_0·C = -1
// ============================================================================

inline Eigen::MatrixXi build_b0Q_hat1(const Eigen::MatrixXi& IF, int hat1_idx) {
    int n = IF.rows();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
    int sig_neg = 0;
    for (int i = 0; i < es.eigenvalues().size(); i++)
        if (es.eigenvalues()(i) < -1e-8) sig_neg++;

    Eigen::MatrixXi m = Eigen::MatrixXi::Zero(n + 1, n + 1);
    m.block(0, 0, n, n) = IF;
    for (int i = 0; i < n; i++) {
        if (i == hat1_idx) {
            // \hat{1}: b_0 · C = -1 (NOT the standard C^2 + 2 = 1)
            m(i, n) = -1;
            m(n, i) = -1;
        } else {
            // Standard: b_0 · C_i = C_i^2 + 2
            m(i, n) = IF(i, i) + 2;
            m(n, i) = IF(i, i) + 2;
        }
    }
    m(n, n) = 9 - sig_neg;
    return m;
}

// ============================================================================
// Blowdown with \hat{1} awareness
// ============================================================================

inline Eigen::MatrixXi blowdown_curve(const Eigen::MatrixXi& IF, int idx) {
    int n = IF.rows();
    Eigen::MatrixXi M = IF;
    for (int j = 0; j < n; j++) {
        if (j == idx) continue;
        for (int k = 0; k < n; k++) {
            if (k == idx) continue;
            M(j, k) += IF(j, idx) * IF(k, idx);
        }
    }
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

struct BlowdownResult {
    Eigen::MatrixXi minimal_IF;
    Eigen::MatrixXi extended_IF;
    int num_blowdowns;
    std::string surface;
    int F_n;
    bool is_P2 = false;
    int b0_sq = 0;
};

BlowdownResult blowdown_hat1(const Eigen::MatrixXi& IF, int hat1_idx) {
    Eigen::MatrixXi cur = build_b0Q_hat1(IF, hat1_idx);
    int bd = 0;
    // Track hat1 index as we remove curves
    int cur_hat1 = hat1_idx;

    while (true) {
        int n = cur.rows();
        int m1 = -1;
        for (int i = 0; i < n - 1; i++) {
            if (cur(i, i) == -1) { m1 = i; break; }
        }
        if (m1 < 0) break;
        cur = blowdown_curve(cur, m1);
        // Update hat1 index
        if (m1 == cur_hat1) cur_hat1 = -1;  // hat1 was blown down
        else if (m1 < cur_hat1) cur_hat1--;
        bd++;
    }

    BlowdownResult res;
    int n = cur.rows();
    int n_curves = n - 1;
    res.num_blowdowns = bd;
    res.F_n = -1;
    res.b0_sq = cur(n - 1, n - 1);
    res.extended_IF = cur;
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
            res.surface = "\\text{UNK}_{2}";
        }
    } else if (n_curves == 0) {
        res.surface = "\\text{pt}";
    } else {
        res.surface = "\\text{UNK}_{" + std::to_string(n_curves) + "}";
    }
    return res;
}

// ============================================================================
// Extended matrix → LaTeX (with b₀ label on last row/col)
// ============================================================================

// Compact: small pmatrix with b₀ row/col
std::string ext_matrix_to_latex(const Eigen::MatrixXi& M) {
    int n = M.rows();
    if (n == 0) return "()";
    // Last index = b₀
    // Format as a small matrix with row/col labels
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

// More readable: show as labeled matrix with b₀
// Rows: C_0, C_1, ..., C_{k-1}, b_0
std::string ext_matrix_labeled_latex(const Eigen::MatrixXi& M, int n_blowdowns) {
    int n = M.rows();
    if (n == 0) return "()";
    int nc = n - 1;  // number of remaining curves
    
    if (nc == 0) {
        // Only b₀ left: just show b₀² 
        return "b_0^2=" + std::to_string(M(0, 0));
    }
    
    // For small matrices (≤ 4×4), show full matrix
    if (n <= 5) {
        return ext_matrix_to_latex(M);
    }
    
    // For larger, just show diagonal + b₀ row
    std::string s = "\\text{diag}=(";
    for (int i = 0; i < nc; i++) {
        if (i) s += ",";
        s += std::to_string(M(i, i));
    }
    s += "),\\; b_0^2=" + std::to_string(M(nc, nc));
    s += ",\\; b_0{\\cdot}C=(";
    for (int i = 0; i < nc; i++) {
        if (i) s += ",";
        s += std::to_string(M(i, nc));
    }
    s += ")";
    return s;
}

// ============================================================================
// T_crit computation for non-unimodular with \hat{1}
// ============================================================================

struct CritInfo {
    int T_sugra;
    double T_crit;
    int T_crit_int;
    bool exact;
    int crit_npos, crit_nneg, crit_nzero;
    BlowdownResult bd_crit;
};

inline CritInfo compute_crit_hat1(const Eigen::MatrixXi& IF, int hat1_idx) {
    CritInfo res;
    int n = IF.rows();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es0(IF.cast<double>());
    res.T_sugra = 0;
    for (int i = 0; i < es0.eigenvalues().size(); i++)
        if (es0.eigenvalues()(i) < -1e-8) res.T_sugra++;

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n + 1, n + 1);
    M.block(0, 0, n, n) = IF.cast<double>();
    for (int i = 0; i < n; i++) {
        double b0C = (i == hat1_idx) ? -1.0 : (IF(i, i) + 2.0);
        M(i, n) = b0C;
        M(n, i) = b0C;
    }

    M(n, n) = 0.0;
    double B = M.determinant();
    M(n, n) = 1.0;
    double A = M.determinant() - B;

    if (std::abs(A) < 1e-10) {
        res.T_crit = 1e9; res.exact = false; res.T_crit_int = -1;
        res.crit_npos = res.crit_nneg = res.crit_nzero = 0;
        return res;
    }

    double b0sq_crit = -B / A;
    res.T_crit = 9.0 - b0sq_crit;
    int tc = (int)std::round(res.T_crit);
    res.exact = (std::abs(res.T_crit - tc) < 1e-6);
    if (res.exact) res.T_crit_int = tc;
    else res.T_crit_int = (int)std::ceil(res.T_crit - 1e-9);

    // Signature at T_crit_int
    int b0sq = 9 - res.T_crit_int;
    M(n, n) = (double)b0sq;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es1(M);
    res.crit_npos = res.crit_nneg = res.crit_nzero = 0;
    for (int i = 0; i < n + 1; i++) {
        double ev = es1.eigenvalues()(i);
        if (std::abs(ev) < 1e-8) res.crit_nzero++;
        else if (ev > 0) res.crit_npos++;
        else res.crit_nneg++;
    }

    // Blowdown extended matrix at T_crit_int (hat1 is NOT blowable: b0·hat1 = -1 ≠ hat1²+2 = 1)
    Eigen::MatrixXi ext = Eigen::MatrixXi::Zero(n + 1, n + 1);
    ext.block(0, 0, n, n) = IF;
    for (int i = 0; i < n; i++) {
        int b0C = (i == hat1_idx) ? -1 : (IF(i, i) + 2);
        ext(i, n) = b0C;
        ext(n, i) = b0C;
    }
    ext(n, n) = b0sq;

    Eigen::MatrixXi cur = ext;
    int bd_count = 0;
    int cur_hat1 = hat1_idx;
    while (true) {
        int nn = cur.rows();
        int m1 = -1;
        for (int i = 0; i < nn - 1; i++) {
            if (i == cur_hat1) continue;  // don't blow down hat1
            if (cur(i, i) == -1) { m1 = i; break; }
        }
        if (m1 < 0) break;
        cur = blowdown_curve(cur, m1);
        if (m1 < cur_hat1) cur_hat1--;
        bd_count++;
    }

    int nc = cur.rows() - 1;
    res.bd_crit.num_blowdowns = bd_count;
    res.bd_crit.extended_IF = cur;
    res.bd_crit.minimal_IF = (nc > 0) ? cur.block(0, 0, nc, nc) : Eigen::MatrixXi();
    res.bd_crit.b0_sq = cur(cur.rows() - 1, cur.cols() - 1);
    res.bd_crit.F_n = -1;
    res.bd_crit.is_P2 = false;
    // Surface identification (skip hat1 for counting)
    int real_curves = nc - 1; // exclude hat1
    if (real_curves == 1 && nc == 2) {
        // One real curve + hat1
        res.bd_crit.surface = "\\text{bd}_{" + std::to_string(bd_count) + "}";
    } else if (real_curves == 0 && nc == 1) {
        res.bd_crit.surface = "\\hat{1}\\text{-only}";
    } else {
        res.bd_crit.surface = "\\text{bd}_{" + std::to_string(bd_count) + "}";
    }
    return res;
}

// ============================================================================
// Gauge algebra LaTeX (from v2)
// ============================================================================

std::string gauge_algebra_latex(int param) {
    switch (param) {
        case 5:  return "\\mathfrak{f}_4";
        case 6:  return "\\mathfrak{e}_6";
        case 7:  return "\\mathfrak{e}_7'";
        case 8:  return "\\mathfrak{e}_7";
        case 12: return "\\mathfrak{e}_8";
        default: return std::to_string(param);
    }
}

std::string interior_link_label(int param) {
    std::string s = std::to_string(param);
    if (s == "331") return "\\overset{3,3}{\\bigcirc}";
    if (s.length() >= 2) {
        std::string digits;
        for (size_t i = 0; i < s.length(); i++) {
            if (i) digits += ",";
            digits += s[i];
        }
        return "\\overset{" + digits + "}{\\otimes}";
    }
    return "\\overset{" + s + "}{\\otimes}";
}

// ============================================================================
// IF → Quiver LaTeX with \hat{1} awareness
// ============================================================================
// \hat{1} curve shown as $\hat{1}$ in blue; other externals would be red

std::string if_to_latex_hat1(const Eigen::MatrixXi& IF, int hat1_idx,
                              const TopoBlockMap& bmap = {}) {
    int n = IF.rows();
    if (n == 0) return "\\emptyset";
    if (n == 1) return std::to_string(-IF(0, 0));

    bool has_blocks = !bmap.blocks.empty();

    // Build adjacency
    std::vector<std::vector<std::pair<int, int>>> adj(n);
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (IF(i, j) != 0) {
                adj[i].push_back({j, IF(i, j)});
                adj[j].push_back({i, IF(j, i)});
            }

    // Start from degree-1 node, prefer not \hat{1}
    int start = -1;
    for (int i = 0; i < n; i++)
        if ((int)adj[i].size() == 1 && i != hat1_idx) { start = i; break; }
    if (start < 0)
        for (int i = 0; i < n; i++)
            if ((int)adj[i].size() == 1) { start = i; break; }
    if (start < 0) start = 0;

    std::vector<bool> visited(n, false);

    // Label for a curve
    auto curve_label = [&](int cur) -> std::string {
        int si = -IF(cur, cur);
        std::string si_str = std::to_string(si);
        if (si >= 10) si_str = "(" + si_str + ")";
        if (cur == hat1_idx) return "\\textcolor{red}{\\hat{1}}";
        return si_str;
    };

    auto skip_block = [&](int cur) -> std::vector<int> {
        auto& blk = bmap.block_of(cur);
        for (int c = blk.start; c < blk.end; c++) visited[c] = true;
        std::vector<int> exits;
        for (int c = blk.start; c < blk.end; c++)
            for (auto& [nb, w] : adj[c])
                if (!visited[nb]) exits.push_back(nb);
        std::sort(exits.begin(), exits.end());
        exits.erase(std::unique(exits.begin(), exits.end()), exits.end());
        return exits;
    };

    std::function<std::string(int, int)> walk = [&](int cur, int from) -> std::string {
        visited[cur] = true;

        // Interior link compression
        if (has_blocks && bmap.is_interior_link(cur)
            && bmap.block_of(cur).end - bmap.block_of(cur).start > 1) {
            std::string label = interior_link_label(bmap.block_of(cur).param);
            auto exits = skip_block(cur);
            if (exits.empty()) return label;
            if (exits.size() == 1) return label + walk(exits[0], cur);
            int main_child = -1;
            std::vector<int> branch_children;
            for (int c : exits) {
                if (main_child < 0 && c != hat1_idx) main_child = c;
                else branch_children.push_back(c);
            }
            if (main_child < 0) {
                main_child = exits[0];
                branch_children.assign(exits.begin() + 1, exits.end());
            }
            std::string branches;
            for (int i = (int)branch_children.size() - 1; i >= 0; i--) {
                if (!branches.empty()) branches += ",";
                branches += walk(branch_children[i], cur);
            }
            std::string main_part = walk(main_child, cur);
            if (!branches.empty())
                return "\\overset{" + branches + "}{" + label + "}" + main_part;
            return label + main_part;
        }

        // Node / standard / \hat{1}
        std::string label;
        if (has_blocks && cur < (int)bmap.curve_to_block.size() && bmap.is_node(cur))
            label = gauge_algebra_latex(bmap.block_of(cur).param);
        else
            label = curve_label(cur);

        std::vector<int> children;
        for (auto& [nb, w] : adj[cur])
            if (!visited[nb]) children.push_back(nb);

        if (children.empty()) return label;
        if (children.size() == 1) return label + walk(children[0], cur);

        // Branch point
        int main_child = -1;
        std::vector<int> branch_children;
        for (int c : children) {
            if (main_child < 0 && c != hat1_idx) main_child = c;
            else branch_children.push_back(c);
        }
        if (main_child < 0) {
            main_child = children[0];
            branch_children.clear();
            for (size_t i = 1; i < children.size(); i++)
                branch_children.push_back(children[i]);
        }

        std::string branches;
        for (int i = (int)branch_children.size() - 1; i >= 0; i--) {
            int bc = branch_children[i];
            // Walk branch as chain
            std::vector<std::string> chain;
            int nd = bc;
            while (nd >= 0) {
                visited[nd] = true;
                chain.push_back(curve_label(nd));
                int next = -1;
                for (auto& [nb, w] : adj[nd])
                    if (!visited[nb]) { next = nb; break; }
                nd = next;
            }
            std::reverse(chain.begin(), chain.end());
            std::string br;
            for (auto& lbl : chain) br += lbl;
            if (!branches.empty()) branches += ",";
            branches += br;
        }

        std::string main_part = walk(main_child, cur);
        if (!branches.empty())
            return "\\overset{" + branches + "}{" + label + "}" + main_part;
        return label + main_part;
    };

    return walk(start, -1);
}

// Multi-target rendering: base quiver with target -1 curves marked + hat1 annotation
std::string if_to_latex_hat1_multi(const Eigen::MatrixXi& base_IF,
                                    const std::vector<int>& targets,
                                    const TopoBlockMap& bmap = {}) {
    int n = base_IF.rows();
    if (n == 0) return "\\emptyset";

    std::set<int> target_set(targets.begin(), targets.end());
    bool has_blocks = !bmap.blocks.empty();

    std::vector<std::vector<std::pair<int, int>>> adj(n);
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (base_IF(i, j) != 0) {
                adj[i].push_back({j, base_IF(i, j)});
                adj[j].push_back({i, base_IF(j, i)});
            }

    int start = 0;
    for (int i = 0; i < n; i++)
        if ((int)adj[i].size() == 1) { start = i; break; }

    std::vector<bool> visited(n, false);

    auto curve_label = [&](int cur) -> std::string {
        int si = -base_IF(cur, cur);
        std::string si_str = std::to_string(si);
        if (si >= 10) si_str = "(" + si_str + ")";
        if (target_set.count(cur))
            return "\\textcolor{blue}{" + si_str + "}";
        return si_str;
    };

    std::function<std::string(int, int)> walk = [&](int cur, int from) -> std::string {
        visited[cur] = true;
        std::string label;
        if (has_blocks && cur < (int)bmap.curve_to_block.size() && bmap.is_node(cur))
            label = gauge_algebra_latex(bmap.block_of(cur).param);
        else
            label = curve_label(cur);

        std::vector<int> children;
        for (auto& [nb, w] : adj[cur])
            if (!visited[nb]) children.push_back(nb);

        if (children.empty()) return label;
        if (children.size() == 1) return label + walk(children[0], cur);

        int main_child = children[0];
        std::vector<int> branch_children(children.begin() + 1, children.end());
        std::string branches;
        for (int i = (int)branch_children.size() - 1; i >= 0; i--) {
            if (!branches.empty()) branches += ",";
            branches += walk(branch_children[i], cur);
        }
        std::string main_part = walk(main_child, cur);
        if (!branches.empty())
            return "\\overset{" + branches + "}{" + label + "}" + main_part;
        return label + main_part;
    };

    std::string quiver = walk(start, -1);
    return quiver + "+\\textcolor{red}{\\hat{1}}";
}

// Dispatcher: picks single or multi rendering
std::string render_hat1_latex(const Hat1Result& r, const TopoBlockMap& bmap = {}) {
    if (r.num_targets() <= 1) {
        return if_to_latex_hat1(r.final_IF, r.hat1_idx, bmap);
    } else {
        return if_to_latex_hat1_multi(r.base_IF, r.all_targets, bmap);
    }
}

// ============================================================================
// \hat{1} generation: attach C^2=-1, b_0C=-1 to -1 and -2 curves
// ============================================================================

// Attach \hat{1} to multiple target curves simultaneously
inline Eigen::MatrixXi attach_hat1_multi(const Eigen::MatrixXi& IF,
                                          const std::vector<int>& targets) {
    int n = IF.rows();
    Eigen::MatrixXi nIF = Eigen::MatrixXi::Zero(n + 1, n + 1);
    nIF.block(0, 0, n, n) = IF;
    nIF(n, n) = -1;  // hat1 self-intersection
    for (int t : targets) {
        nIF(n, t) = 1;
        nIF(t, n) = 1;
    }
    return nIF;
}

std::vector<Hat1Result> generate_hat1(
    const std::vector<CatalogEntry>& catalog,
    const SUGRAConfig& config)
{
    std::vector<Hat1Result> results;

    int processed = 0;
    for (const auto& entry : catalog) {
        if (entry.is_nn() && !config.include_nn) continue;
        if (!entry.is_nn() && !config.include_nk) continue;
        if (entry.T < config.catalog_T_min || entry.T > config.catalog_T_max) continue;

        Eigen::MatrixXi base_IF = reconstruct_IF(entry);
        if (base_IF.rows() == 0) continue;

        auto base_sig = compute_sig(base_IF);
        if (base_sig.sig_neg > config.T_max) continue;

        int n = base_IF.rows();

        // Collect -1 and -2 curves
        std::vector<int> m1_curves, m2_curves;
        for (int i = 0; i < n; i++) {
            if (base_IF(i, i) == -1) m1_curves.push_back(i);
            else if (base_IF(i, i) == -2) m2_curves.push_back(i);
        }

        // (A) Single attachment to -2 curves
        for (int tidx : m2_curves) {
            Eigen::MatrixXi new_IF = attach_curve(base_IF, tidx, -1, 1);
            int hat1_idx = new_IF.rows() - 1;

            auto new_sig = compute_sig(new_IF);
            if (new_sig.sig_neg > config.T_max) continue;

            NHCResult nhc = {true, "", {}, {}, {}};
            if (config.check_nhc) {
                nhc = check_nhc(new_IF, config.cc_budget);
                if (!nhc.passes) continue;
            }
            enhance_hat1_gauge(nhc, new_IF, hat1_idx, -2, config.cc_budget);
            if (!nhc.passes) continue;

            AnomalyResult anom = {0, 0, 0, 0};
            if (nhc.passes) {
                anom = compute_anomaly(new_IF, nhc);
                if (config.check_anomaly && anom.H_neutral < 0) continue;
            }

            Hat1Result r;
            r.catalog_id = entry.id;
            r.catalog_type = entry.type;
            r.base_IF = base_IF;
            r.final_IF = new_IF;
            r.hat1_idx = hat1_idx;
            r.target_idx = tidx;
            r.target_si = -2;
            r.int_num = 1;
            r.all_targets = {tidx};
            r.nhc = nhc;
            r.anomaly = anom;
            r.sig = new_sig;
            r.valid = true;
            results.push_back(r);
        }

        // (B) Subset attachment to -1 curves (all non-empty subsets)
        int nm1 = (int)m1_curves.size();
        if (nm1 > 0) {
            int max_mask = (1 << nm1) - 1;
            for (int mask = 1; mask <= max_mask; mask++) {
                std::vector<int> targets;
                for (int b = 0; b < nm1; b++)
                    if (mask & (1 << b)) targets.push_back(m1_curves[b]);

                Eigen::MatrixXi new_IF = (targets.size() == 1)
                    ? attach_curve(base_IF, targets[0], -1, 1)
                    : attach_hat1_multi(base_IF, targets);
                int hat1_idx = new_IF.rows() - 1;

                auto new_sig = compute_sig(new_IF);
                if (new_sig.sig_neg > config.T_max) continue;

                NHCResult nhc = {true, "", {}, {}, {}};
                if (config.check_nhc) {
                    nhc = check_nhc(new_IF, config.cc_budget);
                    if (!nhc.passes) continue;
                }
                enhance_hat1_gauge(nhc, new_IF, hat1_idx, -1, config.cc_budget);
                if (!nhc.passes) continue;

                AnomalyResult anom = {0, 0, 0, 0};
                if (nhc.passes) {
                    anom = compute_anomaly(new_IF, nhc);
                    if (config.check_anomaly && anom.H_neutral < 0) continue;
                }

                Hat1Result r;
                r.catalog_id = entry.id;
                r.catalog_type = entry.type;
                r.base_IF = base_IF;
                r.final_IF = new_IF;
                r.hat1_idx = hat1_idx;
                r.target_idx = targets[0];
                r.target_si = -1;
                r.int_num = 1;
                r.all_targets = targets;
                r.nhc = nhc;
                r.anomaly = anom;
                r.sig = new_sig;
                r.valid = true;
                results.push_back(r);
            }
        }

        processed++;
        if (config.verbose && processed % 1000 == 0)
            std::cout << processed << " entries → " << results.size() << " \u0125at{1} bases\n";
    }

    // Also try dummy LSTs
    if (config.include_dummy) {
        // Helper: generate hat1 subsets from a base IF
        auto gen_from_dummy = [&](const Eigen::MatrixXi& base_IF, int id, const std::string& type) {
            int n = base_IF.rows();
            std::vector<int> m1_curves, m2_curves;
            for (int i = 0; i < n; i++) {
                if (base_IF(i, i) == -1) m1_curves.push_back(i);
                else if (base_IF(i, i) == -2) m2_curves.push_back(i);
            }
            // Single -2 targets
            for (int tidx : m2_curves) {
                Eigen::MatrixXi new_IF = attach_curve(base_IF, tidx, -1, 1);
                int hat1_idx = new_IF.rows() - 1;
                auto new_sig = compute_sig(new_IF);
                if (new_sig.sig_neg > config.T_max) continue;
                NHCResult nhc = {true, "", {}, {}, {}};
                if (config.check_nhc) {
                    nhc = check_nhc(new_IF, config.cc_budget);
                    if (!nhc.passes) continue;
                }
                enhance_hat1_gauge(nhc, new_IF, hat1_idx, -2, config.cc_budget);
                if (!nhc.passes) continue;
                AnomalyResult anom = {0, 0, 0, 0};
                if (nhc.passes) {
                    anom = compute_anomaly(new_IF, nhc);
                    if (config.check_anomaly && anom.H_neutral < 0) continue;
                }
                Hat1Result r;
                r.catalog_id = id; r.catalog_type = type;
                r.base_IF = base_IF; r.final_IF = new_IF;
                r.hat1_idx = hat1_idx; r.target_idx = tidx;
                r.target_si = -2; r.int_num = 1; r.all_targets = {tidx};
                r.nhc = nhc; r.anomaly = anom; r.sig = new_sig; r.valid = true;
                results.push_back(r);
            }
            // Subsets of -1 targets
            int nm1 = (int)m1_curves.size();
            if (nm1 > 0) {
                int max_mask = (1 << nm1) - 1;
                for (int mask = 1; mask <= max_mask; mask++) {
                    std::vector<int> targets;
                    for (int b = 0; b < nm1; b++)
                        if (mask & (1 << b)) targets.push_back(m1_curves[b]);
                    Eigen::MatrixXi new_IF = (targets.size() == 1)
                        ? attach_curve(base_IF, targets[0], -1, 1)
                        : attach_hat1_multi(base_IF, targets);
                    int hat1_idx = new_IF.rows() - 1;
                    auto new_sig = compute_sig(new_IF);
                    if (new_sig.sig_neg > config.T_max) continue;
                    NHCResult nhc = {true, "", {}, {}, {}};
                    if (config.check_nhc) {
                        nhc = check_nhc(new_IF, config.cc_budget);
                        if (!nhc.passes) continue;
                    }
                    enhance_hat1_gauge(nhc, new_IF, hat1_idx, -1, config.cc_budget);
                    if (!nhc.passes) continue;
                    AnomalyResult anom = {0, 0, 0, 0};
                    if (nhc.passes) {
                        anom = compute_anomaly(new_IF, nhc);
                        if (config.check_anomaly && anom.H_neutral < 0) continue;
                    }
                    Hat1Result r;
                    r.catalog_id = id; r.catalog_type = type;
                    r.base_IF = base_IF; r.final_IF = new_IF;
                    r.hat1_idx = hat1_idx; r.target_idx = targets[0];
                    r.target_si = -1; r.int_num = 1; r.all_targets = targets;
                    r.nhc = nhc; r.anomaly = anom; r.sig = new_sig; r.valid = true;
                    results.push_back(r);
                }
            }
        };

        int a_min = std::max(2, config.catalog_T_min + 1);
        int a_max = config.catalog_T_max + 1;
        for (int nc = a_min; nc <= a_max; nc++)
            gen_from_dummy(build_dummy_A(nc), -nc, "DM:A(" + std::to_string(nc) + ")");

        int d_min = std::max(3, config.catalog_T_min + 1);
        int d_max = config.catalog_T_max + 1;
        for (int nc = d_min; nc <= d_max; nc++)
            gen_from_dummy(build_dummy_D(nc), -(1000 + nc), "DM:D(" + std::to_string(nc) + ")");
    }

    if (config.verbose)
        std::cout << "Total hat{1} bases (before dedup): " << results.size() << "\n";
    return results;
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <catalog> [T_max] [T_min]\n";
        return 1;
    }
    std::string catalog_file = argv[1];
    int T_max = (argc > 2) ? std::atoi(argv[2]) : 193;
    int T_min = (argc > 3) ? std::atoi(argv[3]) : 0;
    std::string suffix = "_T" + std::to_string(T_min) + "_" + std::to_string(T_max);

    auto catalog = load_catalog(catalog_file);
    std::cout << "Catalog: " << catalog.size() << " entries\n";

    // Pre-compute block maps for Nk entries
    std::map<int, TopoBlockMap> catalog_block_maps;
    for (auto& entry : catalog) {
        if (!entry.is_nn()) {
            auto bm = get_block_map(entry);
            if (!bm.blocks.empty())
                catalog_block_maps[entry.id] = std::move(bm);
        }
    }
    std::cout << "Block maps: " << catalog_block_maps.size() << " Nk entries\n";

    // Configuration
    SUGRAConfig config;
    config.T_max = 193;
    config.check_nhc = true;
    config.check_anomaly = true;
    config.check_determinant = false;
    config.catalog_T_min = T_min;
    config.catalog_T_max = T_max;
    config.verbose = true;

    // Generate \hat{1} bases
    auto results = generate_hat1(catalog, config);

    // ── Deduplicate ──
    // Key: eigenvalues + edge structure (including info about which is \hat{1})
    auto dedup_key = [](const Hat1Result& r) -> std::string {
        const auto& IF = r.final_IF;
        int n = IF.rows();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
        auto ev = es.eigenvalues();
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(8);
        ss << "E:";
        for (int i = 0; i < ev.size(); i++) {
            double v = ev(i);
            if (std::abs(v) < 1e-10) v = 0.0;
            ss << v << ";";
        }
        // Edge multiset
        std::vector<std::tuple<int, int, int>> edges;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (IF(i, j) != 0) {
                    int a = IF(i, i), b = IF(j, j);
                    if (a > b) std::swap(a, b);
                    edges.push_back({a, b, IF(i, j)});
                }
        std::sort(edges.begin(), edges.end());
        ss << "|E:";
        for (auto& [a, b, w] : edges) ss << "(" << a << "," << b << "," << w << ")";
        // Also encode which self-int the \hat{1} attaches to and how many
        ss << "|T:" << r.target_si << "|N:" << r.num_targets();
        return ss.str();
    };

    std::set<std::string> seen;
    std::vector<Hat1Result> filtered;
    for (auto& r : results) {
        std::string key = dedup_key(r);
        if (seen.insert(key).second)
            filtered.push_back(std::move(r));
    }
    std::cout << "After dedup: " << filtered.size() << "\n";

    // Sort by T, then target_si, then det
    std::sort(filtered.begin(), filtered.end(),
              [](const Hat1Result& a, const Hat1Result& b) {
                  if (a.anomaly.T != b.anomaly.T) return a.anomaly.T < b.anomaly.T;
                  if (a.target_si != b.target_si) return a.target_si > b.target_si;  // -1 first, then -2
                  if (a.anomaly.H_neutral != b.anomaly.H_neutral)
                      return a.anomaly.H_neutral < b.anomaly.H_neutral;
                  return a.catalog_id < b.catalog_id;
              });

    // Split: unimodular vs non-unimodular
    std::vector<Hat1Result> uni, nonuni;
    for (auto& r : filtered) {
        if (std::abs(r.sig.det) == 1)
            uni.push_back(r);
        else
            nonuni.push_back(r);
    }
    std::cout << "Unimodular: " << uni.size() << "\n";
    std::cout << "Non-unimodular: " << nonuni.size() << "\n";

    // Count by target type
    int n_to_m1 = 0, n_to_m2 = 0;
    for (auto& r : filtered) {
        if (r.target_si == -1) n_to_m1++;
        else if (r.target_si == -2) n_to_m2++;
    }
    std::cout << "  attached to -1: " << n_to_m1 << "\n";
    std::cout << "  attached to -2: " << n_to_m2 << "\n";

    auto get_blocks = [&](int cat_id) -> const TopoBlockMap& {
        static const TopoBlockMap empty;
        auto it = catalog_block_maps.find(cat_id);
        return (it != catalog_block_maps.end()) ? it->second : empty;
    };

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

    auto write_preamble = [&](std::ofstream& tex) {
        tex << "\\documentclass[10pt]{article}\n"
            << "\\usepackage[a4paper,margin=0.6in]{geometry}\n"
            << "\\usepackage{longtable,booktabs,amsmath,amssymb,array}\n"
            << "\\usepackage[dvipsnames]{xcolor}\n"
            << "\\setlength{\\parindent}{0pt}\n"
            << "\\setlength{\\parskip}{2pt}\n"
            << "\\begin{document}\n\n";
    };

    // Sort entries by T, then base spectrum
    auto sort_entries = [&](std::vector<Hat1Result>& v) {
        std::map<int, std::string> keys;
        for (int i = 0; i < (int)v.size(); i++)
            keys[i] = base_spec_key(v[i].base_IF);
        std::vector<int> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            if (v[a].anomaly.T != v[b].anomaly.T) return v[a].anomaly.T < v[b].anomaly.T;
            if (v[a].target_si != v[b].target_si) return v[a].target_si > v[b].target_si;
            return keys[a] < keys[b];
        });
        std::vector<Hat1Result> sorted;
        sorted.reserve(v.size());
        for (int i : idx) sorted.push_back(std::move(v[i]));
        v = std::move(sorted);
    };
    sort_entries(uni);
    sort_entries(nonuni);

    // ── Unimodular: filter by standard b₀ relation ──
    // After blowdown, all remaining curves must satisfy b₀·C_i = C_i²+2.
    // hat{1} always breaks this (b₀·hat1 = -1 ≠ C²+2 = 1), so all are rejected.
    {
        auto is_standard_b0 = [](const BlowdownResult& bd) -> bool {
            auto& M = bd.extended_IF;
            int n = M.rows();
            int b0 = n - 1;
            for (int i = 0; i < b0; i++) {
                if (M(i, b0) != M(i, i) + 2) return false;
            }
            return true;
        };

        int n_std = 0, n_nonstd = 0;
        for (auto& r : uni) {
            auto bd = blowdown_hat1(r.final_IF, r.hat1_idx);
            if (is_standard_b0(bd)) n_std++;
            else n_nonstd++;
        }
        std::cout << "Unimodular: standard b0=" << n_std
                  << ", non-standard b0=" << n_nonstd
                  << " (all rejected)\n";
    }

    // ── Non-unimodular output ──
    {
        // Split by target type
        struct NonuniEntry {
            const Hat1Result* r;
            CritInfo ci;
        };

        std::vector<NonuniEntry> to_m1, to_m2;
        int n_rejected = 0;
        for (auto& r : nonuni) {
            auto ci = compute_crit_hat1(r.final_IF, r.hat1_idx);
            int delta = r.anomaly.H_charged - r.anomaly.V + 29 * r.anomaly.T - 273;
            int delta_crit = delta + 29 * (ci.T_crit_int - r.anomaly.T);
            if (delta_crit > 0) { n_rejected++; continue; }
            NonuniEntry e{&r, ci};
            if (r.target_si == -1) to_m1.push_back(e);
            else to_m2.push_back(e);
        }
        std::cout << "Non-uni: " << nonuni.size() << " total, "
                  << n_rejected << " rejected, "
                  << nonuni.size() - n_rejected << " pass\n";
        std::cout << "  to -1: " << to_m1.size()
                  << ", to -2: " << to_m2.size() << "\n";

        auto write_nonuni = [&](const std::string& fname, const std::string& title,
                                std::vector<NonuniEntry>& entries) {
            std::ofstream tex(fname);
            write_preamble(tex);
            tex << "\\section*{" << title << "}\n\n";

            tex << "$\\hat{1}$ external: $C^2=-1,\\; b_0\\cdot C=-1$.\n\n";

            std::set<std::string> lst_bases;
            for (auto& e : entries) lst_bases.insert(base_spec_key(e.r->base_IF));

            int n_exact = 0, n_frac = 0;
            for (auto& e : entries) {
                if (e.ci.exact) n_exact++; else n_frac++;
            }

            tex << "Total: " << entries.size()
                << " bases (from " << lst_bases.size() << " distinct LST bases, "
                << "$T_{\\text{LST}} \\in [" << T_min << ", " << T_max << "]$).\n\n"
                << "Exact $\\det=0$: " << n_exact
                << ", fractional $T_{\\text{crit}}$: " << n_frac
                << ". Filtered: $\\Delta_{\\text{crit}} \\leq 0$.\n\n"
                << "\\textcolor{red}{$\\hat{1}$} = attached external.\n\n";

            // Statistics
            std::map<int, int> T_count, Tcrit_count;
            for (auto& e : entries) {
                T_count[e.r->anomaly.T]++;
                Tcrit_count[e.ci.T_crit_int]++;
            }

            tex << "\\begin{center}\n";
            tex << "\\begin{tabular}{c|c}\n\\toprule\n$T$ & Count \\\\\n\\midrule\n";
            for (auto& [T, cnt] : T_count) tex << T << " & " << cnt << " \\\\\n";
            tex << "\\bottomrule\n\\end{tabular}\n\\qquad\\qquad\n";
            tex << "\\begin{tabular}{c|c}\n\\toprule\n$T_{\\text{crit}}$ & Count \\\\\n\\midrule\n";
            for (auto& [T, cnt] : Tcrit_count) tex << T << " & " << cnt << " \\\\\n";
            tex << "\\bottomrule\n\\end{tabular}\n";
            tex << "\\end{center}\n\n";

            // |det| distribution
            std::map<int, int> det_count;
            for (auto& e : entries) det_count[std::abs(e.r->sig.det)]++;
            tex << "\\begin{center}\n";
            tex << "\\begin{tabular}{c|c}\n\\toprule\n$|\\det|$ & Count \\\\\n\\midrule\n";
            for (auto& [d, cnt] : det_count) tex << d << " & " << cnt << " \\\\\n";
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
                    tex << "\\smallskip{\\small\\textbf{Base:} $"
                        << if_to_latex_hat1(r.base_IF, -2, get_blocks(r.catalog_id))
                        << "$}\\smallskip\n\n";
                }
                int delta = r.anomaly.H_charged - r.anomaly.V + 29 * r.anomaly.T - 273;
                tex << "$" << render_hat1_latex(r, get_blocks(r.catalog_id)) << "$";
                tex << " \\quad $\\Delta=" << delta
                    << ",\\; |\\det|=" << std::abs(r.sig.det)
                    << ",\\; (n_+,n_-,n_0)=(" << r.sig.sig_pos
                    << "," << r.sig.sig_neg << "," << r.sig.sig_zero << ")$\n\n";
                // T_crit line
                tex << "\\hspace{1em}$T_{\\text{crit}}=";
                if (e.ci.exact) {
                    tex << e.ci.T_crit_int;
                } else {
                    tex << std::fixed << std::setprecision(2) << e.ci.T_crit
                        << " \\to " << e.ci.T_crit_int;
                }
                tex << ",\\; (n_+,n_-,n_0)_{\\text{crit}}=("
                    << e.ci.crit_npos << "," << e.ci.crit_nneg << "," << e.ci.crit_nzero << ")$\n\n";
                // Blowdowned extended matrix on its own line
                auto& em = e.ci.bd_crit.extended_IF;
                int sz = em.rows();
                tex << "\\hspace{2em}$\\left(\\begin{smallmatrix}";
                for (int i = 0; i < sz; i++) {
                    if (i > 0) tex << " \\\\ ";
                    for (int j = 0; j < sz; j++) {
                        if (j > 0) tex << " & ";
                        tex << em(i, j);
                    }
                }
                tex << "\\end{smallmatrix}\\right)_{\\!b_0}$\n\n";
            }
            tex << "\\end{document}\n";
            tex.close();
            std::cout << "Written " << entries.size() << " to " << fname << "\n";
        };

        write_nonuni(("sugra_hat1_nonuni_to_m1" + suffix + ".tex").c_str(),
                     "Non-unimodular $\\hat{1}$ SUGRA Bases: target $= (-1)$", to_m1);
        write_nonuni(("sugra_hat1_nonuni_to_m2" + suffix + ".tex").c_str(),
                     "Non-unimodular $\\hat{1}$ SUGRA Bases: target $= (-2)$", to_m2);
    }

    // Summary
    std::cout << "\n=== Summary ===\n";
    std::cout << "hat{1} bases total: " << filtered.size() << "\n";
    std::cout << "  unimodular: " << uni.size() << "\n";
    std::cout << "  non-unimodular: " << nonuni.size() << "\n";
    std::cout << "  to (-1): " << n_to_m1 << "\n";
    std::cout << "  to (-2): " << n_to_m2 << "\n";

    // T distribution
    std::map<int, int> T_dist;
    for (auto& r : filtered) T_dist[r.anomaly.T]++;
    std::cout << "T distribution:\n";
    for (auto& [T, cnt] : T_dist) std::cout << "  T=" << T << ": " << cnt << "\n";

    return 0;
}
