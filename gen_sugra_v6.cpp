// gen_sugra_v6.cpp
// SUGRA base generator: multi-target external attachment
//
// External properties:
//   C^2 ∈ {-5, -6, -7, -8, -12}
//   b_0·C = C^2 + 2 (standard adjunction)
//   Attaches to: (-1) curves only, intersection number = 1
//   Can attach to MULTIPLE (-1) curves simultaneously
//
// Based on v2 (standard externals) + v5 (multi-target machinery)
//
// Compile:
//   g++ -std=c++17 -O2 -o gen_sugra_v6 gen_sugra_v6.cpp Tensor.C \
//       TopoLineCompact_enhanced.cpp Topology_enhanced.cpp \
//       -I/usr/include/eigen3 -I.
//
// Usage:
//   ./gen_sugra_v6 unified.cat [T_max=193] [T_min=0]

#include "sugra_generator.h"
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <functional>

// ============================================================================
// Block map (same as v2/v5)
// ============================================================================

enum class BlockType { Node, InteriorLink, SideLink, Instanton, External };

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
// v6 Result struct
// ============================================================================

struct ExtResult {
    int catalog_id;
    std::string catalog_type;
    Eigen::MatrixXi base_IF;
    Eigen::MatrixXi final_IF;

    int ext_idx;                     // index of external in final_IF
    int ext_si;                      // self-intersection: -5,-6,-7,-8,-12
    std::vector<int> all_targets;    // indices of (-1) curves in base_IF

    NHCResult nhc;
    AnomalyResult anomaly;
    SigInfo sig;
    bool valid = false;

    int num_targets() const { return (int)all_targets.size(); }
};

// ============================================================================
// Multi-target attachment
// ============================================================================

inline Eigen::MatrixXi attach_ext_multi(const Eigen::MatrixXi& IF,
                                         const std::vector<int>& targets, int ext_si) {
    int n = IF.rows();
    Eigen::MatrixXi nIF = Eigen::MatrixXi::Zero(n + 1, n + 1);
    nIF.block(0, 0, n, n) = IF;
    nIF(n, n) = ext_si;
    for (int t : targets) {
        nIF(n, t) = 1;
        nIF(t, n) = 1;
    }
    return nIF;
}

// ============================================================================
// b₀Q matrix (standard: b₀·C = C²+2 for ALL curves including external)
// ============================================================================

inline Eigen::MatrixXi build_b0Q_matrix(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    int T = 0;
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
        for (int i = 0; i < n; i++)
            if (es.eigenvalues()(i) < -1e-8) T++;
    }
    int b0sq = 9 - T;
    Eigen::MatrixXi ext = Eigen::MatrixXi::Zero(n + 1, n + 1);
    ext.block(0, 0, n, n) = IF;
    for (int i = 0; i < n; i++) {
        int b0C = IF(i, i) + 2;
        ext(i, n) = b0C;
        ext(n, i) = b0C;
    }
    ext(n, n) = b0sq;
    return ext;
}

// ============================================================================
// Blowdown (standard, from v2)
// ============================================================================

inline Eigen::MatrixXi blowdown_curve(const Eigen::MatrixXi& IF, int idx) {
    int n = IF.rows();
    Eigen::MatrixXi nIF = Eigen::MatrixXi::Zero(n - 1, n - 1);
    for (int i = 0, ni = 0; i < n; i++) {
        if (i == idx) continue;
        for (int j = 0, nj = 0; j < n; j++) {
            if (j == idx) continue;
            nIF(ni, nj) = IF(i, j) + IF(i, idx) * IF(idx, j);
            nj++;
        }
        ni++;
    }
    return nIF;
}

struct BlowdownResult {
    Eigen::MatrixXi extended_IF;
    Eigen::MatrixXi minimal_IF;
    int b0_sq;
    int num_blowdowns;
    int F_n;
    bool is_P2;
    std::string surface;
};

inline BlowdownResult blowdown(const Eigen::MatrixXi& IF) {
    BlowdownResult res;
    int n = IF.rows();
    int T = 0;
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
        for (int i = 0; i < n; i++)
            if (es.eigenvalues()(i) < -1e-8) T++;
    }
    int b0sq = 9 - T;

    Eigen::MatrixXi ext = Eigen::MatrixXi::Zero(n + 1, n + 1);
    ext.block(0, 0, n, n) = IF;
    for (int i = 0; i < n; i++) {
        ext(i, n) = IF(i, i) + 2;
        ext(n, i) = IF(i, i) + 2;
    }
    ext(n, n) = b0sq;

    Eigen::MatrixXi cur = ext;
    int bd = 0;
    while (true) {
        int nn = cur.rows();
        int m1 = -1;
        for (int i = 0; i < nn - 1; i++)
            if (cur(i, i) == -1) { m1 = i; break; }
        if (m1 < 0) break;
        cur = blowdown_curve(cur, m1);
        bd++;
    }

    int nc = cur.rows() - 1;
    res.extended_IF = cur;
    res.minimal_IF = (nc > 0) ? cur.block(0, 0, nc, nc) : Eigen::MatrixXi();
    res.b0_sq = cur(cur.rows() - 1, cur.cols() - 1);
    res.num_blowdowns = bd;
    res.F_n = -1;
    res.is_P2 = false;

    if (nc == 0) {
        res.surface = "\\text{pt}";
    } else if (nc == 1) {
        if (cur(0, 0) == 1) {
            res.is_P2 = true;
            res.surface = "\\mathbb{P}^2";
        } else {
            res.surface = "C^2=" + std::to_string(cur(0, 0));
        }
    } else if (nc == 2) {
        int a = cur(0, 0), b = cur(1, 1), off = cur(0, 1);
        if (off == 1) {
            if (a == 0) { res.F_n = -b; res.surface = "\\mathbb{F}_{" + std::to_string(-b) + "}"; }
            else if (b == 0) { res.F_n = -a; res.surface = "\\mathbb{F}_{" + std::to_string(-a) + "}"; }
            else { res.surface = "(" + std::to_string(a) + "," + std::to_string(b) + ")"; }
        } else {
            res.surface = "(" + std::to_string(a) + "," + std::to_string(b) + ")_{" + std::to_string(off) + "}";
        }
    } else {
        res.surface = "\\text{bd}_{" + std::to_string(bd) + "}";
    }
    return res;
}

// ============================================================================
// CritInfo (standard b₀, from v2)
// ============================================================================

struct CritInfo {
    double T_crit;
    int T_crit_int;
    bool exact;
    int crit_npos, crit_nneg, crit_nzero;
    BlowdownResult bd_crit;
};

inline CritInfo compute_crit_info(const Eigen::MatrixXi& IF) {
    CritInfo res;
    int n = IF.rows();

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n + 1, n + 1);
    M.block(0, 0, n, n) = IF.cast<double>();
    for (int i = 0; i < n; i++) {
        double b0C = IF(i, i) + 2.0;
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

    int b0sq = 9 - res.T_crit_int;
    M(n, n) = (double)b0sq;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
    res.crit_npos = res.crit_nneg = res.crit_nzero = 0;
    for (int i = 0; i < n + 1; i++) {
        double ev = es.eigenvalues()(i);
        if (std::abs(ev) < 1e-8) res.crit_nzero++;
        else if (ev > 0) res.crit_npos++;
        else res.crit_nneg++;
    }

    // Blowdown at T_crit
    Eigen::MatrixXi ext = Eigen::MatrixXi::Zero(n + 1, n + 1);
    ext.block(0, 0, n, n) = IF;
    for (int i = 0; i < n; i++) {
        ext(i, n) = IF(i, i) + 2;
        ext(n, i) = IF(i, i) + 2;
    }
    ext(n, n) = b0sq;

    Eigen::MatrixXi cur = ext;
    int bd_count = 0;
    while (true) {
        int nn = cur.rows();
        int m1 = -1;
        for (int i = 0; i < nn - 1; i++)
            if (cur(i, i) == -1) { m1 = i; break; }
        if (m1 < 0) break;
        cur = blowdown_curve(cur, m1);
        bd_count++;
    }

    int nc = cur.rows() - 1;
    res.bd_crit.extended_IF = cur;
    res.bd_crit.minimal_IF = (nc > 0) ? cur.block(0, 0, nc, nc) : Eigen::MatrixXi();
    res.bd_crit.b0_sq = cur(cur.rows() - 1, cur.cols() - 1);
    res.bd_crit.num_blowdowns = bd_count;
    res.bd_crit.F_n = -1;
    res.bd_crit.is_P2 = false;

    if (nc == 0) res.bd_crit.surface = "\\text{pt}";
    else if (nc == 1 && cur(0, 0) == 1) { res.bd_crit.is_P2 = true; res.bd_crit.surface = "\\mathbb{P}^2"; }
    else if (nc == 2) {
        int a = cur(0, 0), b = cur(1, 1);
        if (cur(0, 1) == 1 && (a == 0 || b == 0)) {
            int fn = (a == 0) ? -b : -a;
            res.bd_crit.F_n = fn;
            res.bd_crit.surface = "\\mathbb{F}_{" + std::to_string(fn) + "}";
        } else {
            res.bd_crit.surface = "\\text{bd}_{" + std::to_string(bd_count) + "}";
        }
    } else {
        res.bd_crit.surface = "\\text{bd}_{" + std::to_string(bd_count) + "}";
    }
    return res;
}

// ============================================================================
// LaTeX rendering (from v2, with ext_idx support + multi-target)
// ============================================================================

inline std::string gauge_algebra_latex(int param) {
    switch (param) {
        case 5:  return "\\mathfrak{f}_4";
        case 6:  return "\\mathfrak{e}_6";
        case 7:  return "\\mathfrak{e}_7'";
        case 8:  return "\\mathfrak{e}_7";
        case 12: return "\\mathfrak{e}_8";
        default: return std::to_string(param);
    }
}

inline std::string interior_link_label(int param) {
    switch (param) {
        case 4: return "\\overset{3,3}{\\otimes}";
        case 5: return "\\overset{5,5}{\\otimes}";
        default: return "\\overset{" + std::to_string(param) + "}{\\otimes}";
    }
}

// Single-target rendering (v2-style: external as red number in chain)
std::string if_to_latex(const Eigen::MatrixXi& IF, int ext_idx,
                         const TopoBlockMap& bmap = {}) {
    int n = IF.rows();
    if (n == 0) return "\\emptyset";
    if (n == 1) return std::to_string(-IF(0, 0));

    bool has_blocks = !bmap.blocks.empty();

    std::vector<std::vector<std::pair<int, int>>> adj(n);
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (IF(i, j) != 0) {
                adj[i].push_back({j, IF(i, j)});
                adj[j].push_back({i, IF(j, i)});
            }

    int start = -1;
    for (int i = 0; i < n; i++)
        if ((int)adj[i].size() == 1 && i != ext_idx) { start = i; break; }
    if (start < 0)
        for (int i = 0; i < n; i++)
            if ((int)adj[i].size() == 1) { start = i; break; }
    if (start < 0) start = 0;

    std::vector<bool> visited(n, false);

    auto curve_label = [&](int cur) -> std::string {
        int si = -IF(cur, cur);
        std::string si_str = std::to_string(si);
        if (si >= 10) si_str = "(" + si_str + ")";
        if (cur == ext_idx) return "\\textcolor{red}{" + si_str + "}";
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

    std::function<std::string(int)> walk_branch;

    std::function<std::string(int, int)> walk = [&](int cur, int from) -> std::string {
        visited[cur] = true;

        if (has_blocks && bmap.is_interior_link(cur)
            && bmap.block_of(cur).end - bmap.block_of(cur).start > 1) {
            bool ext_touches_block = false;
            if (ext_idx >= 0) {
                auto& blk = bmap.block_of(cur);
                for (int c = blk.start; c < blk.end && !ext_touches_block; c++)
                    for (auto& [nb, w] : adj[c])
                        if (nb == ext_idx) { ext_touches_block = true; break; }
            }
            if (!ext_touches_block) {
                std::string label = interior_link_label(bmap.block_of(cur).param);
                auto exits = skip_block(cur);
                if (exits.empty()) return label;
                if (exits.size() == 1) return label + walk(exits[0], cur);
                int main_child = -1;
                std::vector<int> branch_children;
                for (int c : exits) {
                    if (main_child < 0 && c != ext_idx) main_child = c;
                    else branch_children.push_back(c);
                }
                if (main_child < 0) {
                    main_child = exits[0];
                    branch_children.assign(exits.begin() + 1, exits.end());
                }
                std::vector<std::string> br_strs;
                for (int i = (int)branch_children.size() - 1; i >= 0; i--)
                    br_strs.push_back(walk_branch(branch_children[i]));
                std::string main_part = walk(main_child, cur);
                if (br_strs.empty()) return label + main_part;
                if (br_strs.size() == 1)
                    return "\\overset{" + br_strs[0] + "}{" + label + "}" + main_part;
                std::string under;
                for (size_t i = 1; i < br_strs.size(); i++) {
                    if (!under.empty()) under += ",";
                    under += br_strs[i];
                }
                return "\\overset{" + br_strs[0] + "}{\\underset{" + under + "}{" + label + "}}" + main_part;
            }
            visited[cur] = true;
        }

        std::string label;
        if (has_blocks && bmap.is_node(cur))
            label = gauge_algebra_latex(bmap.block_of(cur).param);
        else
            label = curve_label(cur);

        std::vector<int> children;
        for (auto& [nb, w] : adj[cur])
            if (!visited[nb]) children.push_back(nb);

        if (children.empty()) return label;
        if (children.size() == 1) return label + walk(children[0], cur);

        int main_child = -1;
        std::vector<int> branch_children;
        for (int c : children) {
            if (main_child < 0 && c != ext_idx) main_child = c;
            else branch_children.push_back(c);
        }
        if (main_child < 0) {
            main_child = children[0];
            branch_children.clear();
            for (size_t i = 1; i < children.size(); i++)
                branch_children.push_back(children[i]);
        }

        std::vector<std::string> br_strs;
        for (int i = (int)branch_children.size() - 1; i >= 0; i--)
            br_strs.push_back(walk_branch(branch_children[i]));

        std::string main_part = walk(main_child, cur);
        if (br_strs.size() == 1)
            return "\\overset{" + br_strs[0] + "}{" + label + "}" + main_part;
        std::string under;
        for (size_t i = 1; i < br_strs.size(); i++) {
            if (!under.empty()) under += ",";
            under += br_strs[i];
        }
        return "\\overset{" + br_strs[0] + "}{\\underset{" + under + "}{" + label + "}}" + main_part;
    };

    walk_branch = [&](int cur) -> std::string {
        visited[cur] = true;
        std::string label;
        if (has_blocks && bmap.is_node(cur))
            label = gauge_algebra_latex(bmap.block_of(cur).param);
        else
            label = curve_label(cur);

        std::vector<int> children;
        for (auto& [nb, w] : adj[cur])
            if (!visited[nb]) children.push_back(nb);

        if (children.empty()) return label;
        if (children.size() == 1) return walk_branch(children[0]) + label;

        int main_child = -1;
        std::vector<int> sub_branches;
        for (int c : children) {
            if (main_child < 0 && c != ext_idx) main_child = c;
            else sub_branches.push_back(c);
        }
        if (main_child < 0) {
            main_child = children[0];
            sub_branches.clear();
            for (size_t i = 1; i < children.size(); i++)
                sub_branches.push_back(children[i]);
        }
        std::vector<std::string> sub_strs;
        for (int i = (int)sub_branches.size() - 1; i >= 0; i--)
            sub_strs.push_back(walk_branch(sub_branches[i]));
        std::string main_part = walk_branch(main_child);
        if (sub_strs.size() == 1)
            return main_part + "\\overset{" + sub_strs[0] + "}{" + label + "}";
        std::string under;
        for (size_t i = 1; i < sub_strs.size(); i++) {
            if (!under.empty()) under += ",";
            under += sub_strs[i];
        }
        return main_part + "\\overset{" + sub_strs[0] + "}{\\underset{" + under + "}{" + label + "}}";
    };

    return walk(start, -1);
}

// Multi-target rendering: base quiver with targets in blue + red external annotation
std::string if_to_latex_multi(const Eigen::MatrixXi& base_IF,
                               const std::vector<int>& targets, int ext_si,
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
        std::vector<std::string> br_strs;
        for (int i = (int)branch_children.size() - 1; i >= 0; i--)
            br_strs.push_back(walk(branch_children[i], cur));
        std::string main_part = walk(main_child, cur);
        if (br_strs.size() == 1)
            return "\\overset{" + br_strs[0] + "}{" + label + "}" + main_part;
        if (br_strs.size() >= 2) {
            std::string under;
            for (size_t i = 1; i < br_strs.size(); i++) {
                if (!under.empty()) under += ",";
                under += br_strs[i];
            }
            return "\\overset{" + br_strs[0] + "}{\\underset{" + under + "}{" + label + "}}" + main_part;
        }
        return label + main_part;
    };

    std::string quiver = walk(start, -1);
    // Annotation: + red external
    int abs_si = -ext_si;
    std::string ext_str = std::to_string(abs_si);
    if (abs_si >= 10) ext_str = "(" + ext_str + ")";
    return quiver + "+\\textcolor{red}{" + ext_str + "}";
}

// Dispatcher
std::string render_ext_latex(const ExtResult& r, const TopoBlockMap& bmap = {}) {
    if (r.num_targets() <= 1) {
        return if_to_latex(r.final_IF, r.ext_idx, bmap);
    } else {
        return if_to_latex_multi(r.base_IF, r.all_targets, r.ext_si, bmap);
    }
}

// ============================================================================
// Generation
// ============================================================================

static const int EXT_SI_LIST[] = {-5, -6, -7, -8, -12};
static const int N_EXT_SI = 5;

std::vector<ExtResult> generate_v6(
    const std::vector<CatalogEntry>& catalog,
    const SUGRAConfig& config)
{
    std::vector<ExtResult> results;

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

        // Collect (-1) curves
        std::vector<int> m1_curves;
        for (int i = 0; i < n; i++)
            if (base_IF(i, i) == -1) m1_curves.push_back(i);

        if (m1_curves.empty()) continue;

        int nm1 = (int)m1_curves.size();
        if (nm1 > 20) nm1 = 20;  // cap to avoid explosion

        // For each external type
        for (int ei = 0; ei < N_EXT_SI; ei++) {
            int ext_si = EXT_SI_LIST[ei];

            // Enumerate all non-empty subsets of (-1) curves
            int max_mask = (1 << nm1) - 1;
            for (int mask = 1; mask <= max_mask; mask++) {
                std::vector<int> targets;
                for (int b = 0; b < nm1; b++)
                    if (mask & (1 << b)) targets.push_back(m1_curves[b]);

                Eigen::MatrixXi new_IF = (targets.size() == 1)
                    ? attach_curve(base_IF, targets[0], ext_si, 1)
                    : attach_ext_multi(base_IF, targets, ext_si);
                int ext_idx = new_IF.rows() - 1;

                auto new_sig = compute_sig(new_IF);
                if (new_sig.sig_neg > config.T_max) continue;

                NHCResult nhc = {true, "", {}, {}, {}};
                if (config.check_nhc) {
                    nhc = check_nhc(new_IF, config.cc_budget);
                    if (!nhc.passes) continue;
                }

                AnomalyResult anom = {0, 0, 0, 0};
                if (nhc.passes) {
                    anom = compute_anomaly(new_IF, nhc);
                    if (config.check_anomaly && anom.H_neutral < 0) continue;
                }

                ExtResult r;
                r.catalog_id = entry.id;
                r.catalog_type = entry.type;
                r.base_IF = base_IF;
                r.final_IF = new_IF;
                r.ext_idx = ext_idx;
                r.ext_si = ext_si;
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
            std::cout << processed << " → " << results.size() << " bases\n";
    }

    // Also dummy LSTs
    if (config.include_dummy) {
        auto gen_from_dummy = [&](const Eigen::MatrixXi& base_IF, int id, const std::string& type) {
            int n = base_IF.rows();
            std::vector<int> m1_curves;
            for (int i = 0; i < n; i++)
                if (base_IF(i, i) == -1) m1_curves.push_back(i);
            if (m1_curves.empty()) return;
            int nm1 = std::min((int)m1_curves.size(), 20);

            for (int ei = 0; ei < N_EXT_SI; ei++) {
                int ext_si = EXT_SI_LIST[ei];
                int max_mask = (1 << nm1) - 1;
                for (int mask = 1; mask <= max_mask; mask++) {
                    std::vector<int> targets;
                    for (int b = 0; b < nm1; b++)
                        if (mask & (1 << b)) targets.push_back(m1_curves[b]);
                    Eigen::MatrixXi new_IF = (targets.size() == 1)
                        ? attach_curve(base_IF, targets[0], ext_si, 1)
                        : attach_ext_multi(base_IF, targets, ext_si);
                    int ext_idx = new_IF.rows() - 1;
                    auto new_sig = compute_sig(new_IF);
                    if (new_sig.sig_neg > config.T_max) continue;
                    NHCResult nhc = {true, "", {}, {}, {}};
                    if (config.check_nhc) {
                        nhc = check_nhc(new_IF, config.cc_budget);
                        if (!nhc.passes) continue;
                    }
                    AnomalyResult anom = {0, 0, 0, 0};
                    if (nhc.passes) {
                        anom = compute_anomaly(new_IF, nhc);
                        if (config.check_anomaly && anom.H_neutral < 0) continue;
                    }
                    ExtResult r;
                    r.catalog_id = id; r.catalog_type = type;
                    r.base_IF = base_IF; r.final_IF = new_IF;
                    r.ext_idx = ext_idx; r.ext_si = ext_si;
                    r.all_targets = targets;
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
        std::cout << "Total v6 bases (before dedup): " << results.size() << "\n";
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

    std::map<int, TopoBlockMap> catalog_block_maps;
    for (auto& entry : catalog) {
        if (!entry.is_nn()) {
            auto bm = get_block_map(entry);
            if (!bm.blocks.empty())
                catalog_block_maps[entry.id] = std::move(bm);
        }
    }

    SUGRAConfig config;
    config.T_max = 193;
    config.check_nhc = true;
    config.check_anomaly = true;
    config.check_determinant = false;
    config.catalog_T_min = T_min;
    config.catalog_T_max = T_max;
    config.verbose = true;

    auto results = generate_v6(catalog, config);

    // ── Deduplicate ──
    auto dedup_key = [](const ExtResult& r) -> std::string {
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
        ss << "|X:" << r.ext_si << "|N:" << r.num_targets();
        return ss.str();
    };

    std::set<std::string> seen;
    std::vector<ExtResult> filtered;
    for (auto& r : results) {
        auto key = dedup_key(r);
        if (seen.insert(key).second)
            filtered.push_back(std::move(r));
    }
    std::cout << "After dedup: " << filtered.size() << " bases\n";

    // Sort by T, then ext_si, then base
    std::sort(filtered.begin(), filtered.end(), [](const ExtResult& a, const ExtResult& b) {
        if (a.anomaly.T != b.anomaly.T) return a.anomaly.T < b.anomaly.T;
        if (a.ext_si != b.ext_si) return a.ext_si > b.ext_si;  // -5 before -12
        return a.sig.det < b.sig.det;
    });

    // ── Block map helper ──
    auto get_blocks = [&](int id) -> TopoBlockMap {
        auto it = catalog_block_maps.find(id);
        if (it != catalog_block_maps.end()) return it->second;
        return {};
    };

    auto base_spec_key = [](const Eigen::MatrixXi& IF) -> std::string {
        std::ostringstream ss;
        int n = IF.rows();
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++)
                ss << IF(i, j) << ",";
            ss << ";";
        }
        return ss.str();
    };

    // ── Split by single vs multi target ──
    std::vector<ExtResult*> single_target, multi_target;
    for (auto& r : filtered) {
        if (r.num_targets() == 1) single_target.push_back(&r);
        else multi_target.push_back(&r);
    }

    std::cout << "\n=== Summary ===\n";
    std::cout << "Total: " << filtered.size() << "\n";
    std::cout << "  single-target: " << single_target.size() << "\n";
    std::cout << "  multi-target: " << multi_target.size() << "\n";

    // By ext_si
    std::map<int, int> ext_count, ext_multi_count;
    for (auto& r : filtered) {
        ext_count[r.ext_si]++;
        if (r.num_targets() > 1) ext_multi_count[r.ext_si]++;
    }
    std::cout << "By ext_si:\n";
    for (auto& [si, cnt] : ext_count)
        std::cout << "  " << si << ": " << cnt << " (" << ext_multi_count[si] << " multi)\n";

    // By T
    std::map<int, int> T_dist;
    for (auto& r : filtered) T_dist[r.anomaly.T]++;
    std::cout << "T distribution:\n";
    for (auto& [T, cnt] : T_dist) std::cout << "  T=" << T << ": " << cnt << "\n";

    // ── standard b₀ filter for unimodular ──
    auto is_standard_b0 = [](const BlowdownResult& bd) -> bool {
        auto& M = bd.extended_IF;
        int n = M.rows();
        int b0 = n - 1;
        for (int i = 0; i < b0; i++)
            if (M(i, b0) != M(i, i) + 2) return false;
        return true;
    };

    // ── Split unimodular / non-unimodular ──
    std::vector<ExtResult*> uni, nonuni;
    for (auto& r : filtered) {
        if (std::abs(r.sig.det) == 1) {
            auto bd = blowdown(r.final_IF);
            if (is_standard_b0(bd)) uni.push_back(&r);
        } else {
            nonuni.push_back(&r);
        }
    }
    std::cout << "Unimodular (std b0): " << uni.size() << "\n";
    std::cout << "Non-unimodular: " << nonuni.size() << "\n";

    // ── LaTeX preamble ──
    auto write_preamble = [](std::ofstream& tex) {
        tex << "\\documentclass[10pt]{article}\n"
            << "\\usepackage[a4paper,margin=0.6in]{geometry}\n"
            << "\\usepackage{longtable}\n\\usepackage{booktabs}\n"
            << "\\usepackage{amsmath,amssymb}\n"
            << "\\usepackage[dvipsnames]{xcolor}\n"
            << "\\usepackage{array}\n"
            << "\\setlength{\\parindent}{0pt}\n"
            << "\\setlength{\\parskip}{2pt}\n"
            << "\\begin{document}\n\n";
    };

    // ── Unimodular output ──
    {
        std::ofstream tex(("sugra_v6_unimodular" + suffix + ".tex").c_str());
        write_preamble(tex);
        tex << "\\section*{6D SUGRA Bases --- Multi-target External ($C^2 \\in \\{-5,\\ldots,-12\\}$, Unimodular)}\n\n";
        tex << "External attaches to $(-1)$ curves only (int.~num.~$= 1$). ";
        tex << "\\textcolor{red}{Red} = external (single-target). ";
        tex << "\\textcolor{blue}{Blue} = target $(-1)$ curves (multi-target).\n\n";
        tex << "Total: " << uni.size() << " bases. ";
        tex << "Filters: NHC, $c(k) \\leq 8$, $H_n \\geq 0$, standard $b_0$.\n\n";

        // T table
        std::map<int, int> T_count;
        for (auto r : uni) T_count[r->anomaly.T]++;
        tex << "\\begin{center}\n\\begin{tabular}{c|c}\n\\toprule\n$T$ & Count \\\\\n\\midrule\n";
        for (auto& [T, cnt] : T_count) tex << T << " & " << cnt << " \\\\\n";
        tex << "\\bottomrule\n\\end{tabular}\n\\end{center}\n\n";

        // Entries
        tex << "\\newpage\n";
        int cur_T = -1;
        std::string cur_base_key = "";
        for (auto rp : uni) {
            auto& r = *rp;
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
                    << if_to_latex(r.base_IF, -2, get_blocks(r.catalog_id))
                    << "$}\\smallskip\n\n";
            }
            int delta = r.anomaly.H_charged - r.anomaly.V + 29 * r.anomaly.T - 273;
            auto bd = blowdown(r.final_IF);
            tex << "$" << render_ext_latex(r, get_blocks(r.catalog_id)) << "$";
            if (r.num_targets() > 1)
                tex << " $\\scriptstyle [" << r.num_targets() << "\\text{-tgt}]$";
            tex << " \\quad $\\Delta=" << delta
                << ",\\; (n_+,n_-,n_0)=(" << r.sig.sig_pos
                << "," << r.sig.sig_neg << "," << r.sig.sig_zero << ")$"
                << " \\quad $\\to " << bd.surface << "$\n\n";
        }
        tex << "\\end{document}\n";
        tex.close();
        std::cout << "Written " << uni.size() << " to sugra_v6_unimodular" << suffix << ".tex\n";
    }

    // ── Non-unimodular output ──
    {
        // Anomaly rejection filter: Δ_crit > 0 → reject
        struct NonUniEntry {
            ExtResult* r;
            CritInfo ci;
        };
        std::vector<NonUniEntry> entries;
        int rejected = 0;
        for (auto r : nonuni) {
            CritInfo ci = compute_crit_info(r->final_IF);
            int delta = r->anomaly.H_charged - r->anomaly.V + 29 * r->anomaly.T - 273;
            int delta_crit = delta + 29 * (ci.T_crit_int - r->anomaly.T);
            if (delta_crit > 0) { rejected++; continue; }
            entries.push_back({r, ci});
        }
        std::cout << "Non-unimodular: " << nonuni.size() << " total, "
                  << rejected << " anomaly-rejected, "
                  << entries.size() << " pass\n";

        std::ofstream tex(("sugra_v6_nonuni" + suffix + ".tex").c_str());
        write_preamble(tex);
        tex << "\\section*{6D SUGRA Bases --- Multi-target External ($C^2 \\in \\{-5,\\ldots,-12\\}$, Non-unimodular)}\n\n";
        tex << "Total: " << entries.size() << " bases (after $\\Delta_{\\text{crit}} \\leq 0$ filter).\n\n";

        std::map<int, int> T_count;
        for (auto& e : entries) T_count[e.r->anomaly.T]++;
        tex << "\\begin{center}\n\\begin{tabular}{c|c}\n\\toprule\n$T$ & Count \\\\\n\\midrule\n";
        for (auto& [T, cnt] : T_count) tex << T << " & " << cnt << " \\\\\n";
        tex << "\\bottomrule\n\\end{tabular}\n\\end{center}\n\n";

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
                    << if_to_latex(r.base_IF, -2, get_blocks(r.catalog_id))
                    << "$}\\smallskip\n\n";
            }
            int delta = r.anomaly.H_charged - r.anomaly.V + 29 * r.anomaly.T - 273;
            tex << "$" << render_ext_latex(r, get_blocks(r.catalog_id)) << "$";
            if (r.num_targets() > 1)
                tex << " $\\scriptstyle [" << r.num_targets() << "\\text{-tgt}]$";
            tex << " \\quad $\\Delta=" << delta
                << ",\\; |\\det|=" << std::abs(r.sig.det)
                << ",\\; (n_+,n_-,n_0)=(" << r.sig.sig_pos
                << "," << r.sig.sig_neg << "," << r.sig.sig_zero << ")$\n\n";

            tex << "\\hspace{1em}$T_{\\text{crit}}=";
            if (e.ci.exact) tex << e.ci.T_crit_int;
            else tex << std::fixed << std::setprecision(2) << e.ci.T_crit
                     << " \\to " << e.ci.T_crit_int;
            tex << ",\\; (n_+,n_-,n_0)_{\\text{crit}}=("
                << e.ci.crit_npos << "," << e.ci.crit_nneg << "," << e.ci.crit_nzero << ")$\n\n";

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
        std::cout << "Written " << entries.size() << " to sugra_v6_nonuni" << suffix << ".tex\n";
    }

    // ── Extended matrix dump ──
    {
        std::string fname = "sugra_v6_ext_matrices" + suffix + ".txt";
        std::ofstream out(fname);
        out << "# Extended intersection matrices for v6 multi-target externals\n";
        out << "# Format: header line then N×N matrix\n";
        out << "# Header: T  det  sig+  sig-  sig0  N  ext_si  num_targets  type\n";
        out << "# Total: " << filtered.size() << " bases\n\n";
        int count = 0;
        for (auto& r : filtered) {
            auto ext = build_b0Q_matrix(r.final_IF);
            int n = ext.rows();
            std::string type = (std::abs(r.sig.det) == 1) ? "U" : "N";
            out << "# " << count++ << "\n";
            out << r.anomaly.T << " " << r.sig.det << " "
                << r.sig.sig_pos << " " << r.sig.sig_neg << " " << r.sig.sig_zero
                << " " << n << " " << r.ext_si << " " << r.num_targets()
                << " " << type << "\n";
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (j > 0) out << " ";
                    out << ext(i, j);
                }
                out << "\n";
            }
            out << "\n";
        }
        out.close();
        std::cout << "Written " << count << " extended matrices to " << fname << "\n";
    }

    return 0;
}
