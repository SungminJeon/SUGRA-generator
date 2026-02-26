// sig_test.cpp
// Extended matrix signature test for non-unimodular SUGRA bases
//
// For each non-unimodular base:
//   - Show quiver (v2 notation), det(IF), initial extended sig
//   - Find T where det(extended) changes sign → n_+ drops 2→1
//
// Output: LaTeX PDF with quiver topology + physics data
//
// Usage:
//   ./sig_test unified.cat [T_max=10] [T_min=0]

#include "sugra_generator.h"
#include <fstream>
#include <cstdlib>
#include <set>
#include <iomanip>
#include <map>
#include <numeric>

// ============================================================================
// Block map (from v2)
// ============================================================================

enum class BlockType { Node, InteriorLink, SideLink, Instanton, External, Unknown };
struct BlockInfo { BlockType type; int param; int start, end; };

struct TopoBlockMap {
    std::vector<BlockInfo> blocks;
    std::vector<int> curve_to_block;
    bool is_node(int c) const {
        return c >= 0 && c < (int)curve_to_block.size()
            && blocks[curve_to_block[c]].type == BlockType::Node;
    }
    bool is_interior_link(int c) const {
        return c >= 0 && c < (int)curve_to_block.size()
            && blocks[curve_to_block[c]].type == BlockType::InteriorLink;
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
// Quiver LaTeX (v2)
// ============================================================================

std::string gauge_algebra_latex(int p) {
    switch (p) {
        case 5:  return "\\mathfrak{f}_4";
        case 6:  return "\\mathfrak{e}_6";
        case 7:  return "\\mathfrak{e}_7'";
        case 8:  return "\\mathfrak{e}_7";
        case 12: return "\\mathfrak{e}_8";
        default: return std::to_string(p);
    }
}

std::string interior_link_label(int p) {
    std::string s = std::to_string(p);
    if (s == "331") return "\\overset{3,3}{\\bigcirc}";
    if (s.length() >= 2) {
        std::string d;
        for (size_t i = 0; i < s.length(); i++) {
            if (i) d += ",";
            d += s[i];
        }
        return "\\overset{" + d + "}{\\otimes}";
    }
    return "\\overset{" + s + "}{\\otimes}";
}

std::string if_to_latex(const Eigen::MatrixXi& IF, int ext_idx = -1,
                        const TopoBlockMap& bmap = {}) {
    int n = IF.rows();
    if (n == 0) return "\\emptyset";
    if (n == 1) {
        if (!bmap.blocks.empty() && bmap.is_node(0))
            return gauge_algebra_latex(bmap.block_of(0).param);
        return std::to_string(-IF(0, 0));
    }
    if (ext_idx == -1) ext_idx = n - 1;
    bool hb = !bmap.blocks.empty();

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

    std::vector<bool> vis(n, false);

    auto curve_label = [&](int c) -> std::string {
        int si = -IF(c, c);
        std::string s = std::to_string(si);
        if (si >= 10) s = "(" + s + ")";
        if (c == ext_idx) return "\\textcolor{red}{" + s + "}";
        return s;
    };

    auto skip_block = [&](int c) -> std::vector<int> {
        auto& blk = bmap.block_of(c);
        for (int i = blk.start; i < blk.end; i++) vis[i] = true;
        std::vector<int> exits;
        for (int i = blk.start; i < blk.end; i++)
            for (auto& [nb, w] : adj[i])
                if (!vis[nb]) exits.push_back(nb);
        std::sort(exits.begin(), exits.end());
        exits.erase(std::unique(exits.begin(), exits.end()), exits.end());
        return exits;
    };

    std::function<std::string(int, int)> walk = [&](int cur, int from) -> std::string {
        vis[cur] = true;

        // Interior link compression
        if (hb && bmap.is_interior_link(cur)
            && bmap.block_of(cur).end - bmap.block_of(cur).start > 1) {
            std::string label = interior_link_label(bmap.block_of(cur).param);
            auto exits = skip_block(cur);
            if (exits.empty()) return label;
            if (exits.size() == 1) return label + walk(exits[0], cur);
            int mc = -1;
            std::vector<int> bc;
            for (int x : exits) {
                if (mc < 0 && x != ext_idx) mc = x;
                else bc.push_back(x);
            }
            if (mc < 0) { mc = exits[0]; bc.assign(exits.begin() + 1, exits.end()); }
            std::string brs;
            for (int i = (int)bc.size() - 1; i >= 0; i--) {
                if (!brs.empty()) brs += ",";
                brs += walk(bc[i], cur);
            }
            std::string mp = walk(mc, cur);
            if (!brs.empty()) return "\\overset{" + brs + "}{" + label + "}" + mp;
            return label + mp;
        }

        // Label
        std::string label;
        if (hb && bmap.is_node(cur))
            label = gauge_algebra_latex(bmap.block_of(cur).param);
        else
            label = curve_label(cur);

        // Children
        std::vector<int> children;
        for (auto& [nb, w] : adj[cur])
            if (!vis[nb]) children.push_back(nb);
        if (children.empty()) return label;
        if (children.size() == 1) return label + walk(children[0], cur);

        // Branch
        int mc = -1;
        std::vector<int> bc;
        for (int x : children) {
            if (mc < 0 && x != ext_idx) mc = x;
            else bc.push_back(x);
        }
        if (mc < 0) {
            mc = children[0];
            bc.clear();
            for (size_t i = 1; i < children.size(); i++) bc.push_back(children[i]);
        }
        std::string brs;
        for (int i = (int)bc.size() - 1; i >= 0; i--) {
            int b = bc[i];
            std::vector<std::string> chain;
            int nd = b;
            while (nd >= 0) {
                vis[nd] = true;
                chain.push_back(curve_label(nd));
                int nx = -1;
                for (auto& [nb, w] : adj[nd])
                    if (!vis[nb]) { nx = nb; break; }
                nd = nx;
            }
            std::reverse(chain.begin(), chain.end());
            std::string br;
            for (auto& l : chain) br += l;
            if (!brs.empty()) brs += ",";
            brs += br;
        }
        std::string mp = walk(mc, cur);
        if (!brs.empty()) return "\\overset{" + brs + "}{" + label + "}" + mp;
        return label + mp;
    };

    return walk(start, -1);
}

// ============================================================================
// Dedup key (v2)
// ============================================================================

std::string dedup_key(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
    auto ev = es.eigenvalues();
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(8) << "E:";
    for (int i = 0; i < ev.size(); i++) {
        double v = ev(i);
        if (std::abs(v) < 1e-10) v = 0.0;
        if (i) ss << ",";
        ss << v;
    }
    struct VD {
        int si, deg;
        std::vector<int> ns;
        bool operator<(const VD& o) const {
            if (si != o.si) return si < o.si;
            if (deg != o.deg) return deg < o.deg;
            return ns < o.ns;
        }
    };
    std::vector<VD> vds(n);
    for (int i = 0; i < n; i++) {
        vds[i].si = IF(i, i);
        vds[i].deg = 0;
        for (int j = 0; j < n; j++)
            if (j != i && IF(i, j) != 0) {
                vds[i].deg++;
                vds[i].ns.push_back(IF(j, j));
            }
        std::sort(vds[i].ns.begin(), vds[i].ns.end());
    }
    std::sort(vds.begin(), vds.end());
    ss << "|V:";
    for (auto& v : vds) {
        ss << "(" << v.si << "," << v.deg << ",{";
        for (int k = 0; k < (int)v.ns.size(); k++) {
            if (k) ss << ",";
            ss << v.ns[k];
        }
        ss << "})";
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
    return ss.str();
}

// ============================================================================
// Extended matrix + T_crit computation
// ============================================================================

struct SigResult {
    // IF data
    int T_sugra;                // sig_neg of curve IF
    int det_IF;                 // det(curve IF)

    // Extended matrix at original T
    int ext_npos, ext_nneg, ext_nzero;
    int det_ext;

    // T_crit: T where det(extended) changes sign
    double T_crit;              // fractional
    int T_crit_floor;           // first integer T where n_+ drops (= ceil or exact)
    bool T_crit_exact;          // T_crit is integer?
};

SigResult compute_sig_result(const Eigen::MatrixXi& IF) {
    SigResult res;
    int n = IF.rows();

    // Curve IF signature
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es0(IF.cast<double>());
    res.T_sugra = 0;
    for (int i = 0; i < es0.eigenvalues().size(); i++)
        if (es0.eigenvalues()(i) < -1e-8) res.T_sugra++;
    res.det_IF = (int)std::round(IF.cast<double>().determinant());

    // Extended matrix at original T
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n + 1, n + 1);
    M.block(0, 0, n, n) = IF.cast<double>();
    for (int i = 0; i < n; i++) {
        M(i, n) = IF(i, i) + 2;
        M(n, i) = IF(i, i) + 2;
    }
    M(n, n) = 9.0 - res.T_sugra;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es1(M);
    res.ext_npos = res.ext_nneg = res.ext_nzero = 0;
    for (int i = 0; i < n + 1; i++) {
        double ev = es1.eigenvalues()(i);
        if (std::abs(ev) < 1e-8) res.ext_nzero++;
        else if (ev > 0) res.ext_npos++;
        else res.ext_nneg++;
    }
    res.det_ext = (int)std::round(M.determinant());

    // T_crit: det(extended) = A * b0^2 + B, b0^2 = 9 - T
    M(n, n) = 0.0;
    double B = M.determinant();
    M(n, n) = 1.0;
    double A = M.determinant() - B;

    if (std::abs(A) > 1e-10) {
        double b0sq_crit = -B / A;
        res.T_crit = 9.0 - b0sq_crit;
        int tc = (int)std::round(res.T_crit);
        res.T_crit_exact = (std::abs(res.T_crit - tc) < 1e-6);
        if (res.T_crit_exact) {
            res.T_crit_floor = tc;
        } else {
            // First integer T where sign changes: ceiling
            res.T_crit_floor = (int)std::ceil(res.T_crit - 1e-9);
        }
    } else {
        res.T_crit = 1e9;
        res.T_crit_exact = false;
        res.T_crit_floor = -1;
    }

    return res;
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv) {
    std::string catalog_file = "unified.cat";
    int T_max = 10, T_min = 0;
    if (argc > 1) catalog_file = argv[1];
    if (argc > 2) T_max = std::atoi(argv[2]);
    if (argc > 3) T_min = std::atoi(argv[3]);

    std::string suffix = "_T" + std::to_string(T_min) + "_" + std::to_string(T_max);

    auto catalog = load_catalog(catalog_file);
    std::cerr << "Catalog: " << catalog.size() << " entries\n";

    // Block maps for Nk
    std::map<int, TopoBlockMap> block_maps;
    for (auto& e : catalog)
        if (!e.is_nn()) {
            auto bm = get_block_map(e);
            if (!bm.blocks.empty()) block_maps[e.id] = std::move(bm);
        }

    // Generate
    SUGRAConfig config;
    config.T_max = 193;
    config.check_nhc = true;
    config.check_anomaly = true;
    config.check_determinant = false;
    config.catalog_T_min = T_min;
    config.catalog_T_max = T_max;
    config.include_dummy = true;
    config.verbose = true;

    auto results = generate_sugra(catalog, config);

    // Dedup
    std::set<std::string> seen;
    std::vector<SUGRAResult> filtered;
    for (auto& r : results) {
        auto k = dedup_key(r.final_IF);
        if (seen.insert(k).second) filtered.push_back(std::move(r));
    }
    std::cerr << "After dedup: " << filtered.size() << "\n";

    // Split
    struct Entry {
        SUGRAResult* r;
        SigResult sig;
    };
    std::vector<Entry> nonuni;
    int n_uni = 0;
    for (auto& r : filtered) {
        if (std::abs(r.sig.det) == 1) {
            n_uni++;
        } else {
            auto sr = compute_sig_result(r.final_IF);
            nonuni.push_back({&r, sr});
        }
    }
    std::cerr << "Uni: " << n_uni << " Non-uni: " << nonuni.size() << "\n";

    // Sort by T, then det
    std::sort(nonuni.begin(), nonuni.end(), [](const Entry& a, const Entry& b) {
        if (a.r->anomaly.T != b.r->anomaly.T) return a.r->anomaly.T < b.r->anomaly.T;
        if (std::abs(a.r->sig.det) != std::abs(b.r->sig.det))
            return std::abs(a.r->sig.det) < std::abs(b.r->sig.det);
        return a.r->catalog_id < b.r->catalog_id;
    });

    auto get_bm = [&](int id) -> const TopoBlockMap& {
        static const TopoBlockMap empty;
        auto it = block_maps.find(id);
        return (it != block_maps.end()) ? it->second : empty;
    };

    // ── Collect statistics ──
    std::map<int, int> npos_dist;                    // ext n_pos → count
    std::map<int, int> gap_dist;                     // T_crit_floor - T_sugra → count
    std::map<int, std::map<int, int>> T_gap;         // T_sugra → gap → count
    std::map<int, std::map<int, int>> ext_gap;       // ext_si → gap → count
    int always1 = 0;

    for (auto& e : nonuni) {
        npos_dist[e.sig.ext_npos]++;
        if (e.sig.ext_npos < 2) { always1++; continue; }
        if (e.sig.T_crit_floor < 0) continue;
        int gap = e.sig.T_crit_floor - e.sig.T_sugra;
        gap_dist[gap]++;
        T_gap[e.sig.T_sugra][gap]++;
        int esi = e.r->externals.empty() ? 0 : e.r->externals[0].ext_si;
        ext_gap[esi][gap]++;
    }

    // ── LaTeX output ──
    std::string outfile = "sig_test" + suffix + ".tex";
    std::ofstream tex(outfile);

    tex << "\\documentclass[9pt]{article}\n"
        << "\\usepackage[a4paper,margin=0.6in]{geometry}\n"
        << "\\usepackage{amsmath,amssymb,booktabs,xcolor,longtable,array,mathtools}\n"
        << "\\setlength{\\parindent}{0pt}\n"
        << "\\setlength{\\parskip}{3pt}\n"
        << "\\begin{document}\n\n";

    // ── Title + explanation ──
    tex << "\\section*{Extended Matrix Signature Test}\n\n"
        << "Non-unimodular SUGRA bases, $T_{\\text{LST}} \\in ["
        << T_min << ",\\," << T_max << "]$ (catalog + dummy A/D-type).\n\n"
        << "Extended matrix: $b_0 \\cdot C_i = C_i^2+2$, $b_0^2 = 9-T$. "
        << "Since curve IF has signature $(1,T)$, Cauchy interlacing gives "
        << "$n_+(\\mathcal{M}_{\\mathrm{ext}}) \\geq 1$ always. "
        << "The $n_+ = 2 \\to 1$ transition occurs at $\\det(\\mathcal{M}_{\\mathrm{ext}})=0$.\n\n"
        << "\\textcolor{red}{Red} = attached external. "
        << "$\\scriptstyle [\\#=n]$ = intersection number.\n\n";

    // ── Summary ──
    tex << "\\subsection*{Summary}\n"
        << "\\begin{center}\\begin{tabular}{lr}\\toprule\n"
        << "Total SUGRA bases & " << filtered.size() << " \\\\\n"
        << "Unimodular ($|\\det|=1$) & " << n_uni << " \\\\\n"
        << "Non-unimodular & " << nonuni.size() << " \\\\\n"
        << "\\quad $n_+(\\text{ext}) = 2$ at original $T$ & " << npos_dist[2] << " \\\\\n"
        << "\\quad $n_+(\\text{ext}) = 1$ at original $T$ & " << always1 << " \\\\\n"
        << "\\bottomrule\\end{tabular}\\end{center}\n\n";

    // ── n_pos distribution ──
    tex << "\\subsection*{$n_+$ Distribution (Extended Matrix at Original $T$)}\n"
        << "\\begin{center}\\begin{tabular}{c|r|r}\\toprule\n"
        << "$n_+$ & Count & \\% \\\\\\midrule\n";
    for (auto& [np, cnt] : npos_dist)
        tex << np << " & " << cnt << " & "
            << std::fixed << std::setprecision(1)
            << (100.0 * cnt / nonuni.size()) << " \\\\\n";
    tex << "\\bottomrule\\end{tabular}\\end{center}\n\n";

    // ── Gap distribution ──
    std::set<int> all_gaps;
    for (auto& [g, c] : gap_dist) all_gaps.insert(g);
    int total_gap = 0;
    for (auto& [g, c] : gap_dist) total_gap += c;

    tex << "\\subsection*{Gap $= T_{\\det=0} - T_{\\text{SUGRA}}$}\n"
        << "\\begin{center}\\begin{tabular}{c|r|r}\\toprule\n"
        << "Gap & Count & \\% \\\\\\midrule\n";
    for (auto& [g, cnt] : gap_dist)
        tex << g << " & " << cnt << " & "
            << std::fixed << std::setprecision(1)
            << (100.0 * cnt / total_gap) << " \\\\\n";
    tex << "\\bottomrule\\end{tabular}\\end{center}\n\n";

    // ── Gap by T ──
    tex << "\\subsection*{Gap by $T_{\\text{SUGRA}}$}\n"
        << "\\begin{center}\n{\\small\\begin{tabular}{c|";
    for (int g : all_gaps) tex << "r";
    tex << "|r}\\toprule\n$T$";
    for (int g : all_gaps) tex << " & " << g;
    tex << " & tot \\\\\\midrule\n";
    for (auto& [T, mp] : T_gap) {
        tex << T;
        int tot = 0;
        for (int g : all_gaps) {
            int c = mp.count(g) ? mp.at(g) : 0;
            tex << " & " << (c ? std::to_string(c) : "");
            tot += c;
        }
        tex << " & " << tot << " \\\\\n";
    }
    tex << "\\bottomrule\\end{tabular}}\\end{center}\n\n";

    // ── Gap by ext type ──
    tex << "\\subsection*{Gap by External Curve}\n"
        << "\\begin{center}\n{\\small\\begin{tabular}{c|";
    for (int g : all_gaps) tex << "r";
    tex << "|r}\\toprule\next";
    for (int g : all_gaps) tex << " & " << g;
    tex << " & tot \\\\\\midrule\n";
    for (auto& [esi, mp] : ext_gap) {
        tex << esi;
        int tot = 0;
        for (int g : all_gaps) {
            int c = mp.count(g) ? mp.at(g) : 0;
            tex << " & " << (c ? std::to_string(c) : "");
            tot += c;
        }
        tex << " & " << tot << " \\\\\n";
    }
    tex << "\\bottomrule\\end{tabular}}\\end{center}\n\n";

    // ── All entries ──
    tex << "\\newpage\n"
        << "\\subsection*{All Non-Unimodular Bases (" << nonuni.size() << ")}\n\n";

    int cur_T = -1;
    for (int idx = 0; idx < (int)nonuni.size(); idx++) {
        auto& e = nonuni[idx];
        auto& r = *e.r;

        if (r.anomaly.T != cur_T) {
            cur_T = r.anomaly.T;
            tex << "\n\\subsubsection*{$T = " << cur_T << "$}\n\n";
        }

        int ext_si = r.externals.empty() ? 0 : r.externals[0].ext_si;
        int int_num = r.externals.empty() ? 0 : r.externals[0].int_num;

        // Quiver
        tex << "$" << if_to_latex(r.final_IF, r.final_IF.rows() - 1,
                                   get_bm(r.catalog_id)) << "$";
        if (int_num > 1) tex << " $\\scriptstyle [\\#=" << int_num << "]$";

        // Data
        tex << " \\quad $\\det=" << e.sig.det_IF << "$";
        tex << ", $(n_+,n_-)=(" << e.sig.ext_npos << "," << e.sig.ext_nneg;
        if (e.sig.ext_nzero) tex << ",0{\\times}" << e.sig.ext_nzero;
        tex << ")$";

        // T_crit
        if (e.sig.ext_npos >= 2 && e.sig.T_crit < 1e8) {
            tex << ", $T_{\\det=0}";
            if (e.sig.T_crit_exact)
                tex << "=" << e.sig.T_crit_floor;
            else
                tex << "\\approx " << std::fixed << std::setprecision(2) << e.sig.T_crit;
            tex << "$";
        } else if (e.sig.ext_npos == 1) {
            tex << ", {\\footnotesize always $n_+=1$}";
        }

        tex << "\n\n";
    }

    tex << "\\end{document}\n";
    tex.close();
    std::cerr << "Written " << outfile << "\n";

    // Compile PDF
    std::string cmd = "pdflatex -interaction=nonstopmode " + outfile + " >/dev/null 2>&1";
    int rc = system(cmd.c_str());
    if (rc == 0)
        std::cerr << "PDF: " << "sig_test" + suffix + ".pdf" << "\n";
    else
        std::cerr << "pdflatex failed (rc=" << rc << "), .tex still available\n";

    return 0;
}
