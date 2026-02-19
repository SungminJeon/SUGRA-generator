// gen_sugra_v3.cpp
// Sidelink-External Pattern Reference for Unimodular SUGRA Bases
//
// For each sidelink type, shows:
//   1. Curve structure (self-intersections)
//   2. Which curves receive externals, and which ext values
//   3. Frequency counts
//
// Compile:
//   g++ -std=c++17 -O2 -o gen_sugra_v3 gen_sugra_v3.cpp Tensor.C \
//       TopoLineCompact_enhanced.cpp Topology_enhanced.cpp \
//       -I/usr/include/eigen3 -I.
//
// Usage:
//   ./gen_sugra_v3 unified.cat [T_max=193] [T_min=0]

#include "sugra_generator.h"
#include <fstream>
#include <cstdlib>
#include <set>
#include <iomanip>
#include <numeric>
#include <map>
#include <algorithm>

// ============================================================================
// Block info (from v2)
// ============================================================================

enum class BlockType { Node, InteriorLink, SideLink, Instanton, External, Unknown };
struct BlockInfo { BlockType type; int param; int start, end; };

struct TopoBlockMap {
    std::vector<BlockInfo> blocks;
    std::vector<int> curve_to_block;
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
        add(b.kind == LKind::g ? BlockType::Node : BlockType::InteriorLink, b.param, sp.kind, sp.param);
    }
    for (auto& s : topo.side_links) add(BlockType::SideLink, s.param, Kind::SideLink, s.param);
    for (auto& i : topo.instantons) add(BlockType::Instanton, i.param, Kind::SideLink, i.param);
    for (auto& e : topo.externals) add(BlockType::External, e.param, Kind::External, e.param);
    result.curve_to_block.resize(offset, -1);
    for (int b = 0; b < (int)result.blocks.size(); b++)
        for (int c = result.blocks[b].start; c < result.blocks[b].end; c++)
            result.curve_to_block[c] = b;
    return result;
}

// ============================================================================
// Quiver rendering for sidelink IF
// ============================================================================

std::string sl_to_quiver(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    if (n == 0) return "\\emptyset";
    if (n == 1) return std::to_string(-IF(0,0));
    
    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            if (IF(i,j) != 0) { adj[i].push_back(j); adj[j].push_back(i); }
    
    int start = 0;
    for (int i = 0; i < n; i++) if ((int)adj[i].size() == 1) { start = i; break; }
    
    std::vector<bool> visited(n, false);
    
    std::function<std::string(int)> walk = [&](int cur) -> std::string {
        visited[cur] = true;
        int si = -IF(cur, cur);
        std::string label = std::to_string(si);
        
        std::vector<int> children;
        for (int nb : adj[cur]) if (!visited[nb]) children.push_back(nb);
        
        if (children.empty()) return label;
        if (children.size() == 1) return label + walk(children[0]);
        
        std::string main_part = walk(children[0]);
        std::string branches;
        for (size_t i = 1; i < children.size(); i++) {
            if (!branches.empty()) branches += ",";
            branches += walk(children[i]);
        }
        return "\\overset{" + branches + "}{" + label + "}" + main_part;
    };
    
    return walk(start);
}

// Render sidelink+external quiver with external curve in red
std::string sl_to_quiver_ext(const Eigen::MatrixXi& IF, int ext_idx) {
    int n = IF.rows();
    if (n == 0) return "\\emptyset";
    if (n == 1) return std::to_string(-IF(0,0));
    
    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            if (IF(i,j) != 0) { adj[i].push_back(j); adj[j].push_back(i); }
    
    // Start from degree-1 node that is not the external
    int start = -1;
    for (int i = 0; i < n; i++)
        if ((int)adj[i].size() == 1 && i != ext_idx) { start = i; break; }
    if (start < 0)
        for (int i = 0; i < n; i++)
            if ((int)adj[i].size() == 1) { start = i; break; }
    if (start < 0) start = 0;
    
    std::vector<bool> visited(n, false);
    
    auto make_label = [&](int cur) -> std::string {
        int si = -IF(cur, cur);
        std::string si_str = std::to_string(si);
        if (cur == ext_idx) return "\\textcolor{red}{" + si_str + "}";
        return si_str;
    };
    
    std::function<std::string(int)> walk = [&](int cur) -> std::string {
        visited[cur] = true;
        std::string label = make_label(cur);
        
        std::vector<int> children;
        for (int nb : adj[cur]) if (!visited[nb]) children.push_back(nb);
        
        if (children.empty()) return label;
        if (children.size() == 1) return label + walk(children[0]);
        
        // Branch: prefer non-external as main
        int main_child = -1;
        std::vector<int> branch_children;
        for (int c : children) {
            if (main_child < 0 && c != ext_idx) main_child = c;
            else branch_children.push_back(c);
        }
        if (main_child < 0) { main_child = children[0]; branch_children.assign(children.begin()+1, children.end()); }
        
        std::string branches;
        for (int i = (int)branch_children.size()-1; i >= 0; i--) {
            if (!branches.empty()) branches += ",";
            branches += walk(branch_children[i]);
        }
        std::string main_part = walk(main_child);
        if (!branches.empty())
            return "\\overset{" + branches + "}{" + label + "}" + main_part;
        return label + main_part;
    };
    
    return walk(start);
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
    
    // Block maps
    std::map<int, TopoBlockMap> bmaps;
    for (auto& e : catalog)
        if (!e.is_nn()) { auto bm = get_block_map(e); if (!bm.blocks.empty()) bmaps[e.id] = std::move(bm); }
    
    // Generate
    SUGRAConfig config;
    config.T_max = 193; config.check_nhc = true; config.check_anomaly = true;
    config.check_determinant = false;
    config.catalog_T_min = T_min; config.catalog_T_max = T_max;
    config.verbose = true;
    
    auto results = generate_sugra(catalog, config);
    std::cout << "Total (before dedup): " << results.size() << "\n";
    
    // Dedup
    auto dk = [](const Eigen::MatrixXi& IF) {
        int n = IF.rows();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
        auto ev = es.eigenvalues();
        std::ostringstream ss; ss << std::fixed << std::setprecision(8) << "E:";
        for (int i = 0; i < ev.size(); i++) { double v = ev(i); if (std::abs(v)<1e-10) v=0; ss << v << ";"; }
        std::vector<std::tuple<int,int,int>> edges;
        for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++)
            if (IF(i,j)!=0) { int a=IF(i,i),b=IF(j,j); if(a>b) std::swap(a,b); edges.push_back({a,b,IF(i,j)}); }
        std::sort(edges.begin(), edges.end());
        ss << "|E:"; for (auto& [a,b,w] : edges) ss << "(" << a << "," << b << "," << w << ")";
        return ss.str();
    };
    std::set<std::string> seen;
    std::vector<SUGRAResult> filtered;
    for (auto& r : results) { auto k = dk(r.final_IF); if (seen.insert(k).second) filtered.push_back(std::move(r)); }
    std::cout << "After dedup: " << filtered.size() << "\n";
    
    // ========================================================================
    // Collect sidelink attachment patterns (unimodular Nk only)
    // ========================================================================
    
    struct AttachKey {
        int local_idx;
        int target_si;
        int ext_si;
        int int_num;
        bool operator<(const AttachKey& o) const {
            return std::tie(local_idx, ext_si, int_num) < std::tie(o.local_idx, o.ext_si, o.int_num);
        }
    };
    
    std::map<int, std::map<AttachKey, int>> sl_patterns;
    std::map<int, std::map<AttachKey, int>> inst_patterns;
    std::map<std::string, int> attach_by_type;
    int uni_nk = 0;
    
    for (auto& r : filtered) {
        if (r.catalog_type == "NN") continue;
        if (std::abs(r.sig.det) != 1) continue;
        uni_nk++;
        
        if (r.externals.empty()) continue;
        auto it = bmaps.find(r.catalog_id);
        if (it == bmaps.end()) continue;
        auto& bmap = it->second;
        
        int tid = r.externals[0].target_idx;
        int bid = bmap.block_id(tid);
        if (bid < 0) continue;
        auto& blk = bmap.blocks[bid];
        
        std::string tname;
        switch (blk.type) {
            case BlockType::Node: tname = "Node"; break;
            case BlockType::InteriorLink: tname = "Link"; break;
            case BlockType::SideLink: tname = "SideLink"; break;
            case BlockType::Instanton: tname = "Instanton"; break;
            case BlockType::External: tname = "SideLink"; break;  // catalog externals are sidelinks
            default: tname = "Other"; break;
        }
        attach_by_type[tname]++;
        
        int local_idx = tid - blk.start;
        int target_si = r.base_IF(tid, tid);
        AttachKey ak{local_idx, target_si, r.externals[0].ext_si, r.externals[0].int_num};
        
        if (blk.type == BlockType::SideLink || blk.type == BlockType::External)
            sl_patterns[blk.param][ak]++;
        else if (blk.type == BlockType::Instanton)
            inst_patterns[blk.param][ak]++;
    }
    
    std::cout << "Unimodular Nk: " << uni_nk << "\n";
    for (auto& [t, c] : attach_by_type) std::cout << "  → " << t << ": " << c << "\n";
    
    // ========================================================================
    // LaTeX output
    // ========================================================================
    
    std::string outfile = "sugra_v3_sidelink_ref" + suffix + ".tex";
    std::ofstream tex(outfile);
    
    tex << "\\documentclass[10pt]{article}\n"
        << "\\usepackage[a4paper,margin=0.6in]{geometry}\n"
        << "\\usepackage{longtable,booktabs,amsmath,amssymb,xcolor,array}\n"
        << "\\setlength{\\parindent}{0pt}\n"
        << "\\setlength{\\parskip}{3pt}\n"
        << "\\begin{document}\n\n"
        << "\\section*{Sidelink--External Attachment Patterns (Unimodular)}\n\n"
        << "From " << uni_nk << " unimodular Nk SUGRA bases"
        << " ($T_{\\text{LST}} \\in [" << T_min << ", " << T_max << "]$).\n\n";
    
    // Summary
    tex << "\\subsection*{Attachment target summary}\n"
        << "\\begin{tabular}{l|r}\n\\toprule\nTarget block type & Count \\\\\n\\midrule\n";
    for (auto& [t, c] : attach_by_type)
        tex << t << " & " << c << " \\\\\n";
    tex << "\\midrule\nTotal & " << uni_nk << " \\\\\n\\bottomrule\n\\end{tabular}\n\n";
    
    // Per sidelink type
    tex << "\\subsection*{Sidelink attachment detail}\n\n"
        << "For each sidelink S($p$): base quiver, then each realized external shown as "
        << "\\textcolor{red}{red} curve attached to the sidelink.\n\n";
    
    for (auto& [sl_param, patterns] : sl_patterns) {
        Tensor t = build_tensor(s(sl_param));
        auto sl_IF = t.GetIntersectionForm();
        bool is_catalog_ext = false;
        if (sl_IF.rows() == 0) {
            t = build_tensor(e(sl_param));
            sl_IF = t.GetIntersectionForm();
            is_catalog_ext = true;
        }
        int ncurves = sl_IF.rows();
        std::string quiver = sl_to_quiver(sl_IF);
        
        tex << "\\medskip\\hrule\\medskip\n"
            << "\\noindent\\textbf{S(" << sl_param << ")}";
        if (is_catalog_ext) tex << " (catalog extra)";
        tex << " --- " << ncurves << " curves: $" << quiver << "$\n\n";
        
        // Collect unique (local_idx, ext_si, int_num) patterns
        struct PatKey {
            int local_idx, ext_si, int_num;
            bool operator<(const PatKey& o) const {
                if (local_idx != o.local_idx) return local_idx < o.local_idx;
                if (ext_si != o.ext_si) return ext_si < o.ext_si;
                return int_num < o.int_num;
            }
        };
        std::map<PatKey, int> unique_patterns;
        for (auto& [ak, cnt] : patterns)
            unique_patterns[{ak.local_idx, ak.ext_si, ak.int_num}] += cnt;
        
        for (auto& [pk, cnt] : unique_patterns) {
            // Build sidelink + external IF
            Eigen::MatrixXi aug_IF = attach_curve(sl_IF, pk.local_idx, pk.ext_si, pk.int_num);
            int ext_idx = aug_IF.rows() - 1;  // last curve is the attached external
            
            // Render with red external
            std::string aug_quiver = sl_to_quiver_ext(aug_IF, ext_idx);
            
            tex << "$" << aug_quiver << "$";
            if (pk.int_num > 1) tex << " $\\scriptstyle [\\#=" << pk.int_num << "]$";
            tex << " \\quad ($\\times " << cnt << "$)\n\n";
        }
        
        int total = 0;
        for (auto& [ak, cnt] : patterns) total += cnt;
        tex << "\\smallskip Total: " << total << "\n\n";
    }
    
    // Instanton section
    if (!inst_patterns.empty()) {
        tex << "\\newpage\n\\subsection*{Instanton attachment detail}\n\n";
        for (auto& [ip, patterns] : inst_patterns) {
            Tensor t = build_tensor(s(ip));
            auto inst_IF = t.GetIntersectionForm();
            int ncurves = inst_IF.rows();
            
            tex << "\\medskip\\hrule\\medskip\n"
                << "\\noindent\\textbf{I(" << ip << ")} --- "
                << ncurves << " curves: $" << sl_to_quiver(inst_IF) << "$\n\n";
            
            struct PatKey {
                int local_idx, ext_si, int_num;
                bool operator<(const PatKey& o) const {
                    if (local_idx != o.local_idx) return local_idx < o.local_idx;
                    if (ext_si != o.ext_si) return ext_si < o.ext_si;
                    return int_num < o.int_num;
                }
            };
            std::map<PatKey, int> unique_patterns;
            for (auto& [ak, cnt] : patterns)
                unique_patterns[{ak.local_idx, ak.ext_si, ak.int_num}] += cnt;
            
            for (auto& [pk, cnt] : unique_patterns) {
                Eigen::MatrixXi aug_IF = attach_curve(inst_IF, pk.local_idx, pk.ext_si, pk.int_num);
                int ext_idx = aug_IF.rows() - 1;
                std::string aug_quiver = sl_to_quiver_ext(aug_IF, ext_idx);
                
                tex << "$" << aug_quiver << "$";
                if (pk.int_num > 1) tex << " $\\scriptstyle [\\#=" << pk.int_num << "]$";
                tex << " \\quad ($\\times " << cnt << "$)\n\n";
            }
            
            int total = 0;
            for (auto& [ak, cnt] : patterns) total += cnt;
            tex << "\\smallskip Total: " << total << "\n\n";
        }
    }
    
    tex << "\\end{document}\n";
    tex.close();
    std::cout << "Written to " << outfile << "\n";
    return 0;
}
