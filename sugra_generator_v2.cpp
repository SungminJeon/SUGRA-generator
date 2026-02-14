// sugra_generator_v2.cpp
// SUGRA base classification with E8 subalgebra gluing rules
// Enhanced version with NHC-aware gauge algebra tracking
// what should be updated for v3...
// 1. renaming the curve we called "external" as "extra"
// 2. the term "external must be used only for additional BPS generator for SUGRA base. 
// 3. with this in mind, we want to estabilish new "external" curve to our class. 
// 4. with this setting, we have to enhance the class Theory_enhanced, TheoryDB_enhanced, Topology_enhanced 

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <thread>
#include <mutex>
#include <atomic>
#include <queue>
#include <cmath>
#include <set>
#include <map>
#include <algorithm>
#include <functional>

#include "Topology_enhanced.h"
#include "TopologyDB_enhanced.hpp"
#include "TopoLineCompact_enhanced.hpp"
#include "Theory_enhanced.h"
#include "Tensor.h"
#include "anomaly_tables.h"  // Anomaly contribution tables
#include <sstream>
#include <unordered_set>

namespace fs = std::filesystem;

// Import types
using ::NodeRef;
using ::AttachmentPoint;
using ::Spec;
using ::TheoryGraph;

// ============================================================================
// GAUGE ALGEBRA DEFINITIONS
// ============================================================================

enum class GaugeAlgebra {
    none = 0,
    su2, su3, su4, su5, su6, su7, su8,
    g2, so7, so8, so9, so10, so11, so12, so13, so14, so16,
    f4, e6, e7, e8
};

int gauge_dim(GaugeAlgebra g) {
    switch (g) {
        case GaugeAlgebra::none: return 0;
        case GaugeAlgebra::su2: return 3;
        case GaugeAlgebra::su3: return 8;
        case GaugeAlgebra::su4: return 15;
        case GaugeAlgebra::su5: return 24;
        case GaugeAlgebra::su6: return 35;
        case GaugeAlgebra::su7: return 48;
        case GaugeAlgebra::su8: return 63;
        case GaugeAlgebra::g2: return 14;
        case GaugeAlgebra::so7: return 21;
        case GaugeAlgebra::so8: return 28;
        case GaugeAlgebra::so9: return 36;
        case GaugeAlgebra::so10: return 45;
        case GaugeAlgebra::so11: return 55;
        case GaugeAlgebra::so12: return 66;
        case GaugeAlgebra::so13: return 78;
        case GaugeAlgebra::so14: return 91;
        case GaugeAlgebra::so16: return 120;
        case GaugeAlgebra::f4: return 52;
        case GaugeAlgebra::e6: return 78;
        case GaugeAlgebra::e7: return 133;
        case GaugeAlgebra::e8: return 248;
        default: return 0;
    }
}

// Central charge for level-1 WZW: c = dim(g) / (1 + h∨)
// This is the E-string budget constraint: sum of c_g must be ≤ 8
double gauge_central_charge(GaugeAlgebra g) {
    switch (g) {
        case GaugeAlgebra::none: return 0.0;
        // SU(n): dim = n²-1, h∨ = n, c = (n²-1)/(n+1) = n-1
        case GaugeAlgebra::su2: return 1.0;    // 3/3 = 1
        case GaugeAlgebra::su3: return 2.0;    // 8/4 = 2
        case GaugeAlgebra::su4: return 3.0;    // 15/5 = 3
        case GaugeAlgebra::su5: return 4.0;    // 24/6 = 4
        case GaugeAlgebra::su6: return 5.0;    // 35/7 = 5
        case GaugeAlgebra::su7: return 6.0;    // 48/8 = 6
        case GaugeAlgebra::su8: return 7.0;    // 63/9 = 7
        // G2: dim = 14, h∨ = 4, c = 14/5 = 2.8
        case GaugeAlgebra::g2: return 14.0/5.0;
        // SO(2n+1) = B_n: dim = n(2n+1), h∨ = 2n-1
        case GaugeAlgebra::so7: return 21.0/6.0;   // B3: 21/6 = 3.5
        case GaugeAlgebra::so9: return 36.0/8.0;   // B4: 36/8 = 4.5
        case GaugeAlgebra::so11: return 55.0/10.0; // B5: 55/10 = 5.5
        case GaugeAlgebra::so13: return 78.0/12.0; // B6: 78/12 = 6.5
        // SO(2n) = D_n: dim = n(2n-1), h∨ = 2n-2
        case GaugeAlgebra::so8: return 28.0/7.0;   // D4: 28/7 = 4
        case GaugeAlgebra::so10: return 45.0/9.0;  // D5: 45/9 = 5
        case GaugeAlgebra::so12: return 66.0/11.0; // D6: 66/11 = 6
        case GaugeAlgebra::so14: return 91.0/13.0; // D7: 91/13 = 7
        case GaugeAlgebra::so16: return 120.0/15.0;// D8: 120/15 = 8
        // Exceptional
        case GaugeAlgebra::f4: return 52.0/10.0;   // 52/10 = 5.2
        case GaugeAlgebra::e6: return 78.0/13.0;   // 78/13 = 6
        case GaugeAlgebra::e7: return 133.0/19.0;  // 133/19 = 7
        case GaugeAlgebra::e8: return 248.0/31.0;  // 248/31 = 8
        default: return 0.0;
    }
}

const char* gauge_name(GaugeAlgebra g) {
    switch (g) {
        case GaugeAlgebra::none: return "none";
        case GaugeAlgebra::su2: return "su2";
        case GaugeAlgebra::su3: return "su3";
        case GaugeAlgebra::su4: return "su4";
        case GaugeAlgebra::su5: return "su5";
        case GaugeAlgebra::su6: return "su6";
        case GaugeAlgebra::su7: return "su7";
        case GaugeAlgebra::su8: return "su8";
        case GaugeAlgebra::g2: return "g2";
        case GaugeAlgebra::so7: return "so7";
        case GaugeAlgebra::so8: return "so8";
        case GaugeAlgebra::so9: return "so9";
        case GaugeAlgebra::so10: return "so10";
        case GaugeAlgebra::so11: return "so11";
        case GaugeAlgebra::so12: return "so12";
        case GaugeAlgebra::so13: return "so13";
        case GaugeAlgebra::so14: return "so14";
        case GaugeAlgebra::so16: return "so16";
        case GaugeAlgebra::f4: return "f4";
        case GaugeAlgebra::e6: return "e6";
        case GaugeAlgebra::e7: return "e7";
        case GaugeAlgebra::e8: return "e8";
        default: return "unknown";
    }
}

// Simple curve to gauge (for isolated curves / externals)
GaugeAlgebra curve_to_gauge(int self_int) {
    switch (std::abs(self_int)) {
        case 1:  return GaugeAlgebra::none;
        case 2:  return GaugeAlgebra::none;  // isolated -2 has no gauge
        case 3:  return GaugeAlgebra::su3;
        case 4:  return GaugeAlgebra::so8;
        case 5:  return GaugeAlgebra::f4;
        case 6:  return GaugeAlgebra::e6;
        case 7:  return GaugeAlgebra::e7;
        case 8:  return GaugeAlgebra::e7;
        case 9:  return GaugeAlgebra::e8;
        case 10: return GaugeAlgebra::e8;
        case 11: return GaugeAlgebra::e8;
        case 12: return GaugeAlgebra::e8;
        default: return GaugeAlgebra::none;
    }
}

// ============================================================================
// PORT-WISE GAUGE ALGEBRA REGISTRY (NHC-aware)
// ============================================================================

struct PortGaugeInfo {
    int self_int;
    GaugeAlgebra gauge;
    int gauge_rank;
};

using PortGaugeTable = std::map<std::pair<LKind, int>, std::vector<PortGaugeInfo>>;

class SpecGaugeRegistry {
public:
    PortGaugeTable table;
    
    SpecGaugeRegistry() {
        initialize();
    }
    
    void initialize() {
        // =====================================================
        // NODES (LKind::g) - Single curves
        // =====================================================
        add_node(1, {{-1, GaugeAlgebra::none, 0}});
        add_node(2, {{-2, GaugeAlgebra::none, 0}});
        add_node(3, {{-3, GaugeAlgebra::su3, 2}});
        add_node(4, {{-4, GaugeAlgebra::so8, 4}});
        add_node(5, {{-5, GaugeAlgebra::f4, 4}});
        add_node(6, {{-6, GaugeAlgebra::e6, 6}});
        add_node(7, {{-7, GaugeAlgebra::e7, 7}});
        add_node(8, {{-8, GaugeAlgebra::e7, 7}});
        add_node(12, {{-12, GaugeAlgebra::e8, 8}});
        
        // =====================================================
        // INTERIOR LINKS (LKind::L)
        // =====================================================
        
        // (1,1): -1
        add_interior(11, {
            {-1, GaugeAlgebra::none, 0}
        });
        
        // (2,2): -1 -3 -1
        add_interior(22, {
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0}
        });
        
        // (3,3): -1 -2 -3 -2 -1
        add_interior(33, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::su2, 1},
            {-3, GaugeAlgebra::so7, 8},
            {-2, GaugeAlgebra::su2, 1},
            {-1, GaugeAlgebra::none, 0}
        });
        
        // (3,3)' with -5: -1 -3 -1 -5 -1 -3 -1
        add_interior(331, {
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0},
            {-5, GaugeAlgebra::f4, 4},
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0}
        });
        
        // (4,4): -1 -2 -3 -1 -5 -1 -3 -2 -1
        add_interior(44, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::su2, 1},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0},
            {-5, GaugeAlgebra::f4, 4},
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-2, GaugeAlgebra::su2, 1},
            {-1, GaugeAlgebra::none, 0}
        });
        
        // (5,5): -1 -2 -2 -3 -1 -5 -1 -3 -2 -2 -1
        add_interior(55, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::sp1, 1},
            {-3, GaugeAlgebra::g2, 2},
            {-1, GaugeAlgebra::none, 0},
            {-5, GaugeAlgebra::f4, 4},
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::g2, 2},
            {-2, GaugeAlgebra::sp1, 1},
            {-2, GaugeAlgebra::none, 0},
            {-1, GaugeAlgebra::none, 0}
        });
        
        // (2,3): -1 -2 -3 -1
        add_interior(23, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::su2, 1},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0}
        });
        add_interior(32, {
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-2, GaugeAlgebra::su2, 1},
            {-1, GaugeAlgebra::none, 0}
        });
        
        // (2,4): -1 -2 -2 -3 -1
        add_interior(24, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::su2, 1},
            {-2, GaugeAlgebra::su2, 1},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0}
        });
        add_interior(42, {
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-2, GaugeAlgebra::su2, 1},
            {-2, GaugeAlgebra::su2, 1},
            {-1, GaugeAlgebra::none, 0}
        });
        
        // (3,4): -1 -3 -1 -5 -1 -3 -2 -1
        add_interior(34, {
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0},
            {-5, GaugeAlgebra::f4, 4},
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-2, GaugeAlgebra::su2, 1},
            {-1, GaugeAlgebra::none, 0}
        });
        add_interior(43, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::su2, 1},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0},
            {-5, GaugeAlgebra::f4, 4},
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0}
        });
        
        // (3,5): -1 -3 -1 -5 -1 -3 -2 -2 -1
        add_interior(35, {
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0},
            {-5, GaugeAlgebra::f4, 4},
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-2, GaugeAlgebra::su2, 1},
            {-2, GaugeAlgebra::none, 0},
            {-1, GaugeAlgebra::none, 0}
        });
        add_interior(53, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::su2, 1},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0},
            {-5, GaugeAlgebra::f4, 4},
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0}
        });
        
        // (4,5): -1 -2 -3 -1 -5 -1 -3 -2 -2 -1
        add_interior(45, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::su2, 1},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0},
            {-5, GaugeAlgebra::f4, 4},
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-2, GaugeAlgebra::su2, 1},
            {-2, GaugeAlgebra::none, 0},
            {-1, GaugeAlgebra::none, 0}
        });
        add_interior(54, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::su2, 1},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0},
            {-5, GaugeAlgebra::f4, 4},
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-2, GaugeAlgebra::su2, 1},
            {-1, GaugeAlgebra::none, 0}
        });
        
        // =====================================================
        // SIDELINKS (LKind::S)
        // TODO: Fill in based on your classification
        // =====================================================
        
        // === Instantons (88x series) - just -1 curves ===
        add_sidelink(1, {{-1, GaugeAlgebra::none, 0}});
        add_sidelink(882, {
            {-1, GaugeAlgebra::none, 0},
            {-1, GaugeAlgebra::none, 0}
        });
        add_sidelink(883, {
            {-1, GaugeAlgebra::none, 0},
            {-1, GaugeAlgebra::none, 0},
            {-1, GaugeAlgebra::none, 0}
        });
        // ... add more 88x
        
        // === Alkali 1 links (9x series) ===
        // 91: -1 -2 -3
        add_sidelink(91, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::su2, 1},
            {-3, GaugeAlgebra::su3, 2}
        });
        
        // 92: -1 -2 -3 -2
        add_sidelink(92, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::su2, 1},
            {-3, GaugeAlgebra::su3, 2},
            {-2, GaugeAlgebra::su2, 1}
        });
        
        // 912: -1 -2 (isolated, no NHC)
        add_sidelink(912, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::none, 0}
        });
        
        // 915: -1 -2 -2 (no -3, no gauge)
        add_sidelink(915, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::none, 0}
        });
        
        // 96: g2 structure -3 -2
        add_sidelink(96, {
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::g2, 2},
            {-2, GaugeAlgebra::none, 0}
        });
        
        // 993: so(7) from -2 -3 -2
        add_sidelink(993, {
            {-1, GaugeAlgebra::none, 0},
            {-2, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::so7, 3},
            {-2, GaugeAlgebra::none, 0}
        });
        
        // === Links with -5 curve ===
        // 918: -1 -5 -1
        add_sidelink(918, {
            {-1, GaugeAlgebra::none, 0},
            {-5, GaugeAlgebra::f4, 4},
            {-1, GaugeAlgebra::none, 0}
        });
        
        // 919: -1 -3 -1 -5
        add_sidelink(919, {
            {-1, GaugeAlgebra::none, 0},
            {-3, GaugeAlgebra::su3, 2},
            {-1, GaugeAlgebra::none, 0},
            {-5, GaugeAlgebra::f4, 4}
        });
        
        // =====================================================
        // EXTERNALS (LKind::E) - single gauge curves
        // =====================================================
        add_external(1, {{-1, GaugeAlgebra::none, 0}});
        add_external(2, {{-2, GaugeAlgebra::none, 0}});
        add_external(3, {{-3, GaugeAlgebra::su3, 2}});
        add_external(4, {{-4, GaugeAlgebra::so8, 4}});
        add_external(5, {{-5, GaugeAlgebra::f4, 4}});
        add_external(6, {{-6, GaugeAlgebra::e6, 6}});
        add_external(7, {{-7, GaugeAlgebra::e7, 7}});
        add_external(8, {{-8, GaugeAlgebra::e7, 7}});
        add_external(12, {{-12, GaugeAlgebra::e8, 8}});
    }
    
    void add_node(int param, std::vector<PortGaugeInfo> ports) {
        table[{LKind::g, param}] = std::move(ports);
    }
    
    void add_interior(int param, std::vector<PortGaugeInfo> ports) {
        table[{LKind::L, param}] = std::move(ports);
    }
    
    void add_sidelink(int param, std::vector<PortGaugeInfo> ports) {
        table[{LKind::S, param}] = std::move(ports);
    }
    
    void add_instanton(int param, std::vector<PortGaugeInfo> ports) {
        table[{LKind::I, param}] = std::move(ports);
    }
    
    void add_external(int param, std::vector<PortGaugeInfo> ports) {
        table[{LKind::E, param}] = std::move(ports);
    }
    
    GaugeAlgebra get_port_gauge(LKind kind, int param, int port_idx) const {
        auto it = table.find({kind, param});
        if (it == table.end()) {
            return GaugeAlgebra::none;
        }
        
        const auto& ports = it->second;
        if (port_idx < 0 || port_idx >= (int)ports.size()) {
            return GaugeAlgebra::none;
        }
        
        return ports[port_idx].gauge;
    }
    
    std::vector<PortGaugeInfo> get_all_ports(LKind kind, int param) const {
        auto it = table.find({kind, param});
        if (it == table.end()) {
            return {};
        }
        return it->second;
    }
    
    bool has_entry(LKind kind, int param) const {
        return table.find({kind, param}) != table.end();
    }
};

// Global registry
static SpecGaugeRegistry g_gauge_registry;

// ============================================================================
// E8 CENTRAL CHARGE CONSTRAINT (Flavor Gluing Rule)
// ============================================================================
// Reference: "Base Classification of 6d N=(1,0) Supergravity" (draft)
//            Section 3.2, eq (3.4)
//
// E-string theory on -1 curve has total left-moving central charge c_L = 8
// Sum of central charges of gauge algebras intersecting -1 must be ≤ 8
//
// For T=1 (Hirzebruch surface): bound relaxes to 20 (heterotic string constraint)
// This reflects eq (2.8) in the draft.

using GaugeCombination = std::multiset<GaugeAlgebra>;

class E8CentralChargeChecker {
public:
    // E-string budget (T > 1)
    static constexpr double E_STRING_BUDGET = 8.0;
    // Heterotic string budget (T = 1)
    static constexpr double H_STRING_BUDGET = 20.0;
    static constexpr double TOLERANCE = 1e-9;
    
    // Check if combination satisfies central charge constraint
    // T = number of tensor multiplets (default T > 1)
    bool is_allowed(const GaugeCombination& combo, int T = 2) const {
        double budget = (T == 1) ? H_STRING_BUDGET : E_STRING_BUDGET;
        double total_c = get_total_central_charge(combo);
        return total_c <= budget + TOLERANCE;
    }
    
    // Get total central charge
    double get_total_central_charge(const GaugeCombination& combo) const {
        double total_c = 0.0;
        for (auto g : combo) {
            if (g != GaugeAlgebra::none) {
                total_c += gauge_central_charge(g);
            }
        }
        return total_c;
    }
    
    // Get the budget based on T
    static double get_budget(int T) {
        return (T == 1) ? H_STRING_BUDGET : E_STRING_BUDGET;
    }
    
    // Also keep dimension check as secondary validation
    bool passes_dim_check(const GaugeCombination& combo) const {
        int total = 0;
        for (auto g : combo) {
            total += gauge_dim(g);
        }
        return total <= 248;
    }
    
    void dump_combo(const GaugeCombination& combo, int T = 2) const {
        double total_c = 0.0;
        int total_dim = 0;
        double budget = get_budget(T);
        std::cout << "Gauge combination (T=" << T << ", budget=" << budget << "): ";
        for (auto g : combo) {
            if (g != GaugeAlgebra::none) {
                std::cout << gauge_name(g) << "(c=" << gauge_central_charge(g) << ") ";
                total_c += gauge_central_charge(g);
                total_dim += gauge_dim(g);
            }
        }
        std::cout << "| total_c=" << total_c << " dim=" << total_dim;
        std::cout << " | " << (total_c <= budget ? "ALLOWED" : "REJECTED") << "\n";
    }
};

// Global E8 checker
static E8CentralChargeChecker g_e8_checker;

// ============================================================================
// CURVE ORIGIN TRACKING
// ============================================================================

struct CurveOrigin {
    enum class Source { Block, SideLink, Instanton, External };
    
    Source source;
    int source_id;
    int port_idx;
    LKind kind;
    int param;
};

struct ComposedIFResult {
    Eigen::MatrixXi IF;
    std::vector<CurveOrigin> origins;
};

// Helper to build Spec from kind and param
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

ComposedIFResult build_IF_with_tracking(const Topology_enhanced& T) {
    ComposedIFResult result;
    
    TheoryGraph G;
    std::vector<NodeRef> blockNodes;
    std::vector<NodeRef> sideNodes;
    std::vector<NodeRef> instNodes;
    std::vector<NodeRef> extNodes;
    
    int global_idx = 0;
    
    // 1. Blocks
    for (size_t i = 0; i < T.block.size(); ++i) {
        const auto& b = T.block[i];
        Spec sp = make_spec(b.kind, b.param);
        
        Tensor t = build_tensor(sp);
        int num_curves = t.GetIntersectionForm().rows();
        
        for (int p = 0; p < num_curves; ++p) {
            result.origins.push_back({
                CurveOrigin::Source::Block,
                (int)i,
                p,
                b.kind,
                b.param
            });
        }
        global_idx += num_curves;
        
        blockNodes.push_back(G.add(sp));
    }
    
    // 2. SideLinks
    for (size_t i = 0; i < T.side_links.size(); ++i) {
        int param = T.side_links[i].param;
        Spec sp = s(param);
        
        Tensor t = build_tensor(sp);
        int num_curves = t.GetIntersectionForm().rows();
        
        for (int p = 0; p < num_curves; ++p) {
            result.origins.push_back({
                CurveOrigin::Source::SideLink,
                (int)i,
                p,
                LKind::S,
                param
            });
        }
        global_idx += num_curves;
        
        sideNodes.push_back(G.add(sp));
    }
    
    // 3. Instantons
    for (size_t i = 0; i < T.instantons.size(); ++i) {
        int param = T.instantons[i].param;
        Spec sp = s(param);
        
        Tensor t = build_tensor(sp);
        int num_curves = t.GetIntersectionForm().rows();
        
        for (int p = 0; p < num_curves; ++p) {
            result.origins.push_back({
                CurveOrigin::Source::Instanton,
                (int)i,
                p,
                LKind::I,
                param
            });
        }
        global_idx += num_curves;
        
        instNodes.push_back(G.add(sp));
    }
    
    // 4. Externals
    for (size_t i = 0; i < T.externals.size(); ++i) {
        int param = T.externals[i].param;
        Spec sp = e(param);
        
        Tensor t = build_tensor(sp);
        int num_curves = t.GetIntersectionForm().rows();
        
        for (int p = 0; p < num_curves; ++p) {
            result.origins.push_back({
                CurveOrigin::Source::External,
                (int)i,
                p,
                LKind::E,
                param
            });
        }
        global_idx += num_curves;
        
        extNodes.push_back(G.add(sp));
    }
    
    // 5. Connections
    for (const auto& conn : T.l_connection) {
        if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
            conn.v >= 0 && conn.v < (int)blockNodes.size()) {
            G.connect(blockNodes[conn.u], blockNodes[conn.v]);
        }
    }
    
    for (const auto& conn : T.s_connection) {
        if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
            conn.v >= 0 && conn.v < (int)sideNodes.size()) {
            G.connect(sideNodes[conn.v], blockNodes[conn.u]);
        }
    }
    
    for (const auto& conn : T.i_connection) {
        if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
            conn.v >= 0 && conn.v < (int)instNodes.size()) {
            G.connect(instNodes[conn.v], blockNodes[conn.u]);
        }
    }
    
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
    
    // 6. Compose IF
    result.IF = G.ComposeIF_Gluing();
    
    return result;
}

// ============================================================================
// E8 GLUING VALIDATION
// ============================================================================

struct E8ValidationResult {
    bool is_valid;
    std::string reason;
    std::vector<GaugeAlgebra> connected_gauges;
    int total_dim;
    int target_global_idx;
};

E8ValidationResult validate_e8_gluing(
    const Topology_enhanced& T_topo,
    int parent_id,
    int parent_type,
    int port_idx,
    int new_ext_curve,
    int T_tensor = 2  // Number of tensor multiplets (default > 1)
) {
    E8ValidationResult result;
    result.is_valid = true;
    result.target_global_idx = -1;
    
    // 1. Build full IF with tracking
    auto composed = build_IF_with_tracking(T_topo);
    const auto& IF = composed.IF;
    const auto& origins = composed.origins;
    
    // 2. Find target -1 curve's global index
    for (size_t gi = 0; gi < origins.size(); ++gi) {
        const auto& orig = origins[gi];
        
        bool match = false;
        switch (parent_type) {
            case 0:
                match = (orig.source == CurveOrigin::Source::Block &&
                         orig.source_id == parent_id &&
                         orig.port_idx == port_idx);
                break;
            case 1:
                match = (orig.source == CurveOrigin::Source::SideLink &&
                         orig.source_id == parent_id &&
                         orig.port_idx == port_idx);
                break;
            case 2:
                match = (orig.source == CurveOrigin::Source::Instanton &&
                         orig.source_id == parent_id &&
                         orig.port_idx == port_idx);
                break;
        }
        
        if (match) {
            result.target_global_idx = (int)gi;
            break;
        }
    }
    
    if (result.target_global_idx < 0) {
        result.is_valid = false;
        result.reason = "Target port not found in composed IF";
        return result;
    }
    
    // 3. Find neighbors in IF
    std::vector<int> neighbor_indices;
    for (int j = 0; j < IF.cols(); ++j) {
        if (j == result.target_global_idx) continue;
        if (IF(result.target_global_idx, j) != 0) {
            neighbor_indices.push_back(j);
        }
    }
    
    // 4. Get gauge algebra for each neighbor from registry
    for (int ni : neighbor_indices) {
        const auto& orig = origins[ni];
        
        GaugeAlgebra g = g_gauge_registry.get_port_gauge(
            orig.kind, orig.param, orig.port_idx
        );
        
        // Fallback: if not in registry, use curve_to_gauge with self-int
        if (g == GaugeAlgebra::none && !g_gauge_registry.has_entry(orig.kind, orig.param)) {
            int self_int = IF(ni, ni);
            g = curve_to_gauge(self_int);
        }
        
        if (g != GaugeAlgebra::none) {
            result.connected_gauges.push_back(g);
        }
    }
    
    // 5. Add new external's gauge
    GaugeAlgebra new_gauge = curve_to_gauge(new_ext_curve);
    if (new_gauge != GaugeAlgebra::none) {
        result.connected_gauges.push_back(new_gauge);
    }
    
    // 6. Central charge check (E-string budget: ≤8 for T>1, ≤20 for T=1)
    GaugeCombination combo(
        result.connected_gauges.begin(),
        result.connected_gauges.end()
    );
    
    double total_c = g_e8_checker.get_total_central_charge(combo);
    double budget = E8CentralChargeChecker::get_budget(T_tensor);
    result.total_dim = 0;  // Also compute dim for reference
    for (auto g : result.connected_gauges) {
        result.total_dim += gauge_dim(g);
    }
    
    if (!g_e8_checker.is_allowed(combo, T_tensor)) {
        result.is_valid = false;
        std::ostringstream oss;
        oss << "Central charge " << total_c << " > " << budget 
            << " (T=" << T_tensor << " budget). Gauges: ";
        for (auto g : result.connected_gauges) {
            if (g != GaugeAlgebra::none) {
                oss << gauge_name(g) << "(c=" << gauge_central_charge(g) << ") ";
            }
        }
        result.reason = oss.str();
        return result;
    }
    
    return result;
}

// ============================================================================
// ANOMALY CONTRIBUTION TABLES
// Now provided by anomaly_tables.h
// ============================================================================

// ============================================================================
// GLUING RULES
// ============================================================================

const std::vector<int> EXTERNAL_CURVES = {-4, -5, -6, -7, -8, -9, -10, -11, -12};
const int TARGET_CURVE = -1;

std::vector<int> get_allowed_intersections(int ext_curve) {
    return {1};
}

// ============================================================================
// CONFIGURATION
// ============================================================================

enum class AttachMode {
    Classify,
    Single,
    VaryInt,
    Multi,
    Exhaustive
};

struct Config {
    std::string input_path;
    std::string output_dir;
    int num_threads = std::thread::hardware_concurrency();
    bool verbose = false;
    bool check_anomaly = true;
    bool allow_sublattice = true;
    bool check_e8_gluing = true;  // NEW: E8 gluing check toggle
    
    AttachMode mode = AttachMode::Single;
    int max_int = 5;
    int max_attach = 2;
};

// ============================================================================
// SUGRA CLASSIFICATION
// ============================================================================

enum class SugraCategory {
    Valid,
    ValidSublattice,
    NonSquareDet,     // det ≠ n², 나중에 external 더 붙여서 재시도
    InvalidSignature,
    InvalidAnomaly,
    InvalidE8Gluing,
    Neither,
    Error
};

const char* category_name(SugraCategory cat) {
    switch (cat) {
        case SugraCategory::Valid: return "Valid";
        case SugraCategory::ValidSublattice: return "ValidSublattice";
        case SugraCategory::NonSquareDet: return "NonSquareDet";
        case SugraCategory::InvalidSignature: return "InvalidSignature";
        case SugraCategory::InvalidAnomaly: return "InvalidAnomaly";
        case SugraCategory::InvalidE8Gluing: return "InvalidE8Gluing";
        case SugraCategory::Neither: return "Neither";
        case SugraCategory::Error: return "Error";
    }
    return "Unknown";
}

struct AnomalyResult {
    int H;
    int V;
    int T;
    int residual;
    
    bool is_satisfied() const { return residual <= 0; }
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

struct SignatureResult {
    int sig_pos;
    int sig_neg;
    int sig_zero;
    int det;
    bool is_unimodular;
    bool is_perfect_square;
    int sublattice_index;
};

bool is_perfect_square(int n) {
    if (n < 0) n = -n;
    if (n == 0) return false;
    int root = static_cast<int>(std::round(std::sqrt(n)));
    return root * root == n;
}

SignatureResult compute_signature(const Eigen::MatrixXd& IF) {
    SignatureResult result;
    result.sig_pos = 0;
    result.sig_neg = 0;
    result.sig_zero = 0;
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(IF);
    auto eigenvalues = solver.eigenvalues();
    
    const double tol = 1e-8;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        double ev = eigenvalues(i);
        if (std::abs(ev) < tol) {
            result.sig_zero++;
        } else if (ev > tol) {
            result.sig_pos++;
        } else {
            result.sig_neg++;
        }
    }
    
    result.det = static_cast<int>(std::round(IF.determinant()));
    result.is_unimodular = (std::abs(result.det) == 1);
    result.is_perfect_square = is_perfect_square(result.det);
    result.sublattice_index = result.is_perfect_square ? 
        static_cast<int>(std::round(std::sqrt(std::abs(result.det)))) : 0;
    
    return result;
}

struct SugraResult {
    SugraCategory category;
    SignatureResult signature;
    AnomalyResult anomaly;
    std::string e8_rejection_reason;
    
    // External gauge information
    std::vector<GaugeAlgebra> external_gauges;
    int external_V;  // Vector multiplets from externals
    int external_H;  // Hypermultiplets from externals
};

SugraResult classify_sugra(const Topology_enhanced& T, const Config& config) {
    SugraResult result;
    result.category = SugraCategory::Error;
    result.external_V = 0;
    result.external_H = 0;
    
    try {
        // Collect external gauge information
        for (const auto& ext : T.externals) {
            GaugeAlgebra g = curve_to_gauge(-ext.param);  // param is positive, curve is negative
            result.external_gauges.push_back(g);
            
            // Get anomaly contribution
            auto contrib = get_external_contrib(ext.param);
            result.external_V += contrib.V;
            result.external_H += contrib.H;
        }
        
        // Build IF with tracking
        auto composed = build_IF_with_tracking(T);
        auto IF = composed.IF;
        
        if (IF.rows() == 0) {
            result.category = SugraCategory::Error;
            return result;
        }
        
        // 1. Compute Signature
        result.signature = compute_signature(IF.cast<double>());
        
        // 2. Check Signature (must have exactly 1 positive eigenvalue)
        if (result.signature.sig_pos != 1) {
            result.category = SugraCategory::InvalidSignature;
            return result;
        }
        
        // 3. Compute Anomaly (always compute for info)
        result.anomaly = compute_anomaly(T, result.signature.sig_neg);
        
        // 4. Check Determinant - categorize based on det type
        if (result.signature.is_unimodular) {
            // |det| = 1: Check anomaly and return Valid if passed
            if (config.check_anomaly && !result.anomaly.is_satisfied()) {
                result.category = SugraCategory::InvalidAnomaly;
                return result;
            }
            result.category = SugraCategory::Valid;
            return result;
        }
        
        if (result.signature.is_perfect_square) {
            // |det| = n²: Sublattice, needs dummy LST attachment later
            if (config.check_anomaly && !result.anomaly.is_satisfied()) {
                result.category = SugraCategory::InvalidAnomaly;
                return result;
            }
            result.category = SugraCategory::ValidSublattice;
            return result;
        }
        
        // |det| ≠ n²: May need more externals
        // Check anomaly - if already violated, more externals won't help
        if (config.check_anomaly && !result.anomaly.is_satisfied()) {
            result.category = SugraCategory::InvalidAnomaly;
            return result;
        }
        result.category = SugraCategory::NonSquareDet;
        return result;
        
    } catch (const std::exception& e) {
        result.category = SugraCategory::Error;
        return result;
    }
}

// ============================================================================
// PORT FINDING
// ============================================================================

struct AttachablePort {
    int parent_id;
    int parent_type;
    int port_idx;
    int self_int;
};

std::vector<int> get_spec_diagonal(const Spec& sp) {
    Tensor t = build_tensor(sp);
    auto IF = t.GetIntersectionForm();
    std::vector<int> diag;
    for (int i = 0; i < IF.rows(); ++i) {
        diag.push_back(IF(i, i));
    }
    return diag;
}

std::vector<AttachablePort> find_attachable_ports(const Topology_enhanced& T, int target_curve) {
    std::vector<AttachablePort> ports;
    
    for (size_t idx = 0; idx < T.block.size(); ++idx) {
        const auto& b = T.block[idx];
        Spec sp = make_spec(b.kind, b.param);
        auto diag = get_spec_diagonal(sp);
        
        for (size_t port_idx = 0; port_idx < diag.size(); ++port_idx) {
            if (diag[port_idx] == target_curve) {
                ports.push_back({(int)idx, 0, (int)port_idx, diag[port_idx]});
            }
        }
    }
    
    for (size_t idx = 0; idx < T.side_links.size(); ++idx) {
        Spec sp = s(T.side_links[idx].param);
        auto diag = get_spec_diagonal(sp);
        
        for (size_t port_idx = 0; port_idx < diag.size(); ++port_idx) {
            if (diag[port_idx] == target_curve) {
                ports.push_back({(int)idx, 1, (int)port_idx, diag[port_idx]});
            }
        }
    }
    
    for (size_t idx = 0; idx < T.instantons.size(); ++idx) {
        Spec sp = s(T.instantons[idx].param);
        auto diag = get_spec_diagonal(sp);
        
        for (size_t port_idx = 0; port_idx < diag.size(); ++port_idx) {
            if (diag[port_idx] == target_curve) {
                ports.push_back({(int)idx, 2, (int)port_idx, diag[port_idx]});
            }
        }
    }
    
    return ports;
}

// ============================================================================
// ATTACHMENT WITH E8 VALIDATION
// ============================================================================

bool attach_external_with_e8_check(
    Topology_enhanced& T,
    const AttachablePort& port,
    int ext_curve,
    const Config& config,
    std::string* rejection_reason = nullptr
) {
    // E8 gluing check for -1 curves
    if (config.check_e8_gluing && port.self_int == -1) {
        auto validation = validate_e8_gluing(
            T, port.parent_id, port.parent_type, port.port_idx, ext_curve
        );
        
        if (!validation.is_valid) {
            if (rejection_reason) {
                *rejection_reason = validation.reason;
            }
            if (config.verbose) {
                std::cerr << "E8 gluing rejected: " << validation.reason << "\n";
                std::cerr << "  Connected gauges: ";
                for (auto g : validation.connected_gauges) {
                    std::cerr << gauge_name(g) << " ";
                }
                std::cerr << "(total dim: " << validation.total_dim << ")\n";
            }
            return false;
        }
    }
    
    // Actual attachment
    int ext_id = T.addExternal(std::abs(ext_curve));
    return T.attachExternal(ext_id, port.parent_id, port.parent_type, port.port_idx);
}

// ============================================================================
// OUTPUT MANAGEMENT
// ============================================================================

struct OutputBuffer {
    std::unordered_map<std::string, std::string> buffers;
    std::mutex mtx;
    
    void append(const std::string& path, const std::string& line) {
        std::lock_guard<std::mutex> lock(mtx);
        buffers[path] += line + '\n';
    }
    
    void flush_to_disk() {
        std::lock_guard<std::mutex> lock(mtx);
        for (const auto& [path, content] : buffers) {
            if (content.empty()) continue;
            
            fs::create_directories(fs::path(path).parent_path());
            std::ofstream out(path, std::ios::app);
            if (out) {
                out << content;
            }
        }
        buffers.clear();
    }
};

std::string get_output_path(const std::string& output_dir, SugraCategory category,
                            const Topology_enhanced& T) {
    std::string cat_str = category_name(category);
    int len = static_cast<int>(T.block.size());
    
    std::string prefix;
    for (int idx = 0; idx < std::min(4, (int)T.block.size()); idx++) {
        switch (T.block[idx].kind) {
            case LKind::g: prefix += 'g'; break;
            case LKind::L: prefix += 'L'; break;
            case LKind::S: prefix += 'S'; break;
            case LKind::I: prefix += 'I'; break;
            case LKind::E: prefix += 'E'; break;
        }
    }
    if (prefix.empty()) prefix = "empty";
    
    std::string dir = output_dir + "/" + cat_str + "/len-" + std::to_string(len);
    fs::create_directories(dir);
    return dir + "/" + prefix + ".txt";
}

// ============================================================================
// STATISTICS
// ============================================================================

struct Stats {
    std::atomic<int> total_input{0};
    std::atomic<int> total_output{0};
    std::atomic<int> valid_count{0};
    std::atomic<int> sublattice_count{0};
    std::atomic<int> nonsquare_det_count{0};  // det ≠ n²
    std::atomic<int> invalid_sig_count{0};
    std::atomic<int> invalid_anomaly_count{0};
    std::atomic<int> invalid_e8_count{0};
    std::atomic<int> neither_count{0};
    std::atomic<int> error_count{0};
    
    void print() const {
        std::cout << "\n=== Statistics ===\n";
        std::cout << "Input topologies:       " << total_input << "\n";
        std::cout << "Output topologies:      " << total_output << "\n";
        std::cout << "  Valid (unimodular):   " << valid_count << "\n";
        std::cout << "  Valid (sublattice):   " << sublattice_count << "\n";
        std::cout << "  NonSquare Det:        " << nonsquare_det_count << "\n";
        std::cout << "  Invalid Signature:    " << invalid_sig_count << "\n";
        std::cout << "  Invalid Anomaly:      " << invalid_anomaly_count << "\n";
        std::cout << "  Invalid E8 Gluing:    " << invalid_e8_count << "\n";
        std::cout << "  Neither:              " << neither_count << "\n";
        std::cout << "  Errors:               " << error_count << "\n";
    }
};

void update_stats(Stats& stats, SugraCategory cat) {
    switch (cat) {
        case SugraCategory::Valid: stats.valid_count++; break;
        case SugraCategory::ValidSublattice: stats.sublattice_count++; break;
        case SugraCategory::NonSquareDet: stats.nonsquare_det_count++; break;
        case SugraCategory::InvalidSignature: stats.invalid_sig_count++; break;
        case SugraCategory::InvalidAnomaly: stats.invalid_anomaly_count++; break;
        case SugraCategory::InvalidE8Gluing: stats.invalid_e8_count++; break;
        case SugraCategory::Neither: stats.neither_count++; break;
        case SugraCategory::Error: stats.error_count++; break;
    }
}

// ============================================================================
// OUTPUT FORMATTING
// ============================================================================

std::string format_output_line(const Topology_enhanced& T, const SugraResult& res, bool verbose) {
    std::string line = TopoLineCompact_enhanced::serialize(T);
    
    // Add external gauge algebras
    if (!res.external_gauges.empty()) {
        line += " | eg=";
        bool first = true;
        for (auto g : res.external_gauges) {
            if (!first) line += ",";
            line += gauge_name(g);
            first = false;
        }
        // Add external anomaly contribution
        line += " eV=" + std::to_string(res.external_V) +
                " eH=" + std::to_string(res.external_H);
    }
    
    line += " | anom=" + std::to_string(res.anomaly.residual + 273) +
            " sig=(" + std::to_string(res.signature.sig_pos) + "," +
            std::to_string(res.signature.sig_neg) + "," +
            std::to_string(res.signature.sig_zero) + ")" +
            " det=" + std::to_string(res.signature.det);
    
    // Add sublattice index for sublattice cases
    if (res.signature.is_perfect_square && !res.signature.is_unimodular) {
        line += " n=" + std::to_string(res.signature.sublattice_index);
    }
    
    if (verbose) {
        line += " # T=" + std::to_string(res.anomaly.T) +
                " H=" + std::to_string(res.anomaly.H) +
                " V=" + std::to_string(res.anomaly.V);
    }
    
    return line;
}

// ============================================================================
// COMBINATION GENERATION
// ============================================================================

void generate_combinations(int n, int k, std::vector<std::vector<int>>& result) {
    if (k > n) return;
    
    std::vector<int> combo(k);
    std::function<void(int, int)> generate = [&](int start, int idx) {
        if (idx == k) {
            result.push_back(combo);
            return;
        }
        for (int i = start; i <= n - k + idx; ++i) {
            combo[idx] = i;
            generate(i + 1, idx + 1);
        }
    };
    generate(0, 0);
}

// ============================================================================
// PROCESSING MODES
// ============================================================================

void process_single_mode(const Topology_enhanced& base, const Config& config,
                         OutputBuffer& output, Stats& stats) {
    auto ports = find_attachable_ports(base, TARGET_CURVE);
    
    if (ports.empty()) return;
    
    for (const auto& port : ports) {
        for (int ext_curve : EXTERNAL_CURVES) {
            Topology_enhanced T = base;
            
            std::string rejection_reason;
            if (!attach_external_with_e8_check(T, port, ext_curve, config, &rejection_reason)) {
                stats.invalid_e8_count++;
                continue;
            }
            
            SugraResult res = classify_sugra(T, config);
            update_stats(stats, res.category);
            
            // Output Valid, ValidSublattice, and NonSquareDet
            if (res.category == SugraCategory::Valid ||
                res.category == SugraCategory::ValidSublattice ||
                res.category == SugraCategory::NonSquareDet) {
                std::string path = get_output_path(config.output_dir, res.category, T);
                std::string line = format_output_line(T, res, config.verbose);
                output.append(path, line);
                stats.total_output++;
            }
        }
    }
}

void process_vary_int_mode(const Topology_enhanced& base, const Config& config,
                           OutputBuffer& output, Stats& stats) {
    auto ports = find_attachable_ports(base, TARGET_CURVE);
    
    if (ports.empty()) return;
    
    for (const auto& port : ports) {
        for (int ext_curve : EXTERNAL_CURVES) {
            auto allowed_ints = get_allowed_intersections(ext_curve);
            
            for (int int_num : allowed_ints) {
                Topology_enhanced T = base;
                
                if (!attach_external_with_e8_check(T, port, ext_curve, config)) {
                    stats.invalid_e8_count++;
                    continue;
                }
                
                SugraResult res = classify_sugra(T, config);
                update_stats(stats, res.category);
                
                if (res.category == SugraCategory::Valid ||
                    res.category == SugraCategory::ValidSublattice ||
                    res.category == SugraCategory::NonSquareDet) {
                    std::string path = get_output_path(config.output_dir, res.category, T);
                    std::string line = format_output_line(T, res, config.verbose);
                    output.append(path, line);
                    stats.total_output++;
                }
            }
        }
    }
}

void process_multi_mode(const Topology_enhanced& base, const Config& config,
                        OutputBuffer& output, Stats& stats) {
    auto ports = find_attachable_ports(base, TARGET_CURVE);
    
    if (ports.empty()) return;
    
    for (int ext_curve : EXTERNAL_CURVES) {
        for (int num_attach = 1; num_attach <= std::min(config.max_attach, (int)ports.size()); ++num_attach) {
            std::vector<std::vector<int>> combos;
            generate_combinations((int)ports.size(), num_attach, combos);
            
            for (const auto& combo : combos) {
                Topology_enhanced T = base;
                bool ok = true;
                
                for (int idx : combo) {
                    if (!attach_external_with_e8_check(T, ports[idx], ext_curve, config)) {
                        ok = false;
                        stats.invalid_e8_count++;
                        break;
                    }
                }
                
                if (!ok) continue;
                
                SugraResult res = classify_sugra(T, config);
                update_stats(stats, res.category);
                
                if (res.category == SugraCategory::Valid ||
                    res.category == SugraCategory::ValidSublattice ||
                    res.category == SugraCategory::NonSquareDet) {
                    std::string path = get_output_path(config.output_dir, res.category, T);
                    std::string line = format_output_line(T, res, config.verbose);
                    output.append(path, line);
                    stats.total_output++;
                }
            }
        }
    }
}

void process_exhaustive_mode(const Topology_enhanced& base, const Config& config,
                             OutputBuffer& output, Stats& stats) {
    auto ports = find_attachable_ports(base, TARGET_CURVE);
    
    if (ports.empty()) return;
    
    for (int num_attach = 1; num_attach <= std::min(config.max_attach, (int)ports.size()); ++num_attach) {
        std::vector<std::vector<int>> port_combos;
        generate_combinations((int)ports.size(), num_attach, port_combos);
        
        for (const auto& port_combo : port_combos) {
            std::vector<std::vector<int>> ext_combos;
            std::function<void(int, std::vector<int>&)> cartesian = [&](int depth, std::vector<int>& current) {
                if (depth == (int)port_combo.size()) {
                    ext_combos.push_back(current);
                    return;
                }
                for (int ext_curve : EXTERNAL_CURVES) {
                    current.push_back(ext_curve);
                    cartesian(depth + 1, current);
                    current.pop_back();
                }
            };
            std::vector<int> temp;
            cartesian(0, temp);
            
            for (const auto& ext_combo : ext_combos) {
                Topology_enhanced T = base;
                bool ok = true;
                
                for (size_t i = 0; i < port_combo.size(); ++i) {
                    if (!attach_external_with_e8_check(T, ports[port_combo[i]], ext_combo[i], config)) {
                        ok = false;
                        stats.invalid_e8_count++;
                        break;
                    }
                }
                
                if (!ok) continue;
                
                SugraResult res = classify_sugra(T, config);
                update_stats(stats, res.category);
                
                if (res.category == SugraCategory::Valid ||
                    res.category == SugraCategory::ValidSublattice ||
                    res.category == SugraCategory::NonSquareDet) {
                    std::string path = get_output_path(config.output_dir, res.category, T);
                    std::string line = format_output_line(T, res, config.verbose);
                    output.append(path, line);
                    stats.total_output++;
                }
            }
        }
    }
}

// ============================================================================
// MAIN PROCESSING
// ============================================================================

void process_topology(const Topology_enhanced& base, const Config& config,
                      OutputBuffer& output, Stats& stats) {
    stats.total_input++;
    
    if (config.mode == AttachMode::Classify) {
        SugraResult res = classify_sugra(base, config);
        update_stats(stats, res.category);
        
        if (res.category == SugraCategory::Error) return;
        
        if (res.category == SugraCategory::Valid ||
            res.category == SugraCategory::ValidSublattice ||
            res.category == SugraCategory::NonSquareDet) {
            std::string path = get_output_path(config.output_dir, res.category, base);
            std::string line = format_output_line(base, res, config.verbose);
            output.append(path, line);
            stats.total_output++;
        }
        return;
    }
    
    switch (config.mode) {
        case AttachMode::Single:
            process_single_mode(base, config, output, stats);
            break;
        case AttachMode::VaryInt:
            process_vary_int_mode(base, config, output, stats);
            break;
        case AttachMode::Multi:
            process_multi_mode(base, config, output, stats);
            break;
        case AttachMode::Exhaustive:
            process_exhaustive_mode(base, config, output, stats);
            break;
        default:
            break;
    }
}

void process_file(const std::string& filepath, const Config& config,
                  OutputBuffer& output, Stats& stats) {
    std::ifstream infile(filepath);
    if (!infile) {
        std::cerr << "Cannot open: " << filepath << "\n";
        return;
    }
    
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        
        Topology_enhanced T;
        if (!TopoLineCompact_enhanced::deserialize(line, T)) {
            continue;
        }
        
        process_topology(T, config, output, stats);
        
        if (config.verbose && stats.total_input % 100 == 0) {
            std::cout << "Processed " << stats.total_input << " topologies...\r" << std::flush;
        }
    }
}

// ============================================================================
// DEBUG UTILITIES
// ============================================================================

void print_e8_validation_detail(
    const Topology_enhanced& T,
    int parent_id,
    int parent_type,
    int port_idx,
    int ext_curve
) {
    auto composed = build_IF_with_tracking(T);
    
    std::cout << "=== E8 Gluing Validation Detail ===\n";
    std::cout << "Target: parent_type=" << parent_type
              << ", parent_id=" << parent_id
              << ", port_idx=" << port_idx << "\n";
    std::cout << "New external: " << ext_curve << "\n\n";
    
    std::cout << "Composed IF (" << composed.IF.rows() << "x" << composed.IF.cols() << "):\n";
    std::cout << composed.IF << "\n\n";
    
    std::cout << "Curve origins:\n";
    for (size_t i = 0; i < composed.origins.size(); ++i) {
        const auto& o = composed.origins[i];
        std::cout << "  [" << i << "] ";
        switch (o.source) {
            case CurveOrigin::Source::Block: std::cout << "Block"; break;
            case CurveOrigin::Source::SideLink: std::cout << "SideLink"; break;
            case CurveOrigin::Source::Instanton: std::cout << "Instanton"; break;
            case CurveOrigin::Source::External: std::cout << "External"; break;
        }
        std::cout << " #" << o.source_id
                  << ", port=" << o.port_idx
                  << ", kind=" << (int)o.kind
                  << ", param=" << o.param;
        std::cout << ", self_int=" << composed.IF(i, i);
        
        GaugeAlgebra g = g_gauge_registry.get_port_gauge(o.kind, o.param, o.port_idx);
        std::cout << ", gauge=" << gauge_name(g);
        std::cout << "\n";
    }
    
    auto result = validate_e8_gluing(T, parent_id, parent_type, port_idx, ext_curve);
    
    std::cout << "\nValidation result: " << (result.is_valid ? "PASS" : "FAIL") << "\n";
    if (!result.is_valid) {
        std::cout << "Reason: " << result.reason << "\n";
    }
    std::cout << "Connected gauges: ";
    for (auto g : result.connected_gauges) {
        std::cout << gauge_name(g) << " ";
    }
    std::cout << "\nTotal dim: " << result.total_dim << "\n";
}

void dump_gauge_registry() {
    std::cout << "=== Spec Gauge Registry ===\n";
    for (const auto& [key, ports] : g_gauge_registry.table) {
        auto [kind, param] = key;
        std::cout << "Kind=" << (int)kind << ", param=" << param << ": ";
        for (size_t i = 0; i < ports.size(); ++i) {
            std::cout << "[" << ports[i].self_int << ":" << gauge_name(ports[i].gauge) << "] ";
        }
        std::cout << "\n";
    }
}

// ============================================================================
// MAIN
// ============================================================================

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n"
              << "\nSUGRA Base Generator v2 (with E8 gluing rules)\n"
              << "\nRequired:\n"
              << "  -i PATH           Input file or directory\n"
              << "  -o DIR            Output directory\n"
              << "\nMode Selection:\n"
              << "  --mode MODE       Attachment mode (default: single)\n"
              << "    classify, single, vary-int, multi, exhaustive\n"
              << "\nAttachment Options:\n"
              << "  --max-attach N    Max simultaneous attachments (default: 2)\n"
              << "\nClassification Options:\n"
              << "  --no-anomaly      Skip anomaly check\n"
              << "  --no-sublattice   Require det = ±1 only\n"
              << "  --no-e8-check     Skip E8 gluing rule check\n"
              << "\nOther:\n"
              << "  -v                Verbose output\n"
              << "  --dump-registry   Print gauge registry and exit\n"
              << "  -h                Show this help\n";
}

AttachMode parse_mode(const std::string& s) {
    if (s == "classify") return AttachMode::Classify;
    if (s == "single") return AttachMode::Single;
    if (s == "vary-int") return AttachMode::VaryInt;
    if (s == "multi") return AttachMode::Multi;
    if (s == "exhaustive") return AttachMode::Exhaustive;
    return AttachMode::Single;
}

const char* mode_name(AttachMode m) {
    switch (m) {
        case AttachMode::Classify: return "classify";
        case AttachMode::Single: return "single";
        case AttachMode::VaryInt: return "vary-int";
        case AttachMode::Multi: return "multi";
        case AttachMode::Exhaustive: return "exhaustive";
    }
    return "unknown";
}

int main(int argc, char* argv[]) {
    Config config;
    bool dump_registry = false;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "-i" && i + 1 < argc) {
            config.input_path = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            config.output_dir = argv[++i];
        } else if (arg == "--mode" && i + 1 < argc) {
            config.mode = parse_mode(argv[++i]);
        } else if (arg == "--max-attach" && i + 1 < argc) {
            config.max_attach = std::stoi(argv[++i]);
        } else if (arg == "--no-anomaly") {
            config.check_anomaly = false;
        } else if (arg == "--no-sublattice") {
            config.allow_sublattice = false;
        } else if (arg == "--no-e8-check") {
            config.check_e8_gluing = false;
        } else if (arg == "-v") {
            config.verbose = true;
        } else if (arg == "--dump-registry") {
            dump_registry = true;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }
    
    if (dump_registry) {
        dump_gauge_registry();
        return 0;
    }
    
    if (config.input_path.empty() || config.output_dir.empty()) {
        std::cerr << "Error: Input and output paths required\n";
        print_usage(argv[0]);
        return 1;
    }
    
    std::cout << "=== SUGRA Base Generator v2 ===\n";
    std::cout << "Input:  " << config.input_path << "\n";
    std::cout << "Output: " << config.output_dir << "\n";
    std::cout << "Mode:   " << mode_name(config.mode) << "\n";
    
    if (config.mode != AttachMode::Classify) {
        std::cout << "Target curve: " << TARGET_CURVE << " (fixed)\n";
        std::cout << "External curves: -4 to -12\n";
        std::cout << "Max attachments: " << config.max_attach << "\n";
    }
    
    std::cout << "Anomaly check: " << (config.check_anomaly ? "ON" : "OFF") << "\n";
    std::cout << "Sublattice allowed: " << (config.allow_sublattice ? "YES" : "NO") << "\n";
    std::cout << "E8 gluing check: " << (config.check_e8_gluing ? "ON" : "OFF") << "\n";
    std::cout << "\n";
    
    OutputBuffer output;
    Stats stats;
    
    if (fs::is_directory(config.input_path)) {
        for (const auto& entry : fs::recursive_directory_iterator(config.input_path)) {
            if (entry.is_regular_file() && entry.path().extension() == ".txt") {
                if (config.verbose) {
                    std::cout << "Processing: " << entry.path().filename() << "\n";
                }
                process_file(entry.path().string(), config, output, stats);
            }
        }
    } else {
        process_file(config.input_path, config, output, stats);
    }
    
    if (config.verbose) {
        std::cout << "\nFlushing output...\n";
    }
    output.flush_to_disk();
    
    stats.print();
    
    std::cout << "\nDone! Output written to: " << config.output_dir << "\n";
    
    return 0;
}
