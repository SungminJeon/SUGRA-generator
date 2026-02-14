// sugra_generator.cpp
// SUGRA base classification with anomaly checking
// Based on external_generator_simple.cpp, adapted for 6D SUGRA

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

#include "Topology_enhanced.h"
#include "TopologyDB_enhanced.hpp"
#include "TopoLineCompact_enhanced.hpp"
#include "Theory_enhanced.h"
#include "Tensor.h"
#include <sstream>
#include <unordered_set>

namespace fs = std::filesystem;

// Import types
using ::NodeRef;
using ::AttachmentPoint;
using ::Spec;
using ::TheoryGraph;

// ============================================================================
// ANOMALY CONTRIBUTION TABLES - EDIT HERE
// ============================================================================

struct AnomalyContrib {
    int H;  // Hypermultiplet contribution
    int V;  // Vector multiplet contribution
};

// Node (LKind::g): param = |self-int| -> (H, V)
// H from nodes is typically 0; V from gauge algebra
AnomalyContrib get_node_contrib(int param) {
    switch (param) {
        case 12: return {0, 248};  // e8
        case 8:  return {0, 133};  // e7
        case 7:  return {0, 133};  // e7
        case 6:  return {0, 78};   // e6
        case 5:  return {0, 52};   // f4
        case 4:  return {0, 28};   // so(8)
        case 3:  return {0, 8};    // su(3)
        case 2:  return {0, 0};    // no gauge
        case 1:  return {0, 0};    // no gauge
        default: return {0, 0};
    }
}

// Interior link (LKind::L): param -> (H, V)
// H from bifundamentals, V from internal gauge curves
AnomalyContrib get_interior_contrib(int param) {
    switch (param) {
        // =====================================================
        // TODO: Fill in from your classification
        // Format: case PARAM: return {H_contribution, V_contribution};
        // =====================================================
        
	// (1,1)
	case 11: return {0, 0};

	case 22: return {0,8};
	case 33: return {16,27};
	case 331: return {0,68};
	case 44: return {16,86};
	case 55: return {16,86};


	// (2,3) type links
	case 23: return {8,17};
	case 32: return {8,17};

        // (2,4) type links
        case 24:  return {8, 17};
	case 42:  return {8, 17};
      
        
        // (3,4) type links
        case 34:  return {8, 77};
	case 43:  return {8, 77};
        
        // (3,5) type links
        case 35:  return {8, 77};
	case 53: return {8,77};

	// (4,5) type links
	case 45: return {16,86};
	case 54: return {16,86};
        
        // Add more as needed...
        
        default:  return {0, 0};
    }
}

// SideLink: param -> (H, V)
// Based on atomic classification appendix
AnomalyContrib get_sidelink_contrib(int param) {
    switch (param) {
        // =====================================================
        // TODO: Fill in from atomic classification
        // Format: case PARAM: return {H_contribution, V_contribution};
        // =====================================================
        
        // Example placeholders (EDIT THESE):
        // case 944: return {4, 0};
        // case 945: return {8, 0};
	
	// SideLink parameters - switch case format
	// case 파라미터: return {곡선들};
	// === instantons (88 notation - blowdown induced) ===
	    case 1: return {0,0};
	    case 882: return {0,0};
	    case 883: return {0,0};
	    case 884: return {0,0};
	    case 885: return {0,0};
	    case 886: return {0,0};
	    case 887: return {0,0};
	    case 8881: return {0,0};
	    case 889: return {0,0};
	    case 8810: return {0,0};
	    case 8811: return {0,0};
	   
      // === interiors ===
      // (1,1)
	    case 11: return {0, 0};

	    case 22: return {0,8};
	    case 33: return {16,27};
	    case 331: return {0,68};
	    case 44: return {16,86};
	    case 55: return {16,86};


		     // (2,3) type links
	    case 23: return {8,17};
	    case 32: return {8,17};

		     // (2,4) type links
	    case 24:  return {8, 17};
	    case 42:  return {8, 17};


		      // (3,4) type links
	    case 34:  return {8, 77};
	    case 43:  return {8, 77};

		      // (3,5) type links
	    case 35:  return {8, 77};
	    case 53: return {8,77};

		     // (4,5) type links
	    case 45: return {16,86};
	    case 54: return {16,86};

      // === alkali 2 links with no -5 ===
	    case 991: return {8,17};
	    case 9920: return {16,27};
	    case 9902: return {16,27};
	    case 993: return {16,27};

      // === alkali 1 links with no -5 ===
	    case 91: return {8,17};
	    case 92: return {16,27};
	    case 93: return {8,17};
	    case 94: return {16,34};
	    case 95: return {16,24};
	    case 96: return {8,25};
	    case 97: return {8,17};
	    case 98: return {16,27};
	    case 99: return {16,34};
	    case 910: return {16,34};
	    case 911: return {8,25};
	    case 912: return {0,8};
	    case 913: return {8,25};
	    case 914: return {8,25};
	    case 915: return {0,16};
	    case 916: return {8,17};
	    case 917: return {8,17};

      // === alkali 3 links with one -5 curve ===
	    case 99910: return {0,60};
	    case 99920: return {8,69};
	    case 99930: return {8,69};

	// === alkali 2 links with one -5 curve ===
	    case 994: return {0,68};
	    case 995: return {8,77};
	    case 996: return {8,77};
	    case 997: return {16,86};
	    case 998: return {16,86};
	    case 999: return {16,86};
	    case 9910: return {8,77};
	    case 9911: return {16,86};
	    case 9912: return {8,69};
	    case 9913: return {8,69};
	    case 9914: return {16,79};

       // === alkali 1 links with one -5 curve ===
	    case 918: return {8,69};
	    case 919: return {16,86};
	    case 920: return {16,86};
	    case 921: return {16,86};
	    case 922: return {8,77};
	    case 923: return {24,96};
	    case 924: return {8,69};
	    case 925: return {16,79};
	    case 926: return {16,86};
	    case 927: return {16,86};
	    case 928: return {16,86};
	    case 929: return {8,77};
	    case 930: return {24,96};
	    case 931: return {24,96};
	    case 932: return {24,96};
	    case 933: return {16,87};
	    case 934: return {0,60};
	    case 935: return {8,77};
	    case 936: return {8,77};
	    case 937: return {8,77};
	    case 938: return {0,68};
	    case 939: return {16,87};
	    case 940: return {8,69};
	    case 941: return {8,69};
	    case 942: return {8,69};
	    case 943: return {0,60};
	    case 944: return {8,69};
	    case 945: return {8,69};

      // === alkali 2 links with two -5 curves ===
	    case 9915: return {8,129};
	    case 9916: return {8,129};
	    case 9917: return {0,120};

       // === alkali 1 links with two -5 curves ===
	    case 946: return {16,138};
	    case 947: return {8,129};
	    case 948: return {16,146};
	    case 949: return {16,146};
	    case 950: return {8,137};
	    case 951: return {16,138};
	    case 952: return {8,129};
	    case 953: return {16,146};
	    case 954: return {8,137};
	    case 955: return {8,129};
	    case 956: return {0,120};
	    case 957: return {0,128};

      // === alkali 1 link with one -5 curve (omitted in appendix) ===
	    case 958: return {0,60};	


	    default:  return {0, 0};
    }
}

// Instanton: param -> (H, V)
// For minimal instantons, typically no contribution
AnomalyContrib get_instanton_contrib(int param) {
    // Minimal instantons have no gauge algebra
    return {0, 0};
}

// External: param -> (H, V)
AnomalyContrib get_external_contrib(int param) {
    switch (param) {
        // =====================================================
        // TODO: Fill in if externals contribute
        // =====================================================
	//
	case 1: return {0,0};
    	case 2: return {0,0};
	case 3: return {0,8};
	case 4: return {0,28};
	case 5: return {0,52};
	case 6: return {0,78};
	case 7: return {0,133};
	case 8: return {0,133};
	case 12: return {0,248};



        default:  return {0, 0};
    }
}

// ============================================================================
// GLUING RULES - EDIT HERE
// ============================================================================

// External curve types to try (gauge curves: -4 to -12)
const std::vector<int> EXTERNAL_CURVES = {-4, -5, -6, -7, -8, -9, -10, -11, -12};

// Target curve is fixed to -1
const int TARGET_CURVE = -1;

// Get allowed intersection numbers for an external curve type
std::vector<int> get_allowed_intersections(int ext_curve) {
    // All gauge curves (-4 to -12) only allow int=1
    return {1};
}

// Check if a curve can have externals attached
bool is_attachable_curve(int self_int, int target_curve) {
    // Only attach to target curve type
    return self_int == target_curve;
}

// Get allowed intersection numbers for a curve
// For -4 to -12: only int=1 allowed
// For -1, -2, -3: can vary
std::vector<int> get_allowed_intersection_numbers(int self_int, int max_int) {
    switch (self_int) {
        // Gauge curves: only int=1
        case -12: case -11: case -10: case -9:
        case -8: case -7: case -6: case -5: case -4:
            return {1};
        
        // -3 curve
        case -3:
            return {1, 2};
        
        // -2 curve
        case -2: {
            std::vector<int> result;
            for (int i = 1; i <= std::min(max_int, 3); ++i) result.push_back(i);
            return result;
        }
        
        // -1 curve (tensor): most flexible
        case -1: {
            std::vector<int> result;
            for (int i = 1; i <= max_int; ++i) result.push_back(i);
            return result;
        }
        
        default:
            return {1};
    }
}

// Legacy function for compatibility
std::vector<int> get_allowed_external_params(int self_int) {
    return get_allowed_intersection_numbers(self_int, 5);
}

// ============================================================================
// Configuration
// ============================================================================

enum class AttachMode {
    Classify,       // Just classify, no attachment
    Single,         // One -1 curve, int=1 only
    VaryInt,        // One -1 curve, vary intersection number
    Multi,          // Multiple -1 curves simultaneously
    Exhaustive      // All combinations
};

struct Config {
    std::string input_path;
    std::string output_dir;
    std::vector<std::string> attachment_specs;
    int num_threads = std::thread::hardware_concurrency();
    bool verbose = false;
    bool classify_only = false;
    bool check_anomaly = true;      // Check gravitational anomaly
    bool allow_sublattice = true;   // Allow det = ±n² (not just ±1)
    
    // Attachment mode settings
    AttachMode mode = AttachMode::Single;
    int max_int = 5;                // Max intersection number
    int max_attach = 2;             // Max simultaneous attachments
    int target_curve = -1;          // Target curve self-int (default: -1)
};

// ============================================================================
// SUGRA Classification
// ============================================================================

enum class SugraCategory {
    Valid,              // All conditions satisfied, det = ±1
    ValidSublattice,    // All conditions satisfied, det = ±n²
    InvalidSignature,   // Signature not (1, -, -, ...)
    InvalidDet,         // det not ±1 or ±n²
    InvalidAnomaly,     // Anomaly cancellation failed
    Neither,            // Other failure
    Error
};

const char* category_name(SugraCategory cat) {
    switch (cat) {
        case SugraCategory::Valid: return "Valid";
        case SugraCategory::ValidSublattice: return "ValidSublattice";
        case SugraCategory::InvalidSignature: return "InvalidSignature";
        case SugraCategory::InvalidDet: return "InvalidDet";
        case SugraCategory::InvalidAnomaly: return "InvalidAnomaly";
        case SugraCategory::Neither: return "Neither";
        case SugraCategory::Error: return "Error";
    }
    return "Unknown";
}

// ============================================================================
// Anomaly Calculation
// ============================================================================

struct AnomalyResult {
    int H;          // Total hypermultiplets
    int V;          // Total vector multiplets
    int T;          // Tensor multiplets (= sig_neg)
    int residual;   // H - V + 29T - 273 (should be 0)
    
    bool is_satisfied() const { return residual <= 0; }
};

AnomalyResult compute_anomaly(const Topology_enhanced& T, int sig_neg) {
    AnomalyResult result;
    result.H = 0;
    result.V = 0;
    result.T = sig_neg;
    
    // 1. Blocks: nodes (g) and interior links (L)
    for (const auto& b : T.block) {
        AnomalyContrib c;
        switch (b.kind) {
            case LKind::g:  // Node
                c = get_node_contrib(b.param);
                break;
            case LKind::L:  // Interior link
                c = get_interior_contrib(b.param);
                break;
            case LKind::S:  // SideLink (if stored as block)
                c = get_sidelink_contrib(b.param);
                break;
            case LKind::I:  // Instanton (if stored as block)
                c = get_instanton_contrib(b.param);
                break;
            case LKind::E:  // External
                c = get_external_contrib(b.param);
                break;
            default:
                c = {0, 0};
                break;
        }
        result.H += c.H;
        result.V += c.V;
    }
    
    // 2. SideLinks (from side_links vector)
    for (const auto& sl : T.side_links) {
        auto c = get_sidelink_contrib(sl.param);
        result.H += c.H;
        result.V += c.V;
    }
    
    // 3. Instantons (minimal → typically no contribution)
    for (const auto& inst : T.instantons) {
        auto c = get_instanton_contrib(inst.param);
        result.H += c.H;
        result.V += c.V;
    }
    
    // 4. Externals
    for (const auto& ext : T.externals) {
        auto c = get_external_contrib(ext.param);
        result.H += c.H;
        result.V += c.V;
    }
    
    // Gravitational anomaly: H - V + 29T - 273 = 0
    result.residual = result.H - result.V + 29 * result.T - 273;
    
    return result;
}

// ============================================================================
// Signature and Determinant Analysis
// ============================================================================

struct SignatureResult {
    int sig_pos;    // Number of positive eigenvalues
    int sig_neg;    // Number of negative eigenvalues
    int sig_zero;   // Number of zero eigenvalues
    int det;        // Determinant of intersection form
    bool is_unimodular;      // |det| = 1
    bool is_perfect_square;  // |det| = n²
    int sublattice_index;    // n if perfect square, else 0
};

bool is_perfect_square(int n) {
    if (n < 0) n = -n;
    if (n == 0) return false;  // 0 is not a valid sublattice index
    int root = static_cast<int>(std::round(std::sqrt(n)));
    return root * root == n;
}

SignatureResult compute_signature(const Eigen::MatrixXd& IF) {
    SignatureResult result;
    result.sig_pos = 0;
    result.sig_neg = 0;
    result.sig_zero = 0;
    
    // Eigenvalues
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
    
    // Determinant
    result.det = static_cast<int>(std::round(IF.determinant()));
    result.is_unimodular = (std::abs(result.det) == 1);
    result.is_perfect_square = is_perfect_square(result.det);
    result.sublattice_index = result.is_perfect_square ? 
        static_cast<int>(std::round(std::sqrt(std::abs(result.det)))) : 0;
    
    return result;
}

// ============================================================================
// Main SUGRA Classification
// ============================================================================

struct SugraResult {
    SugraCategory category;
    SignatureResult signature;
    AnomalyResult anomaly;
};

SugraResult classify_sugra(const Topology_enhanced& T, const Config& config) {
    SugraResult result;
    result.category = SugraCategory::Error;
    
    try {
        // Build TheoryGraph (same as before)
        TheoryGraph G;
        
        std::vector<NodeRef> blockNodes;
        for (const auto& b : T.block) {
            Spec sp;
            switch(b.kind) {
                case LKind::g: sp = n(b.param); break;
                case LKind::L: sp = i(b.param); break;
                case LKind::S: sp = s(b.param); break;
                case LKind::I: sp = s(b.param); break;
                case LKind::E: sp = e(b.param); break;
                default: sp = n(b.param); break;
            }
            blockNodes.push_back(G.add(sp));
        }
        
        std::vector<NodeRef> sideNodes;
        for (const auto& sl : T.side_links) {
            sideNodes.push_back(G.add(s(sl.param)));
        }
        
        std::vector<NodeRef> instNodes;
        for (const auto& inst : T.instantons) {
            instNodes.push_back(G.add(s(inst.param)));
        }
        
        std::vector<NodeRef> extNodes;
        for (const auto& ext : T.externals) {
            extNodes.push_back(G.add(e(ext.param)));
        }
        
        // Connect interior links
        for (const auto& conn : T.l_connection) {
            if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
                conn.v >= 0 && conn.v < (int)blockNodes.size()) {
                G.connect(blockNodes[conn.u], blockNodes[conn.v]);
            }
        }
        
        // Connect sidelinks
        for (const auto& conn : T.s_connection) {
            if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
                conn.v >= 0 && conn.v < (int)sideNodes.size()) {
                G.connect(sideNodes[conn.v], blockNodes[conn.u]);
            }
        }
        
        // Connect instantons
        for (const auto& conn : T.i_connection) {
            if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
                conn.v >= 0 && conn.v < (int)instNodes.size()) {
                G.connect(instNodes[conn.v], blockNodes[conn.u]);
            }
        }
        
        // Connect externals
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
        
        // Get intersection form
        auto IF = G.ComposeIF_Gluing();
        if (IF.rows() == 0) {
            result.category = SugraCategory::Error;
            return result;
        }
        
        // === 1. Compute Signature ===
        result.signature = compute_signature(IF.cast<double>());
        
        // === 2. Check Signature: exactly 1 positive eigenvalue ===
        if (result.signature.sig_pos != 1) {
            result.category = SugraCategory::InvalidSignature;
            return result;
        }
        
        // === 3. Check Determinant ===
        if (!result.signature.is_unimodular) {
            if (!config.allow_sublattice || !result.signature.is_perfect_square) {
                result.category = SugraCategory::InvalidDet;
                return result;
            }
        }
        
        // === 4. Compute and Check Anomaly ===
        result.anomaly = compute_anomaly(T, result.signature.sig_neg);
        
        if (config.check_anomaly && !result.anomaly.is_satisfied()) {
            result.category = SugraCategory::InvalidAnomaly;
            return result;
        }
        
        // === 5. All checks passed ===
        if (result.signature.is_unimodular) {
            result.category = SugraCategory::Valid;
        } else {
            result.category = SugraCategory::ValidSublattice;
        }
        
        return result;
        
    } catch (const std::exception& e) {
        result.category = SugraCategory::Error;
        return result;
    }
}

// ============================================================================
// Attachment Specification Parsing (same as before)
// ============================================================================

struct AttachmentSpec {
    enum Type { Block, SideLink, Instanton } type;
    int index;
    
    std::string describe() const {
        std::string type_str;
        switch (type) {
            case Block: type_str = "b"; break;
            case SideLink: type_str = "s"; break;
            case Instanton: type_str = "i"; break;
        }
        return type_str + "(" + std::to_string(index) + ")";
    }
};

bool parse_attachment_spec(const std::string& spec, AttachmentSpec& out) {
    if (spec.size() < 4) return false;
    
    char type_char = spec[0];
    if (spec[1] != '(' || spec.back() != ')') return false;
    
    std::string num_str = spec.substr(2, spec.size() - 3);
    try {
        int idx = std::stoi(num_str);
        if (idx < 0) return false;
        
        switch (type_char) {
            case 's': case 'S':
                out.type = AttachmentSpec::SideLink;
                out.index = idx;
                return true;
            case 'b': case 'B':
                out.type = AttachmentSpec::Block;
                out.index = idx;
                return true;
            case 'i': case 'I':
                out.type = AttachmentSpec::Instanton;
                out.index = idx;
                return true;
            default:
                return false;
        }
    } catch (...) {
        return false;
    }
}

// ============================================================================
// Port Info and External Attachment (same as before)
// ============================================================================

struct PortInfo {
    int parent_id;
    int parent_type;
    int port_idx;
    int self_int;
    int parent_param;
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

std::vector<PortInfo> get_possible_ports(const Topology_enhanced& T, const AttachmentSpec& spec) {
    std::vector<PortInfo> ports;
    
    switch (spec.type) {
        case AttachmentSpec::Block:
            if (spec.index >= 0 && spec.index < (int)T.block.size()) {
                const auto& block = T.block[spec.index];
                
                Spec sp;
                switch(block.kind) {
                    case LKind::g: sp = n(block.param); break;
                    case LKind::L: sp = i(block.param); break;
                    case LKind::S: sp = s(block.param); break;
                    case LKind::I: sp = s(block.param); break;
                    case LKind::E: sp = e(block.param); break;
                    default: sp = n(block.param); break;
                }
                
                auto diag = get_spec_diagonal(sp);
                
                for (int port_idx = 0; port_idx < (int)diag.size(); ++port_idx) {
                    ports.push_back({
                        spec.index, 0, port_idx, diag[port_idx], block.param
                    });
                }
            }
            break;
            
        case AttachmentSpec::SideLink:
            if (spec.index >= 0 && spec.index < (int)T.side_links.size()) {
                int param = T.side_links[spec.index].param;
                Spec sp = s(param);
                auto diag = get_spec_diagonal(sp);
                
                for (int port_idx = 0; port_idx < (int)diag.size(); ++port_idx) {
                    ports.push_back({
                        spec.index, 1, port_idx, diag[port_idx], param
                    });
                }
            }
            break;
            
        case AttachmentSpec::Instanton:
            if (spec.index >= 0 && spec.index < (int)T.instantons.size()) {
                int param = T.instantons[spec.index].param;
                Spec sp = s(param);
                auto diag = get_spec_diagonal(sp);
                
                for (int port_idx = 0; port_idx < (int)diag.size(); ++port_idx) {
                    ports.push_back({
                        spec.index, 2, port_idx, diag[port_idx], param
                    });
                }
            }
            break;
    }
    
    return ports;
}

bool add_external_at_port(Topology_enhanced& T, const PortInfo& port, int ext_param) {
    int ext_id = T.addExternal(ext_param);
    return T.attachExternal(ext_id, port.parent_id, port.parent_type, port.port_idx);
}

// ============================================================================
// Output Management
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
// Statistics
// ============================================================================

struct Stats {
    std::atomic<int> total_input{0};
    std::atomic<int> total_output{0};
    std::atomic<int> valid_count{0};
    std::atomic<int> sublattice_count{0};
    std::atomic<int> invalid_sig_count{0};
    std::atomic<int> invalid_det_count{0};
    std::atomic<int> invalid_anomaly_count{0};
    std::atomic<int> neither_count{0};
    std::atomic<int> error_count{0};
    
    void print() const {
        std::cout << "\n=== Statistics ===\n";
        std::cout << "Input topologies:       " << total_input << "\n";
        std::cout << "Output topologies:      " << total_output << "\n";
        std::cout << "  Valid (unimodular):   " << valid_count << "\n";
        std::cout << "  Valid (sublattice):   " << sublattice_count << "\n";
        std::cout << "  Invalid Signature:    " << invalid_sig_count << "\n";
        std::cout << "  Invalid Determinant:  " << invalid_det_count << "\n";
        std::cout << "  Invalid Anomaly:      " << invalid_anomaly_count << "\n";
        std::cout << "  Neither:              " << neither_count << "\n";
        std::cout << "  Errors:               " << error_count << "\n";
    }
};

void update_stats(Stats& stats, SugraCategory cat) {
    switch (cat) {
        case SugraCategory::Valid: stats.valid_count++; break;
        case SugraCategory::ValidSublattice: stats.sublattice_count++; break;
        case SugraCategory::InvalidSignature: stats.invalid_sig_count++; break;
        case SugraCategory::InvalidDet: stats.invalid_det_count++; break;
        case SugraCategory::InvalidAnomaly: stats.invalid_anomaly_count++; break;
        case SugraCategory::Neither: stats.neither_count++; break;
        case SugraCategory::Error: stats.error_count++; break;
    }
}

// ============================================================================
// Find Attachable Ports (target curve type only)
// ============================================================================

struct AttachablePort {
    int parent_id;
    int parent_type;    // 0=block, 1=sidelink, 2=instanton
    int port_idx;
    int self_int;
};

std::vector<AttachablePort> find_attachable_ports(const Topology_enhanced& T, int target_curve) {
    std::vector<AttachablePort> ports;
    
    // Search in blocks
    for (size_t idx = 0; idx < T.block.size(); ++idx) {
        const auto& b = T.block[idx];
        Spec sp;
        switch(b.kind) {
            case LKind::g: sp = n(b.param); break;
            case LKind::L: sp = i(b.param); break;
            case LKind::S: sp = s(b.param); break;
            case LKind::I: sp = s(b.param); break;
            case LKind::E: sp = e(b.param); break;
            default: sp = n(b.param); break;
        }
        
        auto diag = get_spec_diagonal(sp);
        for (size_t port_idx = 0; port_idx < diag.size(); ++port_idx) {
            if (diag[port_idx] == target_curve) {
                ports.push_back({(int)idx, 0, (int)port_idx, diag[port_idx]});
            }
        }
    }
    
    // Search in sidelinks
    for (size_t idx = 0; idx < T.side_links.size(); ++idx) {
        Spec sp = s(T.side_links[idx].param);
        auto diag = get_spec_diagonal(sp);
        for (size_t port_idx = 0; port_idx < diag.size(); ++port_idx) {
            if (diag[port_idx] == target_curve) {
                ports.push_back({(int)idx, 1, (int)port_idx, diag[port_idx]});
            }
        }
    }
    
    // Search in instantons
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
// Attachment Helpers
// ============================================================================

bool attach_external(Topology_enhanced& T, const AttachablePort& port, int ext_curve, int intersection_num) {
    // ext_curve is negative (e.g., -1, -2, -3, -5)
    // param stored is |ext_curve| for anomaly calculation
    int ext_id = T.addExternal(std::abs(ext_curve));
    return T.attachExternal(ext_id, port.parent_id, port.parent_type, port.port_idx);
}

// Generate combinations of k elements from n
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
// Mode-based Processing
// ============================================================================

// Helper: format output line with signature and anomaly info
std::string format_output_line(const Topology_enhanced& T, const SugraResult& res, bool verbose) {
    std::string line = TopoLineCompact_enhanced::serialize(T);
    
    // Always add signature and anomaly info
    line += " | anom=" + std::to_string(res.anomaly.residual+273) +
            " sig=(" + std::to_string(res.signature.sig_pos) + "," +
            std::to_string(res.signature.sig_neg) + "," +
            std::to_string(res.signature.sig_zero) + ")" +
            " det=" + std::to_string(res.signature.det);
    
    if (verbose) {
        line += " # T=" + std::to_string(res.anomaly.T) +
                " H=" + std::to_string(res.anomaly.H) +
                " V=" + std::to_string(res.anomaly.V);
    }
    
    return line;
}

void process_single_mode(const Topology_enhanced& base, const Config& config,
                         OutputBuffer& output, Stats& stats) {
    // Find all target curves in base
    auto ports = find_attachable_ports(base, TARGET_CURVE);
    
    if (ports.empty()) {
        if (config.verbose) {
            std::cout << "  No " << TARGET_CURVE << " curves found\n";
        }
        return;
    }
    
    // For each port, try all external curve types with int=1
    for (const auto& port : ports) {
        for (int ext_curve : EXTERNAL_CURVES) {
            Topology_enhanced T = base;
            if (!attach_external(T, port, ext_curve, 1)) continue;
            
            SugraResult res = classify_sugra(T, config);
            update_stats(stats, res.category);
            
            if (res.category == SugraCategory::Valid || 
                res.category == SugraCategory::ValidSublattice) {
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
    
    // For each port, try all external curves with varying intersection numbers
    for (const auto& port : ports) {
        for (int ext_curve : EXTERNAL_CURVES) {
            auto allowed_ints = get_allowed_intersections(ext_curve);
            
            for (int int_num : allowed_ints) {
                Topology_enhanced T = base;
                if (!attach_external(T, port, ext_curve, int_num)) continue;
                
                SugraResult res = classify_sugra(T, config);
                update_stats(stats, res.category);
                
                if (res.category == SugraCategory::Valid || 
                    res.category == SugraCategory::ValidSublattice) {
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
    
    // Try attaching to multiple ports simultaneously with same external type
    for (int ext_curve : EXTERNAL_CURVES) {
        for (int num_attach = 1; num_attach <= std::min(config.max_attach, (int)ports.size()); ++num_attach) {
            std::vector<std::vector<int>> combos;
            generate_combinations((int)ports.size(), num_attach, combos);
            
            for (const auto& combo : combos) {
                Topology_enhanced T = base;
                bool ok = true;
                
                // Attach to each selected port with int=1
                for (int idx : combo) {
                    if (!attach_external(T, ports[idx], ext_curve, 1)) {
                        ok = false;
                        break;
                    }
                }
                
                if (!ok) continue;
                
                SugraResult res = classify_sugra(T, config);
                update_stats(stats, res.category);
                
                if (res.category == SugraCategory::Valid || 
                    res.category == SugraCategory::ValidSublattice) {
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
    
    // Try all combinations: number of attachments × which ports × external curve types
    for (int num_attach = 1; num_attach <= std::min(config.max_attach, (int)ports.size()); ++num_attach) {
        std::vector<std::vector<int>> port_combos;
        generate_combinations((int)ports.size(), num_attach, port_combos);
        
        for (const auto& port_combo : port_combos) {
            // For each port combination, try all external curve type combinations
            // Generate cartesian product of external curve types
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
            
            // Try each combination
            for (const auto& ext_combo : ext_combos) {
                Topology_enhanced T = base;
                bool ok = true;
                
                for (size_t i = 0; i < port_combo.size(); ++i) {
                    if (!attach_external(T, ports[port_combo[i]], ext_combo[i], 1)) {
                        ok = false;
                        break;
                    }
                }
                
                if (!ok) continue;
                
                SugraResult res = classify_sugra(T, config);
                update_stats(stats, res.category);
                
                if (res.category == SugraCategory::Valid || 
                    res.category == SugraCategory::ValidSublattice) {
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
// Processing
// ============================================================================

void process_topology(const Topology_enhanced& base, const Config& config,
                     OutputBuffer& output, Stats& stats) {
    stats.total_input++;
    
    // Classify-only mode
    if (config.mode == AttachMode::Classify) {
        SugraResult res = classify_sugra(base, config);
        update_stats(stats, res.category);
        
        if (res.category == SugraCategory::Error) return;
        
        if (res.category == SugraCategory::Valid || 
            res.category == SugraCategory::ValidSublattice) {
            std::string path = get_output_path(config.output_dir, res.category, base);
            std::string line = format_output_line(base, res, config.verbose);
            
            output.append(path, line);
            stats.total_output++;
        }
        return;
    }
    
    // Attachment modes
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
// Main
// ============================================================================

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n"
              << "\nSUGRA Base Generator\n"
              << "  - Target curve: -1 (fixed)\n"
              << "  - External curves: -4 to -12 (gauge curves)\n"
              << "\nRequired:\n"
              << "  -i PATH           Input file or directory\n"
              << "  -o DIR            Output directory\n"
              << "\nMode Selection:\n"
              << "  --mode MODE       Attachment mode (default: single)\n"
              << "    classify    : Only classify, no external attachment\n"
              << "    single      : One -1 curve, try all external types (-4 to -12)\n"
              << "    vary-int    : Same as single (int=1 fixed for gauge curves)\n"
              << "    multi       : Multiple -1 curves, same external type\n"
              << "    exhaustive  : All combinations of external types\n"
              << "\nAttachment Options:\n"
              << "  --max-attach N    Max simultaneous attachments (default: 2)\n"
              << "\nClassification Options:\n"
              << "  --no-anomaly      Skip anomaly check\n"
              << "  --no-sublattice   Require det = ±1 only\n"
              << "\nOther:\n"
              << "  -v                Verbose output\n"
              << "  -h                Show this help\n"
              << "\nSUGRA Classification Criteria:\n"
              << "  1. Signature: exactly 1 positive eigenvalue\n"
              << "  2. Determinant: ±1 (unimodular) or ±n² (sublattice)\n"
              << "  3. Anomaly: H - V + 29T <= 273\n"
              << "\nExamples:\n"
              << "  # Classify only\n"
              << "  " << prog << " -i input.txt -o out --mode classify\n"
              << "\n"
              << "  # Single external (-4 to -12) on -1 curves\n"
              << "  " << prog << " -i input.txt -o out --mode single\n"
              << "\n"
              << "  # Multiple externals on -1 curves\n"
              << "  " << prog << " -i input.txt -o out --mode multi --max-attach 3\n"
              << "\n"
              << "  # Exhaustive search (different ext types per port)\n"
              << "  " << prog << " -i input.txt -o out --mode exhaustive --max-attach 2\n";
}

AttachMode parse_mode(const std::string& s) {
    if (s == "classify") return AttachMode::Classify;
    if (s == "single") return AttachMode::Single;
    if (s == "vary-int") return AttachMode::VaryInt;
    if (s == "multi") return AttachMode::Multi;
    if (s == "exhaustive") return AttachMode::Exhaustive;
    return AttachMode::Single;  // default
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
        } else if (arg == "--target" && i + 1 < argc) {
            // Target curve is fixed to -1 now
            std::cerr << "Warning: --target is ignored in (fixed to -1)\n";
            ++i;
        } else if (arg == "--max-int" && i + 1 < argc) {
            config.max_int = std::stoi(argv[++i]);
        } else if (arg == "--max-attach" && i + 1 < argc) {
            config.max_attach = std::stoi(argv[++i]);
        } else if (arg == "--no-anomaly") {
            config.check_anomaly = false;
        } else if (arg == "--no-sublattice") {
            config.allow_sublattice = false;
        } else if (arg == "-v") {
            config.verbose = true;
        } else if (arg == "-a" && i + 1 < argc) {
            // Legacy: ignore old attachment specs
            ++i;
        } else if (arg == "--classify-only") {
            // Legacy: map to classify mode
            config.mode = AttachMode::Classify;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }
    
    if (config.input_path.empty() || config.output_dir.empty()) {
        std::cerr << "Error: Input and output paths required\n";
        print_usage(argv[0]);
        return 1;
    }
    
    std::cout << "=== SUGRA Base Generator ===\n";
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
