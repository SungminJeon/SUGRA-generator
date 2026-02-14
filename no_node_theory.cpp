// no_node_theory.cpp
// Generates all no-node LSTs (side links + noble molecules with external)
// Output: LaTeX file with quiver-style notation
// 
// Algorithm: Intersection matrix → Quiver LaTeX directly
// - Follow a_{i,i+1} for linear chain
// - If a_{i,i+1}=0, find j where a_{i,j}!=0, stack curves (i+1 ~ j-1) on j
// - External (last row) attached with red overset
// 
// Deduplication: Compare eigenvalue spectra of intersection forms

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>
#include <set>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "Tensor.h"
#include "Theory_enhanced.h"

// ============================================================================
// All Side Link Params (from Theory_enhanced.h)
// ============================================================================

const std::vector<int> ALL_SIDELINK_PARAMS = {
    // instantons
    1, 882, 883, 884, 885, 886, 887, 8881, 889, 8810, 8811,
       
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

// ============================================================================
// Gluing Rules
// ============================================================================

// Allowed intersection numbers for an external curve (based on external's self-int)
std::vector<int> get_allowed_intersections(int self_int) {
    switch (self_int) {
        case -1: return {1, 2, 3, 4, 5};
        case -2: return {1, 2, 3};
        case -3: return {1, 2};
        default: return {1};  // -4 ~ -12
    }
}

// Possible external curves
const std::vector<int> EXTERNAL_CURVES = {-1, -2, -3, -5};

// Allowed external curves for a base curve (based on base's self-int)
// Rules:
//   -1 base: can attach -1, -2, -3, -5
//   -2 base: can attach -1, -3
//   -3 base: can attach -1, -2
//   -5 base: can attach -1 only
std::vector<int> get_allowed_external_curves(int base_self_int) {
    switch (base_self_int) {
        case -1: return {-1, -2, -3, -5};
        case -2: return {-1, -3};
        case -3: return {-1, -2};
        case -5: return {-1};
        default: return {-1};  // For other curves, only -1 allowed
    }
}

// ============================================================================
// LST Check
// ============================================================================

bool is_LST(const Tensor& t) {
    return t.NullDirection() == 1;
}

// ============================================================================
// P-type Endpoint Check
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
// ============================================================================

// Compute sorted eigenvalue spectrum of an intersection form
std::vector<double> compute_eigenvalue_spectrum(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    if (n == 0) return {};
    
    // Convert to double matrix for eigenvalue computation
    Eigen::MatrixXd IF_double = IF.cast<double>();
    
    // Compute eigenvalues (intersection form is symmetric)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(IF_double);
    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    
    // Convert to vector and sort
    std::vector<double> spectrum(n);
    for (int i = 0; i < n; ++i) {
        spectrum[i] = eigenvalues(i);
    }
    std::sort(spectrum.begin(), spectrum.end());
    
    return spectrum;
}

// Compare two eigenvalue spectra with tolerance
bool are_spectra_equal(const std::vector<double>& s1, const std::vector<double>& s2, double tol = 1e-9) {
    if (s1.size() != s2.size()) return false;
    
    for (size_t i = 0; i < s1.size(); ++i) {
        if (std::abs(s1[i] - s2[i]) > tol) return false;
    }
    return true;
}

// ============================================================================
// LST Result Structure
// ============================================================================

struct LSTResult {
    Eigen::MatrixXi IF;              // Intersection form
    std::string latex;               // LaTeX representation
    std::string description;         // Description (source info)
    std::vector<double> spectrum;    // Eigenvalue spectrum (computed lazily)
    
    void compute_spectrum() {
        spectrum = compute_eigenvalue_spectrum(IF);
    }
};

// Remove duplicates based on eigenvalue spectrum
std::vector<LSTResult> remove_duplicates_by_spectrum(std::vector<LSTResult>& results) {
    // Compute all spectra
    for (auto& r : results) {
        r.compute_spectrum();
    }
    
    std::vector<LSTResult> unique_results;
    
    for (const auto& r : results) {
        bool is_duplicate = false;
        for (const auto& u : unique_results) {
            if (are_spectra_equal(r.spectrum, u.spectrum)) {
                is_duplicate = true;
                break;
            }
        }
        if (!is_duplicate) {
            unique_results.push_back(r);
        }
    }
    
    return unique_results;
}

// ============================================================================
// Matrix to Quiver LaTeX
// ============================================================================

// Helper: external 표기 생성
std::string make_ext_marker(int ext_self_int, int ext_int_num) {
    std::string ext_str = std::to_string(ext_self_int);
    if (ext_int_num == 1) {
        return "\\textcolor{red}{" + ext_str + "}";
    } else {
        return "\\textcolor{red}{" + ext_str + "^{" + std::to_string(ext_int_num) + "}}";
    }
}

std::string matrix_to_quiver_latex(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    if (n == 0) return "\\emptyset";
    if (n == 1) return "{" + std::to_string(-IF(0,0)) + "}";
    
    // === 1. External 분리 (마지막 row) ===
    int ext_attached_to = -1;
    int ext_int_num = 0;
    int ext_self_int = -IF(n-1, n-1);  // External curve의 self-int (양수로)
    int base_n = n - 1;
    
    for (int j = 0; j < base_n; ++j) {
        if (IF(n-1, j) != 0) {
            ext_attached_to = j;
            ext_int_num = IF(n-1, j);
            break;
        }
    }
    
    // === 2. Base topology 순회 ===
    std::vector<std::string> result_parts;
    std::vector<bool> used(base_n, false);
    
    int i = 0;  // 시작점
    
    while (i < base_n) {
        if (used[i]) {
            i++;
            continue;
        }
        
        used[i] = true;
        int diag_i = -IF(i, i);
        
        // a_{i,i+1} 체크
        if (i + 1 < base_n && IF(i, i+1) != 0 && !used[i+1]) {
            // Linear하게 계속
            std::string part = std::to_string(diag_i);
            
            // External이 여기 붙었으면
            if (ext_attached_to == i) {
                part = "\\overset{" + make_ext_marker(ext_self_int, ext_int_num) + "}{" + part + "}";
            }
            
            result_parts.push_back(part);
            i = i + 1;
        }
        else {
            // a_{i,i+1}=0 또는 i+1 이미 사용됨, jump 필요
            int j = -1;
            for (int k = i + 2; k < base_n; ++k) {
                if (IF(i, k) != 0 && !used[k]) {
                    j = k;
                    break;
                }
            }
            
            if (j == -1) {
                // 더 이상 연결 없음 - 마지막 curve
                std::string part = std::to_string(diag_i);
                if (ext_attached_to == i) {
                    part = "\\overset{" + make_ext_marker(ext_self_int, ext_int_num) + "}{" + part + "}";
                }
                result_parts.push_back(part);
                i++;
                continue;
            }
            
            // i 출력
            std::string part_i = std::to_string(diag_i);
            if (ext_attached_to == i) {
                part_i = "\\overset{" + make_ext_marker(ext_self_int, ext_int_num) + "}{" + part_i + "}";
            }
            result_parts.push_back(part_i);
            
            // 중간 curves (i+1 ~ j-1)를 overset chain으로
            // External이 stacked curve에 붙었는지 체크
            bool ext_on_stacked = false;
            int ext_stacked_idx = -1;
            for (int m = i + 1; m < j; ++m) {
                if (ext_attached_to == m) {
                    ext_on_stacked = true;
                    ext_stacked_idx = m;
                    break;
                }
            }
            
            std::string stacked = "";
            for (int m = i + 1; m < j; ++m) {
                used[m] = true;
                std::string mpart = std::to_string(-IF(m, m));
                
                // External이 stacked curve에 붙었으면 overset으로 위에
                if (ext_attached_to == m) {
                    mpart = "\\overset{" + make_ext_marker(ext_self_int, ext_int_num) + "}{" + mpart + "}";
                }
                stacked += mpart;
            }
            
            // j에 stacked를 overset으로 붙임
            int diag_j = -IF(j, j);
            std::string part_j = std::to_string(diag_j);
            
            if (!stacked.empty()) {
                part_j = "\\overset{" + stacked + "}{" + part_j + "}";
            }
            
            // External이 j에 붙었으면
            if (ext_attached_to == j) {
                // 이미 overset이 있으면 underset으로 아래에
                if (!stacked.empty()) {
                    part_j = "\\underset{" + make_ext_marker(ext_self_int, ext_int_num) + "}{" + part_j + "}";
                } else {
                    part_j = "\\overset{" + make_ext_marker(ext_self_int, ext_int_num) + "}{" + part_j + "}";
                }
            }
            
            result_parts.push_back(part_j);
            used[j] = true;
            i = j + 1;
        }
    }
    
    // === 3. 조합 (띄어쓰기 없이 붙여쓰기) ===
    std::ostringstream ss;
    ss << "{";
    for (size_t k = 0; k < result_parts.size(); ++k) {
        ss << result_parts[k];
    }
    ss << "}";
    
    return ss.str();
}

// ============================================================================
// Noble Molecules
// ============================================================================

struct NobleMolecule {
    std::string name;
    std::function<Tensor()> builder;
};

std::vector<NobleMolecule> get_noble_molecules() {
    return {
        // === ADE configurations of -2 curves ===
	// we don't need to arrange all of these type.. since there are too many of them
	/*
        {"A1", []{ Tensor t; t.AT(-2); return t; }},
        {"A2", []{ Tensor t; t.AT(-2); t.AT(-2); return t; }},
        {"A3", []{ Tensor t; for(int i=0;i<3;i++) t.AT(-2); return t; }},
        {"A4", []{ Tensor t; for(int i=0;i<4;i++) t.AT(-2); return t; }},
        {"A5", []{ Tensor t; for(int i=0;i<5;i++) t.AT(-2); return t; }},
        {"A6", []{ Tensor t; for(int i=0;i<6;i++) t.AT(-2); return t; }},
        {"A7", []{ Tensor t; for(int i=0;i<7;i++) t.AT(-2); return t; }},
        {"A8", []{ Tensor t; for(int i=0;i<8;i++) t.AT(-2); return t; }},
        */
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


	// noble molecules with -5 curves (from Appendix D)
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



	//additional noble molecules
	{"1{1,5}13151", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); return t; }},
	{"15132151", []{ Tensor t; t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); return t; }},
	{"21513151", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); return t; }},
	{"1513151315", []{ Tensor t; t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
	{"51232151", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); return t; }},
	




/*
 	// ADE series
	
	{"A", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); return t; }},
	{"D", []{ Tensor t; t.AT(-2); t.ATS(-2,-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); return t; }},
	{"E6", []{ Tensor t; t.AT(-2); t.AT(-2); t.ATS(-2,-2); t.AT(-2); t.AT(-2);  return t; }},
	{"E7", []{ Tensor t; t.AT(-2); t.AT(-2); t.ATS(-2,-2); t.AT(-2); t.AT(-2); t.AT(-2);  return t; }},
	{"E8", []{ Tensor t; t.AT(-2); t.AT(-2); t.ATS(-2,-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); return t; }},

*/

        

        
        // TODO: 논문 Appendix D에서 더 추가
    };
}

// ============================================================================
// Process Side Links (returns results instead of writing directly)
// ============================================================================

std::vector<LSTResult> collect_sidelink_results() {
    std::vector<LSTResult> results;
    
    for (int param : ALL_SIDELINK_PARAMS) {
        Tensor base = build_tensor(s(param));
        if (base.GetT() == 0) continue;
        
        auto IF = base.GetIntersectionForm();
        int n = IF.rows();
        
        // Try attaching external at each port
        for (int port = 1; port <= n; ++port) {  // 1-indexed for ATE
            // Get the self-intersection of the base curve at this port
            int base_self_int = IF(port - 1, port - 1);  // 0-indexed in matrix
            
            // Get allowed external curves for this base curve
            auto allowed_externals = get_allowed_external_curves(base_self_int);
            
            // Try only allowed external curve types
            for (int ext_curve : allowed_externals) {
                auto allowed = get_allowed_intersections(ext_curve);
                
                for (int int_num : allowed) {
                    Tensor t_ext;
                    t_ext.SetIF(IF);
                    t_ext.ATE(ext_curve, port, int_num);
                    
                    if (is_valid_LST(t_ext)) {
                        LSTResult result;
                        result.IF = t_ext.GetIntersectionForm();
                        result.latex = matrix_to_quiver_latex(result.IF);
                        
                        std::ostringstream desc;
                        desc << "sidelink=" << param 
                             << ", port=" << port 
                             << ", ext=" << ext_curve
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
// Process Noble Molecules (returns results instead of writing directly)
// ============================================================================

std::vector<LSTResult> collect_noble_molecule_results() {
    std::vector<LSTResult> results;
    
    for (const auto& nm : get_noble_molecules()) {
        Tensor base = nm.builder();
        if (base.GetT() == 0) continue;
        
        auto IF = base.GetIntersectionForm();
        int n = IF.rows();
        
        for (int port = 1; port <= n; ++port) {
            // Get the self-intersection of the base curve at this port
            int base_self_int = IF(port - 1, port - 1);  // 0-indexed in matrix
            
            // Get allowed external curves for this base curve
            auto allowed_externals = get_allowed_external_curves(base_self_int);
            
            // Try only allowed external curve types
            for (int ext_curve : allowed_externals) {
                auto allowed = get_allowed_intersections(ext_curve);
                
                for (int int_num : allowed) {
                    Tensor t_ext;
                    t_ext.SetIF(IF);
                    t_ext.ATE(ext_curve, port, int_num);
                    
                    if (is_valid_LST(t_ext)) {
                        LSTResult result;
                        result.IF = t_ext.GetIntersectionForm();
                        result.latex = matrix_to_quiver_latex(result.IF);
                        
                        std::ostringstream desc;
                        desc << "noble=" << nm.name
                             << ", port=" << port
                             << ", ext=" << ext_curve
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
// Output Functions
// ============================================================================

// Write intersection form matrix to stream in space-separated format
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

// Write all intersection forms to a file (YH format: matrices separated by blank lines)
void write_intersection_forms_file(const std::string& filename, 
                                   const std::vector<LSTResult>& results) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Cannot open " << filename << " for writing\n";
        return;
    }
    
    for (size_t i = 0; i < results.size(); ++i) {
        write_intersection_form(out, results[i].IF);
        if (i + 1 < results.size()) {
            out << "\n";  // Blank line between matrices
        }
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
        
        // Output eigenvalue spectrum
        data_out << "  Eigenvalues: [";
        for (size_t i = 0; i < r.spectrum.size(); ++i) {
            if (i > 0) data_out << ", ";
            data_out << r.spectrum[i];
        }
        data_out << "]\n";
    }
    
    // Output in align format (3 per row)
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

int main(int argc, char** argv) {
    std::string tex_filename = "no_node_LSTs.tex";
    std::string data_filename = "no_node_LSTs.txt";
    std::string if_filename = "no_node_LSTs_IF.txt";  // Intersection forms file
    
    if (argc >= 2) {
        tex_filename = argv[1];
    }
    if (argc >= 3) {
        data_filename = argv[2];
    }
    if (argc >= 4) {
        if_filename = argv[3];
    }
    
    std::ofstream tex_out(tex_filename);
    std::ofstream data_out(data_filename);
    
    if (!tex_out || !data_out) {
        std::cerr << "Error: Cannot open output files\n";
        return 1;
    }
    
    // ========== Collect all results ==========
    std::cout << "Collecting side link results...\n";
    std::vector<LSTResult> sidelink_results = collect_sidelink_results();
    std::cout << "Found " << sidelink_results.size() << " from side links (before dedup)\n";
    
    std::cout << "Collecting noble molecule results...\n";
    std::vector<LSTResult> noble_results = collect_noble_molecule_results();
    std::cout << "Found " << noble_results.size() << " from noble molecules (before dedup)\n";
    
    // ========== Combine all results ==========
    std::vector<LSTResult> all_results;
    all_results.reserve(sidelink_results.size() + noble_results.size());
    all_results.insert(all_results.end(), sidelink_results.begin(), sidelink_results.end());
    all_results.insert(all_results.end(), noble_results.begin(), noble_results.end());
    
    std::cout << "\nTotal before deduplication: " << all_results.size() << "\n";
    
    // ========== Remove duplicates by eigenvalue spectrum ==========
    std::cout << "Removing duplicates by eigenvalue spectrum...\n";
    std::vector<LSTResult> unique_results = remove_duplicates_by_spectrum(all_results);
    std::cout << "After deduplication: " << unique_results.size() << "\n";
    std::cout << "Removed " << (all_results.size() - unique_results.size()) << " duplicates\n\n";
    
    // ========== Write intersection forms file (YH format) ==========
    write_intersection_forms_file(if_filename, unique_results);
    
    // ========== Write to files ==========
    // LaTeX header (minimal - for input into another document)
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
