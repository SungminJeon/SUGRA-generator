// anomaly_tables.h
// Anomaly contribution tables for 6D SUGRA classification
// Extracted from sugra_generator.cpp
//
// H: Hypermultiplet contribution
// V: Vector multiplet contribution
// Gravitational anomaly condition: H - V + 29T = 273

#pragma once

struct AnomalyContrib {
    int H;  // Hypermultiplet contribution
    int V;  // Vector multiplet contribution
};

// ============================================================================
// NODE CONTRIBUTIONS (LKind::g)
// param = |self-intersection| -> (H, V)
// H from nodes is typically 0; V from gauge algebra dimension
// ============================================================================

inline AnomalyContrib get_node_contrib(int param) {
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

// ============================================================================
// INTERIOR LINK CONTRIBUTIONS (LKind::L)
// param encodes link type (e.g., 33 = (3,3) link)
// H from bifundamentals, V from internal gauge curves
// ============================================================================

inline AnomalyContrib get_interior_contrib(int param) {
    switch (param) {
        // (1,1)
        case 11: return {0, 0};

        // Diagonal links
        case 22: return {0, 8};
        case 33: return {16, 27};
        case 331: return {0, 68};  // (3,3)' with -5
        case 44: return {16, 86};
        case 55: return {16, 86};

        // (2,3) type links
        case 23: return {8, 17};
        case 32: return {8, 17};

        // (2,4) type links
        case 24: return {8, 17};
        case 42: return {8, 17};

        // (3,4) type links
        case 34: return {8, 77};
        case 43: return {8, 77};

        // (3,5) type links
        case 35: return {8, 77};
        case 53: return {8, 77};

        // (4,5) type links
        case 45: return {16, 86};
        case 54: return {16, 86};

        default: return {0, 0};
    }
}

// ============================================================================
// SIDELINK CONTRIBUTIONS (LKind::S)
// Based on atomic classification appendix
// ============================================================================

inline AnomalyContrib get_sidelink_contrib(int param) {
    switch (param) {
        // =============================================================
        // INSTANTONS (88x notation - blowdown induced)
        // =============================================================
        case 1: return {0, 0};
        case 882: return {0, 0};
        case 883: return {0, 0};
        case 884: return {0, 0};
        case 885: return {0, 0};
        case 886: return {0, 0};
        case 887: return {0, 0};
        case 8881: return {0, 0};
        case 889: return {0, 0};
        case 8810: return {0, 0};
        case 8811: return {0, 0};

        // =============================================================
        // INTERIORS (same as interior links, for compatibility)
        // =============================================================
        // (1,1)
        case 11: return {0, 0};

        case 22: return {0, 8};
        case 33: return {16, 27};
        case 331: return {0, 68};
        case 44: return {16, 86};
        case 55: return {16, 86};

        // (2,3) type links
        case 23: return {8, 17};
        case 32: return {8, 17};

        // (2,4) type links
        case 24: return {8, 17};
        case 42: return {8, 17};

        // (3,4) type links
        case 34: return {8, 77};
        case 43: return {8, 77};

        // (3,5) type links
        case 35: return {8, 77};
        case 53: return {8, 77};

        // (4,5) type links
        case 45: return {16, 86};
        case 54: return {16, 86};

        // =============================================================
        // ALKALI 2 LINKS WITH NO -5
        // =============================================================
        case 991: return {8, 17};
        case 9920: return {16, 27};
        case 9902: return {16, 27};
        case 993: return {16, 27};

        // =============================================================
        // ALKALI 1 LINKS WITH NO -5
        // =============================================================
        case 91: return {8, 17};
        case 92: return {16, 27};
        case 93: return {8, 17};
        case 94: return {16, 34};
        case 95: return {16, 24};
        case 96: return {8, 25};
        case 97: return {8, 17};
        case 98: return {16, 27};
        case 99: return {16, 34};
        case 910: return {16, 34};
        case 911: return {8, 25};
        case 912: return {0, 8};
        case 913: return {8, 25};
        case 914: return {8, 25};
        case 915: return {0, 16};
        case 916: return {8, 17};
        case 917: return {8, 17};

        // =============================================================
        // ALKALI 3 LINKS WITH ONE -5 CURVE
        // =============================================================
        case 99910: return {0, 60};
        case 99920: return {8, 69};
        case 99930: return {8, 69};

        // =============================================================
        // ALKALI 2 LINKS WITH ONE -5 CURVE
        // =============================================================
        case 994: return {0, 68};
        case 995: return {8, 77};
        case 996: return {8, 77};
        case 997: return {16, 86};
        case 998: return {16, 86};
        case 999: return {16, 86};
        case 9910: return {8, 77};
        case 9911: return {16, 86};
        case 9912: return {8, 69};
        case 9913: return {8, 69};
        case 9914: return {16, 79};

        // =============================================================
        // ALKALI 1 LINKS WITH ONE -5 CURVE
        // =============================================================
        case 918: return {8, 69};
        case 919: return {16, 86};
        case 920: return {16, 86};
        case 921: return {16, 86};
        case 922: return {8, 77};
        case 923: return {24, 96};
        case 924: return {8, 69};
        case 925: return {16, 79};
        case 926: return {16, 86};
        case 927: return {16, 86};
        case 928: return {16, 86};
        case 929: return {8, 77};
        case 930: return {24, 96};
        case 931: return {24, 96};
        case 932: return {24, 96};
        case 933: return {16, 87};
        case 934: return {0, 60};
        case 935: return {8, 77};
        case 936: return {8, 77};
        case 937: return {8, 77};
        case 938: return {0, 68};
        case 939: return {16, 87};
        case 940: return {8, 69};
        case 941: return {8, 69};
        case 942: return {8, 69};
        case 943: return {0, 60};
        case 944: return {8, 69};
        case 945: return {8, 69};

        // =============================================================
        // ALKALI 2 LINKS WITH TWO -5 CURVES
        // =============================================================
        case 9915: return {8, 129};
        case 9916: return {8, 129};
        case 9917: return {0, 120};

        // =============================================================
        // ALKALI 1 LINKS WITH TWO -5 CURVES
        // =============================================================
        case 946: return {16, 138};
        case 947: return {8, 129};
        case 948: return {16, 146};
        case 949: return {16, 146};
        case 950: return {8, 137};
        case 951: return {16, 138};
        case 952: return {8, 129};
        case 953: return {16, 146};
        case 954: return {8, 137};
        case 955: return {8, 129};
        case 956: return {0, 120};
        case 957: return {0, 128};

        // =============================================================
        // ALKALI 1 LINK WITH ONE -5 CURVE (omitted in appendix)
        // =============================================================
        case 958: return {0, 60};

        default: return {0, 0};
    }
}

// ============================================================================
// INSTANTON CONTRIBUTIONS (LKind::I)
// Minimal instantons typically have no gauge algebra
// ============================================================================

inline AnomalyContrib get_instanton_contrib(int param) {
    (void)param;  // unused
    return {0, 0};
}

// ============================================================================
// EXTERNAL CONTRIBUTIONS (LKind::E)
// param = |self-intersection| of external curve
// ============================================================================

inline AnomalyContrib get_external_contrib(int param) {
    switch (param) {
        case 1:  return {0, 0};
        case 2:  return {0, 0};
        case 3:  return {0, 8};    // su(3)
        case 4:  return {0, 28};   // so(8)
        case 5:  return {0, 52};   // f4
        case 6:  return {0, 78};   // e6
        case 7:  return {0, 133};  // e7
        case 8:  return {0, 133};  // e7
        case 12: return {0, 248};  // e8
        default: return {0, 0};
    }
}
