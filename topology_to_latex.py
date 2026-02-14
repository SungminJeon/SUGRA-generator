#!/usr/bin/env python3
"""
topology_to_latex.py
VERSION 4.0 - 2025-01-28

Converts SUGRA topology line format to LaTeX.
Fully ported from Mathematica version v3.3/v4.0

Usage:
    python topology_to_latex.py input.txt output.tex
    python topology_to_latex.py input.txt output.tex --no-sugra-info
    python topology_to_latex.py --test "LINE"
"""

import sys
import re
import argparse
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Any

# ============================================================================
# INSTANTON POSITION CONFIGURATION
# ============================================================================

INSTANTON_POSITIONS_SINGLE_NODE = {
    1: ["left"],
    2: ["left", "top"],
    3: ["left", "top", "bottom"],
    4: ["left", "top", "bottom", "right"],
}

INSTANTON_POSITIONS_LEFT_END = {
    1: ["top"],
    2: ["top", "bottom"],
    3: ["top", "bottom", "right"],
    4: ["top", "bottom", "right", "left"],
}

INSTANTON_POSITIONS_RIGHT_END = {
    1: ["top"],
    2: ["top", "bottom"],
    3: ["top", "bottom", "left"],
    4: ["top", "bottom", "left", "right"],
}

INSTANTON_POSITIONS_MIDDLE = {
    1: ["top"],
    2: ["top", "bottom"],
    3: ["top", "bottom", "left"],
    4: ["top", "bottom", "left", "right"],
}

# ============================================================================
# GAUGE ALGEBRA MAPPINGS
# ============================================================================

def gauge_algebra(param: int) -> str:
    """Node gauge algebra from param."""
    mapping = {
        4: r"\mathfrak{so}",
        6: r"\mathfrak{e}_6",
        7: r"\mathfrak{e}_7'",
        8: r"\mathfrak{e}_7",
        12: r"\mathfrak{e}_8",
    }
    return mapping.get(param, rf"\mathfrak{{g}}_{{{param}}}")

def external_gauge_to_latex(gauge_name: str) -> str:
    """External gauge name to LaTeX."""
    if not gauge_name:
        return ""
    mapping = {
        "none": "",
        "su2": r"\mathfrak{su}_2",
        "su3": r"\mathfrak{su}_3",
        "su4": r"\mathfrak{su}_4",
        "su5": r"\mathfrak{su}_5",
        "su6": r"\mathfrak{su}_6",
        "su7": r"\mathfrak{su}_7",
        "su8": r"\mathfrak{su}_8",
        "g2": r"\mathfrak{g}_2",
        "so7": r"\mathfrak{so}_7",
        "so8": r"\mathfrak{so}_8",
        "so9": r"\mathfrak{so}_9",
        "so10": r"\mathfrak{so}_{10}",
        "so11": r"\mathfrak{so}_{11}",
        "so12": r"\mathfrak{so}_{12}",
        "so13": r"\mathfrak{so}_{13}",
        "so14": r"\mathfrak{so}_{14}",
        "so16": r"\mathfrak{so}_{16}",
        "f4": r"\mathfrak{f}_4",
        "e6": r"\mathfrak{e}_6",
        "e7": r"\mathfrak{e}_7",
        "e8": r"\mathfrak{e}_8",
    }
    return mapping.get(gauge_name, gauge_name)

# ============================================================================
# CURVE SEQUENCE DEFINITIONS
# ============================================================================

SIDELINK_CURVES = {
    # === Instantons ===
    1: [["1"]],
    882: [["2"], ["1"]],
    883: [["2"], ["2"], ["1"]],
    884: [["2"], ["2"], ["2"], ["1"]],
    885: [["2"], ["2"], ["2"], ["2"], ["1"]],
    886: [["2"], ["2"], ["2"], ["2"], ["2"], ["1"]],
    887: [["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["1"]],
    8881: [["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["1"]],
    889: [["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["1"]],
    8810: [["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["1"]],
    8811: [["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["2"], ["1"]],
    
    # === alkali 1-links with no -5 curve ===
    91: [["3"], ["2", "2"], ["1"]],
    92: [["2"], ["3", "2"], ["1"]],
    93: [["3"], ["2"], ["2"], ["1"]],
    94: [["2"], ["3"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    95: [["2"], ["2"], ["3"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    96: [["3"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    97: [["3"], ["2"], ["1"]],
    98: [["2"], ["3"], ["2"], ["1"]],
    99: [["2"], ["3"], ["1"], ["3"], ["2"], ["1"]],
    910: [["2"], ["2"], ["3"], ["1"], ["3"], ["2"], ["1"]],
    911: [["3"], ["1"], ["3"], ["2"], ["1"]],
    912: [["3"], ["1"]],
    913: [["2"], ["3"], ["1"], ["3"], ["1"]],
    914: [["2"], ["2"], ["3"], ["1"], ["3"], ["1"]],
    915: [["3"], ["1"], ["3"], ["1"]],
    916: [["2"], ["3"], ["1"]],
    917: [["2"], ["2"], ["3"], ["1"]],
    
    # === alkali 2-links with no -5 curve ===
    991: [["2"], ["3", "1"], ["1"]],
    9920: [["1"], ["2"], ["3", "2"], ["1"]],
    9902: [["1"], ["3", "2"], ["2"], ["1"]],
    993: [["2"], ["3", "1"], ["2"], ["1"]],
    
    # === alkali 3-links with one -5 curve ===
    99910: [["1"], ["5", "1"], ["1"], ["3"], ["1"]],
    99920: [["1"], ["5", "1"], ["1"], ["3"], ["2"], ["1"]],
    99930: [["1"], ["5", "1"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    994: [["3"], ["1"], ["5", "1"], ["1"], ["3"], ["1"]],
    995: [["3"], ["1"], ["5", "1"], ["1"], ["3"], ["2"], ["1"]],
    996: [["3"], ["1"], ["5", "1"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    997: [["2"], ["3"], ["1"], ["5", "1"], ["1"], ["3"], ["2"], ["1"]],
    998: [["2"], ["3"], ["1"], ["5", "1"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    999: [["2"], ["2"], ["3"], ["1"], ["5", "1"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    9910: [["2"], ["3"], ["1"], ["5", "1"], ["1"], ["3"], ["1"]],
    9911: [["2"], ["2"], ["3"], ["1"], ["5", "1"], ["1"], ["3"], ["2"], ["1"]],
    9912: [["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    9913: [["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    9914: [["1"], ["5"], ["1"], ["2"], ["3"], ["2"], ["1"]],
    
    # === alkali 1-links with one -5 curve ===
    918: [["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    919: [["3"], ["2"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    920: [["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    921: [["2"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    922: [["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    923: [["2"], ["3"], ["2"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    924: [["5"], ["1"], ["3"], ["2"], ["1"]],
    925: [["5"], ["1"], ["2"], ["3"], ["2"], ["1"]],
    926: [["3"], ["2"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    927: [["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    928: [["2"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    929: [["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    930: [["2"], ["3"], ["2"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    931: [["2"], ["3"], ["1"], ["5"], ["1"], ["2"], ["3"], ["2"], ["1"]],
    932: [["2"], ["2"], ["3"], ["1"], ["5"], ["1"], ["2"], ["3"], ["2"], ["1"]],
    933: [["3"], ["1"], ["5"], ["1"], ["2"], ["3"], ["2"], ["1"]],
    934: [["5"], ["1"], ["3"], ["1"]],
    935: [["3"], ["2"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    936: [["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    937: [["2"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    938: [["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    939: [["2"], ["3"], ["2"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    940: [["5"], ["1"], ["2"], ["3"], ["1"]],
    941: [["1"], ["5"], ["1"], ["2"], ["3"], ["1"]],
    942: [["5"], ["1"], ["2"], ["2"], ["3"], ["1"]],
    943: [["2"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    944: [["2"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    945: [["2"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    
    # === alkali 2-links with two -5 curves ===
    9915: [["1"], ["5"], ["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    9916: [["1"], ["5"], ["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    9917: [["1"], ["5"], ["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    
    # === alkali 1-links with two -5 curves ===
    946: [["5"], ["1"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    947: [["5"], ["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    948: [["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    949: [["2"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    950: [["3"], ["1"], ["5"], ["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    951: [["5"], ["1"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    952: [["5"], ["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    953: [["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    954: [["3"], ["1"], ["5"], ["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    955: [["5"], ["1"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    956: [["5"], ["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    957: [["3"], ["1"], ["5"], ["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    958: [["1"], ["5"], ["1"], ["3"], ["1"]],
    
    # === Basic interior links (also used as sidelinks) ===
    11: [["1"]],
    22: [["1"], ["3"], ["1"]],
    23: [["1"], ["3"], ["2"], ["1"]],
    32: [["1"], ["2"], ["3"], ["1"]],
    33: [["1"], ["2"], ["3"], ["2"], ["1"]],
    24: [["1"], ["3"], ["2"], ["2"], ["1"]],
    42: [["1"], ["2"], ["2"], ["3"], ["1"]],
    34: [["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    43: [["1"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    44: [["1"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    35: [["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    53: [["1"], ["2"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    45: [["1"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    54: [["1"], ["2"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    55: [["1"], ["2"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    
    # === Special 331 connection ===
    331: [["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
}

INTERIOR_CURVES = {
    11: [["1"]],
    22: [["1"], ["3"], ["1"]],
    23: [["1"], ["3"], ["2"], ["1"]],
    32: [["1"], ["2"], ["3"], ["1"]],
    33: [["1"], ["2"], ["3"], ["2"], ["1"]],
    24: [["1"], ["3"], ["2"], ["2"], ["1"]],
    42: [["1"], ["2"], ["2"], ["3"], ["1"]],
    34: [["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    43: [["1"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    44: [["1"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    35: [["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    53: [["1"], ["2"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
    45: [["1"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    54: [["1"], ["2"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["1"]],
    55: [["1"], ["2"], ["2"], ["3"], ["1"], ["5"], ["1"], ["3"], ["2"], ["2"], ["1"]],
    331: [["1"], ["3"], ["1"], ["5"], ["1"], ["3"], ["1"]],
}

INSTANTON_ABBREV = {
    1: "1",
    882: r"I^{\oplus 2}",
    883: r"I^{\oplus 3}",
    884: r"I^{\oplus 4}",
    885: r"I^{\oplus 5}",
    886: r"I^{\oplus 6}",
    887: r"I^{\oplus 7}",
    8881: r"I^{\oplus 8}",
    889: r"I^{\oplus 9}",
    8810: r"I^{\oplus 10}",
    8811: r"I^{\oplus 11}",
}

def get_sidelink_curves(param: int) -> List[List[str]]:
    return SIDELINK_CURVES.get(param, [[str(param)]])

def get_interior_curves(param: int) -> List[List[str]]:
    return INTERIOR_CURVES.get(param, [[str(param)]])

def get_instanton_abbrev(param: int) -> str:
    return INSTANTON_ABBREV.get(param, str(param))

# ============================================================================
# PARSING
# ============================================================================

def parse_line(line: str) -> Optional[Dict[str, Any]]:
    """Parse topology line format."""
    parts = [p.strip() for p in line.split("|")]
    if len(parts) < 6:
        return None
    
    try:
        kinds = [int(x) for x in parts[0].split(",")]
        bparams = [int(x) for x in parts[1].split(",")]
        
        if len(kinds) != len(bparams):
            return None
        
        # Parse S connections
        s_conn = []
        if "(" in parts[2]:
            matches = re.findall(r"\((\d+),(\d+)\)", parts[2])
            s_conn = [(int(a), int(b)) for a, b in matches]
        
        # Parse I connections
        i_conn = []
        if "(" in parts[3]:
            matches = re.findall(r"\((\d+),(\d+)\)", parts[3])
            i_conn = [(int(a), int(b)) for a, b in matches]
        
        # Parse sp
        sp = []
        sp_match = re.search(r"sp=([0-9,]+)", parts[4])
        if sp_match and sp_match.group(1):
            sp = [int(x) for x in sp_match.group(1).split(",") if x]
        
        # Parse ip
        ip = []
        ip_match = re.search(r"ip=([0-9,]+)", parts[5])
        if ip_match and ip_match.group(1):
            ip = [int(x) for x in ip_match.group(1).split(",") if x]
        
        # Parse remaining parts
        e_conn = []
        ep = []
        eg = []
        anom = None
        sig = None
        det = None
        
        for part in parts[6:]:
            if "E=" in part:
                matches = re.findall(r"\((\d+),(\d+),(\d+),(\d+)\)", part)
                e_conn = [(int(a), int(b), int(c), int(d)) for a, b, c, d in matches]
            
            ep_match = re.search(r"ep=([0-9,]+)", part)
            if ep_match and ep_match.group(1):
                ep = [int(x) for x in ep_match.group(1).split(",") if x]
            
            eg_match = re.search(r"eg=([a-z0-9,]+)", part)
            if eg_match:
                eg = eg_match.group(1).split(",")
            
            anom_match = re.search(r"anom=(-?\d+)", part)
            if anom_match:
                anom = int(anom_match.group(1))
            
            sig_match = re.search(r"sig=\((\d+),(\d+),(\d+)\)", part)
            if sig_match:
                sig = (int(sig_match.group(1)), int(sig_match.group(2)), int(sig_match.group(3)))
            
            det_match = re.search(r"det=(-?\d+)", part)
            if det_match:
                det = int(det_match.group(1))
        
        return {
            "kinds": kinds,
            "bparams": bparams,
            "s_conn": s_conn,
            "i_conn": i_conn,
            "sp": sp,
            "ip": ip,
            "e_conn": e_conn,
            "ep": ep,
            "eg": eg,
            "anom": anom,
            "sig": sig,
            "det": det,
        }
    except Exception as e:
        print(f"Parse error: {e}", file=sys.stderr)
        return None

# ============================================================================
# RENDERING
# ============================================================================

def render_position(pos: List[str], ext_bottom: Optional[int], ext_top: Optional[int],
                   gauge_bottom: Optional[str] = None, gauge_top: Optional[str] = None) -> str:
    """Render a single position with optional externals and gauges."""
    base = pos[0]
    top = pos[1] if len(pos) >= 2 else None
    
    result = base
    
    if top is not None:
        result = rf"\overset{{{top}}}{{{result}}}"
    
    if ext_bottom is not None:
        ext_str = str(ext_bottom)
        if gauge_bottom and gauge_bottom != "none" and gauge_bottom != "":
            gauge_latex = external_gauge_to_latex(gauge_bottom)
            ext_str = rf"\overset{{{gauge_latex}}}{{{ext_str}}}"
        result = rf"\overset{{\textcolor{{red}}{{{ext_str}}}}}{{{result}}}"
    
    if ext_top is not None:
        ext_str = str(ext_top)
        if gauge_top and gauge_top != "none" and gauge_top != "":
            gauge_latex = external_gauge_to_latex(gauge_top)
            ext_str = rf"\overset{{{gauge_latex}}}{{{ext_str}}}"
        result = rf"\overset{{\textcolor{{red}}{{{ext_str}}}}}{{{result}}}"
    
    return result

def apply_curves_with_externals(curves: List[List[str]], externals: Dict[int, int],
                                direction: str = "left",
                                gauges: Dict[int, str] = None) -> str:
    """Apply externals to curve sequence."""
    gauges = gauges or {}
    
    # Calculate port ranges
    port_ranges = []
    port_idx = 0
    for pos in curves:
        if len(pos) >= 2:
            port_ranges.append((port_idx, port_idx + 1))
            port_idx += 2
        else:
            port_ranges.append((port_idx,))
            port_idx += 1
    
    # Render each position
    rendered = []
    for i, pos in enumerate(curves):
        ext_bottom = externals.get(port_ranges[i][0])
        gauge_bottom = gauges.get(port_ranges[i][0])
        ext_top = None
        gauge_top = None
        if len(pos) >= 2 and len(port_ranges[i]) >= 2:
            ext_top = externals.get(port_ranges[i][1])
            gauge_top = gauges.get(port_ranges[i][1])
        
        rendered.append(render_position(pos, ext_bottom, ext_top, gauge_bottom, gauge_top))
    
    if direction == "right":
        rendered = list(reversed(rendered))
    
    return "".join(rendered)

def render_interior(param: int, externals: Dict[int, int], gauges: Dict[int, str] = None) -> str:
    """Render interior connection."""
    gauges = gauges or {}
    
    if not externals:
        digits = list(str(param))
        if digits == ['3', '3', '1']:
            return r" \overset{3,3}{\bigcirc} "
        if len(digits) >= 2:
            return rf" \overset{{{','.join(digits)}}}{{\otimes}} "
        return rf" \overset{{{param}}}{{\otimes}} "
    
    curves = get_interior_curves(param)
    content = apply_curves_with_externals(curves, externals, "left", gauges)
    return f" {content} "

def render_side_or_inst(param: int, ext: Dict[int, int], direction: str,
                        gauges: Dict[int, str] = None) -> str:
    """Render sidelink or instanton."""
    gauges = gauges or {}
    
    # param < 0 means instanton
    if param < 0 and not ext:
        return get_instanton_abbrev(-param)
    
    curves = get_sidelink_curves(abs(param))
    return apply_curves_with_externals(curves, ext, direction, gauges)

def build_decorated_node(algebra: str, sidelinks: List, instantons: List,
                        sidelink_externals: Dict, instanton_externals: Dict,
                        node_exts: List, position: str, sp: List, ip: List,
                        sidelink_gauges: Dict = None, instanton_gauges: Dict = None) -> str:
    """Build decorated node with sidelinks and instantons."""
    sidelink_gauges = sidelink_gauges or {}
    instanton_gauges = instanton_gauges or {}
    
    # Build node external string with gauge
    node_ext_str = ""
    if node_exts:
        ext_parts = []
        for ext in node_exts:
            if isinstance(ext, (list, tuple)) and len(ext) >= 2:
                ext_val, ext_gauge = ext[0], ext[1]
            else:
                ext_val, ext_gauge = ext, "none"
            
            if ext_gauge and ext_gauge != "none":
                gauge_latex = external_gauge_to_latex(ext_gauge)
                ext_parts.append(rf"\textcolor{{red}}{{\overset{{{gauge_latex}}}{{{ext_val}}}}}")
            else:
                ext_parts.append(rf"\textcolor{{red}}{{{ext_val}}}")
        node_ext_str = ",".join(ext_parts)
    
    # Check decorations
    has_decorations = len(sidelinks) > 0 or len(instantons) > 0
    
    if not has_decorations:
        if node_ext_str:
            return rf"\underset{{{node_ext_str}}}{{{algebra}}}"
        return algebra
    
    # Sort sidelinks by param descending
    sorted_s = sorted(sidelinks, key=lambda x: x[0], reverse=True)
    
    # Initialize positions
    left_param, right_param, top_param, bottom_param = 0, 0, 0, 0
    left_ext, right_ext, top_ext, bottom_ext = {}, {}, {}, {}
    left_gauge, right_gauge, top_gauge, bottom_gauge = {}, {}, {}, {}
    
    # Select instanton position table
    if position == "single":
        pos_table = INSTANTON_POSITIONS_SINGLE_NODE
    elif position == "left":
        pos_table = INSTANTON_POSITIONS_LEFT_END
    elif position == "right":
        pos_table = INSTANTON_POSITIONS_RIGHT_END
    else:  # middle
        pos_table = INSTANTON_POSITIONS_MIDDLE
    
    # Assign instantons to positions
    inst_count = len(instantons)
    inst_positions = pos_table.get(inst_count, [])
    
    inst_params = [ip[inst_idx] if inst_idx < len(ip) else 1 for inst_idx in instantons]
    inst_exts = [instanton_externals.get(inst_idx, {}) for inst_idx in instantons]
    inst_gauges = [instanton_gauges.get(inst_idx, {}) for inst_idx in instantons]
    
    for i, pos_name in enumerate(inst_positions):
        if i >= len(inst_params):
            break
        if pos_name == "left" and left_param == 0:
            left_param = -inst_params[i]
            left_ext = inst_exts[i]
            left_gauge = inst_gauges[i]
        elif pos_name == "right" and right_param == 0:
            right_param = -inst_params[i]
            right_ext = inst_exts[i]
            right_gauge = inst_gauges[i]
        elif pos_name == "top" and top_param == 0:
            top_param = -inst_params[i]
            top_ext = inst_exts[i]
            top_gauge = inst_gauges[i]
        elif pos_name == "bottom" and bottom_param == 0:
            bottom_param = -inst_params[i]
            bottom_ext = inst_exts[i]
            bottom_gauge = inst_gauges[i]
    
    # Assign sidelinks to remaining positions
    idx = 0
    if position == "right":
        # Right node: right -> top -> bottom
        if idx < len(sorted_s) and right_param == 0:
            sl_param, sl_idx = sorted_s[idx]
            right_param = sl_param
            right_ext = sidelink_externals.get(sl_idx, {})
            right_gauge = sidelink_gauges.get(sl_idx, {})
            idx += 1
        if idx < len(sorted_s) and top_param == 0:
            sl_param, sl_idx = sorted_s[idx]
            top_param = sl_param
            top_ext = sidelink_externals.get(sl_idx, {})
            top_gauge = sidelink_gauges.get(sl_idx, {})
            idx += 1
        if idx < len(sorted_s) and bottom_param == 0:
            sl_param, sl_idx = sorted_s[idx]
            bottom_param = sl_param
            bottom_ext = sidelink_externals.get(sl_idx, {})
            bottom_gauge = sidelink_gauges.get(sl_idx, {})
            idx += 1
    elif position == "middle":
        # Middle node: top -> bottom only
        if idx < len(sorted_s) and top_param == 0:
            sl_param, sl_idx = sorted_s[idx]
            top_param = sl_param
            top_ext = sidelink_externals.get(sl_idx, {})
            top_gauge = sidelink_gauges.get(sl_idx, {})
            idx += 1
        if idx < len(sorted_s) and bottom_param == 0:
            sl_param, sl_idx = sorted_s[idx]
            bottom_param = sl_param
            bottom_ext = sidelink_externals.get(sl_idx, {})
            bottom_gauge = sidelink_gauges.get(sl_idx, {})
            idx += 1
    elif position == "single":
        # Single node: left -> top -> right -> bottom
        if idx < len(sorted_s) and left_param == 0:
            sl_param, sl_idx = sorted_s[idx]
            left_param = sl_param
            left_ext = sidelink_externals.get(sl_idx, {})
            left_gauge = sidelink_gauges.get(sl_idx, {})
            idx += 1
        if idx < len(sorted_s) and top_param == 0:
            sl_param, sl_idx = sorted_s[idx]
            top_param = sl_param
            top_ext = sidelink_externals.get(sl_idx, {})
            top_gauge = sidelink_gauges.get(sl_idx, {})
            idx += 1
        if idx < len(sorted_s) and right_param == 0:
            sl_param, sl_idx = sorted_s[idx]
            right_param = sl_param
            right_ext = sidelink_externals.get(sl_idx, {})
            right_gauge = sidelink_gauges.get(sl_idx, {})
            idx += 1
        if idx < len(sorted_s) and bottom_param == 0:
            sl_param, sl_idx = sorted_s[idx]
            bottom_param = sl_param
            bottom_ext = sidelink_externals.get(sl_idx, {})
            bottom_gauge = sidelink_gauges.get(sl_idx, {})
            idx += 1
    else:  # left end
        # Left node: left -> top -> bottom
        if idx < len(sorted_s) and left_param == 0:
            sl_param, sl_idx = sorted_s[idx]
            left_param = sl_param
            left_ext = sidelink_externals.get(sl_idx, {})
            left_gauge = sidelink_gauges.get(sl_idx, {})
            idx += 1
        if idx < len(sorted_s) and top_param == 0:
            sl_param, sl_idx = sorted_s[idx]
            top_param = sl_param
            top_ext = sidelink_externals.get(sl_idx, {})
            top_gauge = sidelink_gauges.get(sl_idx, {})
            idx += 1
        if idx < len(sorted_s) and bottom_param == 0:
            sl_param, sl_idx = sorted_s[idx]
            bottom_param = sl_param
            bottom_ext = sidelink_externals.get(sl_idx, {})
            bottom_gauge = sidelink_gauges.get(sl_idx, {})
            idx += 1
    
    # Build result
    result = algebra
    
    # Add node external
    if node_ext_str:
        result = rf"\underset{{{node_ext_str}}}{{{result}}}"
    
    # Add top
    if top_param != 0:
        top_content = render_side_or_inst(top_param, top_ext, "left", top_gauge)
        result = rf"\overset{{{top_content}}}{{{result}}}"
    
    # Add bottom
    if bottom_param != 0:
        bottom_content = render_side_or_inst(bottom_param, bottom_ext, "left", bottom_gauge)
        result = rf"\underset{{{bottom_content}}}{{{result}}}"
    
    # Add left
    if left_param != 0:
        left_content = render_side_or_inst(left_param, left_ext, "left", left_gauge)
        result = f"{left_content} {result}"
    
    # Add right (reversed direction)
    if right_param != 0:
        right_content = render_side_or_inst(right_param, right_ext, "right", right_gauge)
        result = f"{result} {right_content}"
    
    return result

def topology_to_latex(parsed: Dict[str, Any]) -> str:
    """Convert parsed topology to LaTeX."""
    kinds = parsed["kinds"]
    bparams = parsed["bparams"]
    s_conn = parsed["s_conn"]
    i_conn = parsed["i_conn"]
    sp = parsed["sp"]
    ip = parsed["ip"]
    e_conn = parsed["e_conn"]
    ep = parsed["ep"]
    eg = parsed["eg"]
    
    # Separate nodes and links
    nodes = []
    links = []
    node_indices = {}
    link_indices = {}
    
    for i, (kind, param) in enumerate(zip(kinds, bparams)):
        if kind == 0:
            node_indices[i] = len(nodes)
            nodes.append(param)
        elif kind == 1:
            link_indices[i] = len(links)
            links.append(param)
    
    # Build sidelinks per node: nodeIdx -> [(param, slIdx), ...]
    sidelinks_per_node = [[] for _ in nodes]
    for block_idx, sl_idx in s_conn:
        if block_idx in node_indices and sl_idx < len(sp):
            node_idx = node_indices[block_idx]
            sidelinks_per_node[node_idx].append((sp[sl_idx], sl_idx))
    
    # Build instantons per node
    instantons_per_node = [[] for _ in nodes]
    for block_idx, inst_idx in i_conn:
        if block_idx in node_indices and inst_idx < len(ip):
            node_idx = node_indices[block_idx]
            instantons_per_node[node_idx].append(inst_idx)
    
    # Build external gauge lookup
    external_gauges = {i: g for i, g in enumerate(eg)}
    
    # Build external connections
    sidelink_externals_per_node = [{} for _ in nodes]
    instanton_externals_per_node = [{} for _ in nodes]
    interior_externals = [{} for _ in links]
    node_externals = [[] for _ in nodes]
    
    sidelink_gauges_per_node = [{} for _ in nodes]
    instanton_gauges_per_node = [{} for _ in nodes]
    interior_gauges = [{} for _ in links]
    
    for parent_id, parent_type, port_idx, ext_id in e_conn:
        ext_val = ep[ext_id] if ext_id < len(ep) else 0
        ext_gauge = external_gauges.get(ext_id, "none")
        
        if parent_type == 0:  # Block/Interior
            if parent_id in link_indices:
                link_idx = link_indices[parent_id]
                interior_externals[link_idx][port_idx] = ext_val
                interior_gauges[link_idx][port_idx] = ext_gauge
            elif parent_id in node_indices:
                node_idx = node_indices[parent_id]
                node_externals[node_idx].append((ext_val, ext_gauge))
        
        elif parent_type == 1:  # SideLink
            for block_idx, sl_idx in s_conn:
                if sl_idx == parent_id and block_idx in node_indices:
                    node_idx = node_indices[block_idx]
                    if sl_idx not in sidelink_externals_per_node[node_idx]:
                        sidelink_externals_per_node[node_idx][sl_idx] = {}
                        sidelink_gauges_per_node[node_idx][sl_idx] = {}
                    sidelink_externals_per_node[node_idx][sl_idx][port_idx] = ext_val
                    sidelink_gauges_per_node[node_idx][sl_idx][port_idx] = ext_gauge
                    break
        
        elif parent_type == 2:  # Instanton
            for block_idx, inst_idx in i_conn:
                if inst_idx == parent_id and block_idx in node_indices:
                    node_idx = node_indices[block_idx]
                    if inst_idx not in instanton_externals_per_node[node_idx]:
                        instanton_externals_per_node[node_idx][inst_idx] = {}
                        instanton_gauges_per_node[node_idx][inst_idx] = {}
                    instanton_externals_per_node[node_idx][inst_idx][port_idx] = ext_val
                    instanton_gauges_per_node[node_idx][inst_idx][port_idx] = ext_gauge
                    break
    
    if not nodes:
        return r"\text{empty}"
    
    # Single node case
    if len(nodes) == 1:
        algebra = gauge_algebra(nodes[0])
        decorated = build_decorated_node(
            algebra, 
            sidelinks_per_node[0], 
            instantons_per_node[0],
            sidelink_externals_per_node[0],
            instanton_externals_per_node[0],
            node_externals[0],
            "single", sp, ip,
            sidelink_gauges_per_node[0],
            instanton_gauges_per_node[0]
        )
        return "{" + decorated + "}"
    
    latex = "{"
    
    # First node
    algebra = gauge_algebra(nodes[0])
    decorated = build_decorated_node(
        algebra,
        sidelinks_per_node[0],
        instantons_per_node[0],
        sidelink_externals_per_node[0],
        instanton_externals_per_node[0],
        node_externals[0],
        "left", sp, ip,
        sidelink_gauges_per_node[0],
        instanton_gauges_per_node[0]
    )
    latex += decorated
    
    # Process remaining nodes
    for i in range(1, len(nodes)):
        link_param = links[i - 1]
        node_param = nodes[i]
        link_exts = interior_externals[i - 1]
        link_gauges = interior_gauges[i - 1]
        
        position = "right" if i == len(nodes) - 1 else "middle"
        
        latex += render_interior(link_param, link_exts, link_gauges)
        
        algebra = gauge_algebra(node_param)
        decorated = build_decorated_node(
            algebra,
            sidelinks_per_node[i],
            instantons_per_node[i],
            sidelink_externals_per_node[i],
            instanton_externals_per_node[i],
            node_externals[i],
            position, sp, ip,
            sidelink_gauges_per_node[i],
            instanton_gauges_per_node[i]
        )
        latex += decorated
    
    latex += "}"
    return latex

# ============================================================================
# FILE PROCESSING
# ============================================================================

def process_file(input_file: str, output_file: str, include_sugra_info: bool = True):
    """Process input file and generate LaTeX output."""
    with open(input_file, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    
    print(f"Processing {len(lines)} topologies...")
    
    results = []
    for i, line in enumerate(lines):
        if (i + 1) % 100 == 0:
            print(f"Processed {i + 1} / {len(lines)}")
        
        parsed = parse_line(line)
        if parsed is None:
            results.append(("% Parse error: " + line[:50], None, []))
            continue
        
        latex = topology_to_latex(parsed)
        
        nodes = [bparams for kind, bparams in zip(parsed["kinds"], parsed["bparams"]) if kind == 0]
        signature = tuple(sorted(nodes))
        
        sugra_info = {
            "anom": parsed["anom"],
            "sig": parsed["sig"],
            "det": parsed["det"],
        }
        
        results.append((latex, sugra_info, signature))
    
    # Group by signature
    groups = defaultdict(list)
    for latex, sugra_info, sig in results:
        groups[sig].append((latex, sugra_info))
    
    sorted_groups = sorted(groups.items(), key=lambda x: -len(x[1]))
    
    # Build output
    output = []
    output.append("% Generated LaTeX file with External curves (in RED)")
    output.append("% Use with amsmath and xcolor packages")
    output.append("% \\usepackage{xcolor}")
    output.append(f"% Total unique: {len(set(r[0] for r in results))}")
    output.append(f"% Total groups: {len(sorted_groups)}")
    if include_sugra_info:
        output.append("% SUGRA info: Delta (anomaly residual), n_+/n_-/n_0 (signature), det")
    output.append("")
    
    for i, (sig, items) in enumerate(sorted_groups, 1):
        seen = set()
        unique_items = []
        for latex, sugra_info in items:
            if latex not in seen:
                seen.add(latex)
                unique_items.append((latex, sugra_info))
        
        output.append(f"% === Group {i}: nodes = {list(sig)} ({len(unique_items)} items) ===")
        output.append(r"\begin{align*}")
        
        for j, (latex, sugra_info) in enumerate(unique_items):
            info_str = ""
            if include_sugra_info and sugra_info.get("anom") is not None:
                info_str = rf" \quad (\Delta={sugra_info['anom']}"
                if sugra_info.get("sig"):
                    s = sugra_info["sig"]
                    info_str += rf", n_{{+}}={s[0]}, n_{{-}}={s[1]}, n_{{0}}={s[2]}"
                if sugra_info.get("det") is not None:
                    info_str += rf", |\det|={abs(sugra_info['det'])}"
                info_str += ")"
            
            prefix = "&" if j == 0 else r" \\" + "\n&"
            output.append(f"{prefix}{latex}{info_str}")
        
        output.append(r"\end{align*}")
        output.append("")
    
    with open(output_file, 'w') as f:
        f.write("\n".join(output))
    
    print(f"Done! Output saved to: {output_file}")
    print(f"Total groups: {len(sorted_groups)}")

def test_single(line: str):
    """Test single line conversion."""
    parsed = parse_line(line)
    if parsed is None:
        print("Parse failed!")
        return
    
    print("Parsed:")
    for k, v in parsed.items():
        print(f"  {k}: {v}")
    
    latex = topology_to_latex(parsed)
    print(f"\nLaTeX:\n{latex}")

# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Convert topology to LaTeX")
    parser.add_argument("input", nargs="?", help="Input file")
    parser.add_argument("output", nargs="?", help="Output file")
    parser.add_argument("--no-sugra-info", action="store_true", help="Exclude SUGRA info")
    parser.add_argument("--test", metavar="LINE", help="Test single line")
    
    args = parser.parse_args()
    
    if args.test:
        test_single(args.test)
    elif args.input and args.output:
        process_file(args.input, args.output, not args.no_sugra_info)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
