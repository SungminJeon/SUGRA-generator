(* ::Package:: *)

(* ================================================================== *)
(* topology_to_latex_external.m *)
(* VERSION 3.3 - 2025-01-05 *)
(* External curves are shown in RED above each curve digit *)
(* Port order: left->right, bottom->top at each position *)
(* ENHANCED: Multiple instanton support with configurable positions *)
(* BUGFIX: Externals are OUTERMOST! Order: base -> top -> extBottom -> extTop *)
(* BUGFIX v3.1: Index-based sidelink assignment (instanton no longer blocks) *)
(* BUGFIX v3.1: Middle node uses top->bottom only (not left) *)
(* NEW v3.2: Repeating pattern compression: node (otimes{link} node)^n *)
(* NEW v3.2: Single node uses "single" position for full 4-position sidelink *)
(* NEW v3.3: SUGRA info support (anom, sig, det) in output *)
(* ================================================================== *)

Print[""];
Print["==============================================="];
Print["  VERSION 3.3 LOADED"];
Print["  Externals OUTERMOST!"];
Print["  Order: base -> top -> extBottom -> extTop"];
Print["  FIXED: Index-based sidelink assignment"];
Print["  FIXED: Middle node top->bottom only"];
Print["  NEW v3.2: Repeating pattern compression"];
Print["  NEW v3.2: Single node full 4-position support"];
Print["  NEW v3.3: SUGRA info (anom, sig, det) support"];
Print["  Date: 2025-01-05"];
Print["==============================================="];
Print[""];


(* ================================================================== *)
(* ===== INSTANTON POSITION CONFIGURATION ===== *)
(* ✨ EDIT HERE to change instanton placement rules *)
(* Positions: "left", "right", "top", "bottom" *)
(* ================================================================== *)

(* Single node (no interior links): instanton positions by count *)
(* Format: count -> {1st position, 2nd position, 3rd position, 4th position} *)
instantonPositionsSingleNode = <|
    1 -> {"left"},
    2 -> {"left", "top"},
    3 -> {"left", "top", "bottom"},
    4 -> {"left", "top", "bottom", "right"}
|>;

(* Left-end node (interior link on right): instanton positions by count *)
(* ✅ v3.1: Instanton goes top first, so sidelink can take left *)
instantonPositionsLeftEnd = <|
    1 -> {"top"},
    2 -> {"top", "bottom"},
    3 -> {"top", "bottom", "right"},
    4 -> {"top", "bottom", "right", "left"}
|>;

(* Right-end node (interior link on left): instanton positions by count *)
(* ✅ v3.1: Instanton goes top first, so sidelink can take right *)
instantonPositionsRightEnd = <|
    1 -> {"top"},
    2 -> {"top", "bottom"},
    3 -> {"top", "bottom", "left"},
    4 -> {"top", "bottom", "left", "right"}
|>;

(* Middle node (interior links on both sides): instanton positions by count *)
instantonPositionsMiddle = <|
    1 -> {"top"},
    2 -> {"top", "bottom"},
    3 -> {"top", "bottom", "left"},
    4 -> {"top", "bottom", "left", "right"}
|>;

(* ===== Gauge Algebra Mapping ===== *)

gaugeAlgebra[param_Integer] := Switch[param,
    4, "\\mathfrak{so}",
    6, "\\mathfrak{e}_6",
    7, "\\mathfrak{e}_7'",
    8, "\\mathfrak{e}_7",
    12, "\\mathfrak{e}_8",
    _, "\\mathfrak{g}_{" <> ToString[param] <> "}"
];

(* ===== Parse Line-Compact Format with External and SUGRA info ===== *)

parseLine[line_String] := Module[{parts, kinds, bparams, sConn, iConn, eConn, sp, ip, ep,
                                   anomVal, sigVal, detVal, sugraInfo},
    parts = StringSplit[line, "|"];
    parts = StringTrim /@ parts;
    
    If[Length[parts] < 6, Return[$Failed]];
    
    (* Parse kinds and bparams *)
    kinds = ToExpression /@ StringSplit[parts[[1]], ","];
    bparams = ToExpression /@ StringSplit[parts[[2]], ","];
    
    If[Length[kinds] != Length[bparams], Return[$Failed]];
    
    (* Parse S connections *)
    sConn = {};
    If[StringContainsQ[parts[[3]], "("],
        sConn = StringCases[parts[[3]], 
            "(" ~~ a:NumberString ~~ "," ~~ b:NumberString ~~ ")" :> 
            {ToExpression[a], ToExpression[b]}];
    ];
    
    (* Parse I connections *)
    iConn = {};
    If[StringContainsQ[parts[[4]], "("],
        iConn = StringCases[parts[[4]], 
            "(" ~~ a:NumberString ~~ "," ~~ b:NumberString ~~ ")" :> 
            {ToExpression[a], ToExpression[b]}];
    ];
    
    (* Parse sp *)
    sp = {};
    If[StringContainsQ[parts[[5]], "="],
        With[{spStr = StringTrim[StringReplace[parts[[5]], "sp=" -> ""]]},
            If[StringLength[spStr] > 0,
                sp = ToExpression /@ StringSplit[spStr, ","];
            ];
        ];
    ];
    
    (* Parse ip *)
    ip = {};
    If[StringContainsQ[parts[[6]], "="],
        With[{ipStr = StringTrim[StringReplace[parts[[6]], "ip=" -> ""]]},
            If[StringLength[ipStr] > 0,
                ip = ToExpression /@ StringSplit[ipStr, ","];
            ];
        ];
    ];
    
    (* Parse E connections: E=(parent_id, parent_type, port_idx, external_id) *)
    eConn = {};
    Do[
        If[StringContainsQ[parts[[i]], "E="],
            eConn = StringCases[parts[[i]], 
                "(" ~~ a:NumberString ~~ "," ~~ b:NumberString ~~ "," ~~ 
                    c:NumberString ~~ "," ~~ d:NumberString ~~ ")" :> 
                {ToExpression[a], ToExpression[b], ToExpression[c], ToExpression[d]}];
        ],
        {i, 7, Length[parts]}
    ];
    
    (* Parse ep *)
    ep = {};
    Do[
        If[StringContainsQ[parts[[i]], "ep="],
            With[{epStr = StringTrim[StringReplace[parts[[i]], "ep=" -> ""]]},
                If[StringLength[epStr] > 0,
                    ep = ToExpression /@ StringSplit[epStr, ","];
                ];
            ];
        ],
        {i, 7, Length[parts]}
    ];
    
    (* Parse SUGRA info: anom=, sig=(+,-,0), det= *)
    anomVal = None;
    sigVal = None;
    detVal = None;
    
    Do[
        (* anom= *)
        If[StringContainsQ[parts[[i]], "anom="],
            With[{anomStr = StringCases[parts[[i]], "anom=" ~~ x:(NumberString | ("-" ~~ NumberString)) :> x]},
                If[Length[anomStr] > 0, anomVal = ToExpression[anomStr[[1]]]];
            ];
        ];
        (* sig=(+,-,0) *)
        If[StringContainsQ[parts[[i]], "sig="],
            With[{sigStr = StringCases[parts[[i]], 
                "sig=(" ~~ a:NumberString ~~ "," ~~ b:NumberString ~~ "," ~~ c:NumberString ~~ ")" :> 
                {ToExpression[a], ToExpression[b], ToExpression[c]}]},
                If[Length[sigStr] > 0, sigVal = sigStr[[1]]];
            ];
        ];
        (* det= *)
        If[StringContainsQ[parts[[i]], "det="],
            With[{detStr = StringCases[parts[[i]], "det=" ~~ x:(NumberString | ("-" ~~ NumberString)) :> x]},
                If[Length[detStr] > 0, detVal = ToExpression[detStr[[1]]]];
            ];
        ],
        {i, 1, Length[parts]}
    ];
    
    sugraInfo = <|"anom" -> anomVal, "sig" -> sigVal, "det" -> detVal|>;
    
    {kinds, bparams, sConn, iConn, sp, ip, eConn, ep, sugraInfo}
];

(* ================================================================== *)
(* ===== CURVE SEQUENCE DEFINITIONS ===== *)
(* Format: { {base}, {base, top}, {base}, ... } *)
(* Each position: {bottom_curve} or {bottom_curve, top_curve} *)
(* Port order: left\[RightArrow]right, bottom\[RightArrow]top at each position *)
(* ================================================================== *)

(* ===== SideLink/Instanton Curve Sequences ===== *)
(* \[FivePointedStar] \:d558\:b4dc\:cf54\:b529: \:c67c\:cabd\[RightArrow]\:c624\:b978\:cabd, \:ac01 \:c704\:ce58\:c5d0\:c11c \:c544\:b798\[RightArrow]\:c704 \[FivePointedStar] *)

getSideLinkCurves[param_Integer] := Switch[param,
    (* === Instantons === *)
    1,    {{"1"}},
    882,  {{"2"}, {"1"}},
    883,  {{"2"}, {"2"}, {"1"}},
    884,  {{"2"}, {"2"}, {"2"}, {"1"}},
    885,  {{"2"}, {"2"}, {"2"}, {"2"}, {"1"}},
    886,  {{"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"1"}},
    887,  {{"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"1"}},
    8881, {{"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"1"}},
    889,  {{"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"1"}},
    8810, {{"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"1"}},
    8811, {{"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"2"}, {"1"}},
    
    (* === alkali 1-links with no -5 curve === *)
    91,  {{"3"}, {"2", "2"}, {"1"}},                    (* 3\overset{2}{2}1 *)
    92,  {{"2"}, {"3", "2"}, {"1"}},                    (* 2\overset{2}{3}1 *)
    93,  {{"3"}, {"2"}, {"2"}, {"1"}},                  (* 3221 *)
    94,  {{"2"}, {"3"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    95,  {{"2"}, {"2"}, {"3"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    96,  {{"3"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    97,  {{"3"}, {"2"}, {"1"}},
    98,  {{"2"}, {"3"}, {"2"}, {"1"}},
    99,  {{"2"}, {"3"}, {"1"}, {"3"}, {"2"}, {"1"}},
    910, {{"2"}, {"2"}, {"3"}, {"1"}, {"3"}, {"2"}, {"1"}},
    911, {{"3"}, {"1"}, {"3"}, {"2"}, {"1"}},
    912, {{"3"}, {"1"}},
    913, {{"2"}, {"3"}, {"1"}, {"3"}, {"1"}},
    914, {{"2"}, {"2"}, {"3"}, {"1"}, {"3"}, {"1"}},
    915, {{"3"}, {"1"}, {"3"}, {"1"}},
    916, {{"2"}, {"3"}, {"1"}},
    917, {{"2"}, {"2"}, {"3"}, {"1"}},
    
    (* === alkali 2-links with no -5 curve === *)
    991,  {{"2"}, {"3", "1"}, {"1"}},                   (* 2\overset{1}{3}1 *)
    9920, {{"1"}, {"2"}, {"3", "2"}, {"1"}},            (* 12\overset{2}{3}1 *)
    9902, {{"1"}, {"3", "2"}, {"2"}, {"1"}},            (* 1\overset{2}{3}21 *)
    993,  {{"2"}, {"3", "1"}, {"2"}, {"1"}},            (* 2\overset{1}{3}21 *)
    
    (* === alkali 3-links with one -5 curve === *)
    99910, {{"1"}, {"5", "1"}, {"1"}, {"3"}, {"1"}},    (* 1\overset{1}{5}131 *)
    99920, {{"1"}, {"5", "1"}, {"1"}, {"3"}, {"2"}, {"1"}},
    99930, {{"1"}, {"5", "1"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    994,   {{"3"}, {"1"}, {"5", "1"}, {"1"}, {"3"}, {"1"}},
    995,   {{"3"}, {"1"}, {"5", "1"}, {"1"}, {"3"}, {"2"}, {"1"}},
    996,   {{"3"}, {"1"}, {"5", "1"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    997,   {{"2"}, {"3"}, {"1"}, {"5", "1"}, {"1"}, {"3"}, {"2"}, {"1"}},
    998,   {{"2"}, {"3"}, {"1"}, {"5", "1"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    999,   {{"2"}, {"2"}, {"3"}, {"1"}, {"5", "1"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    9910,  {{"2"}, {"3"}, {"1"}, {"5", "1"}, {"1"}, {"3"}, {"1"}},
    9911,  {{"2"}, {"2"}, {"3"}, {"1"}, {"5", "1"}, {"1"}, {"3"}, {"2"}, {"1"}},
    9912,  {{"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    9913,  {{"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    9914,  {{"1"}, {"5"}, {"1"}, {"2"}, {"3"}, {"2"}, {"1"}},
    
    (* === alkali 1-links with one -5 curve === *)
    918, {{"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    919, {{"3"}, {"2"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    920, {{"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    921, {{"2"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    922, {{"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    923, {{"2"}, {"3"}, {"2"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    924, {{"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    925, {{"5"}, {"1"}, {"2"}, {"3"}, {"2"}, {"1"}},
    926, {{"3"}, {"2"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    927, {{"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    928, {{"2"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    929, {{"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    930, {{"2"}, {"3"}, {"2"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    931, {{"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"2"}, {"3"}, {"2"}, {"1"}},
    932, {{"2"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"2"}, {"3"}, {"2"}, {"1"}},
    933, {{"3"}, {"1"}, {"5"}, {"1"}, {"2"}, {"3"}, {"2"}, {"1"}},
    934, {{"5"}, {"1"}, {"3"}, {"1"}},
    935, {{"3"}, {"2"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    936, {{"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    937, {{"2"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    938, {{"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    939, {{"2"}, {"3"}, {"2"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    940, {{"5"}, {"1"}, {"2"}, {"3"}, {"1"}},
    941, {{"1"}, {"5"}, {"1"}, {"2"}, {"3"}, {"1"}},
    942, {{"5"}, {"1"}, {"2"}, {"2"}, {"3"}, {"1"}},
    943, {{"2"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    944, {{"2"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    945, {{"2"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    
    (* === alkali 2-links with two -5 curves === *)
    9915, {{"1"}, {"5"}, {"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    9916, {{"1"}, {"5"}, {"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    9917, {{"1"}, {"5"}, {"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    
    (* === alkali 1-links with two -5 curves === *)
    946, {{"5"}, {"1"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    947, {{"5"}, {"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    948, {{"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    949, {{"2"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    950, {{"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    951, {{"5"}, {"1"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    952, {{"5"}, {"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    953, {{"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    954, {{"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    955, {{"5"}, {"1"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    956, {{"5"}, {"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    957, {{"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    958, {{"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    
    11,  {{"1"}},
    22,  {{"1"}, {"3"}, {"1"}},
    23,  {{"1"}, {"3"}, {"2"}, {"1"}},
    32,  {{"1"}, {"2"}, {"3"}, {"1"}},
    33,  {{"1"}, {"2"}, {"3"}, {"2"}, {"1"}},
    24,  {{"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    42,  {{"1"}, {"2"}, {"2"}, {"3"}, {"1"}},
    34,  {{"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    43,  {{"1"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    44,  {{"1"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    35,  {{"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    53,  {{"1"}, {"2"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    45,  {{"1"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    54,  {{"1"}, {"2"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    55,  {{"1"}, {"2"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    
    (* === Special 331 connection === *)
    331, {{"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    (* Default *)
    _, {{ToString[param]}}
];

(* ===== Interior (Link/Connection) Curve Sequences ===== *)
(* \[FivePointedStar] \:d558\:b4dc\:cf54\:b529: \:b178\:b4dc \:c0ac\:c774 \:c5f0\:acb0 \:cee4\:be0c \[FivePointedStar] *)

getInteriorCurves[param_Integer] := Switch[param,
    (* === Basic connections === *)
    11,  {{"1"}},
    22,  {{"1"}, {"3"}, {"1"}},
    23,  {{"1"}, {"3"}, {"2"}, {"1"}},
    32,  {{"1"}, {"2"}, {"3"}, {"1"}},
    33,  {{"1"}, {"2"}, {"3"}, {"2"}, {"1"}},
    24,  {{"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    42,  {{"1"}, {"2"}, {"2"}, {"3"}, {"1"}},
    34,  {{"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    43,  {{"1"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    44,  {{"1"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    35,  {{"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    53,  {{"1"}, {"2"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    45,  {{"1"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    54,  {{"1"}, {"2"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"1"}},
    55,  {{"1"}, {"2"}, {"2"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"2"}, {"2"}, {"1"}},
    
    (* === Special 331 connection === *)
    331, {{"1"}, {"3"}, {"1"}, {"5"}, {"1"}, {"3"}, {"1"}},
    
    (* Default *)
    _, {{ToString[param]}}
];

(* ================================================================== *)
(* ===== RENDERING FUNCTIONS ===== *)
(* ================================================================== *)

(* Render a single position with optional top curve and external *)
(* pos: {base} or {base, top} *)
(* extBottom, extTop: external values or None *)
(* ✅ CORRECT ORDER: base → top → extBottom → extTop (externals OUTERMOST!) *)
renderPosition[pos_List, extBottom_, extTop_] := Module[{base, top, result},
    base = pos[[1]];
    top = If[Length[pos] >= 2, pos[[2]], None];
    
    (* Build: base → top → extBottom → extTop *)
    result = base;
    
    (* Step 1: Add top curve first if exists *)
    If[top =!= None,
        result = "\\overset{" <> top <> "}{" <> result <> "}"
    ];
    
    (* Step 2: Add extBottom (above top!) *)
    If[extBottom =!= None,
        result = "\\overset{\\textcolor{red}{" <> ToString[extBottom] <> "}}{" <> result <> "}"
    ];
    
    (* Step 3: Add extTop (topmost!) *)
    If[extTop =!= None,
        result = "\\overset{\\textcolor{red}{" <> ToString[extTop] <> "}}{" <> result <> "}"
    ];
    
    result
];

(* Apply externals to curve sequence *)
(* curves: {{base}, {base, top}, ...} *)
(* externals: Association[port_idx -> external_value] *)
(* direction: "left" or "right" - right\:c774\:ba74 position \:c21c\:c11c \:b4a4\:c9d1\:c74c *)
applyCurvesWithExternals[curves_List, externals_Association, direction_String:"left"] := Module[
    {result, portIdx, i, pos, extBottom, extTop, 
     portRanges, renderedPositions, renderOrder},
    
    (* 1. \:ba3c\:c800 \:ac01 position\:c758 \:d3ec\:d2b8 \:bc94\:c704 \:acc4\:c0b0 *)
    portRanges = {};
    portIdx = 0;
    For[i = 1, i <= Length[curves], i++,
        pos = curves[[i]];
        If[Length[pos] >= 2,
            AppendTo[portRanges, {portIdx, portIdx + 1}];  (* {bottom, top} *)
            portIdx += 2;
        ,
            AppendTo[portRanges, {portIdx}];  (* {bottom only} *)
            portIdx++;
        ];
    ];
    
    (* 2. \:ac01 position \:b80c\:b354\:b9c1 *)
    renderedPositions = {};
    For[i = 1, i <= Length[curves], i++,
        pos = curves[[i]];
        
        (* \:c774 position\:c758 \:d3ec\:d2b8 \:c778\:b371\:c2a4 *)
        extBottom = If[KeyExistsQ[externals, portRanges[[i, 1]]], 
                       externals[portRanges[[i, 1]]], None];
        extTop = None;
        If[Length[pos] >= 2 && Length[portRanges[[i]]] >= 2,
            extTop = If[KeyExistsQ[externals, portRanges[[i, 2]]], 
                        externals[portRanges[[i, 2]]], None];
        ];
        
        AppendTo[renderedPositions, renderPosition[pos, extBottom, extTop]];
    ];
    
    (* 3. \:bc29\:d5a5\:c5d0 \:b530\:b77c \:c21c\:c11c \:acb0\:c815 *)
    renderOrder = If[direction == "right", 
        Reverse[renderedPositions], 
        renderedPositions
    ];
    
    StringJoin[renderOrder]
];

(* Render interior connection with externals *)
(* External 없으면 축약 표현, 있으면 전개 (otimes 없이) *)
renderInterior[param_Integer, externals_Association] := Module[{curves, content, digits},
    
    (* External이 없으면 기존 축약 표현 사용 *)
    If[Length[externals] == 0,
        digits = IntegerDigits[param];
        If[Length[digits] == 3 && digits == {3, 3, 1},
            (* Special 331 case - 그냥 3,3으로 표시 *)
            Return[" \\overset{3,3}{\\bigcirc} "]
        ];
        If[Length[digits] >= 2,
            Return[" \\overset{" <> StringJoin[Riffle[ToString /@ digits, ","]] <> "}{\\otimes} "],
            Return[" \\overset{" <> ToString[param] <> "}{\\otimes} "]
        ]
    ];
    
    (* External이 있으면 전개 - otimes 없이 *)
    curves = getInteriorCurves[param];
    content = applyCurvesWithExternals[curves, externals, "left"];
    
    " " <> content <> " "
];

(* ================================================================== *)
(* ===== BUILD DECORATED NODE (✨ ENHANCED) ===== *)
(* Supports multiple instantons (up to 4) with configurable positions *)
(* ================================================================== *)

buildDecoratedNode[algebra_String, sidelinks_List, instantons_List, 
                   sideLinkExternals_Association, instantonExternals_Association,
                   nodeExts_List, position_String, sp_List, ip_List] := Module[
    {sortedS, result,
     leftParam, rightParam, topParam, bottomParam,
     leftExt, rightExt, topExt, bottomExt,
     leftCurves, rightCurves, topCurves, bottomCurves,
     nodeExtStr, instCount, instPositions, positionTable,
     instParams, instExts, pos, hasDecorations},
    
    (* Build node external string *)
    nodeExtStr = "";
    If[Length[nodeExts] > 0,
        nodeExtStr = StringJoin[
            Riffle[("\\textcolor{red}{" <> ToString[#] <> "}") & /@ nodeExts, ","]
        ];
    ];
    
    (* Check if there are any decorations *)
    hasDecorations = Length[sidelinks] > 0 || Length[instantons] > 0;
    
    (* If no decorations, return algebra (with nodeExt if exists) *)
    If[!hasDecorations,
        If[nodeExtStr != "",
            Return["\\underset{" <> nodeExtStr <> "}{" <> algebra <> "}"],
            Return[algebra]
        ]
    ];
    
    (* Sort sidelinks by param descending - sidelinks is now {{param, slIdx}, ...} *)
    sortedS = Reverse[SortBy[sidelinks, First]];
    
    (* Initialize all positions *)
    leftParam = 0; rightParam = 0; topParam = 0; bottomParam = 0;
    leftSlIdx = -1; rightSlIdx = -1; topSlIdx = -1; bottomSlIdx = -1;
    leftExt = <||>; rightExt = <||>; topExt = <||>; bottomExt = <||>;
    
    (* ✨ Select instanton position table based on node position *)
    positionTable = Which[
        position == "single", instantonPositionsSingleNode,
        position == "left", instantonPositionsLeftEnd,
        position == "right", instantonPositionsRightEnd,
        position == "middle", instantonPositionsMiddle,
        True, instantonPositionsSingleNode
    ];
    
    (* ✨ Get instanton positions for this count *)
    instCount = Length[instantons];
    instPositions = Lookup[positionTable, instCount, {}];
    
    (* ✨ Build instanton params and externals lists *)
    instParams = Table[ip[[instantons[[i]] + 1]], {i, instCount}];
    instExts = Table[Lookup[instantonExternals, instantons[[i]], <||>], {i, instCount}];
    
    (* ✨ Assign instantons to positions (negative param marks as instanton) *)
    Do[
        pos = instPositions[[i]];
        Which[
            pos == "left" && leftParam == 0,
                leftParam = -instParams[[i]]; leftExt = instExts[[i]],
            pos == "right" && rightParam == 0,
                rightParam = -instParams[[i]]; rightExt = instExts[[i]],
            pos == "top" && topParam == 0,
                topParam = -instParams[[i]]; topExt = instExts[[i]],
            pos == "bottom" && bottomParam == 0,
                bottomParam = -instParams[[i]]; bottomExt = instExts[[i]]
        ],
        {i, instCount}
    ];
    
    (* Assign sidelinks to remaining positions based on node position *)
    (* FIXED v3.2: Use index-based assignment to handle instanton-occupied slots *)
    (* FIXED v3.2: Single node now uses all 4 positions (left -> top -> right -> bottom) *)
    (* FIXED v3.3: Use slIdx for external lookup instead of param *)
    Module[{idx = 1, slParam, slIdx},
        Which[
            position == "right",
                (* Right node: right -> top -> bottom *)
                If[idx <= Length[sortedS] && rightParam == 0, 
                    slParam = sortedS[[idx, 1]]; slIdx = sortedS[[idx, 2]];
                    rightParam = slParam; rightExt = Lookup[sideLinkExternals, slIdx, <||>]; idx++];
                If[idx <= Length[sortedS] && topParam == 0,
                    slParam = sortedS[[idx, 1]]; slIdx = sortedS[[idx, 2]];
                    topParam = slParam; topExt = Lookup[sideLinkExternals, slIdx, <||>]; idx++];
                If[idx <= Length[sortedS] && bottomParam == 0,
                    slParam = sortedS[[idx, 1]]; slIdx = sortedS[[idx, 2]];
                    bottomParam = slParam; bottomExt = Lookup[sideLinkExternals, slIdx, <||>]; idx++],
            position == "middle",
                (* Middle node: top -> bottom only (left/right occupied by interior links) *)
                If[idx <= Length[sortedS] && topParam == 0,
                    slParam = sortedS[[idx, 1]]; slIdx = sortedS[[idx, 2]];
                    topParam = slParam; topExt = Lookup[sideLinkExternals, slIdx, <||>]; idx++];
                If[idx <= Length[sortedS] && bottomParam == 0,
                    slParam = sortedS[[idx, 1]]; slIdx = sortedS[[idx, 2]];
                    bottomParam = slParam; bottomExt = Lookup[sideLinkExternals, slIdx, <||>]; idx++],
            position == "single",
                (* Single node: all 4 positions available - left -> top -> right -> bottom *)
                If[idx <= Length[sortedS] && leftParam == 0,
                    slParam = sortedS[[idx, 1]]; slIdx = sortedS[[idx, 2]];
                    leftParam = slParam; leftExt = Lookup[sideLinkExternals, slIdx, <||>]; idx++];
                If[idx <= Length[sortedS] && topParam == 0,
                    slParam = sortedS[[idx, 1]]; slIdx = sortedS[[idx, 2]];
                    topParam = slParam; topExt = Lookup[sideLinkExternals, slIdx, <||>]; idx++];
                If[idx <= Length[sortedS] && rightParam == 0,
                    slParam = sortedS[[idx, 1]]; slIdx = sortedS[[idx, 2]];
                    rightParam = slParam; rightExt = Lookup[sideLinkExternals, slIdx, <||>]; idx++];
                If[idx <= Length[sortedS] && bottomParam == 0,
                    slParam = sortedS[[idx, 1]]; slIdx = sortedS[[idx, 2]];
                    bottomParam = slParam; bottomExt = Lookup[sideLinkExternals, slIdx, <||>]; idx++],
            True,  (* left end of chain *)
                (* Left node: left -> top -> bottom *)
                If[idx <= Length[sortedS] && leftParam == 0,
                    slParam = sortedS[[idx, 1]]; slIdx = sortedS[[idx, 2]];
                    leftParam = slParam; leftExt = Lookup[sideLinkExternals, slIdx, <||>]; idx++];
                If[idx <= Length[sortedS] && topParam == 0,
                    slParam = sortedS[[idx, 1]]; slIdx = sortedS[[idx, 2]];
                    topParam = slParam; topExt = Lookup[sideLinkExternals, slIdx, <||>]; idx++];
                If[idx <= Length[sortedS] && bottomParam == 0,
                    slParam = sortedS[[idx, 1]]; slIdx = sortedS[[idx, 2]];
                    bottomParam = slParam; bottomExt = Lookup[sideLinkExternals, slIdx, <||>]; idx++]
        ]
    ];
    
    (* Helper to get curves - negative param means instanton *)
    getCurvesForParam[p_] := If[p < 0, getSideLinkCurves[-p], getSideLinkCurves[p]];
    
    (* Helper for instanton abbreviation (when no external) *)
    getInstantonAbbrev[p_Integer] := Switch[p,
        1, "1",
        882, "I^{\\oplus 2}",
        883, "I^{\\oplus 3}",
        884, "I^{\\oplus 4}",
        885, "I^{\\oplus 5}",
        886, "I^{\\oplus 6}",
        887, "I^{\\oplus 7}",
        8881, "I^{\\oplus 8}",
        889, "I^{\\oplus 9}",
        8810, "I^{\\oplus 10}",
        8811, "I^{\\oplus 11}",
        _, ToString[p]  (* fallback *)
    ];
    
    (* Helper to render a position - uses abbreviation for instanton without external *)
    renderSideOrInst[param_, ext_, direction_] := Module[{curves},
        (* param < 0 means instanton *)
        If[param < 0 && Length[ext] == 0,
            (* Instanton without external - use abbreviation *)
            getInstantonAbbrev[-param],
            (* Has external or is sidelink - use curve sequence *)
            curves = getCurvesForParam[param];
            applyCurvesWithExternals[curves, ext, direction]
        ]
    ];
    
    (* Build result *)
    result = algebra;
    
    (* ✨ Add node external as underset - always at bottom under algebra *)
    (* This goes directly under the node, before any sidelinks *)
    If[nodeExtStr != "",
        result = "\\underset{" <> nodeExtStr <> "}{" <> result <> "}"
    ];
    
    (* Add top *)
    If[topParam != 0,
        result = "\\overset{" <> renderSideOrInst[topParam, topExt, "left"] <> "}{" <> result <> "}"
    ];
    
    (* Add bottom *)
    If[bottomParam != 0,
        result = "\\underset{" <> renderSideOrInst[bottomParam, bottomExt, "left"] <> "}{" <> result <> "}"
    ];
    
    (* Add left *)
    If[leftParam != 0,
        result = renderSideOrInst[leftParam, leftExt, "left"] <> " " <> result
    ];
    
    (* Add right - reversed direction *)
    If[rightParam != 0,
        result = result <> " " <> renderSideOrInst[rightParam, rightExt, "right"]
    ];
    
    result
];

(* ================================================================== *)
(* ===== REPEATING PATTERN COMPRESSION (NEW v3.2) ===== *)
(* Finds (link, node) repeats and compresses to (otimes{link} node)^n *)
(* Example: so - 11 - so - 11 - so -> so (otimes{1,1} so)^2 *)
(* ================================================================== *)

(* Helper: render interior link in abbreviated form for compression *)
renderInteriorAbbrev[param_Integer] := Module[{digits},
    digits = IntegerDigits[param];
    If[Length[digits] == 3 && digits == {3, 3, 1},
        (* Special 331 case *)
        " \\overset{3,3}{\\bigcirc} ",
        (* Normal case *)
        If[Length[digits] >= 2,
            " \\overset{" <> StringJoin[Riffle[ToString /@ digits, ","]] <> "}{\\otimes} ",
            " \\overset{" <> ToString[param] <> "}{\\otimes} "
        ]
    ]
];

(* Compression function with external support *)
compressChain[nodes_List, links_List, 
              sidelinksPerNode_List, instantonsPerNode_List,
              sideLinkExternalsPerNode_List, instantonExternalsPerNode_List,
              interiorExternals_List, nodeExternals_List,
              sp_List, ip_List] := Module[
    {latex, i, j, count, nodeAlgebra, position, decorated, linkParam, nodeParam, linkExts},
    
    If[Length[nodes] == 0, Return["\\text{empty}"]];
    
    (* Single node case - use "single" position for full 4-position support *)
    If[Length[nodes] == 1,
        nodeAlgebra = gaugeAlgebra[nodes[[1]]];
        decorated = buildDecoratedNode[
            nodeAlgebra, 
            sidelinksPerNode[[1]], 
            instantonsPerNode[[1]],
            sideLinkExternalsPerNode[[1]],
            instantonExternalsPerNode[[1]],
            nodeExternals[[1]],
            "single", sp, ip
        ];
        Return["{" <> decorated <> "}"]
    ];
    
    latex = "{";
    
    (* First node output *)
    nodeAlgebra = gaugeAlgebra[nodes[[1]]];
    decorated = buildDecoratedNode[
        nodeAlgebra, 
        sidelinksPerNode[[1]], 
        instantonsPerNode[[1]],
        sideLinkExternalsPerNode[[1]],
        instantonExternalsPerNode[[1]],
        nodeExternals[[1]],
        "left", sp, ip
    ];
    latex = latex <> decorated;
    
    (* Process from second node as (link, node) pairs *)
    i = 2;
    While[i <= Length[nodes],
        linkParam = links[[i - 1]];
        nodeParam = nodes[[i]];
        linkExts = interiorExternals[[i - 1]];
        
        (* Check if compression is possible:
           - No decoration on current node
           - No external on interior link *)
        If[Length[sidelinksPerNode[[i]]] > 0 || 
           Length[instantonsPerNode[[i]]] > 0 ||
           Length[nodeExternals[[i]]] > 0 ||
           Length[linkExts] > 0,
            (* Has decoration or external - output directly *)
            position = If[i == Length[nodes], "right", "middle"];
            latex = latex <> renderInterior[linkParam, linkExts];
            nodeAlgebra = gaugeAlgebra[nodeParam];
            decorated = buildDecoratedNode[
                nodeAlgebra, 
                sidelinksPerNode[[i]], 
                instantonsPerNode[[i]],
                sideLinkExternalsPerNode[[i]],
                instantonExternalsPerNode[[i]],
                nodeExternals[[i]],
                position, sp, ip
            ];
            latex = latex <> decorated;
            i++,
            
            (* No decoration - check for repeating pattern *)
            count = 1;
            j = i + 1;
            
            (* Count repeating (link, node) pairs with no decorations *)
            While[j <= Length[nodes] &&
                  links[[j - 1]] == linkParam &&
                  nodes[[j]] == nodeParam &&
                  Length[sidelinksPerNode[[j]]] == 0 &&
                  Length[instantonsPerNode[[j]]] == 0 &&
                  Length[nodeExternals[[j]]] == 0 &&
                  Length[interiorExternals[[j - 1]]] == 0,
                count++;
                j++
            ];
            
            If[count >= 2,
                (* Repeating pattern found - output compressed form *)
                nodeAlgebra = gaugeAlgebra[nodeParam];
                latex = latex <> "(" <> renderInteriorAbbrev[linkParam] <> nodeAlgebra <> ")^{" <> ToString[count] <> "}";
                i = j,
                
                (* No repeat - output directly *)
                position = If[i == Length[nodes], "right", "middle"];
                latex = latex <> renderInterior[linkParam, linkExts];
                nodeAlgebra = gaugeAlgebra[nodeParam];
                decorated = buildDecoratedNode[
                    nodeAlgebra, 
                    sidelinksPerNode[[i]], 
                    instantonsPerNode[[i]],
                    sideLinkExternalsPerNode[[i]],
                    instantonExternalsPerNode[[i]],
                    nodeExternals[[i]],
                    position, sp, ip
                ];
                latex = latex <> decorated;
                i++
            ]
        ]
    ];
    
    latex = latex <> "}";
    latex
];

(* ================================================================== *)
(* ===== CONVERT TOPOLOGY TO LATEX ===== *)
(* ================================================================== *)

(* Returns {latexString, sugraInfo} where sugraInfo is Association with anom, sig, det *)
topologyToLatexWithInfo[line_String] := Module[
    {parsed, kinds, bparams, sConn, iConn, sp, ip, eConn, ep, sugraInfo,
     nodes, links, nodeIndices, linkIndices,
     sidelinksPerNode, instantonsPerNode,
     sideLinkExternalsPerNode, instantonExternalsPerNode, interiorExternals, nodeExternals,
     latex, i, j, k, nodeAlgebra, decorated, linkExts},
    
    parsed = parseLine[line];
    If[parsed === $Failed, Return[{"% Parse error", <||>}]];
    
    {kinds, bparams, sConn, iConn, sp, ip, eConn, ep, sugraInfo} = parsed;
    
    (* Separate nodes and links *)
    nodes = {};
    links = {};
    nodeIndices = {};
    linkIndices = {};
    
    For[i = 1, i <= Length[kinds], i++,
        Which[
            kinds[[i]] == 0, 
                (AppendTo[nodeIndices, i - 1 -> Length[nodes]];
                 AppendTo[nodes, bparams[[i]]]),
            kinds[[i]] == 1, 
                (AppendTo[linkIndices, i - 1 -> Length[links]];
                 AppendTo[links, bparams[[i]]])
        ]
    ];
    nodeIndices = Association[nodeIndices];
    linkIndices = Association[linkIndices];
    
    (* Build sidelinks per node: nodeIdx -> {{param1, slIdx1}, {param2, slIdx2}, ...} *)
    sidelinksPerNode = Table[{}, Length[nodes]];
    For[j = 1, j <= Length[sConn], j++,
        With[{blockIdx = sConn[[j, 1]], slIdx = sConn[[j, 2]]},
            If[KeyExistsQ[nodeIndices, blockIdx] && slIdx + 1 <= Length[sp],
                With[{nodeIdx = nodeIndices[blockIdx] + 1},
                    AppendTo[sidelinksPerNode[[nodeIdx]], {sp[[slIdx + 1]], slIdx}];
                ]
            ]
        ]
    ];
    
    (* Build instantons per node *)
    instantonsPerNode = Table[{}, Length[nodes]];
    For[j = 1, j <= Length[iConn], j++,
        With[{blockIdx = iConn[[j, 1]], instIdx = iConn[[j, 2]]},
            If[KeyExistsQ[nodeIndices, blockIdx] && instIdx + 1 <= Length[ip],
                With[{nodeIdx = nodeIndices[blockIdx] + 1},
                    AppendTo[instantonsPerNode[[nodeIdx]], instIdx];
                ]
            ]
        ]
    ];
    
    (* Build external connections *)
    (* sideLinkExternalsPerNode: nodeIdx -> <|slIdx -> <|port -> extVal|>|> *)
    (* instantonExternalsPerNode: nodeIdx -> <|instIdx -> <|port -> extVal|>|> *)
    (* interiorExternals: linkIdx -> <|port -> extVal|> *)
    (* nodeExternals: nodeIdx -> {extVal1, extVal2, ...} - Block에 직접 붙는 External *)
    sideLinkExternalsPerNode = Table[<||>, Length[nodes]];
    instantonExternalsPerNode = Table[<||>, Length[nodes]];
    interiorExternals = Table[<||>, Length[links]];
    nodeExternals = Table[{}, Length[nodes]];
    
    For[j = 1, j <= Length[eConn], j++,
        With[{parentId = eConn[[j, 1]], parentType = eConn[[j, 2]], 
              portIdx = eConn[[j, 3]], extId = eConn[[j, 4]]},
            With[{extVal = If[extId + 1 <= Length[ep], ep[[extId + 1]], 0]},
                Which[
                    parentType == 0,  (* Block/Interior *)
                        If[KeyExistsQ[linkIndices, parentId],
                            (* It's an interior link *)
                            With[{linkIdx = linkIndices[parentId] + 1},
                                interiorExternals[[linkIdx]] = 
                                    Append[interiorExternals[[linkIdx]], portIdx -> extVal]
                            ],
                            (* It's a Block (node) - External 직접 붙음 *)
                            If[KeyExistsQ[nodeIndices, parentId],
                                With[{nodeIdx = nodeIndices[parentId] + 1},
                                    AppendTo[nodeExternals[[nodeIdx]], extVal]
                                ]
                            ]
                        ],
                    parentType == 1,  (* SideLink *)
                        (* Find which node this sidelink is connected to *)
                        Do[
                            If[sConn[[k, 2]] == parentId,
                                With[{blockIdx = sConn[[k, 1]]},
                                    If[KeyExistsQ[nodeIndices, blockIdx],
                                        With[{nodeIdx = nodeIndices[blockIdx] + 1},
                                            (* Store by slIdx (parentId), not param *)
                                            With[{currentExt = Lookup[sideLinkExternalsPerNode[[nodeIdx]], parentId, <||>]},
                                                sideLinkExternalsPerNode[[nodeIdx]] = 
                                                    Append[
                                                        KeyDrop[sideLinkExternalsPerNode[[nodeIdx]], parentId],
                                                        parentId -> Append[currentExt, portIdx -> extVal]
                                                    ]
                                            ]
                                        ]
                                    ]
                                ];
                                Break[]  (* Only process first match *)
                            ],
                            {k, Length[sConn]}
                        ],
                    parentType == 2,  (* Instanton *)
                        Do[
                            If[iConn[[k, 2]] == parentId && parentId + 1 <= Length[ip],
                                With[{blockIdx = iConn[[k, 1]]},
                                    If[KeyExistsQ[nodeIndices, blockIdx],
                                        With[{nodeIdx = nodeIndices[blockIdx] + 1},
                                            With[{currentExt = Lookup[instantonExternalsPerNode[[nodeIdx]], parentId, <||>]},
                                                instantonExternalsPerNode[[nodeIdx]] = 
                                                    Append[
                                                        KeyDrop[instantonExternalsPerNode[[nodeIdx]], parentId],
                                                        parentId -> Append[currentExt, portIdx -> extVal]
                                                    ]
                                            ]
                                        ]
                                    ]
                                ]
                            ],
                            {k, Length[iConn]}
                        ]
                ]
            ]
        ]
    ];
    
    (* Build LaTeX with compression *)
    latex = compressChain[nodes, links, 
                  sidelinksPerNode, instantonsPerNode,
                  sideLinkExternalsPerNode, instantonExternalsPerNode,
                  interiorExternals, nodeExternals,
                  sp, ip];
    
    {latex, sugraInfo}
];

(* Wrapper for backward compatibility *)
topologyToLatex[line_String] := topologyToLatexWithInfo[line][[1]];

(* ================================================================== *)
(* ===== PROCESS FILE ===== *)
(* ================================================================== *)

extractNodeSignature[line_String] := Module[{parsed, kinds, bparams, nodes},
    parsed = parseLine[line];
    If[parsed === $Failed, Return["error"]];
    {kinds, bparams} = parsed[[1;;2]];
    nodes = {};
    Do[If[kinds[[i]] == 0, AppendTo[nodes, bparams[[i]]]], {i, Length[kinds]}];
    Sort[nodes]
];

processFile[inputFile_String, outputFile_String, includeSugraInfo_:True] := Module[
    {lines, latexLines, sugraInfos, uniqueLatex, output, i, j, validLines, 
     signatures, grouped, sortedGroups, latexWithInfo, infoStr},
    
    If[!FileExistsQ[inputFile],
        Print["Error: Input file not found: ", inputFile];
        Return[$Failed]
    ];
    
    lines = ReadList[inputFile, String];
    validLines = Select[lines, StringLength[#] > 0 &];
    
    Print["Processing ", Length[validLines], " topologies..."];
    
    latexLines = {};
    sugraInfos = {};
    signatures = {};
    For[i = 1, i <= Length[validLines], i++,
        If[Mod[i, 100] == 0, Print["Processed ", i, " / ", Length[validLines]]];
        latexWithInfo = topologyToLatexWithInfo[validLines[[i]]];
        AppendTo[latexLines, latexWithInfo[[1]]];
        AppendTo[sugraInfos, latexWithInfo[[2]]];
        AppendTo[signatures, extractNodeSignature[validLines[[i]]]];
    ];
    
    uniqueLatex = DeleteDuplicates[latexLines];
    Print["After removing exact duplicates: ", Length[uniqueLatex]];
    
    grouped = GatherBy[
        Transpose[{latexLines, signatures, sugraInfos}], 
        #[[2]] &  (* Group by signature *)
    ];
    
    sortedGroups = Reverse[SortBy[grouped, Length]];
    
    output = "% Generated LaTeX file with External curves (in RED)\n";
    output = output <> "% Use with amsmath and xcolor packages\n";
    output = output <> "% \\usepackage{xcolor}\n";
    output = output <> "% Total unique: " <> ToString[Length[uniqueLatex]] <> "\n";
    output = output <> "% Total groups: " <> ToString[Length[sortedGroups]] <> "\n";
    If[includeSugraInfo, 
        output = output <> "% SUGRA info: Delta (anomaly residual), n_+/n_-/n_0 (signature), det\n"
    ];
    output = output <> "\n";
    
    For[i = 1, i <= Length[sortedGroups], i++,
        With[{group = sortedGroups[[i]]},
            With[{sig = group[[1, 2]], 
                  items = DeleteDuplicatesBy[group, First]},  (* Keep first occurrence with its sugraInfo *)
                output = output <> "% === Group " <> ToString[i] <> ": nodes = " <> ToString[sig] <> " (" <> ToString[Length[items]] <> " items) ===\n";
                output = output <> "\\begin{align*}\n";
                For[j = 1, j <= Length[items], j++,
                    With[{latexItem = items[[j, 1]], sugraItem = items[[j, 3]]},
                        (* Build SUGRA info string if available and requested *)
                        infoStr = "";
                        If[includeSugraInfo && sugraItem["anom"] =!= None,
                            infoStr = " \\quad (\\Delta=" <> ToString[sugraItem["anom"]];
                            If[sugraItem["sig"] =!= None,
                                infoStr = infoStr <> ", n_{+}=" <> ToString[sugraItem["sig"][[1]]] <>
                                    ", n_{-}=" <> ToString[sugraItem["sig"][[2]]] <>
                                    ", n_{0}=" <> ToString[sugraItem["sig"][[3]]];
                            ];
                            If[sugraItem["det"] =!= None,
                                infoStr = infoStr <> ", |\\det|=" <> ToString[Abs[sugraItem["det"]]];
                            ];
                            infoStr = infoStr <> ")";
                        ];
                        
                        Which[
                            j == 1,
                                output = output <> "&" <> latexItem <> infoStr,
                            True,
                                output = output <> " \\\\\n&" <> latexItem <> infoStr
                        ]
                    ]
                ];
                output = output <> "\n\\end{align*}\n\n";
            ]
        ]
    ];
    
    WriteString[outputFile, output];
    Close[outputFile];
    
    Print["Done! Output saved to: ", outputFile];
    Print["Total groups: ", Length[sortedGroups]];
    If[includeSugraInfo, Print["SUGRA info: Delta, n_+, n_-, n_0, det"]];
];

(* ===== Usage ===== *)

Print["========================================"];
Print["  Topology to LaTeX Converter v3.3"];
Print["  (WITH EXTERNAL CURVES - RED)"];
Print["  (WITH SUGRA INFO SUPPORT)"];
Print["========================================"];
Print[];
Print["NEW FEATURES in v3.3:"];
Print["  - SUGRA info (anom, sig, det) parsing"];
Print["  - Optional SUGRA info in LaTeX output"];
Print[];
Print["USAGE:"];
Print["  processFile[\"input.txt\", \"output.tex\"]"];
Print["  processFile[\"input.txt\", \"output.tex\", True]   (* with SUGRA info *)"];
Print["  processFile[\"input.txt\", \"output.tex\", False]  (* without SUGRA info *)"];
Print[];
Print["INPUT FORMAT:"];
Print["  <topology> | anom=-50 sig=(1,8,0) det=1"];
Print[];
Print["OUTPUT FORMAT:"];
Print["  {topology} \\quad (\\Delta=-50, n_{+}=1, n_{-}=8, n_{0}=0, \\det=1)"];
Print[];
Print["LaTeX PREAMBLE REQUIRED:"];
Print["  \\usepackage{amsmath}"];
Print["  \\usepackage{amssymb}"];
Print["  \\usepackage{xcolor}"];
Print["  \\allowdisplaybreaks"];
Print[];
Print["========================================"];
