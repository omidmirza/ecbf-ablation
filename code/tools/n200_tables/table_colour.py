#!/usr/bin/env python3
"""Relative colour shading of the generated N200 edition tables, behind a switch.

Every fragment named in COLOUR_TABLES is written with BOTH variants:

    \\ifdefined\\vncolortab\\else\\def\\vncolortab#1{0}\\fi% COLOUR-SWITCH ...
    \\ifnum\\vncolortab{<name>}=1 % COLOUR-SWITCH coloured variant
    <coloured table, then \\vncolourlegend{...} (printed once per document)>
    \\else% COLOUR-SWITCH plain variant
    <the unchanged plain table, byte for byte>
    \\fi% COLOUR-SWITCH end

\\vncolortab{<name>} (release_helper.tex, report/report_helper.tex) is 1 when \\colortablestrue or
\\colortable{<name>,...} was given and 0 otherwise; by default the plain variant is typeset, so the editions are
unchanged. An entry file that defines \\releasecolortables before the helper is read prints the coloured variants
(editions/*_colored.tex). plain_variant() returns the plain text (check_tables.py reads it).

Scheme (one rule in every coloured table):
  * Blocks: a block is what a reader compares in one table block. In the Section VI tables (SCENARIO_TABLES:
    sec6_short, sec6_combined) it is the scenario block: all its setting rows (default setting, S1-S5, S*, at every
    alpha_0) and its baseline box. In the tables with the four variants at one setting it is that setting block;
    in the sweep tables (phi/stress, delta, K, b, checks) it is the block of rows shown together (one setting, one
    L_mu). A baseline box that serves several such blocks is compared with the best reference among the blocks it
    serves at its alpha_0. Nothing outside the table is used.
  * Reference (per block and column with an arrow): the best mean among the CERTIFIED rows of the block that keep
    their risk guarantee (event counts within their levels). Baselines are never the reference (ACI and QT miss
    their levels; c_max is a fixed bound); c=0 and Unfiltered are never coloured. Across S1-S5 every setting is
    green where it is the best of the block and coloured by its gap elsewhere (the tunable trade-off); the
    baselines are coloured against the same reference.
  * Means: soft green when the gap to the reference is within the materiality floor (0.01 m for the clearance
    columns h_min, psi_min, their printing resolution; 1 % of the reference; 1 pp; 0.02 for /c_max) or the value is
    better than the reference; beyond the floor one continuous scale, logarithmic in the gap, light green (2.5 %)
    -> pale yellow (10 %) -> light orange (30 %) -> soft red (>= 100 %); %-columns in percentage points 1, 4, 12,
    >= 40 pp; violation columns (VIOLATION: q_CBF, q_Y, q_unsafe, q_nocert, q_nodep) on a stronger scale,
    1, 2, 5, >= 10 pp.
  * Counts with a level (a column with an arrow): green when the 95 % upper bound meets the level (certified:
    checkmark; baselines: bound <= 5 %), red when it misses it, deeper the further the bound is above the level
    (logarithmic, deepest at 100 %, i.e. 200 of 200); a count of zero events whose bound misses the level only
    because of the sample size (0 of 200 gives 1.5 %, e.g. against delta_1 < 1.5 %) is neutral (no colour).
    Columns without an arrow (m_cert, the excess counts q_A2, q_pw) stay plain, and so do m_tau^use and m_p^max
    (diagnostics) and the counts without a bound (q_CBF and q_unsafe of the checks tables).
  * Collision: the NAME of a method with contact in at least one trajectory (COLLISIONS, set by build_tables.py
    from outcomes.contact of the exports) is purple with a cross; nothing else marks a row.
  * Stress tables (phi_*, stress_*): plain in the coloured editions too; their phi < 1 rows carry no guarantee and
    are never a reference (collect() still reads their counts).
  * Layout: every page block of a coloured table is set in a box of the height of the same block of the plain
    table (plain_height), so the coloured editions keep the page breaks and the page count of the plain ones.
  * The colour key is printed once per edition: inside the key box of Table I (sec6_short, sec6_combined) in the
    paper editions, below Table 4 (default_corridor) in the report; the paper's table_key says the colours are as
    in Table I, the report's points to Table 4.
  * All values black; every coloured box is inset by a white gutter and the rules are drawn in the gutters
    (separator rows): the box of each row group and column and the outer rules, as in the plain variant. The
    values stacked in one box are not separated by a rule: their colours meet halfway between the lines; merged
    boxes are coloured as one unit."""
import math
import re

NMC = 200
# Compact booktabs tables use split mean/SD and event/bound columns rather than w{c} boxes.
COMPACT_TABLES = ("main_maze_arxiv_v8", "main_corridor_arxiv_v8",
                  "main_maze_extended_v8", "main_corridor_extended_v8")
# Fragments that carry a coloured variant (paper and report targets). Not coloured: table_key, design overviews,
# detailed_events_* (counts only), variation_* (temporal statistics), n30_record_* and checks_*_nodef (no edition).
COLOUR_TABLES = ("sec6_combined", "sec6_short", "sec6_maze_arxiv", "star_maze",
                 "default_corridor", "default_maze", "topfive_full_corridor", "topfive_full_maze",
                 "delta_corridor", "delta_maze", "period_corridor", "period_maze", "b_corridor", "b_maze",
                 "phi_corridor", "phi_maze",
                 "checks_corridor_L0005", "checks_corridor_L005", "checks_maze_L0005", "checks_maze_L005",
                 "checks_corridor_arxiv", "checks_maze_arxiv", "stress_corridor_arxiv", "stress_maze_arxiv",
                 "table_key", "table_key_arxiv", "secondary_corridor_arxiv", "secondary_maze_arxiv",
                 "topfive_full_corridor_arxiv_percent", "topfive_full_maze_arxiv_percent",
                 "topfive_full_corridor_onepage", "topfive_full_maze_onepage") + COMPACT_TABLES
KEY_TABLES = ("table_key", "table_key_arxiv")          # no table to shade: the coloured variant adds the colour key
STRESS_PREFIX = ("phi_", "stress_")  # rows with phi < 1 carry no guarantee: plain, never a reference

BASE_LEVEL = 5.0                     # level of the baselines' counts (bound <= 5 %)
MIN_BOUND = 1.5                      # one-sided 95 % bound of 0 of 200 (%)
REL_LO, REL_HI = 2.5, 100.0          # log ramp from 2.5 % (still green) to >= 100 % (soft red)
PP_LO, PP_HI = 1.0, 40.0             # %-columns: same ratios, in percentage points
TIE_REL, TIE_PP, TIE_M, TIE_RATIO = 1.0, 1.0, 0.01, 0.02   # TIE_M: clearance columns (m), printed to 0.01 m
VIOLATION = ('\\rmCBF', 'q_Y', '\\rmunsafe', '\\rmnocert', '\\rmnodep')   # %-columns on the stronger scale
VIOL_PP = (1.0, 1.4, 2.0, 5.0, 10.0)   # pp at the knots of RAMP (light green, yellow-green, yellow, apricot, red)
SCENARIO_TABLES = ("sec6_short", "sec6_combined")   # one block per scenario: all settings and the baseline box
FIXED = {'$c_{\\max}$'}             # coloured, never the reference
PLAIN = {'Unfiltered', '$c=0$'}      # never coloured, never the reference
COLLISIONS = set()                   # (environment, method name as printed): contact in >= 1 trajectory
NONMETHOD = {'$K$', '$L_\\mu$'}     # first-column heads of the checks tables (a row is a check, not a method)
G, T_IN, T_GRP = 0.5, 0.4, 0.6       # gutter, inner rule, group rule (pt)
H_IN = 2 * G + T_IN
INSET = 0.7                          # horizontal inset of a cell colour from the column edge (\\rgcell)
TIE = (0.784, 0.910, 0.769)          # soft green #c8e8c4: as good as the reference, or better


def LT(g):
    return math.log(g / REL_LO) / math.log(REL_HI / REL_LO)


# one continuous scale, starting at the tie green
RAMP = [(0.0, TIE),                          # soft green (gap just beyond the tie threshold, <= 2.5 %)
        (LT(5.0), (0.882, 0.941, 0.761)),   # yellow-green #e1f0c2 (5 %)
        (LT(10.0), (1.000, 0.945, 0.702)),  # pale yellow #fff1b3 (10 %)
        (LT(30.0), (0.996, 0.835, 0.663)),  # light apricot #fed5a9 (30 %)
        (1.0, (0.957, 0.651, 0.651))]       # soft red #f4a6a6 (>= 100 %)
BAD = (0.835, 0.776, 0.937)                 # soft purple #d5c6ef: collision (on the name only)
MISS_LO = (0.973, 0.792, 0.792)             # light red #f8caca: a count just above its level
MISS_HI = (0.890, 0.451, 0.451)             # deep red #e37373: a count bound of 100 % (200 of 200)

MARK_DEF = r"\ifdefined\vncolortab\else\def\vncolortab#1{0}\fi% COLOUR-SWITCH (tools/n200_tables/table_colour.py)"
MARK_IF = r"\ifnum\vncolortab{NAME}=1 % COLOUR-SWITCH coloured variant"
MARK_ELSE = r"\else% COLOUR-SWITCH plain variant"
MARK_END = r"\fi% COLOUR-SWITCH end"


def ramp(t):
    t = min(max(t, 0.0), 1.0)
    for (t0, c0), (t1, c1) in zip(RAMP, RAMP[1:]):
        if t <= t1:
            u = (t - t0) / (t1 - t0)
            return tuple(a + u * (b - a) for a, b in zip(c0, c1))
    return RAMP[-1][1]


def rgb(c):
    return ','.join(f'{x:.3f}' for x in c)


def cellcolor(bg):
    return f'\\rgcell{{{rgb(bg)}}}'


def split_cells(line):
    cells, depth, cur = [], 0, ''
    for ch in line:
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
        if ch == '&' and depth == 0:
            cells.append(cur)
            cur = ''
        else:
            cur += ch
    cells.append(cur)
    return cells


def braced(s, i):
    """s[i] == '{' -> (content, index after the closing brace)"""
    assert s[i] == '{', s[i:i + 20]
    d = 0
    for k in range(i, len(s)):
        if s[k] == '{':
            d += 1
        elif s[k] == '}':
            d -= 1
            if d == 0:
                return s[i + 1:k], k + 1
    raise ValueError(s[i:i + 40])


def preamble_of(line):
    i = line.index('\\begin{tabular}') + len('\\begin{tabular}')
    return braced(line, i)


def parse_preamble(p):
    """-> list of [left material, right material] per column, as TeX code"""
    cols, pend, i = [], [], 0

    def newcol():
        nonlocal pend
        if not cols:
            cols.append([''.join(pend), ''])
        else:
            cols[-1][1] += ''.join(pend)
            cols.append(['', ''])
        pend = []
    while i < len(p):
        ch = p[i]
        if ch == '|':
            pend.append('\\vrule width\\arrayrulewidth')
            i += 1
        elif ch in '!@':
            a, i = braced(p, i + 1)
            pend.append(a)
        elif ch in '><':
            _, i = braced(p, i + 1)
        elif ch in 'wWpmb':
            n = 2 if ch in 'wW' else 1
            i += 1
            for _ in range(n):
                _, i = braced(p, i)
            newcol()
        elif ch in 'lcr':
            i += 1
            newcol()
        elif ch.isspace():
            i += 1
        else:
            raise ValueError(p[i:])
    cols[-1][1] += ''.join(pend)
    return cols


def hkey(cell):
    s = cell.strip()
    if re.match(r'^\\multicolumn\{1\}\{', s):
        _, k = braced(s, len('\\multicolumn{1}'))
        s, _ = braced(s, k)
    m = re.match(r'^\\multirow\{\d+\}\{\*\}', s)
    if m:
        s, _ = braced(s, m.end())
    return re.sub(r'\s+', '', s.replace('\\vnsd', ''))


NOT_COMPARED = ('m_{\\tau^', 'm_{p^{\\max}}')   # m_tau^use, m_p^max: diagnostics, plain like m_cert


def hinfo(key):
    d = -1 if '\\downarrow' in key else (1 if '\\uparrow' in key else 0)
    if any(k in key for k in NOT_COMPARED):
        d = 0
    return d, '(\\%)' in key, '\\rmunsafe' in key


MR = re.compile(r'^(\s*)\\multirow\{(\d+)\}\{\*\}\{(.*)\}(\s*)$', re.S)


def unwrap(cell):
    m = MR.match(cell)
    return (int(m.group(2)), m.group(3)) if m else (1, cell)


COUNT = re.compile(r'^\s*(\d+)(?:--(\d+))?(?:/\d+)?\s*\((\d+(?:\.\d+)?)\)')   # a, a--b or a/m, then (bound)
VAL = re.compile(r'(<?)(-?\d+(?:\.\d+)?)\}?\s*\\vpm\s*(\d+(?:\.\d+)?)')


def parse(cell):
    """('count', checkmark, bound, largest count) | ('mean', mean, sd) | None. A trajectory count without a bound
    (q_CBF, q_unsafe of the checks tables) has no level and is not compared: None (plain)."""
    if '\\na' in cell or not cell.strip() or cell.strip() == '\\vnsd':
        return None
    pct = re.match(r'^(\d+\.\d) \(\d+\.\d\)', cell)
    if pct:
        bound = re.match(r'^\d+\.\d \((\d+\.\d)\)', cell).group(1)
        return ('count', '\\checkmark' in cell, float(bound), float(pct.group(1)))
    m = COUNT.match(cell)
    if m:
        return ('count', '\\checkmark' in cell, float(m.group(3)), int(m.group(2) or m.group(1)))
    m = VAL.search(cell.replace('\\mathbf{', ''))
    if m:
        return ('mean', float(m.group(2)), float(m.group(3)))
    return None


def value_of(p):
    """(value, sd) of a mean, None otherwise"""
    return (p[1], p[2]) if p is not None and p[0] == 'mean' else None


def norm_name(s):
    s = re.sub(r'\s+', '', s)
    m = re.match(r'^\\textcolor\{[^}]*\}\{(.*)\}$', s)
    return m.group(1) if m else s


def floor_of(key, b):
    if '(\\%)' in key:
        return TIE_PP
    if '(m)' in key:
        return TIE_M
    if '/c_{\\max}' in key:
        return TIE_RATIO
    return TIE_REL / 100 * abs(b)


def violation(key):
    """A %-column of violation fractions (q_CBF, q_Y, q_unsafe, q_nocert, q_nodep): the stronger scale."""
    return '(\\%)' in key and any(k in key for k in VIOLATION)


def viol_t(gap):
    """Ramp position of a gap (pp) on the violation scale: log-linear between the knots VIOL_PP."""
    ts = [0.0] + [t for t, _ in RAMP[1:]]
    if gap <= VIOL_PP[0]:
        return 0.0
    for (g0, t0), (g1, t1) in zip(zip(VIOL_PP, ts), zip(VIOL_PP[1:], ts[1:])):
        if gap <= g1:
            return t0 + (t1 - t0) * math.log(gap / g0) / math.log(g1 / g0)
    return 1.0


def grade(v, b, d, pct, key=''):
    """None: within the materiality floor of the reference or better; otherwise the ramp position in [0, 1]"""
    gap = max((v - b) * (-d), 0.0)
    if gap <= floor_of(key, b) + 1e-12:
        return None
    if violation(key):
        return viol_t(gap)
    if pct:
        return max(math.log(max(gap, 1e-9) / PP_LO), 0.0) / math.log(PP_HI / PP_LO)
    if b == 0:
        return 1.0
    return max(math.log(max(100 * gap / abs(b), 1e-9) / REL_LO), 0.0) / math.log(REL_HI / REL_LO)


def bin_of(t):
    if t is None:
        return 'tie'
    if t < LT(5.0):
        return 'green'
    if t < LT(20.0):
        return 'yellow'
    if t < LT(60.0):
        return 'orange'
    return 'red'


def count_ok(p, certified):
    """A count meets its level: certified rows by their checkmark, baselines by a bound <= 5 %."""
    return p[1] if certified else p[2] <= BASE_LEVEL


TUPLE = re.compile(r'\(\s*([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*(\d+)\)')   # (alpha_0, alpha, delta, K)


def count_colour(p, certified, level=None):
    """Green when the count meets its level; red when it misses, deeper the further the bound is above the level
    (log of bound/level, deepest at a bound of 100 %). A certified count whose level is not printed next to it
    (e.g. K = infinity of the checks) is measured from the smallest bound, 0 of 200 (1.5 %). A count of zero events
    that misses its level only because of the sample size (0 of 200: bound 1.5 %) is neutral: (None, True)."""
    if count_ok(p, certified):
        return TIE, True
    if p[3] == 0:
        return None, True
    level = BASE_LEVEL if not certified else (level if level and level < p[2] else MIN_BOUND)
    t = math.log(max(p[2], level) / level) / math.log(100.0 / level)
    t = min(max(t, 0.0), 1.0)
    return tuple(a + t * (b - a) for a, b in zip(MISS_LO, MISS_HI)), False


def verdict_key(b, g, name):
    """Key of a row group's usability verdict, shared by all fragments: baselines by (environment, alpha_0,
    L_mu, name) (no setting: a baseline box); setting rows by (environment, alpha_0, L_mu, setting, name), e.g. a
    phi value of S1 (stress tables print its count, the phi tables do not)."""
    m = re.match(r'^(.*),\$\\alpha_0=([\d.]+)\$$', name)
    a0 = b['a0']
    if m:
        name, a0 = m.group(1), m.group(2)
    name = name.replace('$', '').replace('^\\dagger', '')
    return (b['env'], a0, b['L'] or '0.0005', b['gS'][g], name)


def collect(text, registry, name=''):
    """First pass over all fragments: guarantee verdicts of the row groups that print a count with a direction
    (registry[key] = True when a count misses its level), for the tables that do not print that count."""
    if name in KEY_TABLES or name in COMPACT_TABLES:
        return
    shade(text, registry=registry, collect_only=True, stress=name.startswith(STRESS_PREFIX))


def stress_name(name):
    """True for a phi value below 1 (row name '$0.25$' or '$0.25^\\dagger$')."""
    try:
        return float(name.replace('$', '').replace('^\\dagger', '').strip()) < 1
    except ValueError:
        return False


def nset(s):
    return re.sub(r'\s+', '', s) if s else s


def shade(text, stats=None, registry=None, collect_only=False, stress=False, scenario=False):
    """Coloured copy of a generated fragment (the plain text is not modified); collect_only: verdicts only;
    scenario: one block per scenario (SCENARIO_TABLES): every setting and alpha_0 of the pair shares a reference."""
    stats = stats if stats is not None else {}
    VERDICT = registry if registry is not None else {}
    lines = text.split('\n')
    tabs = []
    cur = None
    pair, env, a0, Lcur = 0, '?', '?', None
    for i, ln in enumerate(lines):
        s = ln.strip()
        m = re.search(r'\\textbf\{(Corridor|Maze)\}', s)
        if m:
            env = m.group(1)
        m = re.search(r'L_\\mu=([\d.]+)', s)
        if m:
            Lcur = m.group(1)
        if '\\begin{tabular}' in s and 'w{c}' in s:
            pre = preamble_of(ln)[0]
            base = '\\hskip\\arrayrulewidth' in pre
            if base:
                pair += 1
            cur = dict(start=i, pre=pre, base=base, pair=pair, headers=[], rows=[], blocks=[],
                       S=None, L=None, k=0, blk=None, pre0=[])
            tabs.append(cur)
            continue
        if cur is None:
            continue
        if '\\end{tabular}' in s:
            cur['end'] = i
            cur = None
            continue
        if s.startswith('\\multicolumn{') and not s.startswith('\\multicolumn{1}'):
            m = re.search(r'alpha_0\s*=\s*([\d.]+)', s) or re.search(r'=\s*\(([\d.]+),', s)
            if m:
                a0 = m.group(1)
            m = re.search(r'L_\\mu=([\d.]+)', s)
            if m:
                cur['L'] = m.group(1)
            m = re.search(r'\\textbf\{(S\d|S\$\^\\ast\$|default setting)\}', s)
            if m:
                cur['S'] = nset(m.group(1))
            m = TUPLE.search(s)
            cur['tuple'] = m.groups() if m else None
            cur['blk'] = None
            cur['k'] = 0
            continue
        if s.startswith('\\Xhline{1.0pt}'):
            cur['blk'] = None
            cur['k'] = 0
            continue
        if s.startswith('\\Xcline') or s.startswith('\\Xhline') or s.startswith('\\hline'):
            cur['k'] = 0
            continue
        if '&' not in s:
            continue
        cells = split_cells(s.rsplit('\\\\', 1)[0])
        keys = [hkey(c) for c in cells]
        if any(hinfo(k)[0] for k in keys):
            cur['headers'].append(keys)
            cur['blk'] = None
            cur['k'] = 0
            continue
        if not cur['headers']:
            cur['pre0'].append(keys[0])   # a head row without arrows (e.g. the L_mu head of the arXiv checks)
            continue
        if cur['blk'] is None:
            setcol = cur['headers'][0][0] == '' and len(cur['headers'][0]) > 1 and \
                cur['headers'][0][1].startswith('Method')
            cur['blk'] = dict(pair=cur['pair'], S=cur['S'], L=cur['L'] or Lcur, base=cur['base'], env=env, a0=a0,
                              tab=len(tabs) - 1, rows=[], setcol=setcol, tuple=cur.get('tuple'),
                              checks=any(k in NONMETHOD for k in [h[0] for h in cur['headers']] + cur['pre0']))
            cur['blocks'].append(cur['blk'])
        H = cur['headers']
        r = cur['k'] % len(H)
        hk = []
        for j in range(len(cells)):
            kk = ''
            for rr in range(r, -1, -1):
                if j < len(H[rr]) and H[rr][j]:
                    kk = H[rr][j]
                    break
            hk.append(kk)
        row = dict(line=i, r=r, cells=cells, keys=hk, tab=len(tabs) - 1, blk=cur['blk'])
        cur['rows'].append(row)
        cur['blk']['rows'].append(row)
        cur['k'] += 1

    blocks = [b for t in tabs for b in t['blocks']]
    # ---- groups (one method / one swept value), setting of each group, verdicts
    for b in blocks:
        certified = any('\\checkmark' in c for row in b['rows'] for c in row['cells'])
        b['cert'] = certified
        groups = []
        for row in b['rows']:
            if row['r'] == 0:
                groups.append([])
            row['g'] = len(groups) - 1
            groups[-1].append(row)
        b['groups'] = groups
        b['nog'] = []       # True: the row group misses its risk guarantee (a count misses its level)
        b['gS'] = []
        b['level'] = []
        b['ga0'] = []       # alpha_0 of the group (None: c_max, which serves every alpha_0)
        for g, rows in enumerate(groups):
            namecol = 1 if b['setcol'] else 0
            name = norm_name(unwrap(rows[0]['cells'][namecol])[1])
            b['gS'].append(nset(norm_name(rows[0]['cells'][0])) if b['setcol'] else b['S'])
            m = TUPLE.search(' '.join(rows[0]['cells'][:namecol + 1]))
            tup = m.groups() if m else b['tuple']
            b['level'].append(100 * float(tup[2]) if tup else None)
            m = re.match(r'^.*,\$\\alpha_0=([\d.]+)\$$', name)
            ga0 = tup[0] if tup else (m.group(1) if m else b['a0'])
            b['ga0'].append(None if name in FIXED or ga0 in (None, '?') else float(ga0))
            miss, has_count = False, False
            for row in rows:
                for j, c in enumerate(row['cells'][1:], 1):
                    p = parse(unwrap(c)[1])
                    if p is not None and p[0] == 'count' and hinfo(row['keys'][j])[0]:
                        has_count = True
                        if not count_ok(p, certified):
                            miss = True
            vkey = verdict_key(b, g, name)
            if has_count:
                VERDICT[vkey] = VERDICT.get(vkey, False) or miss
            b['nog'].append(miss or VERDICT.get(vkey, False))
            b.setdefault('names', []).append(name)
            b.setdefault('stress', []).append(stress and stress_name(name))
    if collect_only:
        return None

    def compat(b1, b2):
        if b1['pair'] != b2['pair'] or b1['env'] != b2['env']:
            return False
        if b1['L'] and b2['L'] and b1['L'] != b2['L']:
            return False
        return True

    def scompat(s1, s2):
        return s1 is None or s2 is None or s1 == s2

    def acompat(a1, a2):
        return a1 is None or a2 is None or abs(a1 - a2) < 1e-9

    cand = {}       # (block, setting, alpha_0, key) -> [(value, sd)]: the certified rows that keep their guarantee
    for bi, b in enumerate(blocks):
        if b['base']:   # baselines are never the reference
            continue
        for row in b['rows']:
            g = row['g']
            if b['nog'][g] or b['stress'][g] or b['names'][g] in FIXED | PLAIN:
                continue
            for j, c in enumerate(row['cells'][1:], 1):
                v = value_of(parse(unwrap(c)[1]))
                if v is not None and hinfo(row['keys'][j])[0]:
                    cand.setdefault((bi, b['gS'][g], b['ga0'][g], row['keys'][j]), []).append(v)

    def best(bi, gS, a0, key, d):
        """reference: best value among the certified rows that keep their guarantee in the same block (same baseline
        box, environment, L_mu; the same setting, or every setting of the scenario in the Section VI tables); a
        baseline row (gS None) gets the best reference of the blocks it serves"""
        pool = [x for (bj, sj, aj, kj), xs in cand.items()
                if kj == key and compat(blocks[bi], blocks[bj])
                and (scenario or (scompat(gS, sj) and acompat(a0, aj))) for x in xs]
        if not pool:
            return None
        return min(pool, key=lambda x: (x[0] * (-d), x[1]))

    # ---- paint plan
    out = list(lines)
    paint = {}      # (line, j) -> (bg, content)
    spans = {}      # (line, j) -> lines covered by a merged box anchored at the first line
    for bi, b in enumerate(blocks):
        for g, rows in enumerate(b['groups']):
            ncols = len(rows[0]['cells'])
            if b['stress'][g]:   # phi < 1: no guarantee, the row stays plain
                stats['stress'] = stats.get('stress', 0) + 1
                continue
            if (b['env'], b['names'][g]) in COLLISIONS:   # purple name with a cross: a collision
                jn = 1 if b['setcol'] else 0
                n, name = unwrap(rows[0]['cells'][jn])
                name = name.strip()
                m = re.match(r'^\\textcolor\{[^}]*\}\{(.*)\}$', name)
                if m:
                    name = m.group(1)
                # name and cross centred over the cell; they may reach into the column padding (1 pt from the
                # rule) and are scaled down only when even that is too narrow
                ws = re.findall(r'w\{c\}\{(\\vnw\{[^}]*\})\}', tabs[b['tab']]['pre'])
                W = f'\\dimexpr{ws[jn]}+2\\tabcolsep-2pt\\relax' if jn < len(ws) else '\\linewidth'
                content = ('\\sbox0{' + name + '\\,$\\scriptstyle\\times$}\\makebox[0pt][c]{'
                           f'\\ifdim\\wd0>{W}\\resizebox{{{W}}}{{!}}{{\\usebox0}}\\else\\usebox0\\fi}}')
                if n > 1:
                    spans[(rows[0]['line'], jn)] = [rw['line'] for rw in rows[:n]]
                paint[(rows[0]['line'], jn)] = (BAD, content)
                stats['collision'] = stats.get('collision', 0) + 1
            if b['names'][g] in PLAIN:   # c=0, Unfiltered: values never coloured
                continue
            level = b['level'][g]
            # merged 2-value count boxes: a column with counts and empty cells in the group
            for j in range(1, ncols):
                cs = [row['cells'][j] if j < len(row['cells']) else '' for row in rows]
                ps = [parse(unwrap(c)[1]) for c in cs]
                d = hinfo(rows[0]['keys'][j])[0] if j < len(rows[0]['keys']) else 0
                if len(rows) > 1 and d and any(p and p[0] == 'count' for p in ps) and \
                        any(not c.strip() or c.strip() == '\\vnsd' for c in cs) and not any(MR.match(c) for c in cs):
                    cols_ok = [count_colour(p, b['cert'], level) for p in ps if p and p[0] == 'count']
                    cols_ok = [(c_, o) for c_, o in cols_ok if c_ is not None]   # neutral zero counts
                    if not cols_ok:
                        stats['neutral'] = stats.get('neutral', 0) + 1
                        continue
                    ok = all(o for _, o in cols_ok)
                    bg = TIE if ok else max((c_ for c_, o in cols_ok if not o), key=lambda c_: -sum(c_))
                    vals = [c.strip() for c, p in zip(cs, ps) if p and p[0] == 'count']
                    content = '\\begin{tabular}{@{}c@{}}' + '\\\\'.join(vals) + '\\end{tabular}'
                    lns = [row['line'] for row in rows]
                    spans[(lns[0], j)] = lns
                    paint[(lns[0], j)] = (bg, content)
                    stats['tie' if ok else 'miss'] = stats.get('tie' if ok else 'miss', 0) + 1
            for idx, row in enumerate(rows):
                for j, c in enumerate(row['cells'][1:], 1):
                    if (row['line'], j) in paint or any(row['line'] in v and jj == j for (l0, jj), v in spans.items()):
                        continue
                    n, inner = unwrap(c)
                    p = parse(inner)
                    if p is None:
                        continue
                    d, pct, uns = hinfo(row['keys'][j])
                    if not d:             # a column without an arrow stays plain
                        continue
                    if p[0] == 'count':   # green when it meets its level, red (deeper = further) otherwise
                        bg, ok = count_colour(p, b['cert'], level)
                        if bg is None:    # zero events, missing only by the sample size: neutral
                            stats['neutral'] = stats.get('neutral', 0) + 1
                            continue
                        stats['tie' if ok else 'miss'] = stats.get('tie' if ok else 'miss', 0) + 1
                    else:
                        ref = best(bi, b['gS'][g], b['ga0'][g], row['keys'][j], d)
                        if ref is None:
                            stats['noref'] = stats.get('noref', 0) + 1
                            continue
                        t = grade(p[1], ref[0], d, pct, row['keys'][j])
                        bg = TIE if t is None else ramp(t)
                        stats[bin_of(t)] = stats.get(bin_of(t), 0) + 1
                        stats['means'] = stats.get('means', 0) + 1
                    if n > 1:
                        spans[(row['line'], j)] = [rw['line'] for rw in rows[idx:idx + n]]
                    paint[(row['line'], j)] = (bg, inner.strip())

    # ---- rewrite the data rows
    covered = {}
    for (l0, j), lns in spans.items():
        for ln_ in lns:
            covered[(ln_, j)] = (l0, lns)
    for t in tabs:
        for row in t['rows']:
            i = row['line']
            body, tail = out[i].rsplit('\\\\', 1)
            cells = split_cells(body)
            for j in range(len(cells)):
                if (i, j) in covered:
                    l0, lns = covered[(i, j)]
                    st = paint.get((l0, j))
                    if st is None:
                        continue
                    bg, content = st
                    if i != lns[-1]:
                        cells[j] = f' {cellcolor(bg)} '
                    else:
                        n = len(lns)
                        fix = (n - 1) * H_IN / 2
                        cells[j] = f' {cellcolor(bg)}\\multirow{{-{n}}}{{*}}[{fix:.2f}pt]{{\\color{{black}}{content}}}'
                elif (i, j) in paint:
                    bg, content = paint[(i, j)]
                    cells[j] = f' {cellcolor(bg)}{{\\color{{black}}{content}}} '
            # positive multirows left (method column, unpainted values): shift by the separators inside the span
            for j in range(len(cells)):
                m = re.match(r'^(\s*)\\multirow\{(\d+)\}\{\*\}\{', cells[j])
                if m:
                    n = int(m.group(2))
                    cells[j] = cells[j][:m.end() - 1] + f'[-{(n - 1) * H_IN / 2:.2f}pt]' + cells[j][m.end() - 1:]
            out[i] = '&'.join(cells) + '\\\\' + tail

    # ---- separator rows (rules drawn in the gutters between the coloured boxes)
    def seprow(cols, fills, hrow=None):
        parts = []
        if hrow is None:
            hrow = 2 * G + T_GRP if 'grp' in fills else (H_IN if any(f == 'in' or isinstance(f, tuple) for f in fills)
                                                        else G)
        for j, ((lm, rm), f) in enumerate(zip(cols, fills)):
            strut = f'\\vrule height{hrow:.2f}pt depth0pt width0pt' if j == 0 else ''
            if f in ('in', 'grp'):
                tt = T_IN if f == 'in' else T_GRP
                fill = f'\\leaders\\hrule height{G + tt:.2f}pt depth-{G}pt\\hfil'
            elif isinstance(f, tuple):
                fill = f'\\hskip{G}pt{{\\color[rgb]{{{rgb(f)}}}\\leaders\\hrule height{hrow:.2f}pt\\hfil}}\\hskip{G}pt'
            elif isinstance(f, list):   # half of the space between two values of one box, in the colour of its value
                fill = f'\\hskip{INSET}pt{{\\color[rgb]{{{rgb(f)}}}\\leaders\\hrule height{hrow:.2f}pt\\hfil}}\\hskip{INSET}pt'
            else:
                fill = '\\hfil'
            parts.append(f'\\omit{lm}{strut}{fill}{rm}')
        return '&'.join(parts) + '\\cr'

    insert_after, remove = {}, set()
    for t in tabs:
        cols = parse_preamble(t['pre'].replace('!{\\hskip\\arrayrulewidth}', '|'))
        nc = len(cols)
        rows = t['rows']
        for k, row in enumerate(rows):
            i = row['line']
            nxt = rows[k + 1] if k + 1 < len(rows) else None
            prev = rows[k - 1] if k else None
            if prev is None or prev['blk'] is not row['blk'] or (prev['line'] != i - 1 and (i - 1) not in remove):
                insert_after.setdefault(i - 1, []).append(seprow(cols, ['gut'] * nc))
            if nxt is not None and nxt['blk'] is row['blk'] and nxt['line'] == i + 1:
                fills = []
                for j in range(nc):
                    cv = covered.get((i, j))
                    if cv and i + 1 in cv[1]:
                        st = paint.get((cv[0], j))
                        fills.append(st[0] if st and st[0] != (1.0, 1.0, 1.0) else 'gut')
                    elif j == 0:
                        n0 = unwrap(next(r_['cells'][0] for r_ in rows if r_['g'] == row['g'] and r_['blk'] is row['blk']))[0]
                        fills.append('gut' if n0 > 1 and nxt['g'] == row['g'] else 'in')
                    elif j >= len(nxt['cells']) or not nxt['cells'][j].replace('\\vnsd', '').strip():
                        fills.append('gut')   # original multirow continuation (unpainted)
                    else:
                        fills.append('in')
                if nxt['g'] != row['g']:
                    fills = ['grp'] * nc
                    insert_after.setdefault(i, []).append(seprow(cols, fills))
                elif 'in' in fills:   # values of one box: no rule, each half in the colour of its own value
                    def half(line, j):
                        st = paint.get((line, j))
                        return list(st[0]) if st and st[0] != (1.0, 1.0, 1.0) else 'gut'
                    up = [half(i, j) if f == 'in' else f for j, f in enumerate(fills)]
                    lo = [half(i + 1, j) if f == 'in' else f for j, f in enumerate(fills)]
                    insert_after.setdefault(i, []).extend([seprow(cols, up, H_IN / 2), seprow(cols, lo, H_IN / 2)])
                else:
                    insert_after.setdefault(i, []).append(seprow(cols, fills))
            elif nxt is not None and nxt['blk'] is row['blk'] and out[i + 1].strip().startswith('\\Xcline'):
                remove.add(i + 1)
                insert_after.setdefault(i, []).append(seprow(cols, ['grp'] * nc))
            else:
                insert_after.setdefault(i, []).append(seprow(cols, ['gut'] * nc))
        if t['base'] and t['rows']:   # vertical rules in the baseline body as in the certified boxes
            out[t['start']] = out[t['start']].replace('!{\\hskip\\arrayrulewidth}', '|')
    res = []
    for i, ln in enumerate(out):
        if i not in remove:
            res.append(ln)
        res.extend(insert_after.get(i, []))
    stats['painted'] = len(paint)
    return '\n'.join(res)


def legend_body():
    r"""The colour key in plain words, two short lines in the IEEE editions (a sample colour stays with its label):
    the scale of the means (green, light green, pale yellow, light orange, soft red), the counts and the collision
    mark."""
    Hh = 6.5

    def sw(c, label=''):
        return (f'\\mbox{{\\raisebox{{-0.5pt}}{{\\textcolor[rgb]{{{rgb(c)}}}{{\\rule{{10pt}}{{{Hh}pt}}}}}}'
                + (f'\\,{label}' if label else '') + '}')
    name = ('\\mbox{{\\setlength{\\fboxsep}{1pt}\\colorbox[rgb]{' + rgb(BAD) + '}{name\\,$\\scriptstyle\\times$}}}')
    # Same swatches, plain wording.
    return ('\\ownerfix{Colors. Means:} ' + sw(TIE, '\\ownerfix{as good as the best certified variant of the same setting, or better,}') + ' '
            + sw(ramp(LT(5.0))) + '\\,' + sw(ramp(LT(10.0))) + '\\,' + sw(ramp(LT(30.0))) + '\\,'
            + sw(ramp(1.0), '\\ownerfix{increasingly worse.}') + ' \\ownerfix{Event counts:} '
            + sw(TIE, '\\ownerfix{95\\% bound within its level ($\\delta$ for certified variants, 5\\% for baselines),}')
            + ' ' + sw(MISS_LO) + '\\,' + sw(MISS_HI, '\\ownerfix{above it, darker for larger bounds.}') + ' '
            + name + '\\ownerfix{: collision in at least one trajectory.}')


def key_legend():
    r"""The colour key, run on after the last sentence of a key box (printed once per document, \vncolourlegend),
    so it adds the fewest lines."""
    return ('\\providecommand{\\vncolourlegend}[1]{#1}'
            '\\vncolourlegend{ ' + legend_body() + '}')


def report_legend():
    r"""The report's colour key, a paragraph below its first coloured table (printed once, \vncolourlegend)."""
    return ('\\providecommand{\\vncolourlegend}[1]{#1}\\vncolourlegend{\\par\\vspace{2pt}{\\raggedright\\noindent '
            + legend_body() + '\\par}}\n')


KEY_END = '\\tabularnewline\\Xhline{\\vnbox}\\end{tabular}}'
# The extended key defines its symbols itself; the coloured variant adds the colour pointer at its end.
KEY_REF = ('N-A: not applicable}}}; bold:', 'N-A: not applicable\\ownerfix{; colors as in Table~\\ref{tab:empirical-performance}}}}}; bold:')


# The report's "Key to all tables" fills its page; its colour key is printed below the first coloured table
# (REPORT_LEGEND_TABLE, Table 4, which ends its page), and the key box points to it within its last line.
REPORT_LEGEND_TABLE = "default_corridor"
REPORT_KEY_POINTER = '\\reponly{ The colors are explained below Table~\\ref{tab:sim-default-corridor}.}'


def with_legend(text, name):
    """Coloured variant of a fragment with a key box: the colour key inside the box (Section VI tables); the
    paper's "Key to all tables" says that the colours are as in Table I, the report's points to the colour key
    below Table 4 (report_legend)."""
    L = text.split('\n')
    k = next(i for i, ln in enumerate(L) if ln.rstrip().endswith('%KEYBOX')) + 1
    assert L[k].count(KEY_END) == 1, name
    if name == "table_key_arxiv":
        legend = key_legend().replace("Event counts:", "Event frequencies:").replace(";", ".")
        L[k] = L[k].replace(KEY_END, legend + KEY_END)
    elif name in KEY_TABLES:
        assert L[k].count(KEY_REF[0]) == 1, name
        L[k] = L[k].replace(KEY_REF[0], KEY_REF[1])
        L[k] = L[k].replace(KEY_END, REPORT_KEY_POINTER + KEY_END)
    else:
        L[k] = L[k].replace(KEY_END, key_legend() + KEY_END)
    return '\n'.join(L)


SEP_PT = 1.36                          # mean height (pt) that the separator rows add to a table row, net of the rules they replace
# the checks tables (one row per group, rules replaced by group separators): calibrated with compile_tests.py
SEP_PT_CHECKS = 0.92
BASELINE = {"paper": 12.0, "report": 13.192}   # table \baselineskip per target (as build_tables.uniform_stretch)
DEFAULT_STRETCH = 0.93                 # \ntablesetup (appendices/n200_assets.tex)


def compensate(text, target, sep=SEP_PT):
    """Keep the coloured table as tall as the plain one: every row loses the height its separator row adds
    (\\arraystretch lowered by SEP_PT / \\baselineskip), so the page blocks of the plain layout still fit."""
    d = sep / BASELINE[target]

    def lower(m):
        return f"\\renewcommand{{\\arraystretch}}{{{float(m.group(1)) - d:.3f}}}"
    text = re.sub(r"\\renewcommand\{\\arraystretch\}\{([\d.]+)\}", lower, text)
    return f"\\renewcommand{{\\arraystretch}}{{{DEFAULT_STRETCH - d:.3f}}}% coloured variant: rows + separators\n" + text


NO_COMPENSATION = ()   # every coloured table keeps the height of its plain variant
# Section VI tables (Table I): the colour key in the key box adds two lines; the coloured variant takes them back
# (rows as tight as their text allows, in pt per row; smaller spaces between the boxes) and plain_height sets it
# at the height of the plain variant.
SEP_PT_LEGEND = {"sec6_combined": 3.3, "sec6_short": 2.95, "sec6_maze_arxiv": 5.5}
SEC6_SPACES = {"2pt": "1pt", "3pt": "1.5pt", "5pt": "1pt", "6pt": "3.5pt"}   # 5pt: sec6_maze_arxiv


def tighten(text):
    """Section VI tables: smaller spaces between the boxes (SEC6_SPACES)."""
    return re.sub(r'\\vspace\{(\d+pt)\}', lambda m: '\\vspace{' + SEC6_SPACES.get(m.group(1), m.group(1)) + '}', text)


BREAK = '\\ntablebreak'
# top-level settings of a page block that the next page blocks of the fragment inherit
CARRY = re.compile(r'^\\(renewcommand\{\\arraystretch\}|normalsize|small|footnotesize|scriptsize|setlength|providecommand)')


def blocks(text):
    """The page blocks of a fragment (split at its \\ntablebreak lines), each with the top-level settings it
    inherits from the blocks before it."""
    out, carry, cur = [], [], []
    for ln in text.rstrip('\n').split('\n'):
        if ln.strip() == BREAK:
            out.append('\n'.join(cur))
            cur = list(carry)
            continue
        if CARRY.match(ln):
            carry.append(ln)
        cur.append(ln)
    out.append('\n'.join(cur))
    return out


def plain_height(name, coloured, plain):
    r"""Every page block of the coloured variant in a box of exactly the height and depth of the same block of the
    plain variant (measured first: height plus depth \vnplainh@<name>@<k>, depth \vnplaind@<name>@<k>; the
    paragraph skip, natural size 0pt, set before the box so that it stretches as in the plain variant), so the
    floats and the running text, and with them every page break of the edition, are those of the plain layout. A
    coloured block taller than its plain block keeps its own height and is reported in the log (VN-COLOUR-TALL)."""
    P, C = blocks(plain), blocks(coloured)
    assert len(P) == len(C), name
    head = ['\\ifdefined\\vnplainbox\\else\\newbox\\vnplainbox\\newbox\\vncolourbox\\fi\\par']
    for k, b in enumerate(P, 1):
        head.append('\\setbox\\vnplainbox\\vbox{\\parskip=0pt\\relax%\n' + b + '\n}\\expandafter\\xdef\\csname vnplainh@' + name
                    + f'@{k}' + '\\endcsname{\\the\\dimexpr\\ht\\vnplainbox+\\dp\\vnplainbox\\relax}'
                    + '\\expandafter\\xdef\\csname vnplaind@' + name + f'@{k}' + '\\endcsname{\\the\\dp\\vnplainbox}%')
    body = []
    for k, b in enumerate(C, 1):
        h = '\\csname vnplainh@' + name + f'@{k}' + '\\endcsname'
        d = '\\csname vnplaind@' + name + f'@{k}' + '\\endcsname'
        body.append('\\setbox\\vncolourbox\\vbox{\\parskip=0pt\\relax%\n' + b + '\n}%\n'
                    + '\\ifdim\\dimexpr\\ht\\vncolourbox+\\dp\\vncolourbox\\relax>' + h
                    + '\\typeout{VN-COLOUR-TALL ' + name + f' {k} ' + '\\the\\dimexpr\\ht\\vncolourbox+\\dp\\vncolourbox-'
                    + h + '\\relax}\\else\\setbox\\vncolourbox\\vbox to' + h + '{\\unvbox\\vncolourbox\\vss}'
                    + '\\ht\\vncolourbox=\\dimexpr' + h + '-' + d + '\\relax\\dp\\vncolourbox=' + d + '\\fi'
                    + '\\vskip\\parskip\\box\\vncolourbox\\par')
    return '\n'.join(head) + '\n' + ('\n' + BREAK + '\n').join(body) + '\n'


LEGEND_TABLES = ("sec6_combined", "sec6_short", "sec6_maze_arxiv") + KEY_TABLES   # the key boxes that carry the colour key


def compact_shade(name, plain, records, stats):
    """Shades the compact split columns without changing the plain table or numeric records.

    Only the mean half is shaded. SDs, peak-risk diagnostics and undefined values remain neutral. Event halves
    share the existing bound-versus-level color. The comparison pool contains only certified rows retaining
    their guarantee, within each setting; baseline rows use compatible settings at their alpha_0.
    """
    import json
    grouped, rows = {}, []
    for e in records:
        if e['fragment'] != name or e['target'] != 'paper':
            continue
        key = (e['table'], e['method'], json.dumps(e['row'], sort_keys=True))
        if key not in grouped:
            grouped[key] = []
            rows.append(grouped[key])
        grouped[key].append(e)
    assert rows and all(len(es) == 6 for es in rows), name
    specs = ('m_c', 'T_mis', 'q_CBF', 'pk', 'qdep', 'qcert')
    assert all(tuple(e['spec'] for e in es) == specs for es in rows), name
    keys = {'m_c': r'\widehat m_c/c_{\max}\downarrow',
            'T_mis': r'\widehat T_{\rm mis}\downarrow(s)',
            'q_CBF': r'\widehat q^{+}_{\rmCBF}\downarrow(\%)'}

    def event(e):
        if e['spec'] not in ('qdep', 'qcert') or '---' in e['printed']:
            return None
        a, u = e['printed'].split('&')
        bound = float(re.search(r'\(([\d.]+)\)', u).group(1))
        count = e.get('count', float(a.split('--')[-1]))
        certified = not e['row'].get('baseline', False)
        level = 100 * e['row']['delta'] if certified else BASE_LEVEL
        return ('count', bound <= level, bound, count), certified, level

    def eligible(es):
        return not es[0]['row'].get('baseline', False) and all(
            count_ok(x[0], x[1]) for e in es if (x := event(e)) is not None)

    def ref(e):
        row, values = e['row'], []
        for es in rows:
            if not eligible(es):
                continue
            other = es[0]['row']
            if row.get('baseline'):
                if e['method'] != 'Conservative' and abs(row['alpha0'] - other['alpha0']) > 1e-9:
                    continue
            elif row.get('setting') != other.get('setting'):
                continue
            candidate = next(x for x in es if x['spec'] == e['spec'])
            if 'mean_str' in candidate:
                values.append((float(candidate['mean_str']), float(candidate['sd_str'])))
        return min(values) if values else None

    replacements = []
    def painted(part, bg):
        # Explicit zero overhang keeps the adjacent SD half neutral and preserves the natural table width.
        new = rf'\vncompactcell{{{rgb(bg)}}}{{\color{{black}}{part}}}'
        replacements.append((new, part))
        return new

    output, ri = [], 0
    for line in plain.splitlines():
        chunks = line[:-2].split(' & ') if line.endswith('\\\\') else []
        if ri >= len(rows) or len(chunks) != 7 or chunks[1:] != [e['printed'] for e in rows[ri]]:
            output.append(line)
            continue
        es = rows[ri]
        method, row = es[0]['method'], es[0]['row']
        scenario = 'Maze' if row['scenario'] == 'maze' else 'Corridor'
        if (scenario, norm_name(chunks[0])) in COLLISIONS:
            original = chunks[0]
            # Zero-width cross to the left of the name preserves the compact table's width.
            cross = r'\llap{$\scriptstyle\times$\,}' + original
            replacements.append((cross, original))
            chunks[0] = painted(cross, BAD)
            stats['collision'] = stats.get('collision', 0) + 1
        if method not in ('Unfiltered', 'Zero margin'):
            for j, e in enumerate(es, 1):
                if '---' in e['printed'] or e['spec'] == 'pk':
                    continue
                halves = e['printed'].split('&')
                assert len(halves) == 2, e
                ev = event(e)
                if ev is not None:
                    bg, ok = count_colour(*ev)
                    if bg is None:
                        continue
                    chunks[j] = '&'.join(painted(part, bg) for part in halves)
                    stats['tie' if ok else 'miss'] = stats.get('tie' if ok else 'miss', 0) + 1
                elif e['spec'] in keys and 'mean_str' in e:
                    best = ref(e)
                    if best is None:
                        continue
                    k = keys[e['spec']]
                    t = grade(float(e['mean_str']), best[0], -1, e['spec'] == 'q_CBF', k)
                    chunks[j] = painted(halves[0], TIE if t is None else ramp(t)) + '&' + halves[1]
                    stats['means'] = stats.get('means', 0) + 1
                else:
                    continue
                stats['painted'] = stats.get('painted', 0) + 1
        output.append(' & '.join(chunks) + r'\\')
        ri += 1
    assert ri == len(rows), (name, ri, len(rows))
    coloured = '\n'.join(output) + '\n'
    restored = coloured
    for new, old in reversed(replacements):
        restored = restored.replace(new, old)
    assert restored == plain, (name, 'colour adapter changed content or layout')
    assert stats.get('painted', 0) > 0, name
    setup = (r'\makeatletter\providecommand{\vncompactcell}[1]{\gdef\CT@cell@color{'
             r'\@tempdimb0pt\@tempdimc0pt\CT@color[rgb]{#1}\global\let\CT@cell@color\relax}}'
             r'\makeatother' + '\n')
    return setup + coloured


def wrap(name, plain, stats=None, target="paper", registry=None, records=None):
    """Fragment body with both variants; the plain variant is `plain` unchanged. `registry`: verdicts from
    collect() over all fragments (both targets)."""
    stats = stats if stats is not None else {}
    if name in COMPACT_TABLES:
        coloured = compact_shade(name, plain, records or [], stats)
    elif name in KEY_TABLES or name.startswith(STRESS_PREFIX):   # phi stress tables: plain in every edition
        coloured = plain
    else:
        coloured = shade(plain, stats, registry=dict(registry or {}), stress=name.startswith(STRESS_PREFIX),
                         scenario=name in SCENARIO_TABLES)
        if coloured == plain:
            return plain
        if name not in NO_COMPENSATION:
            coloured = compensate(coloured, target, SEP_PT_LEGEND.get(name, SEP_PT_CHECKS if name.startswith('checks_')
                                                                     else SEP_PT))
    if name in LEGEND_TABLES:
        coloured = with_legend(coloured, name)
    if name in SEP_PT_LEGEND:
        coloured = tighten(coloured)
    if coloured != plain and name not in COMPACT_TABLES:
        coloured = plain_height(name, coloured, plain)
    if name in ("main_maze_arxiv_v8", "main_maze_extended_v8"):
        legend = legend_body()
        if name == "main_maze_arxiv_v8":
            legend = legend.replace("Event counts:", "Event frequencies:")
        coloured += (r"\par\vspace{2pt}{\footnotesize\raggedright\noindent "
                     + r"\vncolourlegend{" + legend + r"}\par}" + "\n")
    if target == "report" and name == REPORT_LEGEND_TABLE:
        coloured += report_legend()
    return '\n'.join([MARK_DEF, MARK_IF.replace('NAME', name), coloured.rstrip('\n'),
                      MARK_ELSE, plain.rstrip('\n'), MARK_END]) + '\n'


def plain_variant(text):
    """The plain (uncoloured) text of a generated fragment: the lines before the switch plus its else branch.
    Identity for a fragment without the switch."""
    if MARK_DEF not in text:
        return text
    L = text.split('\n')
    a = L.index(MARK_DEF)
    e = L.index(MARK_ELSE)
    f = len(L) - 1 - L[::-1].index(MARK_END)
    return '\n'.join(L[:a] + L[e + 1:f]) + '\n'
