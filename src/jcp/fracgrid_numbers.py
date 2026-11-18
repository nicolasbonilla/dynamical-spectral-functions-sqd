# -*- coding: utf-8 -*-
"""fracgrid_numbers.py -- the numbers of Table tab:fracgrid (Sec. S7) and of every sentence that quotes it, computed
from data/eta_sweep.json, and the passages as printed.  Read by src/verify_jcp.py (check_fracgrid,
check_directed_rounding), which requires every passage to be printed verbatim.

Deposit 3.5.0: numbers(), cell() and edits() of the one-off rewriting script of the JCP package (2026-11-08, pass 7
of the number sweep), copied verbatim; only the data path is the deposit's.  That script also rewrote the sources
and appended the replaced text to the record of withdrawn numbers (docs/JCP_WITHDRAWN_NUMBERS.md); that part is not
deposited, so nothing here writes a file.  Each entry of edits() is (file, the text before 2026-11-08, the text
printed now).
"""
import io
import math
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))                      # <repo>/src/jcp
DATA = os.path.normpath(os.path.join(HERE, '..', '..', 'data'))
ETAS = ['0.05', '0.075', '0.1', '0.15', '0.18', '0.25', '0.35', '0.5']
FIT_ETAS = ['0.05', '0.1', '0.15', '0.25', '0.5']


def numbers():
    """Every number the rewritten passages print, computed from data/eta_sweep.json."""
    d = json.load(io.open(os.path.join(DATA, 'eta_sweep.json'), encoding='utf-8'))
    pts = {p['L']: p for p in d['points']}
    F = d['fits']
    n = {}
    n['rows'] = {L: [pts[L]['frac'][e] for e in ETAS] for L in sorted(pts)}
    n['bf'] = [F[e]['all_L']['fraction']['b'] for e in FIT_ETAS]
    n['r2f'] = [F[e]['all_L']['fraction']['R2'] for e in FIT_ETAS]
    n['bs'] = [F[e]['all_L']['absolute_support']['b'] for e in FIT_ETAS]
    n['r2s'] = [F[e]['all_L']['absolute_support']['R2'] for e in FIT_ETAS]
    n['bs6'] = [F[e]['L_ge_6']['absolute_support']['b'] for e in FIT_ETAS]
    n['shift'] = n['bs'][0] - n['bs'][-1]
    n['shift_f'] = n['bf'][0] - n['bf'][-1]
    n['factor'] = n['bf'][-1] / n['bf'][0]
    n['depth'] = sorted(set(p['n_lanczos'] for p in d['points']))
    n['pub15'] = [pts[L]['frac']['0.15'] for L in sorted(pts)]
    n['L12_05'] = pts[12]['frac']['0.5']
    n['L14'] = d['published_sequence_L4_to_L14']['values'][-1]
    return n


def cell(x, fmt):
    return '$' + (fmt % x) + '$'


def edits(n):
    bf, bs, bs6, r2f, r2s = n['bf'], n['bs'], n['bs6'], n['r2f'], n['r2s']
    sh = n['shift']
    assert abs(n['shift'] - n['shift_f']) < 1e-9, 'support and fraction shifts must be equal (S = frac*dim)'
    assert n['depth'] == [260], n['depth']
    pub = ', '.join('$%.4f$' % x for x in n['pub15'])
    rows = '\n'.join('$%d$%s & %s\\\\' % (L, ' ' * (3 - len(str(L))), ' & '.join(cell(x, '%.3f') for x in v))
                     for L, v in n['rows'].items())
    blk = lambda v, f: ' & '.join([cell(v[0], f), '---', cell(v[1], f), cell(v[2], f), '---', cell(v[3], f), '---',
                                   cell(v[4], f)])
    E = []
    # ---- sm_prices.tex ------------------------------------------------------------------------
    E.append(('sm_prices.tex',
              'applies to the object fitted: $1.082$ is the fit to the published points, while the $1.090$ in the\n'
              'bottom block of Table~\\ref{tab:fracgrid} is the fit to that table\'s own coarser replication grid at the\n'
              'same resolution. They are two fits to two datasets and the difference is theirs, not a rounding.',
              'applies to the object fitted: $1.082$ is the fit to the published points, and it is also the\n'
              '$\\eta=0.15\\,t$ entry of the bottom block of Table~\\ref{tab:fracgrid}, whose $\\eta=0.15\\,t$ column\n'
              'reproduces those points.'))
    E.append(('sm_prices.tex',
              'across the decade $\\eta=0.05\\to0.50\\,t$ both shift by $0.148$ in absolute value',
              'across the decade $\\eta=0.05\\to0.50\\,t$ both shift by $%.3f$ in absolute value' % sh))
    E.append(('sm_prices.tex',
              '$0.148$ is three quarters of the quantity itself: $e^{-0.13\\,L}$ becomes $e^{-0.27\\,L}$, and the',
              '$%.3f$ is more than half of the quantity itself: $e^{-%.2f\\,L}$ becomes $e^{-%.2f\\,L}$, and the'
              % (sh, -bf[0], -bf[-1])))
    E.append(('sm_prices.tex',
              '$0.148$ carries it from $1.165$ to $1.017$---a seventh of its own value, $13\\%$ measured against the\n'
              'larger end and $15\\%$ against the smaller---',
              '$%.3f$ carries it from $%.3f$ to $%.3f$---a seventh of its own value, $%d\\%%$ measured against the\n'
              'larger end and $%d\\%%$ against the smaller---' % (sh, bs[0], bs[-1], round(100 * sh / bs[0]),
                                                                  int(100 * sh / bs[-1]))))
    E.append(('sm_prices.tex',
              'snapshots; \\texttt{scaling\\_lanczos\\_mf.py}), same threshold, same depth $n_\\ell=150$; only $\\eta$\n'
              'changes. The scan over subspace size in this replication is coarser than the bisection that produced\n'
              'the published points, and returns the smallest scanned fraction meeting the threshold, so the\n'
              '$\\eta=0.15\\,t$ column is an upper bound on them: against the published $0.9167$, $0.8200$, $0.5561$,\n'
              '$0.3625$, $0.1700$ it reads $0.917$, $0.820$, $0.560$, $0.380$, $0.180$, agreeing at the two smallest\n'
              'sizes and standing $0.7\\%$, $4.8\\%$ and $5.9\\%$ high at the three largest. The same granularity is why\n'
              'the fraction exponent of this column, $-0.2012$, is not identical to the $-0.2093$ obtained by fitting\n'
              'the published points themselves over $L=4$--$12$. Bottom block: least-squares exponents.',
              'snapshots), same threshold; only $\\eta$ changes (\\texttt{eta\\_sweep.py},\n'
              '\\texttt{data/eta\\_sweep.json}: Lanczos depth $%d$, a coarse scan over the subspace size followed by\n'
              'bisection to $\\dim/200$). The $\\eta=0.15\\,t$ column reproduces the published points, %s, and its\n'
              'fraction exponent is the $%.4f$ of those points over $L=4$--$12$. Bottom block: least-squares\n'
              'exponents.' % (n['depth'][0], pub, bf[2])))
    E.append(('sm_prices.tex',
              'criticizes exactly that elsewhere in this section. The fraction exponent moves by a factor $2.2$ over\n'
              'the decade, from $-0.1266$ to $-0.2743$, while the support exponent moves from $1.165$ to $1.017$; in\n'
              'absolute value both move by $0.148$, which is the whole content of the\n'
              'comparison.}',
              'criticizes exactly that elsewhere in this section. The fraction exponent moves by a factor $%.1f$ over\n'
              'the decade, from $%.4f$ to $%.4f$, while the support exponent moves from $%.3f$ to $%.3f$; in\n'
              'absolute value both move by $%.3f$, which is the whole content of the\n'
              'comparison.}' % (n['factor'], bf[0], bf[-1], bs[0], bs[-1], sh)))
    old_rows = ('$4$  & $0.917$ & $0.917$ & $0.917$ & $0.917$ & $0.917$ & $0.917$ & $0.917$ & $0.875$\\\\\n'
                '$6$  & $0.900$ & $0.890$ & $0.870$ & $0.820$ & $0.800$ & $0.730$ & $0.680$ & $0.630$\\\\\n'
                '$8$  & $0.740$ & $0.680$ & $0.640$ & $0.560$ & $0.540$ & $0.480$ & $0.420$ & $0.380$\\\\\n'
                '$10$ & $0.520$ & $0.480$ & $0.440$ & $0.380$ & $0.340$ & $0.260$ & $0.220$ & $0.200$\\\\\n'
                '$12$ & $0.340$ & $0.280$ & $0.240$ & $0.180$ & $0.160$ & $0.140$ & $0.120$ & $0.100$\\\\')
    E.append(('sm_prices.tex', old_rows, rows))
    old_fit = ('fraction exponent, $L{=}4$--$12$ & $-0.1266$ & --- & $-0.1681$ & $-0.2012$ & --- & $-0.2395$ & --- & $-0.2743$\\\\\n'
               '$R^{2}$                          & $0.904$   & --- & $0.918$   & $0.929$   & --- & $0.970$   & --- & $0.982$\\\\\n'
               'support exponent, $L{=}4$--$12$  & $1.165$   & --- & $1.123$   & $1.090$   & --- & $1.052$   & --- & $1.017$\\\\\n'
               '$R^{2}$                          & $0.999$   & --- & $0.999$   & $0.998$   & --- & $0.999$   & --- & $0.999$\\\\\n'
               'support exponent, $L{=}6$--$12$  & $1.136$   & --- & $1.088$   & $1.053$   & --- & $1.022$   & --- & $0.992$\\\\')
    new_fit = ('fraction exponent, $L{=}4$--$12$ & %s\\\\\n' % blk(bf, '%.4f') +
               '$R^{2}$                          & %s\\\\\n' % blk(r2f, '%.3f') +
               'support exponent, $L{=}4$--$12$  & %s\\\\\n' % blk(bs, '%.3f') +
               '$R^{2}$                          & %s\\\\\n' % blk(r2s, '%.3f') +
               'support exponent, $L{=}6$--$12$  & %s\\\\' % blk(bs6, '%.3f'))
    E.append(('sm_prices.tex', old_fit, new_fit))
    E.append(('sm_prices.tex',
              'is claimed here. Over the decade of $\\eta$ the fraction exponent moves by a factor $2.2$, from\n'
              '$-0.1266$ at $\\eta=0.05\\,t$ to $-0.2743$ at $\\eta=0.50\\,t$. At $\\eta=0.50\\,t$ the fraction at $L=12$\n'
              'is already $0.10$, a quarter above the $L=14$ value at the resolution used throughout this work and\n'
              'reached two sizes earlier:',
              'is claimed here. Over the decade of $\\eta$ the fraction exponent moves by a factor $%.1f$, from\n'
              '$%.4f$ at $\\eta=0.05\\,t$ to $%.4f$ at $\\eta=0.50\\,t$. At $\\eta=0.50\\,t$ the fraction at $L=12$\n'
              'is already $%.3f$, within $%d\\%%$ of the $L=14$ value at the resolution used throughout this work and\n'
              'reached two sizes earlier:' % (n['factor'], bf[0], bf[-1], n['L12_05'],
                                              math.ceil(100 * (n["L12_05"] / n["L14"] - 1)))))
    # ---- sec_7.tex ----------------------------------------------------------------------------
    E.append(('sec_7.tex',
              'support exponents differ by the constant $b_{\\dim}=1.291$ and move together, by $0.148$ over the decade\n'
              '$\\eta=0.05\\to0.50\\,t$ (Table~\\ref{tab:fracgrid}): a seventh of the support exponent, three quarters of\n'
              'the fraction exponent.',
              'support exponents differ by the constant $b_{\\dim}=1.291$ and move together, by $%.3f$ over the decade\n'
              '$\\eta=0.05\\to0.50\\,t$ (Table~\\ref{tab:fracgrid}): a seventh of the support exponent, more than half\n'
              'of the fraction exponent.' % sh))
    E.append(('sec_7.tex',
              '$2.2$, from $-0.1266$ to $-0.2743$ (Table~\\ref{tab:fracgrid}), and even at fixed $\\eta$ the local',
              '$%.1f$, from $%.4f$ to $%.4f$ (Table~\\ref{tab:fracgrid}), and even at fixed $\\eta$ the local'
              % (n['factor'], bf[0], bf[-1])))
    # ---- caption of Fig. 12 -------------------------------------------------------------------
    E.append(('figs/caption_fig12.tex',
              '$\\lvert\\mathcal{S}\\rvert$ moves from $1.165$ to $1.017$ and that of the fraction in panel~(b) by a\n'
              'factor $2.2$, from $-0.1266$ to $-0.2743$ (Table~\\ref{tab:fracgrid}, $L=4$--$12$). In absolute value\n'
              'both move by the same $0.148$:',
              '$\\lvert\\mathcal{S}\\rvert$ moves from $%.3f$ to $%.3f$ and that of the fraction in panel~(b) by a\n'
              'factor $%.1f$, from $%.4f$ to $%.4f$ (Table~\\ref{tab:fracgrid}, $L=4$--$12$). In absolute value\n'
              'both move by the same $%.3f$:' % (bs[0], bs[-1], n['factor'], bf[0], bf[-1], sh)))
    # ---- sec_3_1.tex (main text) --------------------------------------------------------------
    E.append(('sec_3_1.tex',
              'at which the resource scan of Sec.~\\ref{sec:prices} operates (Fig.~\\ref{fig:scaling}; the coarser\n'
              'replication grid of Table~\\ref{tab:fracgrid} reads $0.180$ at $L=12$ for the same operating point).',
              'at which the resource scan of Sec.~\\ref{sec:prices} operates (Fig.~\\ref{fig:scaling}).'))
    return E
