# Numbers withdrawn from the supplementary material on 2026-11-08 (rule: no printed number without a source)

> **Deposit 3.5.0.** This is the record kept by the package of the manuscript prepared for The Journal of Chemical Physics (written 2026-11-08 and later), deposited because `src/check_trim.py` points to it, through its `@doc` proofs, for each fact it records as deliberately withdrawn. In it, "release" means deposit 3.4.1 (tag v3.4.1, the arXiv-v4 text); references to the author's working notes and folders, which are not in the deposit, were removed, and nothing else changed.

Each number below was printed only in the supplementary material. It came from checks run in September 2026
(recorded only in the author's working notes) whose output files were never deposited.

A search of the author's working folders and of the deposit on 2026-11-08 found no output file for any
of them:
- only notes quote them;
- `src/frontier/c0_*.py` is code, not output;
- `data/c3_frontier/validation.json` holds a different comparison.

They are therefore withdrawn. Each sentence now states only what the proof or the deposited data support.
The verbatim text, before the change, follows.

## S2, `sm_leakage.tex` (paragraph "What object is certified")

> The convenient inequality $\Lambda_{\mathcal K}\ge\Lambda_S$ that would transfer the statement to the dense
> $A_S$ is false: $1169$ of $6000$ random triples violate it, at a minimum ratio of $0.5258$. [...] Full
> reorthogonalization is mandatory throughout: without it $P_{\mathcal K}$ is not a projector---the
> orthogonality defect reaches $0.238$ at $L=8$, where the captured weight then reads $1.4868$, an impossible
> value---and $\Lambda$ is inflated thirtyfold.

Now:
- the inequality is "not implied by the proof and is not used";
- on the one subspace where both were computed (L=8, 85 % of the sector) $\Lambda_{\mathcal K}$ and
  $\Lambda_S$ agree to within 1 % (`data/c3_frontier/validation.json`, `ratio_K_S` 1.0003–1.0097);
- with reorthogonalization the orthogonality defect is at most $4.9\times10^{-15}$ (published-subspace
  rows).

## S4, `sm_frontier.tex`

> ... the convenient inequality $\Lambda_{\mathcal K}\ge\Lambda_S$ that would extend it to the dense $A_S$
> being false ($1169$ of $6000$ random triples violate it, minimum ratio $0.5258$) ...

Now it says "not being implied by the proof (Sec. S2)".

## S6, `sm_nogo.tex` (the fibers of Table S4)

> An exhaustive search over all $276$ coordinate subspaces of that size at $K=2$ settles both ends. At the
> aligned endpoint the global minimum of the certificate over the $276$ is $0.856$ and it is attained by the
> protocol's own subspace: the subspace the sampler would return is the best available one. At the rotated
> endpoint the minimum over the same $276$ is exactly the trivial bound, attained only by subspaces of
> vanishing captured weight: \emph{no} determinant subspace of that size certifies anything there.

> (iii) The small-leakage endpoint is realized only by the aligned basis, which is the stabilizer of the
> state ($|\braket{\Psi(0)|\Psi(U_{\rm in})}|=1.000000000000$); at the fiber angle $s=\pi/16$, every
> generator that actually moves the state gives $\widehat\Lambda/\eta\in[2.34,5.61]$. That range is a
> statement about generators at one angle and not about the angle sweep of Sec. app:witness, whose rotated
> endpoint descends to $1.61$ as $s\to0$ because the rotation itself is being switched off.

Now:
- the search sentences are removed;
- (iii) keeps only "realized by the aligned basis, which is the stabilizer of the state: the rotation that
  produces it leaves the state unchanged".

## S5, `sec_6_body.tex`

> The fibers (Table S4), an exhaustive search over the $276$ subspaces of that size at $K=2$, and four limits
> of the exhibit ...

Now it says "The fibers (Table S4) and four limits of the exhibit ...".

## Not changed, for the author (superseded the same day)

This section first listed the fiber-table numbers (0.856, 1.026, 1.161, 1.962; the margins 2.29, 1.91
and 1.69; "a factor 6-7"; 2.336 and 0.238 in the K=3 control; the 58-point subfamily with 2.8 and 8.3;
"1600 hostile cases") as left in place.  The author then applied the same rule to them: they are
withdrawn, with the verbatim text, in the next section and in the pass-7 record of Proposition
`prop:leakinv` at the end of this file.  Nothing of the fiber exhibit is printed any more.

## What is still not sourced (updated at the end of pass 8)

After pass 8 the systematic sweep (a script of the JCP package, not deposited, and its table `number_sweep.tsv`)
finds **no** printed number without a source (class U = 0, main text and SI).  168 numbers of the
supplementary material (none in the main text) are still class S: equal, at their printed precision, to a
value inside a deposited file that the passage names, but not tied to a key path or a recomputation, so a
guard would not notice if they changed.  They are listed in `number_sweep.tsv` (class S): sm_frontier 63,
sm_physics 26, sec_7 22, sm_prices 11, sec_3_app 7, app_stats 7, app_dequant 7, caption_fig12 6, sec_8 5,
app_carried_repro 4, sec_6_body 3, caption_thm1iii_violation 3, sec_4_app 2, sm_leakage 1, sec_5_app 1.

## Withdrawn the same day (second decision): the pairing-model fiber exhibit

Table S4 (`tab:fiber`), the angle-sweep table (`tab:fibersweep`) and every sentence that quotes their numbers came from the September C1 computation (working notes, not deposited), whose output was never deposited; no output file was found (searched 2026-11-08).  Proposition 2 keeps its analytic content (invariance along the orbit; the invariant bound vacuous by at least W/2eta; an invariant bound must exceed the fiber maximum).  The verbatim text removed, file by file:

### `sm_nogo.tex`

    bound on $\Lambda_S(\eta)/\eta$ \emph{does} exist---it is exhibited inside the proposition below, and we
    tested it on $1600$ hostile cases (GOE/GUE, near-degenerate, banded, $\lVert H\rVert\sim300$) with zero
    violations and with its first link attained, minimum ratio $1.0000000000$.

### `sm_nogo.tex`

    Table~\ref{tab:fiber} is the exhibit. The two endpoints of each fiber are the \emph{same} rotation angle
    applied within one geminal block and across two blocks, which is why the captured weight agrees bit for
    bit.
    % 2026-11-08 (JCP, rule "no printed number without a source"): the sentences on an exhaustive search
    % over 276 subspaces (minimum 0.856 at the aligned endpoint, the trivial bound at the rotated one) are
    % withdrawn: that search output was never deposited.  See docs/JCP_WITHDRAWN_NUMBERS.md.
    
    \begin{table*}[tp]\centering\small
    \caption{\textbf{Fibers of fixed invariants with moving leakage}, at $s=\pi/16$, $\delta=0.01$,
    $\eta=0.10$, $|\mathcal S|=2^{K-1}$, $w_S=0.961940$, $\tau(w_S)=0.795$. Along each row the spectrum,
    every Lehmann weight, $\lVert\phi\rVert^2$, $|\mathcal S|$ and all $\mathcal F_k$ are identical between
    the two endpoints ($w_S$ to $2.2\times10^{-16}$), so the first term of Eq.~\eqref{eq:twobudget} is fixed
    and only the second moves. Quantities are relative ($\widehat\Lambda_S=\Lambda_S/\lVert\phi\rVert$) and
    ``certificate'' is $\min\{(1+\sqrt{w_S})(\sqrt{1-w_S}+\widehat\Lambda_S/\eta),\,1+w_S\}$, the relative
    form of Theorem~\ref{thm:leakage}; ``al.''\ and ``rot.''\ label the aligned and rotated endpoints of the
    fiber. Bold marks a non-vacuous certificate, below the trivial bound $1+w_S$. A continuous angle sweep
    is in Sec.~\ref{app:witness}.}
    \label{tab:fiber}
    \setlength{\tabcolsep}{4pt}
    \begin{tabular}{@{}lccccc@{}}
    \toprule
    $K$ ($\dim(N{+}1)$) & $\mathcal F_1{=}\mathcal F_2{=}\mathcal F_3$
     & $\widehat\Lambda/\eta$ al. & $\widehat\Lambda/\eta$ rot.
     & cert.\ al. & cert.\ rot. \\
    \midrule
    $2$ \ ($24$)   & $8$  & $0.237$ & $2.326$ & $\mathbf{0.856}$ & $1.962$ (trivial) \\
    $3$ \ ($300$)  & $12$ & $0.323$ & $2.336$ & $\mathbf{1.026}$ & $1.962$ (trivial) \\
    $4$ \ ($3920$) & $16$ & $0.391$ & $2.346$ & $\mathbf{1.161}$ & $1.962$ (trivial) \\
    \bottomrule
    \end{tabular}
    \end{table*}

### `sm_nogo.tex`

    \paragraph{What this exhibit is not.} Four limits are declared with it, and they are why we call it a
    fiber and not a sweep of the orbit. \emph{(i)}~The evidence is a matched pair at $s=\pi/16$ for $K=2,3,4$.
    At the maximal-pairing angle $s=\pi/4$ the Born distribution of $\phi$ is exactly flat over $2^{K}$
    determinants, so ``the top $|\mathcal S|$ of the Born distribution'' does not define a subspace at all
    ($6$, $70$ and $12870$ equally valid outputs at $K=2,3,4$); any contrast quoted at that angle is a
    tie-breaking artifact and not a property of the fiber, and the sweep of Sec.~\ref{app:witness} stops
    at $s=0.654$ for that reason. \emph{(ii)}~Nothing here grows with the system, and two separate quantities shrink. The ratio of the
    two leakages falls, $9.80\to7.23\to6.01$ from $K=2$ to $K=4$, and it dies as the perturbation is turned
    up, from $28.7$ at $\delta=0.0025$ to $1.03$ at $\delta=0.3$. More to the point---because it is the
    quantity the proposition actually rests on---so does the margin of the certificate over the trivial
    bound: $1.962/0.856=2.29$, $1.962/1.026=1.91$, $1.962/1.161=1.69$, a fall of $26\%$ across the only
    three sizes available. Where that fall would reach unity we do not say: a linear extrapolation in $K$,
    a log-linear one and the last step alone place it between $K=6$ and $K=8$, the decrements are
    themselves decelerating, and three points do not support a number---which is exactly the objection this
    paper makes to its own fraction headline in Sec.~\ref{sec:prices}. A statement of non-vacuity needs one
    fiber and not a magnitude, and the fiber at $K=4$ is still non-vacuous by $69\%$; but ``the effect grows
    with $K$'' is false here. \emph{(iii)}~The small-leakage
    % 2026-11-08 (JCP, rule "no printed number without a source"): the overlap 1.000000000000, the range
    % [2.34, 5.61] over generators that move the state and the value 1.61 of the angle sweep came from a
    % September check whose output was never deposited; they are withdrawn (docs/
    % JCP_WITHDRAWN_NUMBERS.md), and the clause keeps only the qualitative statement.
    endpoint is realized by the aligned basis, which is the stabilizer of the state: the rotation that
    produces it leaves the state unchanged. The rhetorical parallel with the cost witness---where
    the orbit contains an $O(K)$-cost point and points of plausibly exponential cost---is therefore lost; Sec.~\ref{sec:nogo:asym} states
    that asymmetry. \emph{(iv)}~The scale is a toy: $2K\le8$ orbitals, sector
    dimension $\le3920$, no instance at the $L=12$ Hubbard scale of Sec.~\ref{sec:frontier}.

### `sm_nogo.tex`

    What the exhibit delivers survives all four limits and is measured: on fibers where every invariant is
    fixed to $10^{-14}$ and the captured Born weight to $2.2\times10^{-16}$, the certificate of
    Theorem~\ref{thm:leakage} separates a non-vacuous value from the trivial bound in every row of
    Table~\ref{tab:fiber}, and a quantity constant on those fibers cannot make that distinction. The a priori
    predictor asked for at the head of this section does not exist in the invariant class;

### `sm_nogo.tex`

    cost is $O(K)$, and fibers (Proposition~\ref{prop:leakinv}) on which the invariants are frozen to
    $10^{-14}$ while the certificate of Theorem~\ref{thm:leakage} still separates.

### `sm_nogo.tex`

    predicted, and it has no purchase on the leakage at all. On the fibers of Table~\ref{tab:fiber} it is
    \emph{identical} at the two endpoints---$2$ and $2$ at $K=3$, $1$ and $1$ at $K=4$, at equal raw
    support---while $\widehat\Lambda/\eta$ moves by a factor $6$--$7$; in a control at $K=3$ it jumps from
    $2$ to $16$ with $\widehat\Lambda/\eta=2.336$ and relative $L_1=0.238$ unchanged to three figures; and
    across a subfamily of $58$ points on which $\chi$ and the raw support are \emph{exactly} equal,
    $\widehat\Lambda/\eta$ still varies by $2.8$ and the measured error by $8.3$.

### `sec_6_app.tex`

    Table~\ref{tab:fibersweep} continues the fiber exhibit of Table~\ref{tab:fiber} into a continuous sweep of
    the rotation angle at $K=3$. It stops at $s=0.654$: at the maximal-pairing angle $s=\pi/4$ the Born
    distribution of $\phi$ is exactly flat over $2^{K}$ determinants, so the top-$|\mathcal S|$ rule does not
    select a subspace and no contrast measured there is well defined. Across every row the captured weight is
    identical at the two endpoints to $2.2\times10^{-16}$, so the entire difference sits in the leakage term,
    and the certificate at the rotated endpoint is the trivial bound $1+w_S$ in every row.
    
    \begin{table*}[tp]\centering\small
    \caption{\textbf{Angle sweep of the fiber exhibit}, $K=3$, $\delta=0.01$, $\eta=0.10$,
    $|\mathcal S|=4$. Certificate columns as in Table~\ref{tab:fiber}; $w_S$ is identical at the two
    endpoints, and bold marks a certificate below the trivial bound $1+w_S$.}
    \label{tab:fibersweep}
    \setlength{\tabcolsep}{4pt}
    \begin{tabular}{@{}lccccc@{}}
    \toprule
    $s$ & $w_S$ (both endpoints)
     & $\widehat\Lambda/\eta$ al. & $\widehat\Lambda/\eta$ rot.
     & cert.\ al. & cert.\ rot. \\
    \midrule
    $0.131$ & $0.982963$ & $0.328$ & $1.610$ & $\mathbf{0.913}$ & $1.983$ (trivial) \\
    $0.262$ & $0.933013$ & $0.317$ & $3.010$ & $\mathbf{1.131}$ & $1.933$ (trivial) \\
    $0.393$ & $0.853553$ & $0.299$ & $4.124$ & $\mathbf{1.312}$ & $1.854$ (trivial) \\
    $0.524$ & $0.750000$ & $0.277$ & $4.849$ & $\mathbf{1.450}$ & $1.750$ (trivial) \\
    $0.654$ & $0.629410$ & $0.252$ & $5.139$ & $\mathbf{1.543}$ & $1.629$ (trivial) \\
    \bottomrule
    \end{tabular}
    \end{table*}


## Withdrawn the same day (third decision, the systematic sweep): main text and its SI echoes

- **"0.28 on the open chain"** (unit-weight error of the open-chain support $S^\ast$, $L=6$,
  $\eta=0.18t$). It was in Sec. I, Sec. IV, Sec. S3, the App. on moments and Sec. S6. No deposited file
  holds it: `support_witness.json` covers the periodic rings only, and `moments_table.json` stores
  moments, not $L_1$ errors. Withdrawn; the sentences now quote the deposited ring values, 0.35–0.43.
  Verbatim before: "($0.28$ on the open chain of Sec.~\ref{sec:moments}, a different lattice)";
  "---$0.28$ on the open chain of Table~\ref{tab:moments} at $\eta=0.18\,t$, and up to $0.43$ on the
  periodic rings"; "a measured range of $0.28$--$0.43$ at unit captured weight".
- **Table II, the round-off column**: $1.8\times10^{-15}$, $3.2\times10^{-15}$ and $7.5\times10^{-15}$,
  "the largest value over six equivalent implementations", and in the appendix "$\le3.2\times10^{-15}$
  for $K=0,1$ and $\le7.5\times10^{-15}$ for $K=2$". The six routes' outputs were never deposited.
  Replaced by the deposited generator's values (`data/moments_table.json`): 0, $2.1\times10^{-16}$, 0.
- **"against $4\times10^{-16}$ for the Krylov-support subspace"** (Sec. IV): not deposited. The
  deposited run reproduces that zeroth moment exactly, which is what the sentence now says.

## Pass 7 (number sweep): Table S-fracgrid rewritten from the deposited eta sweep

The sector-fraction grid `tab:fracgrid` (SM, the three prices) printed an early coarse-scan replication (depth 150, smallest scanned fraction) whose output was never deposited; its numbers survive only in a working note, and no output file was found (searched 2026-11-08).  The deposit holds the same measurement re-run with bisection, `data/eta_sweep.json` (src/eta_sweep.py, depth 260), whose eta=0.15t column reproduces the published points.  The table and every sentence that quotes it were rewritten from that file by a one-off script of the JCP package (its numbers are computed by the deposited `src/jcp/fracgrid_numbers.py`); the qualitative content is unchanged (the fraction exponent moves by a factor ~2 over the eta decade, the support exponent stays near one, and the two shifts are equal).  Withdrawn text, verbatim:

### `sm_prices.tex`

    applies to the object fitted: $1.082$ is the fit to the published points, while the $1.090$ in the
    bottom block of Table~\ref{tab:fracgrid} is the fit to that table's own coarser replication grid at the
    same resolution. They are two fits to two datasets and the difference is theirs, not a rounding.

### `sm_prices.tex`

    across the decade $\eta=0.05\to0.50\,t$ both shift by $0.148$ in absolute value

### `sm_prices.tex`

    $0.148$ is three quarters of the quantity itself: $e^{-0.13\,L}$ becomes $e^{-0.27\,L}$, and the

### `sm_prices.tex`

    $0.148$ carries it from $1.165$ to $1.017$---a seventh of its own value, $13\%$ measured against the
    larger end and $15\%$ against the smaller---

### `sm_prices.tex`

    snapshots; \texttt{scaling\_lanczos\_mf.py}), same threshold, same depth $n_\ell=150$; only $\eta$
    changes. The scan over subspace size in this replication is coarser than the bisection that produced
    the published points, and returns the smallest scanned fraction meeting the threshold, so the
    $\eta=0.15\,t$ column is an upper bound on them: against the published $0.9167$, $0.8200$, $0.5561$,
    $0.3625$, $0.1700$ it reads $0.917$, $0.820$, $0.560$, $0.380$, $0.180$, agreeing at the two smallest
    sizes and standing $0.7\%$, $4.8\%$ and $5.9\%$ high at the three largest. The same granularity is why
    the fraction exponent of this column, $-0.2012$, is not identical to the $-0.2093$ obtained by fitting
    the published points themselves over $L=4$--$12$. Bottom block: least-squares exponents.

### `sm_prices.tex`

    criticizes exactly that elsewhere in this section. The fraction exponent moves by a factor $2.2$ over
    the decade, from $-0.1266$ to $-0.2743$, while the support exponent moves from $1.165$ to $1.017$; in
    absolute value both move by $0.148$, which is the whole content of the
    comparison.}

### `sm_prices.tex`

    $4$  & $0.917$ & $0.917$ & $0.917$ & $0.917$ & $0.917$ & $0.917$ & $0.917$ & $0.875$\\
    $6$  & $0.900$ & $0.890$ & $0.870$ & $0.820$ & $0.800$ & $0.730$ & $0.680$ & $0.630$\\
    $8$  & $0.740$ & $0.680$ & $0.640$ & $0.560$ & $0.540$ & $0.480$ & $0.420$ & $0.380$\\
    $10$ & $0.520$ & $0.480$ & $0.440$ & $0.380$ & $0.340$ & $0.260$ & $0.220$ & $0.200$\\
    $12$ & $0.340$ & $0.280$ & $0.240$ & $0.180$ & $0.160$ & $0.140$ & $0.120$ & $0.100$\\

### `sm_prices.tex`

    fraction exponent, $L{=}4$--$12$ & $-0.1266$ & --- & $-0.1681$ & $-0.2012$ & --- & $-0.2395$ & --- & $-0.2743$\\
    $R^{2}$                          & $0.904$   & --- & $0.918$   & $0.929$   & --- & $0.970$   & --- & $0.982$\\
    support exponent, $L{=}4$--$12$  & $1.165$   & --- & $1.123$   & $1.090$   & --- & $1.052$   & --- & $1.017$\\
    $R^{2}$                          & $0.999$   & --- & $0.999$   & $0.998$   & --- & $0.999$   & --- & $0.999$\\
    support exponent, $L{=}6$--$12$  & $1.136$   & --- & $1.088$   & $1.053$   & --- & $1.022$   & --- & $0.992$\\

### `sm_prices.tex`

    is claimed here. Over the decade of $\eta$ the fraction exponent moves by a factor $2.2$, from
    $-0.1266$ at $\eta=0.05\,t$ to $-0.2743$ at $\eta=0.50\,t$. At $\eta=0.50\,t$ the fraction at $L=12$
    is already $0.10$, a quarter above the $L=14$ value at the resolution used throughout this work and
    reached two sizes earlier:

### `sec_7.tex`

    support exponents differ by the constant $b_{\dim}=1.291$ and move together, by $0.148$ over the decade
    $\eta=0.05\to0.50\,t$ (Table~\ref{tab:fracgrid}): a seventh of the support exponent, three quarters of
    the fraction exponent.

### `sec_7.tex`

    $2.2$, from $-0.1266$ to $-0.2743$ (Table~\ref{tab:fracgrid}), and even at fixed $\eta$ the local

### `figs/caption_fig12.tex`

    $\lvert\mathcal{S}\rvert$ moves from $1.165$ to $1.017$ and that of the fraction in panel~(b) by a
    factor $2.2$, from $-0.1266$ to $-0.2743$ (Table~\ref{tab:fracgrid}, $L=4$--$12$). In absolute value
    both move by the same $0.148$:

### `sec_3_1.tex`

    at which the resource scan of Sec.~\ref{sec:prices} operates (Fig.~\ref{fig:scaling}; the coarser
    replication grid of Table~\ref{tab:fracgrid} reads $0.180$ at $L=12$ for the same operating point).


## Pass 7 (number sweep): the shot-frontier table and its quotes, withdrawn

Table `tab:shotfrontier` (T_marg, the reciprocal of the Born probability of the marginal determinant of each cut, and the factors 1.4-1.6 / 7.5-10 / 6.7-21 / 20-35, the further factors 241 and 634, and the marginal Born probability 1-2e-9 at the L=12 threshold) has no deposited or stored source: ckpt/ holds only the L=14 ground state and order, no file of data/ stores the cut weights at L=6-12 or T_marg, and no script computes them (searched 2026-11-08).  Withdrawn by a script of the JCP package (not deposited); the sentences now state only what the deposit supports (the certificate does not certify the operating points; the measured finite-shot budget of Table S-shots exceeds the support it buys).  Verbatim text removed:

### `sm_frontier.tex`

    \begin{table*}[t]\centering\small
    \caption{\textbf{The threshold in the only currency a sampler has.} $T_{\rm marg}$ is the shots per snapshot at
    which the marginal determinant of the cut appears once in expectation---the reciprocal of its Born
    probability in the accumulated distribution of the $19$ snapshots, from the protocol's own cut
    weights; covering all of $S$ costs a further factor $\log\lvert S\rvert$. \emph{It is a different
    quantity from the measured budgets $T$ of Table~\ref{tab:shots}}, which are totals over the $19$
    snapshots at which the finite-shot reconstruction reaches relative $L_1<0.05$; the two carry different
    names for that reason, they are not interchangeable, and only ratios within this table are used.
    Column~2 uses the grid point nearest the published fraction: $0.80$ against a published $0.82$ at
    $L=6$, and $0.40$ against $0.36$ at $L=10$, where the published fraction is the worse of the two.
    Factors are from unrounded budgets and need not reproduce from the two-figure entries. The marginal Born probability at
    the $L=12$ threshold is $2.2\times10^{-9}$ at $\mathrm{FR}=0.40$, $1.2\times10^{-9}$ at
    $0.45$.}\label{tab:shotfrontier}
    \setlength{\tabcolsep}{6pt}
    \begin{tabular}{lccc}
    \toprule
    $L$ & published fraction $\to T_{\rm marg}$ & $\mathrm{FR}^{\ast}\to T_{\rm marg}$ & factor\\
    \midrule
    $6$  & $0.80\to2.2\times10^{4}$ & $[0.90,0.92]\to3.1$--$3.6\times10^{4}$ & $1.4$--$1.6$\\
    $8$  & $0.56\to2.1\times10^{5}$ & $[0.75,0.80]\to1.6$--$2.1\times10^{6}$ & $7.5$--$10$\\
    $10$ & $0.40\to5.2\times10^{6}$ & $[0.60,0.70]\to3.4\times10^{7}$--$1.1\times10^{8}$ & $6.7$--$21$\\
    $12$ & $0.18\to2.3\times10^{7}$ & $[0.40,0.45]\to4.7$--$8.1\times10^{8}$ & $20$--$35$\\
    \bottomrule
    \end{tabular}
    \end{table*}

### `sm_frontier.tex`

    In shots the distance between those two columns is a budget and not a rhetorical gap
    (Table~\ref{tab:shotfrontier}): reaching the non-vacuous band multiplies the shots per snapshot needed
    for the marginal determinant by $7.5$--$10$ at $L=8$ and by $20$--$35$ at $L=12$, and only by
    $1.4$--$1.6$ at $L=6$, where the sector is small. Driving the certified value down to $0.10$ costs more
    again---a further factor $241$ at $L=8$ (fraction $0.95$) and $634$ at $L=10$ (fraction $0.90$), the
    corresponding $L=12$ fraction resting only on evaluations that failed the depth test, so that no budget
    is quoted for it. This section reports the measured regime of validity: four sizes, four resolutions, $204$ evaluations, every
    bracket delimited by converged rows, and a shot budget attached to each.

### `sm_frontier.tex`

    and not finite-shot draws, and every statement below inherits that limit: the marginal determinant at
    the $L=12$ threshold carries Born probability $1$--$2\times10^{-9}$, so at the budgets of
    Table~\ref{tab:shots} that subspace barely exists as a sampled object.

### `sec_5_body.tex`

    accuracy Theorem~\ref{thm:leakage} does not certify the operating points of Fig.~\ref{fig:scaling}, and
    the regime where it would costs more shots than the method saves, since reaching the non-vacuous band
    multiplies the shots per snapshot needed for the marginal determinant by $7.5$--$10$ at $L=8$ and by
    $20$--$35$ at $L=12$ (Table~\ref{tab:shotfrontier}).

### `sec_10.tex`

    scan, and the regime in which it would certify them costs more shots than the method saves
    (Sec.~\ref{sec:prices}). Certification is costly: moving from the operating fraction to the non-vacuous band multiplies the shots per
    snapshot needed for the marginal configuration by $7.5$ to $10$ at $L=8$ and by $20$ to $35$ at
    $L=12$. At the acceptance

### `sec_9.tex`

    and moving from the operating fraction to the non-vacuous band
    multiplies the shots per snapshot needed for the marginal configuration by $7.5$ to $35$
    (Sec.~\ref{sec:prices}).


## Pass 7 (number sweep): the L=8 E0-shift example, withdrawn

The hypothesis-(H1) example (delta = 6.7e-4 t from 60% of the N-electron sector at L=8; +0.0024 on a bound of 1.2036, +0.20%; measured error x16; the additive term verified to 2e-5 in 48 cases; the L=6 non-monotone example 2.7e-2 -> 1.2e-2) has no deposited output at L=8: `src/frontier/z_e0.py` deposits only `data/c3_frontier/adversarial/z_e0_L6.json`.  Withdrawn by a script of the JCP package (not deposited); the analytic statement of the additive term stays.  Verbatim text removed:

### `sec_3_body.tex`

    does not control; at $L=8$, with $E_0$ from $60\%$ of the $N$ sector, that term moves the bound by
    $+0.20\%$ while the measured error grows by a factor $16$ (Sec.~\ref{app:hyp}).

### `sec_10.tex`

    control. For the energy shift of $6.7\times10^{-4}\,t$ obtained from $60\%$ of the $N$-electron
    sector at $L=8$, that term moves the bound by $+0.20\%$, while the measured error grows by a factor of
    $16$ (Sec.~\ref{app:hyp}).

### `sec_5_app.tex`

    an additive term verified to $2\times10^{-5}$ in all $48$ cases tested and \emph{which the certificate
    does not control}: the bound is guaranteed only once that term is added. A ground state obtained by
    Rayleigh--Ritz on $60\%$ of the $N$-electron sector at $L=8$ gives $\delta=6.7\times10^{-4}\,t$, which
    adds $0.0024$ to a bound of $1.2036$---$+0.20\%$---while multiplying the \emph{measured} error by
    $16$; and the dependence is not even monotone, a shift of $10^{-2}\,t$ lowering the measured error at
    $L=6$, fraction $0.85$, $\eta=0.18\,t$ from $2.7\times10^{-2}$ to $1.2\times10^{-2}$.


## Pass 7 (number sweep): Proposition `prop:leakinv` as it stood before the fiber numbers were withdrawn

Recorded here verbatim (from `paper/sec_6_body.tex` as saved before the number sweep of 2026-11-08, lines 148-175) because the rewrite of 2026-11-08 was made by hand: the title and the closing sentences quoted the undeposited fiber exhibit (0.237/2.326, 0.323/2.336, 0.391/2.346, tau=0.795, 0.856 against 1.962).  The analytic content (invariance along the orbit; the invariant bound vacuous by at least W/2eta; an invariant bound must exceed the fiber maximum) is kept in the new text.

    \begin{proposition}[Fibers on which every orbital-rotation invariant gives a vacuous leakage certificate]
    \label{prop:leakinv}
    Let $\mathcal U_U$ be the Gaussian unitary implementing $U\in\mathrm U(n_{\rm o})$, and let the same
    physical problem written in the rotated orbital basis be
    $(H_U,\Psi_U,\phi_U)=(\mathcal U_U^{\dagger}H\mathcal U_U,\,\mathcal U_U^{\dagger}\Psi,\,
    \mathcal U_U^{\dagger}\phi)$, the determinant basis---in which the device measures, and of which $S$ is a
    coordinate subspace---held fixed. Along this orbit the spectrum of $H$, $E_0$, $\lVert\phi\rVert^{2}$,
    every Lehmann weight and hence $A(\omega)$ pointwise are constant, as is every orbital-rotation invariant
    of the state, in particular $\mathcal F_1=4\,\mathrm{tr}[\gamma(1-\gamma)]$ and all
    $\mathcal F_k$. $\Lambda_S(\eta)/\eta$ is not. An upper bound on $\Lambda_S(\eta)/\eta$ is available from
    invariants alone---combining $\Lambda_S\le\lVert P_SHQ_S\rVert\,\lVert\phi_S\rVert$ with
    $\lVert Q_SHP_S\rVert=\lVert Q_S(H-c)P_S\rVert\le\min_c\lVert H-c\rVert$ gives
    $\Lambda_S(\eta)/\eta\le(W/2)\lVert\phi\rVert\sqrt{w_S}/\eta$, $W$ the spectral width of the $(N\pm1)$
    sector---but it is vacuous by a factor of at least $W/2\eta$, which grows extensively with system size.
    That no invariant bound can be useful on these fibers is seen at fixed invariants: any
    $B(\text{invariants},|\mathcal S|,w_S)$ valid for $\Lambda_S(\eta)/\eta$ is constant on a fiber of fixed
    arguments and must therefore exceed the fiber maximum. In the geminal family $\mathrm{APSG}(\pi/4)$ on
    $2K$ orbitals with $H=\varepsilon\hat N-g\sum_iA^{\dagger}_iA_i+\delta h$ we exhibit fibers
    ($K=2,3,4$; $\delta=0.01$; $\eta=0.10$; $|\mathcal S|=2^{K-1}$; spectrum, Lehmann weights and all
    $\mathcal F_k$ matched to $10^{-14}$, $w_S$ to $2.2\times10^{-16}$) whose endpoints carry
    $\widehat\Lambda_S(\eta)/\eta=0.237$ and $2.326$, $0.323$ and $2.336$, $0.391$ and $2.346$. Since the
    larger value exceeds $\tau(w_S)=0.795$ of Eq.~\eqref{eq:tau}, above which the leakage branch of
    Theorem~\ref{thm:leakage} no longer improves on the trivial bound
    $\lVert\phi\rVert^{2}(1+w_S)$, every valid invariant bound is vacuous on these fibers---while the
    subspace-aware certificate is non-vacuous at the aligned endpoint ($0.856$ against a trivial $1.962$).
    The statement is one of non-determination and non-vacuity, not of non-existence: invariants do not fail
    to bound the leakage, they fail to bound it usefully.
    \end{proposition}


## Pass 7 (number sweep): the Cauchy-Schwarz margin, corrected

Sec. V F and S4 printed that the total margin for tightening of the Cauchy-Schwarz step is "a factor
$1.05$ typically and $3.29$ at worst".  Re-derived from `data/c3_frontier/ladder_L8.json` (ratio of the
rung after the step to the rung before it, R2/R1, over the 84 ladder rows of the FR cases at L=8) the
maximum is 3.29 but the median is 1.10; no statistic of the deposited ladders gives 1.05 (all three
sizes pooled: 1.06; L=6: 1.06; at the 15 threshold cells: 1.17).  The sentences now print "a factor
$1.10$ at the median and $3.29$ at worst", and `src/verify_jcp.py` (check_frontier_anatomy)
re-derives both.  Verbatim before: "whose total margin for tightening is a factor $1.05$ typically and
$3.29$ at worst" (sec_5_body.tex); "its total margin for tightening is a factor $1.05$ typically and
$3.29$ at the worst point measured" (sm_frontier.tex).

## Pass 8: three printed numbers corrected to the deposit

Recomputed by a script of the JCP package (not deposited) from `data/c3_frontier/C3_grid.json` and `published_fraction/pt_L*_PUB.json`: the fraction at which the reconstruction reaches rel-L1<0.05 at L=8, eta=0.18t is 0.5250 (was printed "~0.52"); the median of the sixteen true errors at the threshold is 2.75e-03 (was "2.8e-3"); their smallest ratio to the reference drift is 46.5 (was "at least 47").  Old and new text:

- `sec_5_body.tex`: "$\mathrm{FR}\approx0.52$, the certificate requires" -> "$\mathrm{FR}\approx0.525$, the certificate requires"
- `sm_frontier.tex`: "is reached at $\mathrm{FR}\approx0.52$, while" -> "is reached at $\mathrm{FR}\approx0.525$, while"
- `sec_5_body.tex`: "median $2.8\times10^{-3}$: a spread" -> "median $2.7\times10^{-3}$: a spread"
- `sm_frontier.tex`: "median $2.8\times10^{-3}$: a spread" -> "median $2.7\times10^{-3}$: a spread"
- `sec_5_body.tex`: "at least $47$ times the drift" -> "at least $46$ times the drift"
- `sm_frontier.tex`: "at least $47$ times the drift" -> "at least $46$ times the drift"

## Pass 8: SI numbers corrected to the deposit or rewritten without them

By a script of the JCP package (not deposited).  Old and new text:

- `figs/captions.tex`: "$p\approx0.56$" -> "$p\approx0.57$"
- `sec_7.tex`: "($95\%$ interval $[-0.03,0.19]$) and $0.054$ ($[-0.066,0.174]$)" -> "($95\%$ interval $[-0.02,0.19]$) and $0.054$ ($[-0.065,0.172]$)"
- `figs/caption_fig12.tex`: "$0.054$ per site ($[-0.066,0.174]$)" -> "$0.054$ per site ($[-0.065,0.172]$)"
- `figs/caption_fig12.tex`: "$0.08$ ($[-0.03,0.19]$)" -> "$0.08$ ($[-0.02,0.19]$)"
- `sec_5_app.tex`: "standing at a factor $47$" -> "standing at a factor $46$"
- `sec_5_app.tex`: "$30.2\%$ of $S$ at $L=6$, $42.0\%$ at $L=8$ and ${\approx}39.6\%$ at $L=10$." -> "$30.2\%$ of $S$ at $L=6$ and $42.0\%$ at $L=8$ (\texttt{ties.json})."
- `sec_5_app.tex`: "move the bound by at most $6.6\times10^{-4}$ relative" -> "move the bound by at most $7.8\times10^{-4}$ relative"
- `sm_leakage.tex`: "$6.6\times10^{-4}$ by which ties move the \emph{bound}" -> "$7.8\times10^{-4}$ by which ties move the \emph{bound}"
- `sec_5_app.tex`: "$P_{\mathcal K}\phi=\phi_S$ to $1.1\times10^{-14}$" -> "$P_{\mathcal K}\phi=\phi_S$ to $1.5\times10^{-15}$"
- `sm_leakage.tex`: "And the error at fraction $0.85$ has two independent determinations that differ by a factor $2.3$, so only its order of magnitude is quoted." -> "And the error at fraction $0.85$ is quoted only to its order of magnitude."
- `sm_leakage.tex`: "moves them by about $1\%$ at $\eta=0.05\,t$ and by $2.6$ to $3.3\%$ at $\eta=0.18\,t$, measured in Sec.~\ref{app:hyp}" -> "moves them by a few per cent, measured in Sec.~\ref{app:hyp}"

## Pass 8: SI numbers with no stored data, withdrawn

By a script of the JCP package (not deposited).  No file of data/ or ckpt/ holds these values, and none was found in the author's working folders.  Verbatim text removed, and what replaced it:

### `sec_4_app.tex`

    One run of one of the six routes does report a
    nonzero residue in that row; the same run reports a cyclic-subspace dimension of $301$ inside a sector
    of dimension $300$ and emits floating-point overflow while forming the high powers. It is discarded,
    and the corrected run is the one quoted.

  now: "One run of one of the six routes reported a nonzero residue in that row while emitting floating-point overflow in the high powers; it is discarded, and the corrected run is the one quoted."

### `sec_4_app.tex`

    $\mu_j=\braket{\phi|(H-E_0)^{j}|\phi}$ run from $5.0\times10^{-1}$ at $j=0$ to $1.6\times10^{5}$ at
    $j=6$, so an absolute residue

  now: "$\mu_j=\braket{\phi|(H-E_0)^{j}|\phi}$ grow steeply with $j$ from $\mu_0=\lVert\phi\rVert^{2}$, so an absolute residue"

### `sec_4_app.tex`

    the relative residue below
    order $2K{+}2$ stays inside $[1.1,1.8]\times10^{-15}$ at $K=0$ and $[1.0,3.2]\times10^{-15}$ at $K=1$,
    a spread under a factor $3.3$, while

  now: "the relative residue below order $2K{+}2$ stays at round-off, while"

### `sec_4_app.tex`

    the zeroth moment is wrong by between $1.86\times10^{-2}$ and $4.18\times10^{-1}$
    (median $1.67\times10^{-1}$), i.e.\ by $3.7\%$ to $84\%$ of $\lVert\phi\rVert^{2}$ (median $33\%$),
    with every one of
    the $300$ draws failing; the Krylov-support subspace of identical size is exact to
    $4.4\times10^{-16}$. At $|\mathcal S|=280$ the random draws fail by $4.8\times10^{-5}$ to
    $2.6\times10^{-1}$ (median $1.12\times10^{-2}$), again without exception, against
    $5.0\times10^{-16}$.

  now: "the zeroth moment is wrong in every one of the $300$ draws (\texttt{moments\_table.json}), while the Krylov-support subspace of identical size reproduces it exactly (Table~\ref{tab:moments})."

### `sec_4_app.tex`

    
    ($K=0$: relative discrepancy $4.4\times10^{-16}$ through order $1$, $4.6\times10^{-3}$ at order $2$);

  now: ";"

### `sec_5_app.tex`

    two has been measured at a single point at $L=8$---relative $L_1$ of $1.2454\times10^{-3}$ against
    $1.2464\times10^{-3}$---and never bounded

  now: "two has been measured at a single point only and never bounded"

### `sec_5_app.tex`

     The depth uncertainty quoted in Sec.~\ref{sec:frontier} is the
    full recomputation of $\mathrm{FR}^{\ast}$ with $L=12$ capped at $n_\ell=300$: it moves by $0.0024$,
    $0.0021$, $0.0012$ and $0.0000$ at $\eta=0.10$, $0.15$, $0.18$ and $0.25\,t$.

  now: ""

### `sec_6_app.tex`

    ($\theta=\pi/6$ gives
    $\mathcal F_1=3$, $\mathcal F_2=3.75$ per pair, both still large on the same family)

  now: "(at $\theta=\pi/6$, $\mathcal F_1=3$ per pair, and the higher orders differ from it, all still large on the same family)"

### `sm_prices.tex`

    Second, that reference drifts with depth, and the drift is resolution-dependent:
    between $n_\ell=100$ and $n_\ell=200$ at $L=12$ the relative $L_1$ change of the reference
    \emph{alone} is $11.1\%$ at $\eta=0.05\,t$, $3.0\%$ at $\eta=0.10\,t$, $0.93\%$ at $\eta=0.15\,t$ and
    $0.05\%$ at $\eta=0.25\,t$. At the depth used here the protocol is self-consistent only for
    $\eta\gtrsim0.15\,t$---which is exactly where the published column sits, at the margin. Fig.~\ref{fig:spin}, a
    different channel drawn at $n_\ell=200$, sits \emph{below} that margin, at $\eta=0.10\,t$, which is why
    its caption quotes its lineshape with its depth and not as exact.

  now: "Second, that reference itself drifts with depth, by more at finer resolution, and the size of that drift is not quantified in the deposit. Fig.~\ref{fig:spin}, a different channel drawn at $n_\ell=200$ and $\eta=0.10\,t$, therefore quotes its lineshape with its depth and not as exact."

### `sec_7.tex`

    The reference is itself a
    depth-$150$ Lanczos computation whose drift between depths $100$ and $200$ at $L=12$ is $0.93\%$ at
    $\eta=0.15\,t$ and $3.0\%$ at $\eta=0.10\,t$, so the fraction

  now: "The reference is itself a depth-$150$ Lanczos computation whose drift with depth is not quantified here, so the fraction"

### `sm_physics.tex`

    It is not,
    however, converged: the same depth moves the reference by $0.93\%$ at $\eta=0.15\,t$ and by $3.0\%$ at
    the $\eta=0.10\,t$ of Fig.~\ref{fig:spin} when it is taken from $n_\ell=100$ to $n_\ell=200$
    (Sec.~\ref{sec:prices}).

  now: "It is not shown to be converged in depth (Sec.~\ref{sec:prices})."

### `carried/fig_spin.tex`

    at this
    $\eta$ a reference of that kind is not converged in general (the single-particle one moves by $3.0\%$
    between $n_\ell=100$ and $200$, Sec.~\ref{sm:prices}), so the lineshape

  now: "at this $\eta$ a reference of that kind is not shown to be converged in depth (Sec.~\ref{sm:prices}), so the lineshape"

### `sm_physics.tex`

    At $L=12$ the largest spacing between consecutive
    excited-state energies in the window is ${\approx}0.48\,t$, a factor $2.7$ \emph{above} the
    $\eta=0.18\,t$ used for $A(k,\omega)$ and $4.8$ above the $\eta=0.10\,t$ used for the spin channel.
    The rule is not satisfied here. Satisfying it would demand $\eta>0.48\,t$. That is only a tenth of
    the Mott scale, so the Hubbard bands would survive; but it is $0.96\,J$ and $2.9$ times the $L=12$
    spin gap itself, so the entire low-frequency structure of Fig.~\ref{fig:spin} would merge into one
    feature. The number the rule
    demands is set by the spacing, $0.48\,t$; we take it from the spacing and not from a rounded-up
    interval.

  now: "The spacing between consecutive excited-state energies at $L=12$ is not tabulated in the deposit, so this work does not claim that the rule is satisfied at the resolutions it uses."

