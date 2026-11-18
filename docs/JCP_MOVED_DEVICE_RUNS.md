# Text that left the JCP build on 2026-10-28 (verbatim, from the release sources)

> **Deposit 3.5.0.** This is the record kept by the package of the manuscript prepared for The Journal of Chemical Physics (written 2026-10-28 and later), deposited because `src/check_trim.py` points to it, through its `@doc` proofs, for each fact it records as deliberately withdrawn. In it, "release" means deposit 3.4.1 (tag v3.4.1, the arXiv-v4 text); references to the author's working notes and folders, which are not in the deposit, were removed, and nothing else changed.

Why: on the author's instruction (JCP package of 2026-10-28), the two device runs leave the
abstract and the main text and are stated in ONE line of Sec. S8 ("Device runs."), because
they add no evidence to the work: the ibm_fez run retained the whole (N+1) sector from
bitstrings read in reversed qubit order, so its A(omega) is the exact one by coverage, and the
ibm_marrakesh run has no random-configuration control.  Nothing is deleted from the deposit:
data/heron_spectral.json, data/hw_lucj_n2_result.json, data/hw_bitorder_check.json,
src/hw_bitorder_check.py, notebooks/, paper/carried/fig_heron.tex and
paper/figs/fig_hardware_hero_frag.tex all stay.  The abstract was rewritten results-first the
same day; its old text is below too.  This file is the "@doc" proof that the JCP
registrations of src/check_trim.py point to for each withdrawn fact.
Deposit 3.4.1 (tag v3.4.1, the arXiv-v4 text) keeps the full account.

## main.tex

### release lines 2-22 (replaced)

```tex
Rayleigh--Ritz on a subspace from sample-based quantum diagonalization yields the spectral function
of a probe state; we ask what bounds its error. The probe's captured Born weight $w$ does not do so
alone: an analytic counterexample forces any bound indexed on it to its trivial value. On Hubbard
rings at broadening $\eta=0.18\,t$, the probe's whole support ($w=1$) still gives a relative error
of $0.35$--$0.43$: Rayleigh--Ritz displaces the poles a subspace keeps. The error is bounded instead
by the missed weight and the Hamiltonian coupling across the subspace boundary. Given the exact
ground state, we prove for any orthogonal projector that the relative $L_1$ error lies between $1-w$
and $\min\{1+w,(1+\sqrt{w})(\sqrt{1-w}+\widehat\Lambda(\eta)/\eta)\}$, where $\widehat\Lambda$,
that coupling at resolution $\eta$, is computed from the Ritz pairs and the retained probe. We test
it on $L$-site rings with classical Born draws of an exactly time-evolved probe, or their
infinite-shot ranking; our device runs are execution records only. On the ranked subspaces the
bound is vacuous at every operating fraction of the resource scan, at all five sizes up to a sector
of $10\,306\,296$ determinants: its leakage branch exceeds the trivial bound $1+w$ by factors of
$1.66$ to $8.61$. Only well above those fractions is it informative, and there it is calibrated: at
$L=6$--$12$ it first beats the trivial bound at a true error of $(1.2$--$7.3)\times10^{-3}$, though
an unproven fit in $1-w$ is tighter on 11 of 14 splits. For the Born-ranked $A(k,\omega)$ at $85\%$
of its sector, the bound is $350$--$900$ times the true error. In the determinant basis, zero
leakage at full weight needs every symmetry-allowed determinant: the probe's measured Krylov
support at $L=6$--$14$. Plane-wave orbitals reduce a momentum probe's support to one block, $1/L$ of
an exponentially growing sector. At full weight a second-order certificate follows; whether
one exists below full weight, where every operating fraction lies, is open.
```

now reads:

```tex
Sample-based quantum diagonalization returns sampled configurations, and Rayleigh--Ritz in their
span yields the spectral function of a probe state but, by itself, no bound on its error. The
probe's captured Born weight $w$, the natural diagnostic, does not supply one alone: an analytic
counterexample forces any bound indexed on it to its trivial value, and on Hubbard rings the probe's
own support, where $w=1$, still carries an error of order one. Given the exact ground state, we
prove for any orthogonal projector that the relative $L_1$ error lies between $1-w$ and
$\min\{1+w,(1+\sqrt{w})(\sqrt{1-w}+\widehat\Lambda(\eta)/\eta)\}$, where $\widehat\Lambda$, the
Hamiltonian coupling across the subspace boundary at resolution $\eta$, is computed from the Ritz
pairs and the retained probe. We test the bound on $L$-site rings with classical Born draws of an
exactly time-evolved probe or their infinite-shot ranking. On the ranked subspaces, at every
operating fraction of the resource scan (the sector fraction needed for the accuracy target) and at
all five sizes up to a sector of $10\,306\,296$ determinants, its leakage branch exceeds the trivial
bound $1+w$ by factors of $1.66$ to $8.61$. The bound is informative only well above those
fractions; for the Born-ranked $A(k,\omega)$ at $85\%$ of its sector, it is $350$--$900$ times the
true error. In the determinant basis, zero leakage at full weight needs every symmetry-allowed
determinant, the probe's Krylov support measured at the same five sizes. At full weight a
second-order certificate follows; whether one exists below full weight, where every operating
fraction lies, is open.
```

## sec_1.tex

### release lines 101-101 (replaced)

```tex
\section{The operational question}
```

now reads:

```tex
\section{Introduction}
```

### release lines 210-210 (replaced)

```tex
been measured on devices far larger than the $12$-qubit spectral run of this work. They include a
```

now reads:

```tex
been measured on quantum devices. They include a
```

### release lines 216-219 (replaced)

```tex
device and compared with neutron data for copper sulfate~\cite{quantinuum_dsf2026}. The two device
runs of this work are execution records only: the spectral one retained its entire sector, so its
agreement with exact diagonalization is exact by coverage rather than a test of the method
(Sec.~\ref{sec:scope:hw}). \emph{(iv)}~Computable Gaussian-invariant quantifiers of fermionic
```

now reads:

```tex
device and compared with neutron data for copper sulfate~\cite{quantinuum_dsf2026}. \emph{(iv)}~Computable Gaussian-invariant quantifiers of fermionic
```

### release lines 242-242 (replaced)

```tex
the device runs, the molecular suite and the relation to classical methods. The systems are the half-filled
```

now reads:

```tex
the noise sweeps, the molecular suite and the relation to classical methods. The systems are the half-filled
```

### release lines 248-248 (replaced)

```tex
discussions are in the Supplemental Material~\cite{supplemental}, cited as Sec.~S$n$.
```

now reads:

```tex
discussions are in the supplementary material~\cite{supplemental}, cited as Sec.~S$n$.
```

## sec_2.tex

### release lines 102-108 (deleted)

```tex
\raggedright $A(\omega)$; hardware, \texttt{ibm\_fez}, $L{=}6$, $U/t{=}4$, $\eta{=}0.15\,t$; $\dim=300$
 & \raggedright Fig.~\ref{fig:heron}(a)
 & \raggedright post-selected from device bitstrings read in reversed bit order, so every retained shot is a
   device error (Sec.~\ref{sec:scope:hw}); $|\mathcal S|=300$ of $300$, the entire sector, so the reconstruction
   equals the exact $A(\omega)$ by coverage
 & \raggedright $3.5\times10^{5}$: $5\times10^{4}$ on each of $7$ circuits, $k{=}0$ (the $t{=}0$
   reference) to $k{=}6$\tabularnewline
```

## sec_9.tex

### release lines 154-154 (replaced)

```tex
\section{Hardware runs, classical reach, and scope}
```

now reads:

```tex
\section{Noise sweeps, classical reach, and scope}
```

### release lines 157-157 (replaced)

```tex
\subsection{The two hardware runs, and what ``exact by coverage'' means}
```

now reads:

```tex
\subsection{The noise sweeps and the molecular suite}
```

### release lines 160-187 (replaced)

```tex
Two device runs were made, and they are the only quantum hardware used in this work. The first, on the IBM~Heron processor \texttt{ibm\_fez}, sampled the $L=6$ periodic Hubbard ring at $U/t=4$
with $12$ qubits and $3.5\times10^{5}$ shots, and post-selection on particle number and $S_z$, with no
configuration recovery, retained all
$300$ determinants of the $(N{+}1)$ sector, so the reconstructed $A_{\rm hw}$ equals the exact $A$ \emph{by coverage}:
Rayleigh--Ritz on the whole sector is its diagonalization, and the vanishing error is an identity rather
than a measurement of device fidelity (Fig.~\ref{fig:heron}(a)). The circuits prepared a
computational-basis $(N{+}1)$ determinant rather than $\hat c^{\dagger}_p\lvert\Psi_0\rangle$ and applied
$R_{xx}R_{yy}$ hopping without Jordan--Wigner strings, so they are not the primitive of
Fig.~\ref{fig:circ}. The post-selection, moreover, read each bitstring in reversed qubit order, which
reflects the ring and exchanges the two spin species: under that reading none of the $3.5\times10^{5}$
shots of a noiseless simulation of the same seven circuits is retained, against all of them under the
correct reading, so every one of the $300$ retained determinants arose from a device error
(\texttt{src/hw\_bitorder\_check.py}). The run is therefore an execution record: the steps of the pipeline execute end to end on a
processor, but the run does not show that the device sampled the intended distribution, and it
establishes nothing about accuracy at complete or incomplete coverage. The second, stretched $\mathrm{N_2}$ at $24$ qubits on
\texttt{ibm\_marrakesh}, reaches $0.59\pm0.14$~mHa of the exact active-space energy after recovery,
below the $29.5$~mHa of the noiseless simulation of the same circuit (Fig.~\ref{fig:heron}(b)). The
record does not store the recovered subspace sizes on which that comparison turns, so it cannot say
whether device noise helps, although a study under a simple noise model reports the mechanism~\cite{vaquero2026}.
That study also finds that recovery from randomly generated configurations alone builds accurate CI
spaces, and no such control was run here. What is
controlled is recovery against no recovery on the same corrupted bitstrings, and there the gain grows
with the noise (Figs.~\ref{fig:noise} and~\ref{fig:noiserec}). Neither run is a spectral record in the sense of
Refs.~\cite{vilchez2025,granet2026akw,quantinuum_dsf2026}, whose spectra are determined by device
data; the full account, the job identifiers and the shot arithmetic are in
Sec.~\ref{sm:scope}.

\input{carried/fig_heron.tex}
```

now reads:

```tex
The noise results of this work are controlled classical simulations: under a synthetic per-qubit
bit-flip channel, recovery is compared with no recovery on the same corrupted bitstrings, and the gain
grows with the noise (Figs.~\ref{fig:noise} and~\ref{fig:noiserec}). Where the recovered error falls
below the noiseless value, that comparison is not used here: it is decided by the size of the
recovered subspace, which the sweeps do not record (Sec.~\ref{sm:scope}). A study under a simple
noise model reports the mechanism~\cite{vaquero2026}. That study also finds that recovery from
randomly generated configurations alone builds accurate CI spaces, and no such control was run here.
```

### release lines 191-191 (replaced)

```tex
\paragraph{The molecular suite.} On the seven Hamiltonians of the nineteen-molecule suite (Sec.~\ref{sec:scope:mol}, Table~\ref{tab:molscope}) whose sectors have dimension $20$ to $300$, the same primitive reproduces the orbital-projected $A(\omega)$ to a relative $L_1$ between $2.2\times10^{-6}$ and $1.0\times10^{-4}$ from $21\%$ to $50\%$ of the sector, a consistency check of the primitive outside the Hubbard family and not a scaling statement. The other twelve reconstructions are exact by construction, because their sampled subspaces, of at most twelve determinants, already contain the probe's Krylov space, so the suite cannot show that the reconstruction is efficient across chemistry.
```

now reads:

```tex
\paragraph{The molecular suite.} On the seven Hamiltonians of the nineteen-molecule suite (Sec.~\ref{sec:scope:mol}, Table~\ref{tab:molscope}) whose sectors have dimension $20$ to $300$, the sampling primitive of Sec.~\ref{sec:method} reproduces the orbital-projected $A(\omega)$ to a relative $L_1$ between $2.2\times10^{-6}$ and $1.0\times10^{-4}$ from $21\%$ to $50\%$ of the sector, a consistency check of the primitive outside the Hubbard family and not a scaling statement. The other twelve reconstructions are exact by construction, because their sampled subspaces, of at most twelve determinants, already contain the probe's Krylov space, so the suite cannot show that the reconstruction is efficient across chemistry.
```

### release lines 408-417 (replaced)

```tex
The quantum-hardware results are deposited as the arrays that enter the figures, together with the
backend name, the shot count and the IBM~Quantum job identifier for each of the two runs
(\texttt{data/heron\_spectral.json}, \texttt{data/hw\_lucj\_n2\_result.json}). \emph{The raw measured
bitstrings are not part of the deposit}: the counts can be obtained only from the IBM~Quantum
service, by retrieving the recorded job identifier (\texttt{service.job(...)}) in place of the
submission cell of the notebooks in \texttt{notebooks/}, which requires an account and a job whose
results are still retained, so the classical post-processing of the raw counts cannot be re-run from the deposited
files alone. The hardware runs themselves cannot be regenerated in any case---the jobs are closed, and
a device executed today has a different calibration, so new counts would be new data rather than a
reproduction. The dequantization control of Sec.~\ref{app:dequant} is regenerated seed by seed by \texttt{src/gflow\_dequant.py}, which needs \texttt{pyscf} and \texttt{torch} and writes \texttt{data/gflow.json}. The $L=14$ ground-state checkpoint ($94$~MB), on which the largest vacuity factor and the $L=14$ support count rest, is deposited in \texttt{ckpt/} with its SHA-256 checksums. Known gaps between the deposited code and the
```

now reads:

```tex
The dequantization control of Sec.~\ref{app:dequant} is regenerated seed by seed by \texttt{src/gflow\_dequant.py}, which needs \texttt{pyscf} and \texttt{torch} and writes \texttt{data/gflow.json}. The $L=14$ ground-state checkpoint ($94$~MB), on which the largest vacuity factor and the $L=14$ support count rest, is deposited in \texttt{ckpt/} with its SHA-256 checksums. Known gaps between the deposited code and the
```

## sec_10.tex

### release lines 221-222 (replaced)

```tex
We acknowledge the use of IBM~Quantum services for this work; the two device runs reported in
Sec.~\ref{sec:scope:hw} were executed on \texttt{ibm\_fez} and \texttt{ibm\_marrakesh}. The views
```

now reads:

```tex
We acknowledge the use of IBM~Quantum services for this work; the two device runs noted in
Sec.~\ref{sm:scope} were executed on \texttt{ibm\_fez} and \texttt{ibm\_marrakesh}. The views
```

## sm_scope.tex

### release lines 1-1 (replaced)

```tex
\section{Hardware and scope: full account}
```

now reads:

```tex
\section{Scope: full account}
```

### release lines 10-28 (replaced)

```tex
\paragraph{The two hardware runs.} Both device runs are reported in full and they are the entirety of the quantum hardware used in this
work. The first was executed on the IBM~Heron processor \texttt{ibm\_fez} (job
\texttt{d9s16avpemts73ct6g8g}): the $L=6$ periodic Hubbard ring at $U/t=4$, $\eta=0.15\,t$, $12$ qubits,
$3.5\times10^{5}$ computational-basis shots---$5\times10^{4}$ on each of seven circuits
$k=0,\dots,6$, of which $k=0$ is the $t=0$ reference---pooled into one subspace,
under measurement twirling, Pauli twirling and dynamical decoupling, with no readout-error
mitigation applied to the counts and no zero-noise extrapolation
(Fig.~\ref{fig:heron}(a) is the result; the circuits prepared the computational-basis $(N{+}1)$ determinant with orbitals $0$--$3\!\uparrow$ and $0$--$2\!\downarrow$ occupied rather than $\hat c^{\dagger}_p\lvert\Psi_0\rangle$, and their hopping layer applies $R_{xx}R_{yy}$ rotations without Jordan--Wigner strings, so they are not the primitive of Fig.~\ref{fig:circ}; because the retained subspace is the whole sector, neither affects the reported $A_{\rm hw}$). No configuration recovery was applied. The notebook post-selected on particle number and $S_z$ after
reading each bitstring in reversed qubit order (\texttt{int(bs[::-1],2)}), which maps qubit $q$ to $11-q$:
site $i$ to $5-i$ and spin $\uparrow$ to $\downarrow$. The circuits conserve both quantities and start
from four $\uparrow$ and three $\downarrow$ electrons, so under that reading their own output falls in
the $S_z=-\tfrac12$ sector and is discarded: in a noiseless statevector simulation of the same seven circuits,
$0$ of $3.5\times10^{5}$ shots are retained, against all $3.5\times10^{5}$ and the full sector of $300$
under the correct reading (\texttt{src/hw\_bitorder\_check.py}, \texttt{data/hw\_bitorder\_check.json}).
The $|\mathcal S|=300$ determinants retained from the device therefore all came from device errors that
changed the spin balance; the device counts needed to redo the post-selection are not part of the
deposit (Sec.~\ref{sec:scope:data}). The recorded relative $L_1$ is $0.0$: the reconstructed $A_{\rm hw}$ agrees with the exact $A$
to $7.9\times10^{-15}$ at every frequency, not byte for byte, because the retained subspace spans the
whole sector.
```

now reads:

```tex
\paragraph{Device runs.} Two exploratory runs on IBM~Heron processors add no evidence to this work,
and no result of it rests on them: on \texttt{ibm\_fez}, circuits that are not the primitive of
Fig.~\ref{fig:circ} sampled the $L=6$ ring at $U/t=4$ ($12$ qubits, $3.5\times10^{5}$ shots), and the
post-selection read every bitstring in reversed qubit order, under which the circuits' own output is
discarded, so each of the $300$ retained determinants arose from a device error and the retained
subspace, the whole $(N{+}1)$ sector, returns the exact $A(\omega)$ by coverage; on
\texttt{ibm\_marrakesh}, stretched $\mathrm{N_2}$ ($24$ qubits) was run with no random-configuration
control and without recording the recovered subspace sizes; their arrays, shot counts and job
identifiers are deposited (\texttt{data/heron\_spectral.json}, \texttt{data/hw\_lucj\_n2\_result.json};
bit-order check \texttt{src/hw\_bitorder\_check.py}, \texttt{data/hw\_bitorder\_check.json}), their raw
bitstrings are not, and neither run can be regenerated.
```

### release lines 30-64 (replaced)

```tex
That number is exact \emph{by coverage}, and the distinction is not rhetorical: at
$|\mathcal S|=\dim$, Rayleigh--Ritz on $\operatorname{span}(\mathcal S)$ \emph{is} diagonalization of
the full sector, $P_{\mathcal S}HP_{\mathcal S}=H$, and a vanishing error is an identity rather than a
measurement of device fidelity. The run establishes only that the steps of the pipeline execute end
to end on a processor---transpiled circuits, device counts, post-selection to a particle-number and
$S_z$ sector, and the classical Lehmann assembly of Sec.~\ref{sec:method}. The post-selection itself was
wrong: it read every bitstring in reversed bit order, so every determinant it kept came from a device
error, and the run is not evidence that the device sampled the intended distribution. It establishes
nothing about accuracy, neither at the complete coverage it reached, where agreement is an identity, nor
at incomplete coverage, the regime the rest of this paper is about.
The runs are therefore described as \emph{executed on} the named processors, not as demonstrating the method on them, and the coverage qualification applies wherever they are mentioned.

\paragraph{The second run, and a comparison the record cannot decide.} Stretched $\mathrm{N_2}$ ($R=2.0$~\AA,
CAS(10e,12o)/cc-pVDZ, $24$ qubits) was executed on \texttt{ibm\_marrakesh} (job
\texttt{da125f2ein7c73bcsqs0}, $10^{4}$ shots, CCSD-seeded local unitary cluster-Jastrow (LUCJ) circuit). Averaged over $8$
recovery seeds on the cached device counts it reaches $0.59\pm0.14$~mHa of the exact active-space value
(best seed $0.40$~mHa), while the noiseless statevector simulation of the same circuit stops at
$29.5$~mHa (Fig.~\ref{fig:heron}(b)). That contrast invites the reading that device noise \emph{helps}, and the record cannot support it, because the mechanism named for it is a statement about a quantity the
record does not contain: recovery is said to turn corrupted bitstrings into valid configurations and
thereby \emph{enlarge} the variational subspace---a valid-shot economy, not a fidelity gain---so the
comparison is decided by the two recovered subspace sizes, and neither is stored. The same gap runs
through the controlled sweeps, whose scripts save at most the size of the naively post-selected subspace
and not of the recovered one (the spectral script of Fig.~\ref{fig:noise} saves neither), so the places
where an error \emph{below} the noiseless value is reported beyond the seed bands---$\mathrm{NH_3}$ from
$2.53$ to $1.24$~mHa at $\varepsilon=10\%$, CO from $1.73$ to $0.63$~mHa at $20\%$ and $\mathrm{C_2}$
from $0.184$ to $0.026$~mHa at $30\%$ in Fig.~\ref{fig:noiserec}(a), and the spectral relative $L_1$ of
Fig.~\ref{fig:noise}, from $0.105$ at $\varepsilon=0$ to $0.045$ at $13\%$---rest on the same unrecorded
variable, at the price of one line in each script, both of which already build the subspace.%
% 2026-09-30 (final pass), layout only: \looseness=-1 sets this paragraph one line shorter
% (its last line held two words).  The final-review rewording of Sec. S5 moved the page
% breaks that follow, and p. 55 came out 4.97 pt overfull again at the \pagebreak below (the
% footnote of this page booked at half height, see that comment).  No word changes; the
% paragraph has no \cite, so its lines are the same in main.pdf and sm_pra.pdf.
\looseness=-1\relax

```

now reads:

```tex
\paragraph{Noise and recovery.} Recovery is said to turn corrupted bitstrings into valid
configurations and thereby \emph{enlarge} the variational subspace---a valid-shot economy, not a
fidelity gain---so a comparison with the noiseless value is decided by the size of the recovered
subspace. The scripts of the controlled sweeps save at most the size of the naively post-selected
subspace and not of the recovered one (the spectral script of Fig.~\ref{fig:noise} saves neither), so
the places where an error \emph{below} the noiseless value is reported beyond the seed
bands---$\mathrm{NH_3}$ from $2.53$ to $1.24$~mHa at $\varepsilon=10\%$, CO from $1.73$ to $0.63$~mHa
at $20\%$ and $\mathrm{C_2}$ from $0.184$ to $0.026$~mHa at $30\%$ in Fig.~\ref{fig:noiserec}(a), and
the spectral relative $L_1$ of Fig.~\ref{fig:noise}, from $0.105$ at $\varepsilon=0$ to $0.045$ at
$13\%$---rest on that unrecorded variable, at the price of one line in each script, both of which
already build the subspace.
```

### release lines 66-66 (replaced)

```tex
$\mathrm{N_2}$ in a large active space, with the same LUCJ ansatz used for the run above,
```

now reads:

```tex
$\mathrm{N_2}$ in a large active space, with a local unitary cluster-Jastrow (LUCJ) ansatz,
```

### release lines 70-72 (replaced)

```tex
improved energies''~\cite{vaquero2026}: the same molecule, the same ansatz and the same mechanism, and there it is carried by a controlled comparison under a simple noise model that the work introduces. The same work also reports that ``recovery alone (starting from randomly generated configurations) can efficiently construct accurate CI spaces''~\cite{vaquero2026}. Our record contains no control of that kind (uniformly random bitstrings passed through the same recovery with the same eight seeds), so it cannot show that the device counts contributed anything to the $0.59$~mHa beyond what recovery supplies on its own. What this work cannot claim is its own measurement of the effect, not the effect: our two arms differ in
the one quantity neither of them recorded, so they cannot be evidence either way, whatever the
literature independently establishes.
```

now reads:

```tex
improved energies''~\cite{vaquero2026}, carried there by a controlled comparison under a simple noise model that the work introduces. The same work also reports that ``recovery alone (starting from randomly generated configurations) can efficiently construct accurate CI spaces''~\cite{vaquero2026}. No such control was run in this work.
```

### release lines 84-92 (replaced)

```tex
\paragraph{Scale.} These runs are $12$ and $24$ qubits, and the same class of observable has been
executed elsewhere at larger scale---retarded spin Green's functions of the one-dimensional
Fermi--Hubbard model on IBM~Heron~r2 at up to $30$ qubits~\cite{vilchez2025}, an ARPES-type
$A(k,\omega)$ of a free-fermion $27$-site chain at $54$ qubits~\cite{granet2026akw}---and, at a
comparable $20$ sites on trapped ions, dynamical structure factors against scattering
data~\cite{quantinuum_dsf2026}.
% Paragraph break (2026-09-28): the three citations above print as [5]-[7] in main.pdf and as
% two-digit numbers in sm_pra.pdf; ending the paragraph here keeps that width difference from
% reflowing the rest of it, which pra_split.py's page comparison would report as a difference.
```

now reads:

```tex
\paragraph{Shot cost of the momentum-resolved observable.} At $L=8$ the momentum-resolved
observable is sixteen channels, $8.3\times10^{7}$ shots as measured for Fig.~\ref{fig:akwsampled} to
reach its relative $L_1$ of $2.16\times10^{-3}$. That is not its cost at the $0.05$ threshold of the
resource scan, which was not measured: half that budget already gives $6.16\times10^{-3}$.
```

### release lines 94-118 (deleted)

```tex
In those works the device data
determine the reported spectrum, as they do in the four-qubit trapped-ion Hubbard-dimer reconstruction
of Ref.~\cite{greenediniz2024}; in our spectral run they do not: the retained subspace is the whole
$300$-dimensional $(N{+}1)$ sector and every retained determinant came from a device error, so the
spectral function it returns is the exact-diagonalization one, the run is not an entry in that class,
and it is not offered as a scale result. The run is also at a size where the coverage question is
degenerate, and Sec.~\ref{sec:prices} says how degenerate in
shots: the $3.5\times10^{5}$ spent---seven circuits at $5\times10^{4}$ shots each, pooled into one
subspace---are $14.4$ times the $(2.43\pm0.34)\times10^{4}$ that Table~\ref{tab:shots} measures for one
orbital seed and one branch at $L=6$, so complete coverage of a $300$-dimensional sector is what the
budget predicts for a correctly read run, as the noiseless simulation above confirms; the coverage
this run recorded came from device errors and does not test that prediction.\footnote{Table~\ref{tab:shots} replicates the repository protocol, $K=18$ with $19$ snapshots,
at $U/t=8$ rather than the $U/t=4$ of the device run; this is a consistency statement between two
budgets, not an identity.} The first size at which a device run would test the reconstruction rather
than only the execution of the pipeline is $L=8$: sector $3\,920$, at $(2.70\pm0.13)\times10^{5}$ shots per orbital seed and
branch---a factor $1.3$ \emph{below} what was already spent here, so what confined this run to $L=6$
was not the shot budget. The momentum-resolved observable at that size is a different matter: sixteen
channels, $8.3\times10^{7}$ shots as measured for Fig.~\ref{fig:akwsampled} to reach its relative $L_1$ of
$2.16\times10^{-3}$, over two orders of magnitude above this job. That is not its cost at the $0.05$
threshold used above, which was not measured: half that budget already gives $6.16\times10^{-3}$.
% 2026-09-30, layout only: the page ends here.  Without the break the page builder ran one
% line further, counting the two-line footnote above at half its height (revtex4-2 one-column
% \count\footins=500), and the page came out 4.97 pt overfull (p. 55).  The break is vertical
% only, so no line of any paragraph is reset, in main.pdf or in sm_pra.pdf.  No text changes.
\pagebreak
```

## app_carried_repro.tex

### release lines 83-83 (replaced)

```tex
suite is $300$, the dimension of the $L=6$ Hubbard sector of Sec.~\ref{sec:scope:hw}, against
```

now reads:

```tex
suite is $300$, the dimension of the $L=6$ Hubbard sector, against
```

## carried/fig_noise.tex

### release lines 11-11 (replaced)

```tex
second comparison is not used here (Sec.~\ref{sec:scope:hw}): it is decided by the size of the recovered
```

now reads:

```tex
second comparison is not used here (Sec.~\ref{sm:scope}): it is decided by the size of the recovered
```

## carried/fig_noiserec.tex

### release lines 10-10 (replaced)

```tex
that second comparison is not used here (Sec.~\ref{sec:scope:hw}): it is decided by the size of
```

now reads:

```tex
that second comparison is not used here (Sec.~\ref{sm:scope}): it is decided by the size of
```

## carried/fig_lattice.tex

### release lines 20-21 (replaced)

```tex
$\omega$ is measured from the chemical potential $\mu=U/2$, whereas Figs.~\ref{fig:method},
\ref{fig:akwsampled} and~\ref{fig:heron} measure it from $E_0$),
```

now reads:

```tex
$\omega$ is measured from the chemical potential $\mu=U/2$, whereas Figs.~\ref{fig:method}
and~\ref{fig:akwsampled} measure it from $E_0$),
```

## carried/fig_heron.tex (no longer \input; the file itself is unchanged)

```tex
% carried verbatim from the v2 tree (hardware.tex); NOT part of this session's deliverable.
% CENTRING, measured 2026-09-19 on the composed page at 600 dpi.  The tikz bounding
% box of this two-panel fragment is not symmetric about its ink: the ink is 421.6 pt
% wide but the box TeX centres is 404.2 pt, because panel (b)'s title overhangs the
% box on the right while panel (a)'s reserved label space on the left is not filled.
% \centering therefore put the graphic 17.3 pt right of the text-block centre, while
% the other sixteen figures sit between -2.9 and +1.3 pt.  \leftskip = -2d plus 1fil
% biases a centred line by exactly -d (here d = 17.3 pt) and changes nothing else;
% measured again after the fix: +0.01 pt.  A \hspace* pair does NOT work here: it
% starts the paragraph before the fragment is read, which turns the end-of-line
% spaces of the fragment's three \definecolor lines into real interword glue and
% delivers only 12.7 pt of the 17.3 pt asked for.  The fragment itself is untouched
% (it is one of the eight byte-identical artefacts check_figures.py regenerates).
\begin{figure*}[tp]\centering
\leftskip=-34.6pt plus 1fil\relax
\input{figs/fig_hardware_hero_frag.tex}
\caption{\textbf{The complete quantum-hardware record of this work: two runs on IBM~Heron processors.}
\textbf{(a)}~Electron-addition branch of the local spectral function, $A^{+}_{0\uparrow}(\omega)$, of the
$L=6$ periodic Hubbard ring at $U/t=4$,
$\eta=0.15\,t$ ($12$ qubits; the classical validation of Fig.~\ref{fig:method} is a different instance,
$U/t=8$), computed on the subspace post-selected from $350{,}000$ computational-basis bitstrings---$50{,}000$ from each of seven
circuits---measured on \texttt{ibm\_fez} (measurement and Pauli
twirling and dynamical decoupling; no readout-error mitigation and no configuration recovery). The retained shots
covered the complete $(N{+}1)$ sector ($|\mathcal S|\!=\!300/300$), so the reconstruction is exact by
coverage; they were post-selected in reversed bit order, so every retained shot is a device error
(Sec.~\ref{sec:scope:hw}). The panel is an execution record, and neither a fidelity measurement nor a
claim of advantage.
\textbf{(b)}~Ground-state energy error of stretched $\mathrm{N_2}$ ($R=2.0$~\AA, CAS(10e,12o)/cc-pVDZ,
$24$ qubits) versus self-consistent configuration-recovery step, on \texttt{ibm\_marrakesh} ($10{,}000$ shots,
CCSD-seeded local unitary cluster-Jastrow (LUCJ) circuit). Recovery on the device counts (solid) reaches
$0.59\pm0.14$~mHa of the exact FCI value---well below chemical accuracy---whereas the noiseless statevector
simulation of the \emph{same} circuit (dashed) is stranded at $29.5$~mHa (its recorded recovery history
ends at step~$3$, unchanged since step~$2$). That contrast is not read here as device noise \emph{helping}: the
mechanism such a reading needs is decided by the two recovered subspace sizes, and the record
stores neither (Sec.~\ref{sec:scope:hw}). The
solid curve is the mean over $8$ recovery seeds on the cached device counts; the shaded envelope and the
capped whiskers at every step both show the $\pm1$ population s.d.\ (normalized by $8$, not $7$) across those $8$ seeds (widest at
step~$0$, $31.1\pm4.0$~mHa, where it starts together with the noiseless curve, and tapering as recovery
converges to $0.59\pm0.14$~mHa; best seed $0.40$~mHa). No random-bitstring control was run, so the
panel cannot separate what the device counts contributed from what recovery alone supplies; the
controlled recovery-versus-no-recovery comparison, under a simulated bit-flip model, is
Fig.~\ref{fig:noiserec}. These two runs are the entirety of the quantum
hardware used here; every other result in this work is an exact classical simulation.}\label{fig:heron}\end{figure*}

```
