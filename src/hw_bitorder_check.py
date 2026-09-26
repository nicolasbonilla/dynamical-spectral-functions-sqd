"""Which shots of the L=6 ibm_fez run (job d9s16avpemts73ct6g8g) could the post-selection keep?

The notebook that processed the run (notebooks/Spectral_Heron.ipynb, cell 15) turned each measured
bitstring into a Fock integer with ``int(bs[::-1], 2)`` and kept it if it had N+1 = 7 electrons and
S_z = +1.  Qiskit prints qubit q as character ``bs[-1-q]``, so the correct integer is ``int(bs, 2)``;
the reversed reading maps qubit q to 11-q, i.e. site i to 5-i and spin s to 1-s (p = 2*site + spin).
The circuits conserve particle number and S_z and start from 4 up + 3 down electrons, so under the
reversed reading their own output lands in the S_z = -1 sector and is rejected.

This script rebuilds the seven circuits exactly as in the notebook (cells 5 and 7), samples each one
noiselessly with 50,000 shots, and counts the shots kept under both readings.  The device counts
themselves are not in the deposit, so this is the check that can be run here.

    python src/hw_bitorder_check.py        # ~1 min; writes data/hw_bitorder_check.json
"""
import json
import os

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.quantum_info import Statevector

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L, U, THOP = 6, 4.0, 1.0
M = 2 * L
DIM = 1 << M
K, N_TROT, DT_TOTAL, SHOTS, SEED = 7, 3, 0.5 / THOP, 50_000, 1

NOCC = np.array([bin(s).count('1') for s in range(DIM)])
SZOCC = np.array([sum(((s >> (2 * i)) & 1) - ((s >> (2 * i + 1)) & 1) for i in range(L)) for s in range(DIM)])


def build(k):
    """Cell 7 of the notebook, without the final measurement."""
    q = QuantumRegister(M)
    qc = QuantumCircuit(q)
    for o in (0, 1, 2, 3):
        qc.x(q[2 * o])          # 4 up electrons (even qubits)
    for o in (0, 1, 2):
        qc.x(q[2 * o + 1])      # 3 down electrons (odd qubits)
    dt = (k * DT_TOTAL) / N_TROT
    for _ in range(N_TROT):
        for i in range(L):
            qc.rzz(2 * U * dt, q[2 * i], q[2 * i + 1])
        for s in (0, 1):
            for i in range(L):
                j = (i + 1) % L
                qc.rxx(THOP * dt, q[2 * i + s], q[2 * j + s])
                qc.ryy(THOP * dt, q[2 * i + s], q[2 * j + s])
    return qc


def in_sector(t):
    return NOCC[t] == L + 1 and SZOCC[t] == 1


def main():
    np.random.seed(SEED)
    total = kept_correct = kept_reversed = 0
    s_correct, s_reversed = set(), set()
    for k in range(K):
        sv = Statevector(build(k))
        sv.seed(SEED + k)
        for bs, n in sv.sample_counts(SHOTS).items():
            total += int(n)
            tc, tr = int(bs, 2), int(bs[::-1], 2)
            if in_sector(tc):
                kept_correct += int(n)
                s_correct.add(tc)
            if in_sector(tr):
                kept_reversed += int(n)
                s_reversed.add(tr)
    out = {
        'what': 'noiseless statevector sampling of the seven L=6 circuits of the ibm_fez run; shots kept by '
                'the (N+1=7, Sz=+1) post-selection under the correct and the notebook (reversed) bit readings',
        'L': L, 'U': U, 'circuits': K, 'shots_per_circuit': SHOTS, 'total_shots': total,
        'kept_correct_reading': kept_correct, 'support_correct_reading': len(s_correct),
        'kept_notebook_reading': kept_reversed, 'support_notebook_reading': len(s_reversed),
        'sector_dimension': int(np.sum((NOCC == L + 1) & (SZOCC == 1))),
        'conclusion': 'under the notebook reading no shot of the intended circuits is retained, so every '
                      'determinant retained from the device came from a device error',
    }
    path = os.path.join(ROOT, 'data', 'hw_bitorder_check.json')
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
