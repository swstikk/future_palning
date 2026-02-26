# ⚛️ Quantum Bioinformatics — Deep Chapter-wise Syllabus PART 3
## Phase 2: QC Theory + Qiskit (Weeks 11-18) + Phase 3: Algorithms (Weeks 19-26)

---

# PHASE 2: QUANTUM COMPUTING THEORY + QISKIT (Weeks 11-18)

---

## Module QC1.1: Qubits & Quantum States

> **PREREQUISITES: Part 1 (Math) + Part 2 (Physics) ALL gates passed.**
> You must know: complex numbers, vectors in ℂ², inner product, Born rule,
> |0⟩=[1,0]ᵀ, |1⟩=[0,1]ᵀ, |+⟩, |-⟩, matrix-vector multiplication.
> From Part 2: measurement postulate, eigenvalues, Pauli matrices.

```
QC1.1.1  The Qubit — Physical vs Mathematical
├── Classical bit vs Qubit (THE fundamental difference):
│   ┌──────────────────────┬──────────────────────────────────┐
│   │ Classical Bit         │ Qubit                            │
│   ├──────────────────────┼──────────────────────────────────┤
│   │ Value: 0 OR 1         │ Value: α|0⟩ + β|1⟩ (BOTH!)      │
│   │ State: bit b∈{0,1}    │ State: vector in ℂ²              │
│   │ Deterministic          │ Probabilistic until measured     │
│   │ Copy freely            │ NO-CLONING theorem              │
│   │ Read without changing  │ Measurement DESTROYS superposition│
│   │ n bits → n values      │ n qubits → 2ⁿ amplitudes       │
│   └──────────────────────┴──────────────────────────────────┘
│
├── General single-qubit state:
│   |ψ⟩ = α|0⟩ + β|1⟩     where α,β ∈ ℂ
│   Constraint: |α|² + |β|² = 1  (probabilities sum to 1)
│
│   WORKED EXAMPLES:
│   |ψ₁⟩ = |0⟩ → α=1, β=0. P(0)=1, P(1)=0.  (definitely |0⟩)
│   |ψ₂⟩ = |+⟩ = (1/√2)|0⟩+(1/√2)|1⟩ → P(0)=P(1)=1/2 (coin flip)
│   |ψ₃⟩ = (√3/2)|0⟩+(1/2)|1⟩ → P(0)=3/4, P(1)=1/4  (biased)
│   |ψ₄⟩ = (1/√2)|0⟩+(i/√2)|1⟩ → P(0)=P(1)=1/2 but DIFFERENT from |+⟩!
│   (Same probabilities, different PHASE → different physics)
│
├── Physical realizations (how real qubits work):
│   ┌─────────────────────┬────────────────────────────────────┐
│   │ Technology           │ What is |0⟩ and |1⟩                │
│   ├─────────────────────┼────────────────────────────────────┤
│   │ Superconducting (IBM)│ Current flowing ↻ or ↺ in loop    │
│   │ Trapped ion (IonQ)   │ Electron in ground vs excited state│
│   │ Photonic (Xanadu)    │ Polarization: horizontal/vertical  │
│   │ Spin qubit           │ Electron spin up ↑ or down ↓      │
│   └─────────────────────┴────────────────────────────────────┘
│   IBM Eagle processor: 127 superconducting qubits (2023)
│   IBM Heron: 133 qubits, 2-qubit error <1% (2024)
│
├── Code:
│   import numpy as np
│   from qiskit.quantum_info import Statevector
│   # Create various qubit states:
│   psi_0 = Statevector([1, 0])     # |0⟩
│   psi_plus = Statevector([1/np.sqrt(2), 1/np.sqrt(2)])  # |+⟩
│   psi_biased = Statevector([np.sqrt(3)/2, 1/2])         # 75/25 state
│   psi_phase = Statevector([1/np.sqrt(2), 1j/np.sqrt(2)]) # same probs, diff phase
│   for name, sv in [('|0⟩',psi_0),('|+⟩',psi_plus),('biased',psi_biased),('phase',psi_phase)]:
│       print(f"{name}: valid={sv.is_valid()}, probs={sv.probabilities()}")
│
└── Exit check:
    1. Is [0.6, 0.8] a valid qubit state? |0.6|²+|0.8|²=0.36+0.64=1 ✓
    2. Is [0.5, 0.5] valid? |0.5|²+|0.5|²=0.5 ≠ 1 ✗ (not normalized!)
    3. What is |α|² if α=(1+i)/2? |(1+i)/2|² = (1²+1²)/4 = 2/4 = 0.5

QC1.1.2  Bloch Sphere — The Qubit Visualization Tool
├── ANY single qubit state can be written as:
│   |ψ⟩ = cos(θ/2)|0⟩ + e^(iφ)·sin(θ/2)|1⟩
│   θ ∈ [0, π]:  polar angle (how far from north pole)
│   φ ∈ [0, 2π): azimuthal angle (which direction on equator)
│   → Every qubit maps to a point on a SPHERE (the Bloch sphere)
│
├── KEY STATES on the Bloch sphere (MEMORIZE):
│   ┌────────────────┬──────────────┬──────────────────────────┐
│   │ State          │ (θ, φ)       │ Bloch vector [x,y,z]     │
│   ├────────────────┼──────────────┼──────────────────────────┤
│   │ |0⟩            │ (0, -)       │ [0, 0, +1]  North pole   │
│   │ |1⟩            │ (π, -)       │ [0, 0, -1]  South pole   │
│   │ |+⟩            │ (π/2, 0)     │ [+1, 0, 0]  +x axis      │
│   │ |-⟩            │ (π/2, π)     │ [-1, 0, 0]  -x axis      │
│   │ |+i⟩=(|0⟩+i|1⟩)/√2│ (π/2, π/2)│ [0, +1, 0]  +y axis   │
│   │ |-i⟩=(|0⟩-i|1⟩)/√2│ (π/2, 3π/2)│[0, -1, 0]  -y axis   │
│   └────────────────┴──────────────┴──────────────────────────┘
│
├── WORKED EXAMPLE — find Bloch coordinates:
│   |ψ⟩ = (√3/2)|0⟩ + (1/2)|1⟩
│   cos(θ/2) = √3/2 → θ/2 = π/6 → θ = π/3
│   e^(iφ)·sin(θ/2) = 1/2 → e^(iφ)·(1/2) = 1/2 → e^(iφ) = 1 → φ = 0
│   Bloch: x=sin(π/3)cos(0)=√3/2, y=sin(π/3)sin(0)=0, z=cos(π/3)=1/2
│   → [√3/2, 0, 1/2] — between north pole and +x axis, upper hemisphere
│
├── Gates AS rotations on the Bloch sphere:
│   X gate = π rotation about x-axis: |0⟩ ↔ |1⟩ (north↔south)
│   Y gate = π rotation about y-axis: |0⟩ → i|1⟩
│   Z gate = π rotation about z-axis: |+⟩ ↔ |-⟩ (swaps +x and -x)
│   H gate = π rotation about (x+z)/√2 axis: |0⟩ ↔ |+⟩
│   Rx(θ) = rotation by θ about x-axis
│   Ry(θ) = rotation by θ about y-axis (moves from pole toward equator)
│   Rz(θ) = rotation by θ about z-axis (changes φ only)
│
├── Code:
│   from qiskit.visualization import plot_bloch_multivector
│   from qiskit.quantum_info import Statevector
│   # Visualize several states:
│   for state_label, sv in [('|0⟩', [1,0]), ('|+⟩', [1/np.sqrt(2), 1/np.sqrt(2)]),
│                            ('|+i⟩', [1/np.sqrt(2), 1j/np.sqrt(2)])]:
│       psi = Statevector(sv)
│       fig = plot_bloch_multivector(psi)
│       fig.suptitle(state_label)
│       fig.savefig(f'bloch_{state_label}.png')
│
└── Exit check:
    1. Where is Ry(π/2)|0⟩ on the Bloch sphere?
       Ry(π/2)|0⟩ = cos(π/4)|0⟩+sin(π/4)|1⟩ = |+⟩ → +x axis ✓
    2. Where is Rz(π/2)|+⟩?
       Rz(π/2) rotates by π/2 about z → |+⟩ goes to |+i⟩ (+y axis) ✓
    3. Plot all 6 key states on Bloch sphere using Qiskit.

QC1.1.3  Phase: Global vs Relative — The Subtlety That Makes Quantum Work
├── Global phase: e^(iα)|ψ⟩ is PHYSICALLY IDENTICAL to |ψ⟩
│   WHY? Born rule: P = |⟨φ|e^(iα)ψ⟩|² = |e^(iα)|²|⟨φ|ψ⟩|² = |⟨φ|ψ⟩|²
│   The e^(iα) cancels! → undetectable by any measurement
│   Example: |0⟩ and i|0⟩ and -|0⟩ are the SAME physical state
│
├── Relative phase: PHYSICALLY OBSERVABLE and CRUCIAL
│   |+⟩ = (1/√2)(|0⟩ + |1⟩)  → relative phase between |0⟩,|1⟩ is 0
│   |-⟩ = (1/√2)(|0⟩ - |1⟩)  → relative phase is π (the minus sign)
│   Same P(0)=P(1)=1/2 for BOTH states!
│   But ⟨X⟩ = +1 for |+⟩, ⟨X⟩ = -1 for |-⟩ → DIFFERENT measurements!
│
├── How relative phase arises:
│   Z gate: Z|+⟩ = Z(1/√2)(|0⟩+|1⟩) = (1/√2)(|0⟩-|1⟩) = |-⟩
│   Z didn't change probabilities but added π relative phase!
│   This is INVISIBLE in Z-measurement but VISIBLE in X-measurement
│
├── INTERFERENCE — why phase matters:
│   |ψ₁⟩ = (1/√2)(|0⟩+|1⟩) → apply H → H|+⟩ = |0⟩ (constructive)
│   |ψ₂⟩ = (1/√2)(|0⟩-|1⟩) → apply H → H|-⟩ = |1⟩ (destructive)
│   SAME probabilities before H, DIFFERENT outcomes after H!
│   This is quantum INTERFERENCE — the engine of all quantum speedups
│
├── VQE/Grover connection:
│   Grover oracle: marks target with -1 phase (relative phase change)
│   Diffuser: amplifies marked state via interference
│   VQE: rotation gate angles θ control relative phase → change ⟨H⟩
│
├── Code — demonstrating phase matters:
│   from qiskit import QuantumCircuit
│   from qiskit.quantum_info import Statevector
│   # State 1: |+⟩, then H → should give |0⟩
│   qc1 = QuantumCircuit(1); qc1.h(0); qc1.h(0)
│   print(Statevector(qc1))  # [1, 0] = |0⟩ ✓
│   # State 2: |-⟩, then H → should give |1⟩
│   qc2 = QuantumCircuit(1); qc2.x(0); qc2.h(0); qc2.h(0)
│   print(Statevector(qc2))  # [0, 1] = |1⟩ ✓
│   # Same probs (50/50) before final H, different outcomes!
│
└── Exit check:
    1. Are e^(iπ/4)|+⟩ and |+⟩ physically different? NO (global phase)
    2. Are |+⟩ and |-⟩ physically different? YES (relative phase)
    3. Apply H to (1/√2)(|0⟩+i|1⟩). What's the result?
       H[(1/√2)(|0⟩+i|1⟩)] = (1/2)[(1+i)|0⟩+(1-i)|1⟩]
       P(0) = |1+i|²/4 = 2/4 = 1/2, P(1) = 1/2 (equal, but different from |+⟩ case!)
    4. Verify in Qiskit code.

═══════════════════════════════════════════
 GATE TO QC1.2 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know: qubit = α|0⟩+β|1⟩ with |α|²+|β|²=1
 □ Can check if a vector is a valid qubit state
 □ Know 6 key Bloch sphere states (table above) from memory
 □ Can find θ,φ for a given state (e.g. (√3/2)|0⟩+(1/2)|1⟩ → θ=π/3, φ=0)
 □ Know: gates = rotations on Bloch sphere (X=x-rot, Z=z-rot, H=diagonal)
 □ Know: global phase undetectable; relative phase → interference
 □ Demonstrated interference: H|+⟩=|0⟩ vs H|-⟩=|1⟩ in code
 □ Created Statevector objects and used plot_bloch_multivector
═══════════════════════════════════════════
```

---

## Module QC1.2: Quantum Gates — Complete Reference

> **PREREQUISITES: QC1.1 gate passed.**
> Need: matrix-vector multiplication (M2.2), Pauli matrices (M2.3),
> Euler's formula (M1.1), unitary = U†U=I (M2.3).

```
QC1.2.1  Single-Qubit Gates — Complete Reference
├── PAULI GATES (memorize):
│   X = [[0,1],[1,0]]   "bit flip / NOT"
│     X|0⟩=|1⟩,  X|1⟩=|0⟩
│   Y = [[0,-i],[i,0]]  "bit+phase flip"
│     Y|0⟩=i|1⟩, Y|1⟩=-i|0⟩
│   Z = [[1,0],[0,-1]]  "phase flip"
│     Z|0⟩=|0⟩,  Z|1⟩=-|1⟩
│
├── PHASE GATES:
│   S = [[1,0],[0,i]]  (π/2 phase on |1⟩);  S²=Z
│   T = [[1,0],[0,e^(iπ/4)]]  (π/4 phase);  T²=S
│   Phase chain: T→T=S→S=Z→Z=I
│   WORKED: S|+⟩ = (1/√2)(|0⟩+i|1⟩) = |+i⟩ (rotated to +y)
│
├── HADAMARD:
│   H = (1/√2)[[1,1],[1,-1]]
│   H|0⟩=|+⟩, H|1⟩=|-⟩, H|+⟩=|0⟩, H|-⟩=|1⟩.  H²=I.
│   Key: H switches Z-basis ↔ X-basis
│
├── ROTATION GATES (VQE parameters live here):
│   Rx(θ) = [[cos(θ/2), -i·sin(θ/2)], [-i·sin(θ/2), cos(θ/2)]]
│   Ry(θ) = [[cos(θ/2), -sin(θ/2)], [sin(θ/2), cos(θ/2)]]
│   Rz(θ) = [[e^(-iθ/2), 0], [0, e^(iθ/2)]]
│   Ry(0)=I, Ry(π/3)|0⟩=(√3/2)|0⟩+(1/2)|1⟩, Ry(π)=-iY
│   VQE: optimizer adjusts θ₁,θ₂,... in Ry gates → minimize ⟨H⟩
│
├── UNIVERSAL SET: {H, T, CNOT} → ANY unitary (Solovay-Kitaev thm)
│
├── Code:
│   from qiskit import QuantumCircuit
│   from qiskit.quantum_info import Statevector, Operator
│   import numpy as np
│   qc = QuantumCircuit(1); qc.h(0); qc.s(0)
│   print(Statevector(qc))  # [0.707, 0.707j] = |+i⟩ ✓
│   T = Operator.from_label('T'); S = Operator.from_label('S')
│   print(np.allclose((T@T).data, S.data))  # True: T²=S ✓
│
└── Exit check:
    1. S|+⟩=? (answer: |+i⟩); T|+⟩=? (answer: (1/√2)(|0⟩+e^(iπ/4)|1⟩))
    2. Ry(π/3)|0⟩=? (answer: (√3/2)|0⟩+(1/2)|1⟩)
    3. Rz(π)|+⟩=? (answer: -i|-⟩ = |-⟩ up to global phase)

QC1.2.2  Two-Qubit Gates — Entanglement Creators
├── CNOT (CX) — MOST important 2-qubit gate:
│   4×4 matrix (basis |00⟩,|01⟩,|10⟩,|11⟩):
│   [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]
│   Rule: target FLIPPED only if control=|1⟩
│   Truth table:
│   |00⟩→|00⟩, |01⟩→|01⟩, |10⟩→|11⟩, |11⟩→|10⟩
│
│   WORKED — Bell state |Φ+⟩ step by step:
│   Start:  |00⟩
│   H on q0: (1/√2)(|0⟩+|1⟩)⊗|0⟩ = (1/√2)(|00⟩+|10⟩)
│   CNOT:    (1/√2)(|00⟩+|11⟩) = |Φ+⟩  ✓
│   (|10⟩→|11⟩ because ctrl=1 flips target)
│
├── CZ: minus only on |11⟩ → CZ|11⟩=-|11⟩, rest unchanged
├── SWAP: exchanges qubits → SWAP|01⟩=|10⟩; decomp = 3 CNOTs
│
├── Hardware error reality:
│   ┌──────────────┬──────────┬───────────┐
│   │ Gate type     │ Error    │ Time      │
│   ├──────────────┼──────────┼───────────┤
│   │ 1-qubit (Rz) │ ~0.01%   │ ~20 ns    │
│   │ 2-qubit (CX) │ ~0.5-1%  │ ~200 ns   │
│   │ Measurement  │ ~1-2%    │ ~500 ns   │
│   └──────────────┴──────────┴───────────┘
│   → VQE circuit design = MINIMIZE 2-qubit gate count!
│
├── Code:
│   qc = QuantumCircuit(2); qc.h(0); qc.cx(0,1)
│   bell = Statevector(qc)
│   print(bell.data)           # [0.707, 0, 0, 0.707]
│   print(bell.probabilities())# [0.5, 0, 0, 0.5] ✓
│
└── Exit check:
    1. CNOT|+0⟩ = ? → |Φ+⟩ (Bell state)
    2. CNOT|-0⟩ = ? → |Φ-⟩ = (1/√2)(|00⟩-|11⟩)
    3. Verify CZ = (H⊗I)·CNOT·(H⊗I) using np.kron

QC1.2.3  Circuit Building in Qiskit
├── Reading circuit diagrams:
│   Left→right = time; lines = qubits; boxes = gates
│   Dot+⊕ = CNOT (dot=control, ⊕=target)
│
├── Qiskit API:
│   qc = QuantumCircuit(n_qubits, n_classical_bits)
│   qc.h(0); qc.x(1); qc.ry(theta, 0); qc.rz(phi, 1)
│   qc.s(0); qc.t(0)           # phase gates
│   qc.cx(0,1); qc.cz(0,1)     # 2-qubit
│   qc.measure(0, 0)            # qubit→classical
│   qc.measure_all()            # all qubits
│   qc.draw('mpl')              # matplotlib
│
├── Code — H⊗H superposition:
│   qc = QuantumCircuit(2); qc.h([0,1])
│   sv = Statevector(qc)
│   print(sv.probabilities())   # [0.25, 0.25, 0.25, 0.25] ✓
│
└── Exit check:
    Build circuits for:
    1. |00⟩ → equal superposition of all 4 basis states (H⊗H)
    2. |00⟩ → |Φ+⟩ Bell state (H, CNOT)
    3. |00⟩ → |Ψ+⟩ = (1/√2)(|01⟩+|10⟩) (X(0), H(0), CNOT(0,1))
    Draw all 3 with qc.draw('mpl').

═══════════════════════════════════════════
 GATE TO QC1.3 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Write X,Y,Z,H,S,T matrices from memory
 □ Know: S²=Z, T²=S, H²=I
 □ Know Ry/Rz forms; can compute Ry(π/3)|0⟩ by hand
 □ Universal set: {H, T, CNOT} can build any gate
 □ CNOT truth table: |00⟩→|00⟩, |01⟩→|01⟩, |10⟩→|11⟩, |11⟩→|10⟩
 □ Created Bell state: H(0)→CNOT(0,1) gives |Φ+⟩
 □ Know hardware errors: 2-qubit ~10x worse → minimize in VQE
 □ Built and drew 3+ circuits in Qiskit
═══════════════════════════════════════════
```

## Module QC1.3: Quantum Measurement

```
QC1.3.1  Projective Measurement — Z Basis
├── Default measurement = computational (Z) basis
├── For state |ψ⟩ = α|0⟩ + β|1⟩:
│   P(outcome 0) = |α|²   → post-measurement state: |0⟩
│   P(outcome 1) = |β|²   → post-measurement state: |1⟩
│
├── Key: measurement DESTROYS superposition (irreversible!)
├── Multiple shots needed: 1 measurement gives 0 or 1, NOT probability
│   Need N shots → count(0)/N ≈ |α|²  (gets better as N→∞)
│
├── Code:
│   qc = QuantumCircuit(1, 1)
│   qc.h(0)           # |+⟩ state
│   qc.measure(0, 0)
│   from qiskit_aer import AerSimulator
│   sim = AerSimulator()
│   result = sim.run(qc, shots=10000).result()
│   counts = result.get_counts()
│   print(counts)  # {'0': ~5000, '1': ~5000}
│
└── Exit check: Verify P(|0⟩)=0.75 for |ψ⟩=(√3/2)|0⟩+(1/2)|1⟩ using 10000 shots.

QC1.3.2  Measuring in Different Bases
├── To measure in X-basis ({|+⟩,|-⟩}):
│   Apply H THEN measure in Z-basis
│   ⟨X⟩ measurement circuit: H → measure
│
├── To measure in Y-basis ({|i⟩,|-i⟩}):
│   Apply S†-H THEN measure in Z-basis
│   ⟨Y⟩ circuit: Sdg → H → measure
│
├── General: rotate to eigenbasis of observable, then Z-measure
│
└── VQE link: H_mol = Σₖ cₖ Pₖ (Pauli strings)
    Each Pₖ needs its own measurement circuit with basis rotations
    ⟨ZZX⟩: no rotation for Z qubits, H for X qubit, then measure all

QC1.3.3  Expectation Values — The VQE Cost Function
├── Estimator primitive (Qiskit 1.0+):
│   from qiskit.primitives import StatevectorEstimator
│   estimator = StatevectorEstimator()
│   job = estimator.run([(qc, observable)])
│   result = job.result()
│   energy = result[0].data.evs  # ⟨ψ|H|ψ⟩
│
├── For shot-based (real hardware):
│   from qiskit_ibm_runtime import EstimatorV2
│   # Same API, just different backend
│
└── Exit check:
    Circuit: Ry(π/3)|0⟩. Observable: Z.
    Compute ⟨Z⟩ analytically: cos(π/3) = 0.5.
    Verify with StatevectorEstimator.

═══════════════════════════════════════════
 GATE TO QC1.4 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know: Z-basis measurement gives |0⟩ or |1⟩ with P=|α|²,|β|²
 □ Know: measurement DESTROYS superposition (irreversible)
 □ Ran 10000-shot simulation, verified P≈0.75 for biased state
 □ Know: X-basis = H then Z-measure; Y-basis = Sdg→H then Z-measure  
 □ Know: VQE measures each Pauli string with basis rotation
 □ Used StatevectorEstimator to compute ⟨Z⟩ for Ry(π/3)|0⟩ = 0.5
═══════════════════════════════════════════
```

---

## Module QC1.4: Entanglement

```
QC1.4.1  Product States vs Entangled States
├── Product state: |ψ⟩ = |ψ_A⟩ ⊗ |ψ_B⟩ (can factor)
│   |00⟩ = |0⟩⊗|0⟩  ← product state
│
├── Entangled state: CANNOT factor
│   |Φ+⟩ = (1/√2)(|00⟩+|11⟩) ← testing:
│   If |Φ+⟩ = (α|0⟩+β|1⟩)⊗(γ|0⟩+δ|1⟩) = αγ|00⟩+αδ|01⟩+βγ|10⟩+βδ|11⟩
│   Need αδ=0, βγ=0, αγ=βδ=1/√2 → contradiction → entangled!
│
└── Schmidt decomposition: |ψ_AB⟩ = Σₖ √λₖ |uₖ⟩|vₖ⟩
    # of non-zero λₖ > 1 → entangled
    Entanglement entropy: S = -Σλₖ log λₖ

QC1.4.2  The Four Bell States
├── |Φ+⟩ = (1/√2)(|00⟩ + |11⟩)  circuit: H(0), CNOT(0,1)
├── |Φ-⟩ = (1/√2)(|00⟩ - |11⟩)  circuit: H(0), CNOT(0,1), Z(0)
├── |Ψ+⟩ = (1/√2)(|01⟩ + |10⟩)  circuit: H(0), CNOT(0,1), X(0)
├── |Ψ-⟩ = (1/√2)(|01⟩ - |10⟩)  circuit: H(0), CNOT(0,1), XZ(0)
│
└── All four: orthonormal, maximally entangled, form a complete basis

QC1.4.3  Why VQE Needs Entanglement
├── H₂ molecule: 2 electrons, each in LCAO orbital
│   Classical product state wavefunction = Hartree-Fock (HF) approximation
│   HF misses electron CORRELATION energy (entanglement effect)
│   HF error for H₂: ~0.054 Hartree (way above chemical accuracy)
│
├── VQE UCCSD ansatz explicitly includes correlated (entangled) terms
│   These model excitations: |0011⟩ ↔ |1100⟩ (orbital swaps)
│   This is exactly what entangling 2-qubit gates create!
│
└── BIO link:
    Aromatic rings (benzene, DNA bases): π-electron system is entangled
    Classical HF can't capture ring current, resonance stabilization properly
    Quantum simulation → correct reaction energies for DNA interstrand crosslinks


QC1.4.4  Quantum Teleportation — Entanglement in Action
├── Goal: transfer |ψ⟩ from Alice to Bob using 1 shared Bell pair + 2 classical bits
│   NO physical qubit is transported — only the STATE is transferred!
│
├── Protocol step-by-step:
│   Step 0: Alice and Bob share |Φ+⟩ = (1/√2)(|00⟩+|11⟩)
│   Step 1: Alice has unknown |ψ⟩ = α|0⟩+β|1⟩. Total 3-qubit state:
│           |ψ⟩⊗|Φ+⟩ = (α|0⟩+β|1⟩)⊗(1/√2)(|00⟩+|11⟩)
│   Step 2: Alice applies CNOT(her ψ qubit → her Bell qubit)
│   Step 3: Alice applies H to her ψ qubit
│   Step 4: Alice measures both her qubits → gets 2 classical bits (m₁m₂)
│   Step 5: Alice sends m₁m₂ to Bob (classical channel)
│   Step 6: Bob applies corrections:
│           If m₂=1: apply X.  If m₁=1: apply Z.
│           Bob now has |ψ⟩ = α|0⟩+β|1⟩ exactly! ✓
│
├── Why this matters:
│   1) Proves entanglement is a RESOURCE (consumed in teleportation)
│   2) Foundation of quantum networks and quantum repeaters
│   3) No information travels faster than light (needs classical bits!)
│
├── Code:
│   qc = QuantumCircuit(3, 2)
│   # Create Bell pair between q1 and q2 (Alice's and Bob's)
│   qc.h(1); qc.cx(1, 2)
│   # Alice's unknown state on q0 (e.g., Ry(π/3)|0⟩)
│   qc.ry(np.pi/3, 0)
│   qc.barrier()
│   # Teleportation protocol
│   qc.cx(0, 1)   # CNOT: q0→q1
│   qc.h(0)        # H on q0
│   qc.measure(0, 0); qc.measure(1, 1)  # Alice measures
│   # Bob's corrections (classically controlled)
│   qc.x(2).c_if(1, 1)   # if m₂=1: X on Bob
│   qc.z(2).c_if(0, 1)   # if m₁=1: Z on Bob
│   # q2 now holds |ψ⟩!
│
└── Exit check: teleport Ry(π/4)|0⟩ and verify Bob's qubit matches.

QC1.4.5  Superdense Coding — Sending 2 Classical Bits via 1 Qubit
├── Reverse of teleportation: send 2 classical bits using 1 shared Bell pair
│   Alice encodes: 00→I, 01→X, 10→Z, 11→XZ on her qubit
│   Sends her qubit to Bob
│   Bob applies CNOT+H and measures → gets 2 bits
│
├── Code:
│   qc = QuantumCircuit(2, 2)
│   qc.h(0); qc.cx(0, 1)   # Bell pair
│   # Alice wants to send "10": applies Z
│   qc.z(0)
│   # Bob decodes
│   qc.cx(0, 1); qc.h(0)
│   qc.measure([0,1], [0,1])
│   # Should get '10' with 100% probability
│
└── Exit check: encode and decode all 4 messages (00,01,10,11). Verify each.

═══════════════════════════════════════════
 GATE TO C2.x — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know: product state = factorable; entangled = cannot factor
 □ Can prove |Φ+⟩ entangled (contradiction argument)
 □ Know all 4 Bell states from memory with circuits
 □ VQE needs entanglement for electron correlation
 □ Hartree-Fock misses correlation energy → VQE fixes this
 □ Teleportation: can explain 6-step protocol and code it
 □ Superdense coding: encode/decode all 4 messages
═══════════════════════════════════════════
```

---

---

## Module QC1.5: NISQ Noise & Error Mitigation

> **PREREQUISITES: QC1.4 gate passed.**
> This is what makes real quantum computing HARD. Simulators are noise-free;
> real hardware is noisy. Every VQE run on real IBM hardware needs these techniques.

```
QC1.5.1  What is NISQ?
├── NISQ = Noisy Intermediate-Scale Quantum
│   "Noisy": every gate has ~0.1-1% error; qubits decohere in ~100μs
│   "Intermediate-Scale": 50-1000 qubits (not enough for error correction)
│   "Quantum": still quantum — can do things classical can't (maybe)
│
├── Error sources:
│   ┌─────────────────┬──────────────────────────────────────┐
│   │ Error type       │ What happens                         │
│   ├─────────────────┼──────────────────────────────────────┤
│   │ Gate error       │ Rotation angle slightly wrong         │
│   │ Decoherence (T1) │ Qubit spontaneously decays |1⟩→|0⟩  │
│   │ Dephasing (T2)   │ Phase information lost randomly       │
│   │ Readout error    │ Measure |0⟩ but get "1" (~1-2%)      │
│   │ Crosstalk        │ Gate on qubit A affects qubit B       │
│   └─────────────────┴──────────────────────────────────────┘
│
├── Impact on VQE:
│   Noisy ⟨H⟩ → optimizer converges to WRONG energy
│   More gates = more noise accumulated → shallow circuits preferred
│   IBM Heron (2024): T1≈300μs, 2-qubit error≈0.5%, max ~50 useful qubits
│
└── Rule: circuit depth × error_rate < 1, otherwise results are garbage

QC1.5.2  Error Mitigation Techniques (NOT error correction!)
├── Readout Error Mitigation:
│   Problem: measure |0⟩ but hardware says "1" sometimes
│   Fix: calibrate confusion matrix M, then apply M⁻¹ to raw counts
│   Code:
│   from qiskit_ibm_runtime import QiskitRuntimeService
│   # Calibration circuits measure |00...0⟩ and |11...1⟩
│   # Build M from calibration data, invert, apply to results
│
├── Zero-Noise Extrapolation (ZNE):
│   Idea: run circuit at noise levels 1x, 2x, 3x → extrapolate to 0x noise
│   How to increase noise: insert identity pairs (CNOT·CNOT=I) → same logic, more noise
│   Extrapolate: fit polynomial to (noise_level, energy) → evaluate at noise=0
│   Best for: VQE energy estimation
│
├── Pauli Twirling:
│   Converts coherent errors → stochastic errors (easier to handle)
│   Insert random Pauli gates before/after each CNOT, undo on other side
│   Average over many random choices → symmetric noise
│
├── Dynamical Decoupling:
│   Insert X pulses during idle times to refocus dephasing
│   Like a quantum echo — cancels slow noise
│   qc.add_calibration('dd', ...)  # Qiskit transpiler can auto-add
│
└── Exit check:
    Run VQE for H=Z on FakeWashingtonV2 (noisy sim):
    1. Without mitigation → energy ≈ -0.85 (should be -1.0)
    2. With readout mitigation → energy ≈ -0.95
    3. Observe the improvement.

═══════════════════════════════════════════
 GATE TO C2.x — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know: NISQ = noisy, 50-1000 qubits, no error correction
 □ Know 5 error sources: gate, T1, T2, readout, crosstalk
 □ Know: circuit depth × error_rate must be < 1
 □ Know readout mitigation: calibrate confusion matrix → invert
 □ Know ZNE: run at 1x,2x,3x noise → extrapolate to 0x
 □ Know Pauli twirling and dynamical decoupling concepts
 □ Ran noisy VQE with and without mitigation → saw improvement
═══════════════════════════════════════════
```

## Module C2.1-C2.3: Qiskit — From Basics to Parameterized Circuits

```
C2.1  Qiskit Environment Setup
├── Install:
│   pip install qiskit qiskit-aer qiskit-ibm-runtime qiskit-nature pyscf
│
├── Basic imports:
│   from qiskit import QuantumCircuit, transpile
│   from qiskit.quantum_info import Statevector, Operator
│   from qiskit_aer import AerSimulator
│
├── First circuit — verify setup:
│   qc = QuantumCircuit(2)
│   qc.h(0); qc.cx(0,1)
│   sv = Statevector(qc)
│   print(sv)  # Statevector([0.707+0j, 0, 0, 0.707+0j])
│
└── If import fails: check Python 3.9+ and pip install --upgrade qiskit

C2.2  Qiskit Simulators — When to Use Each
├── StatevectorSimulator:
│   Stores full 2ⁿ-dim state vector
│   Exact (no noise, no sampling), max ≈ 30 qubits (RAM limited)
│   USE: VQE development and debugging
│
├── QASMSimulator / AerSimulator (shot-based):
│   Samples outcomes N times (like real hardware)
│   Output: counts dict {'00': 512, '11': 488}
│   USE: Realistic simulation, checking measurement circuits
│
├── Noise model simulation:
│   from qiskit_aer.noise import NoiseModel
│   from qiskit_ibm_runtime.fake_provider import FakeWashingtonV2
│   backend = FakeWashingtonV2()
│   noise_model = NoiseModel.from_backend(backend)
│   USE: Approximate real hardware NISQ behavior
│
└── Code — compare exact vs shot:
    qc = QuantumCircuit(1, 1); qc.h(0); qc.measure(0, 0)
    # Exact
    sv = Statevector(QuantumCircuit(1).h(0))
    print(sv.probabilities())  # [0.5, 0.5]
    # Shot-based
    result = AerSimulator().run(qc, shots=1000).result()
    print(result.get_counts())  # {'0': ~500, '1': ~500}

C2.3  Parameterized Circuits — The Core of VQE
├── Why parameterized:
│   VQE ansatz |ψ(θ)⟩ depends on continuous parameters θ = (θ₁,...,θₙ)
│   We need to OPTIMIZE these parameters → need to evaluate for many θ
│   Parameterized circuit = template, fill in θ at runtime
│
├── Building parameterized circuits:
│   from qiskit.circuit import ParameterVector, Parameter
│   theta = ParameterVector('θ', 4)  # 4 parameters
│   qc = QuantumCircuit(2)
│   qc.ry(theta[0], 0)
│   qc.ry(theta[1], 1)
│   qc.cx(0, 1)
│   qc.ry(theta[2], 0)
│   qc.ry(theta[3], 1)
│
├── Assigning parameters:
│   bound_qc = qc.assign_parameters({theta: [0.1, 0.2, 0.3, 0.4]})
│
├── Parameter-shift gradient rule:
│   ∂⟨E⟩/∂θₖ = [⟨E(θₖ + π/2)⟩ - ⟨E(θₖ - π/2)⟩] / 2
│   This gives EXACT gradient on real hardware (no finite difference needed!)
│
└── MINI VQE (Toy Problem):
    # Find ground state of H = Z (eigenvalue -1 at |1⟩)
    theta = Parameter('θ')
    qc = QuantumCircuit(1)
    qc.ry(theta, 0)   # ansatz: Ry(θ)|0⟩
    # Cost: ⟨Z⟩ = ⟨0|Ry(-θ)ZRy(θ)|0⟩ = cos(θ)
    # Minimize cos(θ) → θ=π → |ψ(π)⟩=|1⟩, ⟨Z⟩=-1 ✓
    # THIS IS VQE LOGIC IN ITS SIMPLEST FORM


C2.4  Real IBM Hardware Access — Running on Actual Quantum Computers
├── IBM Quantum account setup:
│   1. Create account: quantum.ibm.com
│   2. Get API token from dashboard
│   3. Save token:
│      from qiskit_ibm_runtime import QiskitRuntimeService
│      QiskitRuntimeService.save_account(channel="ibm_quantum", token="YOUR_TOKEN")
│
├── Choosing a backend:
│   service = QiskitRuntimeService()
│   backends = service.backends()
│   # Filter by qubit count and queue length
│   backend = service.least_busy(min_num_qubits=4)
│   print(backend.name, backend.num_qubits)
│
├── Transpilation — matching circuit to hardware:
│   from qiskit import transpile
│   # Hardware has limited connectivity (not all qubits connect)
│   qc_transpiled = transpile(qc, backend=backend, optimization_level=3)
│   # optimization_level: 0=none, 1=light, 2=medium, 3=heavy (best for VQE)
│   print(f"Depth: {qc_transpiled.depth()}, CNOTs: {qc_transpiled.count_ops().get('cx',0)}")
│
├── Running on real hardware:
│   from qiskit_ibm_runtime import EstimatorV2, SamplerV2
│   estimator = EstimatorV2(backend)
│   job = estimator.run([(qc_transpiled, hamiltonian)])
│   result = job.result()
│   # WARNING: job may queue for minutes/hours depending on traffic
│
├── Key differences from simulator:
│   ┌──────────────────┬────────────────┬─────────────────┐
│   │                   │ Simulator       │ Real Hardware    │
│   ├──────────────────┼────────────────┼─────────────────┤
│   │ Speed             │ Instant         │ Mins-hours queue │
│   │ Noise             │ None (default)  │ Always present   │
│   │ Max qubits        │ ~30 (RAM)       │ 127-1000+        │
│   │ Gate set          │ Any             │ Native only      │
│   │ Connectivity      │ All-to-all      │ Limited topology │
│   └──────────────────┴────────────────┴─────────────────┘
│
└── Exit check:
    Run Bell state circuit on real IBM hardware.
    Compare counts {'00':~50%, '11':~50%} with error ~2-5%.
    Observe noise vs perfect simulator result.

═══════════════════════════════════════════
 GATE TO QC2.1 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Qiskit installed; Bell state circuit runs
 □ Know 3 simulators: Statevector, Aer (shots), Noise model
 □ Built parameterized circuit with ParameterVector
 □ Used assign_parameters() to bind values
 □ Know parameter-shift: ∂E/∂θ = [E(θ+π/2)-E(θ-π/2)]/2
 □ Mini VQE: H=Z, ansatz=Ry(θ), θ=π gives ⟨Z⟩=-1
 □ Set up IBM Quantum account + saved API token
 □ Transpiled circuit for real backend (optimization_level=3)
 □ Ran Bell state on real hardware and observed noise
═══════════════════════════════════════════
```

---

# PHASE 3: QC ALGORITHMS (Weeks 19-26)

---

## Module QC2.1: Variational Principle — The Foundation of VQE

```
QC2.1.1  Rayleigh-Ritz Variational Principle
├── Statement:
│   For ANY normalized state |ψ⟩: ⟨ψ|H|ψ⟩ ≥ E₀
│   where E₀ is the TRUE ground state energy
│
├── Proof:
│   Expand |ψ⟩ in energy eigenstates: |ψ⟩ = Σₙ cₙ|Eₙ⟩
│   ⟨ψ|H|ψ⟩ = Σₙ |cₙ|² Eₙ   (H acts on each eigenstate)
│   Since E₀ ≤ Eₙ for all n:
│   Σₙ |cₙ|² Eₙ ≥ Σₙ |cₙ|² E₀ = E₀ · Σₙ|cₙ|² = E₀  ✓
│   (Used: Σₙ|cₙ|² = 1 by normalization)
│
├── Key consequence:
│   Minimize ⟨ψ(θ)|H|ψ(θ)⟩ over θ → approach E₀ from above
│   More expressive ansatz → closer to true E₀
│
└── Analogy:
    You're guessing what the lowest point on a landscape is.
    You can check height at any point.
    Keep picking points until you find the lowest you can reach.
    VQE is this process on quantum hardware.

QC2.1.2  Ansatz Design — What Makes a Good Ansatz
├── Requirements:
│   1. Expressibility: can reach states close to true ground state
│   2. Efficiency: few gates, low circuit depth (NISQ constraint)
│   3. Trainability: gradients don't vanish (avoid barren plateaus)
│
├── Hardware-Efficient Ansatz (HEA):
│   Layer structure: [Ry]⊗n → CNOT_entangling → [Ry]⊗n → ...
│   Number of parameters: O(n × depth)
│   Pros: Native gate set, low depth
│   Cons: May not capture true ground state chemistry
│
├── UCCSD (Unitary Coupled Cluster Singles Doubles):
│   Chemistry-motivated: models single/double electron excitations
│   For H₂ (4 spin-orbitals): 1 parameter, 1 CNOT
│   For larger molecules: exponentially many terms (classically intractable)
│   Qiskit: from qiskit_nature.second_q.circuit.library import UCCSD
│
└── Code (EfficientSU2 built-in):
    from qiskit.circuit.library import EfficientSU2
    ansatz = EfficientSU2(num_qubits=4, reps=2)
    print(f"Parameters: {ansatz.num_parameters}")  # shows count
    ansatz.decompose().draw('mpl')

QC2.1.3  Classical Optimizers in VQE
├── COBYLA (Constrained Optimization By Linear Approximations):
│   Derivative-free, robust to shot noise
│   Best for: small circuits, hardware runs
│   Starting choice for VQE
│   from scipy.optimize import minimize
│   minimize(cost_fn, x0, method='COBYLA')
│
├── SPSA (Simultaneous Perturbation Stochastic Approximation):
│   Uses only 2 cost evaluations per gradient step (vs 2n for param-shift)
│   Good for hardware with many parameters
│
├── Adam (gradient-based):
│   Uses param-shift for gradients
│   Faster convergence than COBYLA for smooth landscapes
│   from qiskit_algorithms.optimizers import ADAM
│
├── Barren plateau problem:
│   For random initial parameters on n-qubit circuits:
│   ∂E/∂θₖ → 0 exponentially fast as n increases
│   Gradient is exponentially small → training fails!
│   Solutions: layerwise training, problem-specific init, shallow circuits
│
└── Exit check (VQE Gate Exam):
    H = 0.5·(Z⊗Z) + 0.5·(X⊗I)  [create as SparsePauliOp]
    Ansatz: Ry(θ₁)⊗Ry(θ₂) → CNOT
    Find ground state energy using COBYLA. Compare to exact eigvals.

═══════════════════════════════════════════
 GATE TO QC3.1 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know variational principle: ⟨ψ|H|ψ⟩ ≥ E₀ for all |ψ⟩
 □ Know 3 ansatz needs: expressibility, efficiency, trainability
 □ Know HEA vs UCCSD tradeoffs
 □ Know COBYLA vs SPSA vs Adam
 □ Know barren plateau problem and mitigations
 □ Ran full VQE on 2-qubit Hamiltonian
═══════════════════════════════════════════
```

---

## Module QC3.0: Deutsch-Jozsa Algorithm — First Quantum Speedup

> **PREREQUISITES: QC2.1 gate passed.**
> This is historically the FIRST algorithm proving quantum > classical.
> Understand this before Grover — it introduces oracle + interference pattern.

```
QC3.0.1  The Problem — Constant vs Balanced
├── Given: a function f:{0,1}ⁿ→{0,1}
│   CONSTANT: f(x)=0 for ALL x, or f(x)=1 for ALL x
│   BALANCED: f(x)=0 for exactly HALF of inputs, 1 for the other half
│   PROMISE: f is either constant OR balanced (nothing else)
│   TASK: determine which one
│
├── Classical: worst case need 2ⁿ/2+1 queries (check >half the inputs)
│   For n=10: up to 513 evaluations
│   Quantum: Deutsch-Jozsa needs EXACTLY 1 query!  (exponential speedup)
│
└── This is the canonical "quantum computers CAN be faster" proof

QC3.0.2  The Algorithm — Step by Step
├── Circuit:
│   |0⟩⊗ⁿ|1⟩ → H⊗(n+1) → Uf (oracle) → H⊗ⁿ (on first n qubits) → Measure
│
├── WORKED for n=2:
│   Step 1: Start with |00⟩|1⟩ = |001⟩
│   Step 2: Apply H⊗³ → uniform superposition ⊗ |-⟩ on ancilla
│   Step 3: Apply oracle Uf:
│           Uf|x⟩|-⟩ = (-1)^f(x)|x⟩|-⟩  (phase kickback!)
│           If CONSTANT (f=0): all phases same → H brings back to |00⟩
│           If BALANCED: phases cancel → H gives non-|00⟩ → some qubit = 1
│   Step 4: Measure first n qubits
│           ALL zeros → CONSTANT
│           ANY non-zero → BALANCED
│
├── Why it works: INTERFERENCE
│   Constant f → all amplitudes interfere constructively at |0...0⟩
│   Balanced f → amplitudes at |0...0⟩ cancel to zero (destructive)
│
├── Code (2-qubit, balanced oracle f(x)=x₁ XOR x₂):
│   qc = QuantumCircuit(3, 2)  # 2 input qubits + 1 ancilla
│   qc.x(2); qc.barrier()     # ancilla to |1⟩
│   qc.h([0,1,2])             # Step 2: superposition
│   qc.barrier()
│   # Oracle for f(x)=x₁⊕x₂ (balanced):
│   qc.cx(0, 2); qc.cx(1, 2)
│   qc.barrier()
│   qc.h([0,1])               # Step 4: H on input qubits only
│   qc.measure([0,1], [0,1])
│   # Result: NOT '00' → balanced ✓
│
└── Exit check:
    1. Build constant oracle (f=0): just identity. Should measure '00'.
    2. Build balanced oracle (f=x₀): CNOT(0,2). Should measure non-'00'.
    3. Verify both in Qiskit with 1000 shots.

═══════════════════════════════════════════
 GATE TO QC3.1 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know the constant vs balanced promise problem
 □ Classical needs 2ⁿ/2+1 queries; Deutsch-Jozsa needs 1
 □ Know phase kickback: Uf|x⟩|-⟩ = (-1)^f(x)|x⟩|-⟩
 □ Constant → all-zeros measurement; balanced → non-zero
 □ Coded both constant and balanced oracles, verified
 □ Understand: interference (constructive/destructive) is the engine
═══════════════════════════════════════════
```

---

## Module QC3.1: Grover's Algorithm

```
QC3.1.1  Unstructured Search — Classical vs Quantum
├── Problem: N items in unsorted database, 1 marked item (f(x)=1)
│   Classical: check items one by one → O(N) queries
│   Quantum: Grover → O(√N) queries → quadratic speedup
│
├── For N=10⁶: classical=10⁶ checks, Grover=1000 checks
│   For human genome (N=3×10⁹ bp k-mers): 3×10⁹ vs 55,000 queries
│
└── This is THE proven quantum speedup with widest applicability

QC3.1.2  The Oracle — Marking Target States
├── Phase oracle: Uf|x⟩ = (-1)^f(x)|x⟩
│   f(x)=0 (not target): state unchanged
│   f(x)=1 (target x*): amplitude gets multiplied by -1
│
├── Oracle for |11⟩ = CZ gate (marks |11⟩ with -1 phase):
│   from qiskit import QuantumCircuit
│   oracle = QuantumCircuit(2)
│   oracle.cz(0, 1)   # flips phase of |11⟩ only
│
├── Phase kickback trick:
│   Target qubit in |-⟩ = (1/√2)(|0⟩-|1⟩)
│   CNOT from control: |c⟩|-⟩ → (-1)^c|c⟩|-⟩
│   Phase "kicks back" to control qubit!
│
└── Generalizing: oracle for string |w₁w₂...wₙ⟩:
    Apply X to qubits where wᵢ=0, then n-controlled-Z, then X again

QC3.1.3  The Diffuser — Inversion About Mean
├── W = 2|ψ⟩⟨ψ| - I  where |ψ⟩ = H⊗ⁿ|0...0⟩ (uniform superposition)
├── Effect: reflects about mean amplitude
│   Small amplitudes pushed more negative
│   Large amplitudes boosted
│   Iterating: target amplitude grows, others shrink
│
├── Circuit for n qubits:
│   H⊗ⁿ → X⊗ⁿ → n-controlled-Z → X⊗ⁿ → H⊗ⁿ
│
└── Code (Grover on 2 qubits, target |11⟩):
    qc = QuantumCircuit(2)
    # Superposition
    qc.h([0,1])
    # Oracle (CZ marks |11⟩)
    qc.cz(0,1)
    # Diffuser
    qc.h([0,1]); qc.x([0,1])
    qc.h(1); qc.cx(0,1); qc.h(1)   # 2-qubit controlled-Z
    qc.x([0,1]); qc.h([0,1])
    # Measure
    qc.measure_all()
    # For N=4, 1 iteration. Should get |11⟩ with ~100% probability

QC3.1.4  BIO Application — Genomic k-mer Search
├── Problem: find specific DNA k-mer in genome database
│   k=20 bp: N=4²⁰ ≈ 10¹² possible k-mers → massive search space
│
├── Quantum approach:
│   Encode each k-mer as a quantum state using 2k qubits (2 bits/nucleotide)
│   A=00, T=11, C=01, G=10
│   Oracle: phase-mark the target k-mer
│   Grover iterations: O(√N) ≈ O(√4^k) = O(2^k) — better than classical 4^k
│
└── Active research: major challenge = QRAM (quantum random access memory)
    Loading the database into superposition is itself an O(N) operation
    Quantum advantage only realized if QRAM loading is efficient


---

## Module QC3.2: Quantum Fourier Transform (QFT)

> **PREREQUISITES: QC3.1 gate passed.**
> QFT is the quantum version of the classical Discrete Fourier Transform.
> It's the ENGINE inside Shor's algorithm and QPE.

```
QC3.2.1  What QFT Does
├── Classical DFT: converts time-domain → frequency-domain  (O(N log N) = FFT)
│   Quantum QFT: same transformation on quantum amplitudes  (O(n²) gates, n=log₂N)
│
├── Definition: QFT|j⟩ = (1/√N) Σₖ e^(2πijk/N) |k⟩
│   Maps computational basis → Fourier basis
│   N=2ⁿ amplitudes transformed using only O(n²) gates!
│
├── For n=1: QFT = H (Hadamard IS the 1-qubit Fourier transform!)
│   For n=2: QFT = H⊗I · controlled-phase(π/2) · I⊗H · SWAP
│
├── Circuit structure (n qubits):
│   For each qubit j (from top to bottom):
│     1. Apply H to qubit j
│     2. Apply controlled-Rk gates from qubits j+1,...,n-1
│        where Rk = phase gate with angle 2π/2^k
│   3. Reverse qubit order (SWAP gates)
│
└── Code:
    from qiskit.circuit.library import QFT
    qft = QFT(num_qubits=3)
    qft.decompose().draw('mpl')

QC3.2.2  Why QFT Matters
├── Inside Shor's algorithm: QFT finds the PERIOD of modular exponentiation
│   Period finding → factoring → breaks RSA encryption
│   (We won't implement Shor fully, but understanding QFT is essential)
│
├── Inside QPE: QFT converts phase information → binary representation
│   QPE uses INVERSE QFT at the end to read out the phase
│
└── Exit check:
    Build 3-qubit QFT circuit manually (H + controlled-phase gates + SWAPs).
    Apply to |101⟩ and verify output matches Qiskit's QFT library.

═══════════════════════════════════════════
 GATE TO QC3.3 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know: QFT = quantum version of DFT; O(n²) gates vs O(N log N) classical
 □ Know: 1-qubit QFT = Hadamard
 □ Can describe QFT circuit: H + controlled-Rk + SWAPs
 □ Know: QFT is inside Shor's and QPE
 □ Built 3-qubit QFT circuit in Qiskit
═══════════════════════════════════════════
```

---

## Module QC3.3: Quantum Phase Estimation (QPE)

> **PREREQUISITES: QC3.2 gate passed.**
> QPE estimates eigenvalues of unitary operators.
> For quantum chemistry: QPE gives EXACT ground state energy (unlike VQE which approximates).
> VQE is for NISQ (noisy, now); QPE is for fault-tolerant (clean, future).

```
QC3.3.1  The Problem QPE Solves
├── Given: unitary U with eigenvector |u⟩ and eigenvalue e^(2πiφ)
│   TASK: estimate the phase φ (a number between 0 and 1)
│
├── Why this matters for chemistry:
│   Hamiltonian H has ground state |E₀⟩ with energy E₀
│   U = e^(-iHt) is unitary with eigenvalue e^(-iE₀t)
│   QPE extracts E₀ from the phase! → exact ground state energy
│   VQE: approximate E₀ (good for NISQ)
│   QPE: exact E₀ (needs fault-tolerant QC → future)
│
└── Precision: t ancilla qubits → phase accuracy 1/2^t

QC3.3.2  The Algorithm
├── Circuit:
│   t ancilla qubits (initialized |0⟩) + eigenstate register |u⟩
│   Step 1: H⊗t on ancillas (create superposition)
│   Step 2: Controlled-U^(2^k) from ancilla k to eigenstate
│   Step 3: Inverse QFT on ancillas
│   Step 4: Measure ancillas → binary representation of φ
│
├── WORKED for 1-ancilla:
│   U = Z (eigenvalues ±1 → phases 0 and 0.5)
│   |u⟩ = |1⟩ (eigenvalue -1 = e^(2πi·0.5) → phase = 0.5)
│   QPE should return |1⟩ on ancilla (binary 0.1₂ = 0.5)
│
├── Code:
│   from qiskit.circuit.library import PhaseEstimation
│   import numpy as np
│   # Estimate phase of Z gate acting on |1⟩
│   Z_gate = QuantumCircuit(1); Z_gate.z(0)
│   qpe = PhaseEstimation(num_evaluation_qubits=3, unitary=Z_gate.to_gate())
│   # Prepend: set eigenstate to |1⟩
│   full = QuantumCircuit(4, 3)
│   full.x(3)  # eigenstate = |1⟩
│   full.compose(qpe, inplace=True)
│   full.measure(range(3), range(3))
│   # Should measure '100' (binary 0.100 = 0.5 → phase = 0.5)
│
└── Exit check:
    1. QPE on Z with |1⟩: should get phase=0.5 (energy eigenvalue = -1)
    2. QPE on T gate (phase=1/8): should get binary 0.001 with 3 ancillas
    3. Understand: more ancilla qubits → higher precision

═══════════════════════════════════════════
 GATE COMPLETE — QC3.3 checked:
═══════════════════════════════════════════
 □ Know: QPE estimates eigenvalue phase of unitary U
 □ Know: chemistry use: U=e^(-iHt) → extracts E₀
 □ Know: VQE = NISQ approximate; QPE = fault-tolerant exact
 □ Know circuit: H ancillas → controlled-U^2k → inverse QFT → measure
 □ Ran QPE on Z gate → phase=0.5 verified
 □ Know: t ancilla qubits → precision 1/2^t
═══════════════════════════════════════════
```

═══════════════════════════════════════════
 ⭐ MASTER QC GATE — PART 3 COMPLETE ⭐
═══════════════════════════════════════════
 □ QC1.1: Qubit states, Bloch sphere, phase
 □ QC1.2: All gates, universal set, CNOT
 □ QC1.3: Measurement, basis rotation, expectation
 □ QC1.4: Entanglement, Bell states, teleportation, superdense coding
 □ QC1.5: NISQ noise, error mitigation (ZNE, readout, twirling)
 □ C2.x: Qiskit, simulators, parameterized circuits, real hardware
 □ QC2.1: Variational principle, ansatz, optimizers
 □ QC3.0: Deutsch-Jozsa (constant vs balanced, 1 query)
 □ QC3.1: Grover oracle, diffuser, bio app
 □ QC3.2: Quantum Fourier Transform
 □ QC3.3: Quantum Phase Estimation (exact eigenvalues)
═══════════════════════════════════════════
```

---

# ✅ COMPLETE TO-DO LIST — PART 3 (QC THEORY + ALGORITHMS)

## QC1.1 Qubits
- [ ] Classical bit vs qubit table
- [ ] |ψ⟩=α|0⟩+β|1⟩; normalization check
- [ ] Physical realizations (4 types)
- [ ] Statevector in Qiskit
- [ ] Bloch sphere 6 states memorized
- [ ] Find θ,φ for given state by hand
- [ ] Gates = rotations on Bloch
- [ ] Global vs relative phase
- [ ] Interference demo: H|+⟩=|0⟩ vs H|-⟩=|1⟩
- [ ] QC1.1 GATE ✓

## QC1.2 Gates
- [ ] X,Y,Z,H,S,T matrices from memory
- [ ] S²=Z, T²=S, H²=I
- [ ] S|+⟩=|+i⟩ by hand
- [ ] Ry(θ), Rz(θ) forms; Ry(π/3)|0⟩ by hand
- [ ] Universal set {H,T,CNOT}
- [ ] CNOT truth table; 4×4 matrix
- [ ] Bell state step-by-step
- [ ] CZ, SWAP gates
- [ ] Hardware error rates table
- [ ] 3 circuits built and drawn
- [ ] QC1.2 GATE ✓

## QC1.3 Measurement
- [ ] Z-basis: P(0)=|α|², destroys superposition
- [ ] 10000-shot verification
- [ ] X-basis = H→Z; Y-basis = Sdg→H→Z
- [ ] VQE Pauli string measurement
- [ ] StatevectorEstimator: ⟨Z⟩ for Ry(π/3)|0⟩
- [ ] QC1.3 GATE ✓

## QC1.4 Entanglement
- [ ] Product vs entangled; prove |Φ+⟩ entangled
- [ ] Schmidt decomposition concept
- [ ] 4 Bell states from memory + circuits
- [ ] VQE needs entanglement (correlation energy)
- [ ] BIO: DNA π-electrons
- [ ] QC1.4 GATE ✓

## C2.x Qiskit
- [ ] Install + first circuit
- [ ] 3 simulator types
- [ ] Parameterized circuit with ParameterVector
- [ ] assign_parameters()
- [ ] Parameter-shift gradient rule
- [ ] Mini VQE: H=Z, Ry(θ), θ=π →⟨Z⟩=-1
- [ ] C2.x GATE ✓

## QC2.1 VQE
- [ ] Variational principle proof
- [ ] HEA vs UCCSD ansatz
- [ ] EfficientSU2 from library
- [ ] COBYLA, SPSA, Adam optimizers
- [ ] Barren plateau problem
- [ ] Full 2-qubit VQE run
- [ ] QC2.1 GATE ✓

## QC3.1 Grover
- [ ] Classical O(N) vs quantum O(√N)
- [ ] Phase oracle construction
- [ ] Phase kickback
- [ ] Diffuser circuit
- [ ] 2-qubit Grover code
- [ ] Optimal iterations ≈(π/4)√N
- [ ] BIO: k-mer search + QRAM challenge
- [ ] QC3.1 GATE ✓

---


## QC1.4 (Extended) Teleportation + Superdense
- [ ] Teleportation 6-step protocol — explain and code
- [ ] Superdense coding: encode/decode 00,01,10,11
- [ ] QC1.4 Extended GATE ✓

## QC1.5 NISQ Noise
- [ ] NISQ definition: noisy, 50-1000 qubits
- [ ] 5 error sources: gate, T1, T2, readout, crosstalk
- [ ] depth × error_rate < 1 rule
- [ ] Readout error mitigation (confusion matrix)
- [ ] Zero-noise extrapolation (ZNE)
- [ ] Pauli twirling concept
- [ ] Dynamical decoupling concept
- [ ] Noisy VQE run with vs without mitigation
- [ ] QC1.5 GATE ✓

## C2.4 Real Hardware
- [ ] IBM Quantum account + API token
- [ ] Transpilation with optimization_level=3
- [ ] Run Bell state on real hardware
- [ ] Compare noisy real vs perfect simulator
- [ ] C2.4 GATE ✓

## QC3.0 Deutsch-Jozsa
- [ ] Constant vs balanced problem
- [ ] Classical 2^(n-1)+1 vs quantum 1 query
- [ ] Phase kickback: Uf|x⟩|-⟩ = (-1)^f(x)|x⟩|-⟩
- [ ] All-zeros = constant; non-zero = balanced
- [ ] Code constant + balanced oracles
- [ ] QC3.0 GATE ✓

## QC3.2 QFT
- [ ] QFT = quantum DFT; O(n²) gates
- [ ] 1-qubit QFT = Hadamard
- [ ] Circuit: H + controlled-Rk + SWAPs
- [ ] QFT inside Shor's and QPE
- [ ] Built 3-qubit QFT
- [ ] QC3.2 GATE ✓

## QC3.3 QPE
- [ ] QPE estimates eigenvalue phase of U
- [ ] Chemistry: U=e^(-iHt) → extracts E₀
- [ ] VQE = NISQ approx; QPE = exact (future)
- [ ] Circuit: H→controlled-U^2k→inv QFT→measure
- [ ] QPE on Z gate → phase=0.5
- [ ] Precision: t ancillas → 1/2^t
- [ ] QC3.3 GATE ✓

## ⭐ MASTER SIGN-OFF — PART 3

- [ ] All 11 module gates passed (QC1.1→QC3.3)
- [ ] Can build circuits in Qiskit + ran on real IBM hardware
- [ ] Mini VQE + full VQE completed
- [ ] Deutsch-Jozsa + Grover coded and tested
- [ ] QFT + QPE understood and coded
- [ ] Know NISQ noise + error mitigation
- [ ] **READY FOR PART 4 — QUANTUM CHEMISTRY 🚀**

---

## 🌳 Part 3 Module Dependency Tree

```mermaid
graph TD
    P1["Part 1: Math ✅"] --> P2["Part 2: Physics ✅"]
    P2 --> QC11["QC1.1: Qubits"]
    QC11 --> QC12["QC1.2: Gates"]
    QC12 --> QC13["QC1.3: Measurement"]
    QC13 --> QC14["QC1.4: Entanglement"]
    QC14 --> C2["C2.x: Qiskit"]
    C2 --> QC21["QC2.1: VQE"]
    QC14 --> QC15["QC1.5: NISQ Noise"]
    QC15 --> C2["C2.x: Qiskit + Hardware"]
    C2 --> QC21["QC2.1: VQE"]
    QC21 --> QC30["QC3.0: Deutsch-Jozsa"]
    QC30 --> QC31["QC3.1: Grover"]
    QC31 --> QC32["QC3.2: QFT"]
    QC32 --> QC33["QC3.3: QPE"]
    
    QC12 -.->|"CNOT creates"| QC14
    QC15 -.->|"noise model"| C2
    QC14 -.->|"correlation"| QC21
    QC32 -.->|"engine for"| QC33
    
    QC33 --> P4["Part 4: Quantum Chemistry"]
    QC31 --> BIO["BIO: k-mer Search"]
    
    style P1 fill:#2d6a4f,color:#fff
    style P2 fill:#2d6a4f,color:#fff
    style QC11 fill:#1d3557,color:#fff
    style QC12 fill:#1d3557,color:#fff
    style QC13 fill:#457b9d,color:#fff
    style QC14 fill:#457b9d,color:#fff
    style C2 fill:#e76f51,color:#fff
    style QC21 fill:#f4a261,color:#000
    style QC15 fill:#457b9d,color:#fff
    style QC30 fill:#e9c46a,color:#000
    style QC31 fill:#e9c46a,color:#000
    style QC32 fill:#e9c46a,color:#000
    style QC33 fill:#e9c46a,color:#000
    style P4 fill:#264653,color:#fff
    style BIO fill:#2a9d8f,color:#fff
```

