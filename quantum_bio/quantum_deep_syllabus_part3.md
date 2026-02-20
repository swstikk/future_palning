# ⚛️ Quantum Bioinformatics — Deep Chapter-wise Syllabus PART 3
## Phase 2: QC Theory + Qiskit (Weeks 11-18) + Phase 3: Algorithms (Weeks 19-26)

---

# PHASE 2: QUANTUM COMPUTING THEORY + QISKIT (Weeks 11-18)

---

## Module QC1.1: Qubits & Quantum States

```
QC1.1.1  The Qubit — Physical vs Mathematical
├── What to learn:
│   ├── Classical bit: 0 OR 1 (definite)
│   ├── Qubit: α|0⟩ + β|1⟩ (BOTH, with complex weights)
│   │   Before measurement: genuinely in superposition
│   │   After measurement: collapses to |0⟩ or |1⟩
│   ├── Physical realizations:
│   │   Superconducting qubits (IBM, Google): Josephson junction
│   │   Trapped ions (IonQ): single ionized atoms
│   │   Photonic qubits (PsiQuantum): photon polarization
│   └── Mathematical: state lives in ℂ² with |α|²+|β|²=1
│
├── Normalization constraint:
│   |α|² + |β|² = 1
│   This is NOT optional — it ensures probabilities sum to 1
│   Geometrically: state lives on the surface of unit sphere in ℂ²
│
└── Code:
    from qiskit.quantum_info import Statevector
    psi = Statevector([1/np.sqrt(2), 1j/np.sqrt(2)])
    print(psi.is_valid())  # True (normalized)
    print(psi.probabilities())  # [0.5, 0.5]

QC1.1.2  Bloch Sphere — The Qubit Visualization Tool
├── Any single qubit state:
│   |ψ⟩ = cos(θ/2)|0⟩ + e^(iφ)·sin(θ/2)|1⟩
│   θ ∈ [0, π]: polar angle (latitude)
│   φ ∈ [0, 2π): azimuthal angle (longitude)
│
├── Key locations on sphere:
│   θ=0:    |0⟩ = north pole
│   θ=π:    |1⟩ = south pole
│   θ=π/2, φ=0:    |+⟩ = +x axis
│   θ=π/2, φ=π:    |-⟩ = -x axis
│   θ=π/2, φ=π/2:  |i⟩ = +y axis
│   θ=π/2, φ=3π/2: |-i⟩= -y axis
│
├── Gates as rotations:
│   X gate = π rotation about x-axis (|0⟩ ↔ |1⟩)
│   Z gate = π rotation about z-axis (|+⟩ ↔ |-⟩)
│   H gate = π rotation about x+z axis diagonal
│   Rx(θ): rotation by θ about x-axis
│
└── Code:
    from qiskit.visualization import plot_bloch_vector
    # |+⟩ state = equator, φ=0
    plot_bloch_vector([1, 0, 0])  # [x, y, z] Bloch coordinates

QC1.1.3  Phase: Global vs Relative
├── Global phase: e^(iα)|ψ⟩ is physically identical to |ψ⟩
│   You can never detect global phase by any measurement
│
├── Relative phase: PHYSICALLY OBSERVABLE
│   |+⟩ = (1/√2)(|0⟩ + |1⟩)
│   |-⟩ = (1/√2)(|0⟩ - |1⟩)
│   Same |α|²=|β|²=½, but DIFFERENT states — they interfere differently!
│
├── Relative phase matters for:
│   Quantum interference (Grover uses this!)
│   Gate action: Z gate introduces relative phase (|0⟩→|0⟩, |1⟩→-|1⟩)
│
└── Exit check:
    Draw these states on Bloch sphere: |0⟩, |1⟩, |+⟩, |-⟩, (1/√2)(|0⟩+i|1⟩).
    Which have same probabilities but different physical states?
```

---

## Module QC1.2: Quantum Gates — Complete Reference

```
QC1.2.1  Single-Qubit Gates (Matrix Reference)
├── Pauli gates (must memorize — NO calculator):
│   X = [[0,1],[1,0]]     Y = [[0,-i],[i,0]]     Z = [[1,0],[0,-1]]
│
├── Phase gates:
│   S = [[1,0],[0,i]]  (T-gate in two applications: S²=Z)
│   T = [[1,0],[0,e^(iπ/4)]]  (π/8 gate)
│
├── Hadamard:
│   H = (1/√2)[[1,1],[1,-1]]
│
├── Rotation gates (continuous, key for VQE):
│   Rx(θ) = [[cos(θ/2), -i·sin(θ/2)], [-i·sin(θ/2), cos(θ/2)]]
│   Ry(θ) = [[cos(θ/2), -sin(θ/2)],   [sin(θ/2),    cos(θ/2)]]
│   Rz(θ) = [[e^(-iθ/2),   0],         [0,       e^(iθ/2)]]
│
│   Ry(0)=I, Ry(π)=iY, Ry(π/2) = H·phase  ← nearly Hadamard
│
└── Universal set: {H, T, CNOT} — can approximate ANY unitary gate!

QC1.2.2  Two-Qubit Gates
├── CNOT (CX gate):
│   Matrix (4×4, computational basis |00⟩,|01⟩,|10⟩,|11⟩):
│   [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]
│   Rule: if control=|0⟩ → target unchanged, if control=|1⟩ → target flipped
│   CNOT |00⟩ = |00⟩,  CNOT |10⟩ = |11⟩,  CNOT |11⟩ = |10⟩
│
├── CZ gate:
│   Applies Z to target ONLY if control = |1⟩
│   Symmetric: CZ|11⟩ = -|11⟩, all others unchanged
│
├── SWAP gate:
│   Exchanges two qubits: SWAP|01⟩ = |10⟩
│   Equivalent to 3 CNOTs: SWAP = CNOT₀₁·CNOT₁₀·CNOT₀₁
│
├── Why 2-qubit gates are expensive:
│   On real hardware: 2-qubit gates have 10x higher error than 1-qubit
│   Gate fidelity: 1-qubit ≈ 99.9%, 2-qubit ≈ 99.0-99.5% (IBM 2024)
│   VQE circuit design = minimize 2-qubit gate count!
│
└── Code (Bell state from |00⟩):
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(2)
    qc.h(0)          # H on qubit 0: |00⟩ → |+0⟩
    qc.cx(0, 1)      # CNOT: |+0⟩ → |Φ+⟩ = (|00⟩+|11⟩)/√2
    print(Statevector(qc).data)  # [1/√2, 0, 0, 1/√2]

QC1.2.3  Circuit Notation & Building Circuits in Qiskit
├── Reading circuit diagrams:
│   Left to right = time
│   Horizontal lines = qubit "wires"
│   Boxes = single-qubit gates
│   Vertical lines with dot/circle = controlled gates
│
├── Qiskit circuit building:
│   qc = QuantumCircuit(n_qubits, n_classical_bits)
│   qc.h(0)          # Hadamard on qubit 0
│   qc.rx(theta, 1)  # Rx rotation on qubit 1
│   qc.cx(0, 1)      # CNOT: control=0, target=1
│   qc.measure(0, 0) # measure qubit 0 → classical bit 0
│   qc.draw('mpl')   # visualize
│
└── Exit check:
    Build circuit: |00⟩ → (1/2)(|00⟩+|01⟩+|10⟩+|11⟩) using minimum gates.
    Solution: H⊗H. Verify statevector shows equal amplitudes.
```

---

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
```

---

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
```
