# ⚛️ Quantum Bioinformatics — Deep Chapter-wise Syllabus PART 2
## Phase 1 Physics Track (Ph1.1 – Ph2.3)

---

# PHASE 1 — PHYSICS TRACK (Weeks 3-10, Parallel with Math)

---

## Module Ph1.1: Classical Mechanics — Energy & Hamiltonian Intuition

```
Ph1.1.1  Work, Energy, Conservative Forces
├── What to learn:
│   ├── Work: W = F · d (force times displacement, dot product)
│   ├── Kinetic energy: KE = ½mv² = p²/2m  ← memorize both forms
│   ├── Potential energy: PE = mgh (gravity), PE = ½kx² (spring)
│   ├── Total mechanical energy: E = KE + PE = constant (conservation)
│   └── Conservative force: work done = -(change in PE), path-independent
│
├── Momentum form (critical for QM):
│   KE = p²/2m where p = mv is linear momentum
│   In quantum: p becomes an OPERATOR p̂ = -iℏ(∂/∂x)
│   KE then becomes: KE → p̂²/2m = (-ℏ²/2m)(∂²/∂x²)
│   This IS the kinetic energy term in Schrödinger equation!
│
└── Exit check: Write H = T + V for a diatomic molecule (two masses + spring).
    Answer: H = p₁²/2m₁ + p₂²/2m₂ + ½k(x₁-x₂)²

Ph1.1.2  Units & Constants You Must Know
├── Reduced Planck constant: ℏ = h/2π = 1.055 × 10⁻³⁴ J·s
├── Electron mass: mₑ = 9.109 × 10⁻³¹ kg
├── Elementary charge: e = 1.602 × 10⁻¹⁹ C
├── Boltzmann: kB = 1.381 × 10⁻²³ J/K
├── Useful: 1 Hartree = 27.21 eV = 627.5 kcal/mol (quantum chemistry unit)
└── VQE energies are in Hartree. |1 mHa error| = "chemical accuracy"

Ph1.1.3  Harmonic Oscillator (Prototype QM System)
├── Classical: mass m on spring with constant k
│   F = -kx,  V(x) = ½kx²,  H = p²/2m + ½kx²
│   Natural frequency: ω = √(k/m)
│
├── Quantum: energy is QUANTIZED (not continuous!)
│   Eₙ = ℏω(n + ½),  n = 0,1,2,...
│   Ground state (n=0): E₀ = ℏω/2  ← zero-point energy, even at T=0!
│   This "zero-point energy" is how molecules vibrate even at absolute zero
│
├── BIO link:
│   Molecular vibrations (IR spectroscopy) = quantum harmonic oscillator
│   C-H bond stretches at ~3000 cm⁻¹ → discrete frequency levels
│   DNA base pair bending, protein backbone oscillations: same physics
│
└── Exit check: H₂ has ω ≈ 1.32×10¹⁴ rad/s. Compute E₀ and E₁. E₁-E₀ = ℏω.
```

---

## Module Ph1.2: Hamiltonian Mechanics

```
Ph1.2.1  Lagrangian and Conjugate Momentum
├── What to learn:
│   ├── Lagrangian: L = T - V (kinetic minus potential)
│   ├── Compare: Hamiltonian H = T + V (kinetic PLUS potential)
│   ├── Conjugate momentum: pᵢ = ∂L/∂q̇ᵢ  (derivative wrt velocity)
│   └── For simple particle: L = ½mv² - V → p = ∂L/∂v = mv ✓
│
└── Why two formulations?
    Lagrangian: useful for deriving equations of motion
    Hamiltonian: useful for quantization ("promote p,q to operators")

Ph1.2.2  Hamilton's Equations
├── Two 1st-order equations replace one 2nd-order Newton's law:
│   dqᵢ/dt = +∂H/∂pᵢ   (position changes w/ momentum)
│   dpᵢ/dt = -∂H/∂qᵢ   (momentum changes w/ force)
│
├── Example for particle: H = p²/2m + V(x)
│   dx/dt = ∂H/∂p = p/m = v ✓
│   dp/dt = -∂H/∂x = -dV/dx = F ✓ (Newton's 2nd law recovered)
│
└── Quantum transition:
    Classical: {A,B} Poisson bracket
    Quantum: [Â,B̂]/iℏ commutator
    This is the canonical quantization recipe!

Ph1.2.3  Hamiltonian as Generator of Evolution
├── H generates time evolution of any observable A:
│   dA/dt = {A, H}  (Poisson bracket form)
│   dÂ/dt = [Â,Ĥ]/iℏ + ∂Â/∂t  (quantum Heisenberg equation)
│
├── Key insight:
│   The Hamiltonian is not just "total energy"
│   It is the mathematical machine that DRIVES time evolution
│   Both classically AND quantum mechanically
│
├── BIO link:
│   AMBER, CHARMM, GROMACS force fields = classical Hamiltonian mechanics
│   V(r) = Σbonds k(r-r₀)² + Σangles kθ(θ-θ₀)² + Σtorsions + VLJZ + Velec
│   These simulate protein folding, drug-receptor binding
│
└── Exit check: H = p²/2m + ½kx². Write Hamilton's equations. Solve to get x(t)=A·cos(ωt+φ).
```

---

## Module Ph1.3: Wave Mechanics

```
Ph1.3.1  Classical Wave Equation
├── What to learn:
│   ├── Wave equation: ∂²y/∂t² = v²·(∂²y/∂x²)
│   ├── General solution: y(x,t) = A·sin(kx - ωt + φ)
│   │   k = 2π/λ (wavenumber),  ω = 2πf (angular frequency)
│   │   v = ω/k = λf (wave speed)
│   ├── Amplitude A: maximum displacement
│   └── Phase φ: initial offset
│
├── Code:
│   import numpy as np, matplotlib.pyplot as plt
│   x = np.linspace(0, 10, 500)
│   t_vals = [0, 0.5, 1.0]
│   for t in t_vals:
│       y = np.sin(2*np.pi*(x - t))  # v=1, λ=1
│       plt.plot(x, y, label=f't={t}')
│   plt.legend(); plt.show()
│
└── This IS how you visualize quantum wavefunctions!

Ph1.3.2  Superposition of Waves
├── Two waves with same wavelength, different phase:
│   y₁ = A·sin(kx - ωt)
│   y₂ = A·sin(kx - ωt + φ)
│   y_total = 2A·cos(φ/2)·sin(kx - ωt + φ/2)
│
├── φ=0: constructive interference (amplitudes add)
├── φ=π: destructive interference (amplitudes cancel)
├── Quantum superposition = WAVE superposition:
│   |ψ⟩ = α|0⟩ + β|1⟩ looks like two waves superposed
│   Interference of amplitudes determines measurement probability
│
└── Quantum algorithm design = engineering interference:
    Make CORRECT answers interfere constructively
    Make WRONG answers interfere destructively
    Example: Grover's algorithm works exactly this way!

Ph1.3.3  Standing Waves = Quantization
├── Opposite waves: y = A·sin(kx)·cos(ωt)
│   Fixed boundary: y(0,t) = y(L,t) = 0
│   Discrete allowed wavelengths: λₙ = 2L/n, n=1,2,3,...
│   → Only certain frequencies allowed! This IS quantization.
│
├── Same logic for electrons:
│   Electron in atom = standing wave around nucleus
│   Only discrete "fitting" wavelengths allowed
│   → Only discrete energy levels (orbitals)
│
├── de Broglie: λ = h/p (matter wavelength)
│   For electron at v = 10⁶ m/s:
│   λ = 6.63×10⁻³⁴ / (9.1×10⁻³¹ × 10⁶) ≈ 0.73 nm
│   This is comparable to atomic spacing → wave nature matters!
│
└── BIO link:
    DNA absorbs UV at 260 nm via π-electron resonance (standing wave modes)
    Protein β-sheets, α-helices = specific vibrational standing wave patterns
```

---

## Module Ph2.1: Schrödinger Equation ⛔ BLOCKER

```
Ph2.1.1  Time-Dependent Schrödinger Equation (TDSE)
├── The equation:
│   iℏ ∂|ψ⟩/∂t = Ĥ|ψ⟩
│
├── Every term explained:
│   i   = imaginary unit (complex number)
│   ℏ   = reduced Planck constant
│   ∂/∂t = partial time derivative
│   |ψ⟩  = quantum state (wavefunction)
│   Ĥ   = Hamiltonian operator (total energy operator)
│
├── For single particle in 1D:
│   Ĥ = -(ℏ²/2m)(∂²/∂x²) + V(x)
│   iℏ ∂ψ/∂t = [-(ℏ²/2m)(∂²/∂x²) + V(x)] ψ
│
├── Compare to classical: F=ma → d²x/dt² = -dV/dx
│   In QM: wavefunction ψ evolves, not position x
│
└── What ψ(x,t) means:
    ψ is complex-valued. You don't measure ψ directly.
    Born rule: probability density = |ψ(x,t)|²
    P(finding particle in [a,b]) = ∫[a to b] |ψ(x,t)|² dx

Ph2.1.2  Time-Independent Schrödinger Equation (TISE) — The VQE Target
├── When V doesn't depend on time → separate ψ(x,t) = ψ(x)·T(t)
│   Time part: T(t) = e^(-iEt/ℏ)  (pure phase rotation)
│   Space part: Ĥψ(x) = Eψ(x)  ← THIS IS TISE
│
├── TISE is an EIGENVALUE EQUATION:
│   Ĥ is the operator (matrix in finite QC)
│   ψ are eigenvectors = stationary states ("orbitals")
│   E are eigenvalues = allowed energy levels
│
├── Code (particle in box, numerical):
│   L, n_points = 1.0, 1000
│   x = np.linspace(0, L, n_points)
│   hbar, m = 1.0, 1.0  # natural units
│   # Build 2nd derivative matrix (finite difference)
│   dx = x[1]-x[0]
│   diag = -2*np.ones(n_points)
│   off  = np.ones(n_points-1)
│   D2 = (np.diag(diag) + np.diag(off,1) + np.diag(off,-1)) / dx**2
│   H = -(hbar**2/(2*m)) * D2   # particle in box V=0
│   E, psi = np.linalg.eigh(H)  # Hermitian eigenvalue solver
│   # Ground state energy should ≈ π²/2 in natural units
│   print(f"E₀ = {E[0]:.4f}, exact = {np.pi**2/2:.4f}")
│
└── BIO link:
    Molecular orbital theory: electrons in molecule = particle in 3D box
    HOMO-LUMO gap (frontier orbitals) determines:
    - Whether UV/visible light is absorbed (DNA photolesions → mutations)
    - Reactivity (drug-receptor binding sites)
    - Protein redox behavior (electron transfer in cellular respiration)

Ph2.1.3  Particle in a Box — Full Analytical Solution
├── Setup: V(x)=0 for 0<x<L, V=∞ at walls
│
├── Step 1: Solve TISE with V=0
│   -(ℏ²/2m)(d²ψ/dx²) = Eψ
│   d²ψ/dx² = -(2mE/ℏ²)ψ = -k²ψ where k=√(2mE)/ℏ
│   General: ψ(x) = A·sin(kx) + B·cos(kx)
│
├── Step 2: Apply boundary conditions
│   ψ(0)=0 → B=0, so ψ(x) = A·sin(kx)
│   ψ(L)=0 → sin(kL)=0 → kL = nπ → k = nπ/L
│
├── Step 3: Normalize ψ
│   ∫₀ᴸ |ψ|² dx = 1 → A = √(2/L)
│   ψₙ(x) = √(2/L) · sin(nπx/L), n=1,2,3,...
│
├── Step 4: Energy levels
│   E = ℏ²k²/2m = n²π²ℏ²/(2mL²) ∝ n²
│   E₁ : E₂ : E₃ = 1 : 4 : 9 (integers squared)
│
└── Exit check:
    Verify ψ₁ is normalized: ∫₀¹|√2·sin(πx)|²dx = 1.
    Plot ψ₁, ψ₂, ψ₃ and |ψ|² using matplotlib.
    Identify nodes. ψₙ has (n-1) interior nodes.
```

---

## Module Ph2.2: Quantum Postulates

```
Ph2.2.1  The Six Postulates (Commit to Memory)
├── P1 STATE:
│   A quantum system is described by |ψ⟩ in Hilbert space ℋ.
│   All physical information is encoded in |ψ⟩.
│
├── P2 OBSERVABLE:
│   Every measurable quantity is a Hermitian operator Â on ℋ.
│   Examples: energy Ĥ, position x̂, momentum p̂=-iℏ∂/∂x
│
├── P3 MEASUREMENT OUTCOME:
│   A measurement of Â can only yield an eigenvalue aₙ (where Â|aₙ⟩=aₙ|aₙ⟩).
│   This is WHY eigenvalues = possible energy levels!
│
├── P4 STATE COLLAPSE:
│   After measuring aₙ, state immediately collapses to |aₙ⟩.
│   This is irreversible and information-destroying.
│
├── P5 BORN RULE:
│   P(obtaining aₙ from state |ψ⟩) = |⟨aₙ|ψ⟩|²
│   This is the fundamental probabilistic rule of QM.
│
└── P6 TIME EVOLUTION:
    |ψ(t)⟩ evolves by TDSE: iℏ d|ψ⟩/dt = Ĥ|ψ⟩
    Solution: |ψ(t)⟩ = e^(-iĤt/ℏ)|ψ(0)⟩

Ph2.2.2  Commutators & Uncertainty
├── Commutator: [Â,B̂] = ÂB̂ - B̂Â
│   [Â,B̂] = 0 → Â and B̂ can be measured simultaneously
│   [Â,B̂] ≠ 0 → cannot know BOTH values precisely
│
├── Pauli commutators (memorize):
│   [X,Y] = 2iZ
│   [Y,Z] = 2iX
│   [Z,X] = 2iY
│   [X,Z] = -2iY  (implies X and Z CANNOT be measured simultaneously)
│
├── Heisenberg uncertainty:
│   ΔA · ΔB ≥ ½|⟨[Â,B̂]⟩|
│   Δx · Δp ≥ ℏ/2  (position-momentum uncertainty)
│
├── VQE link:
│   Measuring ⟨X₀Z₁⟩ requires Z-basis measurement after basis rotation.
│   The non-commutativity of Pauli terms means H must be measured
│   GROUPWISE (commuting groups), requiring multiple circuit executions.
│
└── Exit check:
    Compute [X,Z] = XZ - ZX using matrix multiplication.
    Answer: -2iY. Verify numerically.

Ph2.2.3  Expectation Values and the VQE Energy
├── General expectation:
│   ⟨Â⟩ = ⟨ψ|Â|ψ⟩ = Σₙ aₙ |⟨aₙ|ψ⟩|² (weighted average of eigenvalues)
│
├── For density matrix (mixed state):
│   ⟨Â⟩ = Tr(ρÂ)
│
├── VQE cost function:
│   ⟨E(θ)⟩ = ⟨ψ(θ)|H|ψ(θ)⟩
│   H = Σₖ cₖ Pₖ  (sum of Pauli strings)
│   ⟨H⟩ = Σₖ cₖ ⟨Pₖ⟩  (linearity of expectation)
│   Each ⟨Pₖ⟩ measured separately on quantum hardware
│
└── Exit check:
    State |ψ⟩ = (2/3)|0⟩ + (i√5/3)|1⟩.
    Compute ⟨Z⟩ = ⟨ψ|Z|ψ⟩. Check: real number? Yes (Z is Hermitian).
    Compute ΔZ = √(⟨Z²⟩ - ⟨Z⟩²). Check ΔZ ≥ 0.
```

---

## Module Ph2.3: Dirac Notation ⛔ BLOCKER for All QC

```
Ph2.3.1  Ket, Bra, and Their Relationship
├── Ket |ψ⟩: the quantum state — a column vector in Hilbert space
│   |0⟩ = [[1],[0]],  |1⟩ = [[0],[1]],  |+⟩ = (1/√2)[[1],[1]]
│
├── Bra ⟨ψ|: the dual — a row vector (conjugate transpose of ket)
│   ⟨0| = [1, 0],  ⟨1| = [0, 1],  ⟨+| = (1/√2)[1, 1]
│   Note: ⟨ψ| = |ψ⟩† (dagger operation)
│
├── Bra-Ket (inner product):
│   ⟨φ|ψ⟩: bra × ket = scalar
│   ⟨0|1⟩ = [1,0]·[0,1]ᵀ = 0  (orthogonal)
│   ⟨+|+⟩ = 1 (normalized)
│
└── Ket-Bra (outer product):
    |ψ⟩⟨φ|: ket × bra = matrix (operator!)
    |0⟩⟨0| = [[1,0],[0,0]]  (projection onto |0⟩)
    |1⟩⟨1| = [[0,0],[0,1]]  (projection onto |1⟩)
    |0⟩⟨0| + |1⟩⟨1| = I    (completeness!)

Ph2.3.2  Operators in Dirac Notation
├── Action of operator on ket:  Â|ψ⟩ = |new state⟩
├── Sandwich (matrix element):  ⟨φ|Â|ψ⟩ = scalar
│   This evaluates to: bra⟨φ| times operator Â times ket|ψ⟩
│   ⟨0|X|0⟩ = [1,0]·[[0,1],[1,0]]·[1,0]ᵀ = [1,0]·[0,1]ᵀ = 0
│   ⟨1|X|0⟩ = [0,1]·[[0,1],[1,0]]·[1,0]ᵀ = [0,1]·[0,1]ᵀ = 1
│
├── Expectation value notation:
│   ⟨Â⟩_ψ = ⟨ψ|Â|ψ⟩  ← THE VQE COST FUNCTION
│
└── Spectral decomposition (any Hermitian Â):
    Â = Σₙ aₙ|aₙ⟩⟨aₙ|
    Z = (+1)|0⟩⟨0| + (-1)|1⟩⟨1|  ← write Z in Dirac form!

Ph2.3.3  Projection Operators & Measurement
├── Projector onto state |n⟩: P̂ₙ = |n⟩⟨n|
│   Properties: P̂ₙ² = P̂ₙ (idempotent), P̂ₙ is Hermitian
│
├── Measurement of observable Â with eigenstates {|aₙ⟩}:
│   Before: |ψ⟩ = Σₙ cₙ|aₙ⟩
│   Outcome aₙ occurs with probability P(aₙ) = |cₙ|² = |⟨aₙ|ψ⟩|²
│   After measuring aₙ: state collapses to |aₙ⟩
│
├── Code translation:
│   # Everything above in NumPy:
│   ket_0 = np.array([[1],[0]])  # column vector
│   bra_0 = ket_0.conj().T      # row vector
│   proj_0 = ket_0 @ bra_0       # outer product = [[1,0],[0,0]]
│
└── Exit check (THE FINAL DIRAC EXAM):
    Given |ψ⟩ = (√3/2)|0⟩ + (1/2)(1+i)/√2|1⟩  [normalize first!]
    1. Is it normalized? If not, normalize it.
    2. Compute ⟨Z⟩ = ⟨ψ|Z|ψ⟩ by hand.
    3. Compute P(|0⟩) and P(|1⟩) using Born rule.
    4. Write Z in spectral form: Z = Σₙ zₙ|zₙ⟩⟨zₙ|.
    5. Read this without confusion: ⟨ψ(θ)|Ĥ|ψ(θ)⟩ — explain every symbol.

Ph2.3.4  Multi-qubit States in Dirac Notation
├── |00⟩ = |0⟩⊗|0⟩ = [1,0,0,0]ᵀ
├── |01⟩ = |0⟩⊗|1⟩ = [0,1,0,0]ᵀ
├── Entangled state (CANNOT factor):
│   |Φ+⟩ = (1/√2)(|00⟩ + |11⟩)
│   Proof it can't factor: assume |Φ+⟩ = (α|0⟩+β|1⟩)⊗(γ|0⟩+δ|1⟩)
│   = αγ|00⟩ + αδ|01⟩ + βγ|10⟩ + βδ|11⟩
│   Need: αδ=0 AND βγ=0 AND αγ=βδ=1/√2 → contradiction!
│
└── Exit check:
    Write the 4 Bell states in both ket notation and vector form.
    Verify all are orthonormal using NumPy inner products.
```
