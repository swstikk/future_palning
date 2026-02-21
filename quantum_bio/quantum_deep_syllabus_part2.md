# ⚛️ Quantum Bioinformatics — Deep Chapter-wise Syllabus PART 2
## Phase 1 Physics Track (Ph1.1 – Ph2.3)

---

# PHASE 1 — PHYSICS TRACK (Weeks 3-10, Parallel with Math)

---

## Module Ph1.1: Classical Mechanics — Energy & Hamiltonian Intuition

> **PREREQUISITES from Class 11 ISC/CBSE you MUST already know:**
> Newton's laws (F=ma), scalar vs vector, distance/displacement, speed/velocity.
> If shaky → revise HC Verma Ch 1-4 before starting.

```
Ph1.1.1  Newton's Laws → Energy Approach (WHY we shift to energy)
├── What you already know (Class 11 recap):
│   ├── Newton 1: No net force → constant velocity (inertia)
│   ├── Newton 2: F = ma → F = m(dv/dt) → F = dp/dt
│   ├── Newton 3: Action = -Reaction
│   └── Problem: F=ma works for 1 object. For 10²³ atoms? IMPOSSIBLE.
│
├── Why energy formulation is better:
│   ├── Energy is a SCALAR (not vector) → much simpler math
│   ├── Conservation of energy: E_total = constant (if no friction)
│   ├── Works for ANY number of particles (molecules, proteins, DNA)
│   └── Quantum mechanics is ENTIRELY built on the energy formulation
│       (Schrödinger equation = energy equation, NOT force equation)
│
├── Key terminology (memorize before moving forward):
│   ├── System: the object(s) we are studying
│   ├── State: everything needed to describe the system (position + velocity)
│   ├── Observable: anything we can measure (position, energy, momentum)
│   └── Conservation: a quantity that doesn't change with time
│
└── Self-check: Can you explain why F=ma is NOT used in quantum mechanics?
    (Because at atomic scale, you can't track exact position + velocity simultaneously.
     Energy formulation avoids this by working with probabilities.)

Ph1.1.2  Work — The Definition That Connects Force to Energy
├── What to learn:
│   ├── Work done by constant force: W = F·d·cosθ
│   │   F = force magnitude, d = displacement, θ = angle between F and d
│   │   θ=0 → W = Fd (force along motion, max work)
│   │   θ=90° → W = 0 (force perpendicular, e.g. circular motion)
│   │   θ=180° → W = -Fd (force opposes motion, e.g. friction)
│   │
│   ├── Work done by variable force: W = ∫F(x)dx
│   │   (This is integration — you learned in Math M1.2.4)
│   │   Example: spring → F = -kx → W = ∫₀ˣ(-kx)dx = -½kx²
│   │
│   ├── Work-Energy Theorem:
│   │   Net work done on object = Change in kinetic energy
│   │   W_net = ΔKE = ½mv²_final - ½mv²_initial
│   │   This is WHY work and energy are connected!
│   │
│   └── Units: Work is in Joules (J) = kg·m²/s²
│
├── Worked example (DO THIS BY HAND):
│   A 2kg ball is pushed with F=10N over d=3m along the floor (θ=0):
│   W = 10 × 3 × cos(0) = 30 J
│   If ball started from rest: 30 = ½(2)(v²) → v = √30 = 5.48 m/s
│
└── Why this matters for QM:
    The Hamiltonian H = T + V stores all energy information.
    In quantum: "work" becomes "operator acting on state."
    You won't use W=Fd in QM, but the CONCEPT of energy transfer remains.

Ph1.1.3  Kinetic Energy — Both Forms You Must Know
├── Standard form: KE = ½mv²
│   m = mass (kg), v = velocity (m/s)
│   Example: electron (m=9.1×10⁻³¹ kg) at v=10⁶ m/s
│   KE = ½(9.1×10⁻³¹)(10⁶)² = 4.55×10⁻¹⁹ J
│
├── Momentum form: KE = p²/(2m)   ← THIS IS THE QUANTUM FORM
│   Derive: p = mv → v = p/m → KE = ½m(p/m)² = p²/2m  ✓
│   Why preferred in QM:
│   │   In quantum mechanics, momentum p becomes an OPERATOR:
│   │   p̂ = -iℏ(d/dx)   (the hat means "operator", not just a number)
│   │   KE operator = p̂²/(2m) = [-iℏ(d/dx)]²/(2m)
│   │                = (-ℏ²/2m)(d²/dx²)
│   │   THIS exact expression appears in Schrödinger equation!
│   └── You don't need to "understand" operators yet — but memorize this form.
│
├── Code verification:
│   m = 9.109e-31   # electron mass (kg)
│   v = 1e6          # 10⁶ m/s
│   p = m * v        # momentum
│   KE_v = 0.5 * m * v**2       # from velocity
│   KE_p = p**2 / (2 * m)       # from momentum
│   print(f"KE(v) = {KE_v:.3e}, KE(p) = {KE_p:.3e}")  # same!
│
└── Gate: you MUST be able to switch between ½mv² and p²/2m instantly.

Ph1.1.4  Potential Energy — Types You Need
├── 1. Gravitational PE: V = mgh
│   ├── h = height from reference point
│   ├── Near Earth's surface (constant g)
│   └── Example: 1kg at 10m → V = 1×9.8×10 = 98 J
│
├── 2. Spring/Elastic PE: V = ½kx²
│   ├── k = spring constant (N/m), x = displacement from equilibrium
│   ├── Key: V ∝ x² → parabolic shape → minimum at x=0
│   ├── This is the HARMONIC OSCILLATOR potential (used everywhere in QM)
│   └── Example: k=200 N/m, x=0.05 m → V = ½(200)(0.05)² = 0.25 J
│
├── 3. Coulomb PE (electrostatic): V = kₑq₁q₂/r   ← CRITICAL FOR QM
│   ├── kₑ = 8.99×10⁹ N·m²/C²
│   ├── q₁, q₂ = charges, r = distance between them
│   ├── Negative → attraction (opposite charges), Positive → repulsion
│   ├── In atoms: V(r) = -kₑZe²/r (electron-nucleus attraction)
│   ├── In molecules: Σ electron-nuclear + Σ electron-electron + Σ nuclear-nuclear
│   └── THIS is what VQE actually solves: find E for this V
│
├── Force from potential:
│   F = -dV/dx  (force is negative gradient of potential)
│   Spring: F = -d(½kx²)/dx = -kx  ← Hooke's law recovered!
│   Gravity: F = -d(mgh)/dh = -mg  ← gravitational force recovered!
│   THIS relationship between F and V is fundamental. Memorize: F = -dV/dx.
│
└── Self-check:
    Given V(x) = 3x² + 2x, find F(x).
    Answer: F = -dV/dx = -(6x + 2) = -6x - 2.

Ph1.1.5  Total Energy, Conservation, and the Hamiltonian
├── Total mechanical energy:
│   E = KE + PE = T + V = ½mv² + V(x) = constant (if no friction)
│
├── The Hamiltonian H:
│   ├── Definition: H = T + V = p²/2m + V(x)
│   ├── This is just "total energy written with momentum instead of velocity"
│   ├── Classical: H is a number (e.g., H = 3.5 Joules)
│   ├── Quantum: H becomes an OPERATOR Ĥ = -ℏ²/2m·d²/dx² + V(x)
│   └── VQE goal: find minimum eigenvalue of Ĥ → ground state energy
│
├── Conservation in practice:
│   Ball at height h, drops to ground:
│   Initial: KE=0, PE=mgh → E = mgh
│   Final: PE=0, KE=½mv² → E = ½mv²
│   mgh = ½mv² → v = √(2gh)   ← no need for F=ma!
│
├── Worked example for H₂-like molecule:
│   Two H atoms, masses m₁=m₂=mₚ, connected by "bond" (like spring)
│   H = p₁²/(2mₚ) + p₂²/(2mₚ) + V(|x₁-x₂|)
│   ↑ KE of atom 1   ↑ KE of atom 2   ↑ interaction potential
│   In real H₂: V includes electron-electron, electron-nuclear, nuclear-nuclear
│   VQE computes the quantum version of this exact H
│
└── Exit check:
    1. Write H = T + V for a particle in a gravitational field
       Answer: H = p²/2m + mgx
    2. Write H for a spring system
       Answer: H = p²/2m + ½kx²
    3. Is energy conserved if there's friction? NO — only conservative systems.

Ph1.1.6  Units & Constants You Must Memorize
├── SI units:
│   ├── Force: Newton (N) = kg·m/s²
│   ├── Energy: Joule (J) = kg·m²/s² = N·m
│   ├── Momentum: kg·m/s
│   └── Angular momentum: kg·m²/s = J·s
│
├── Quantum constants (MEMORIZE — will use daily):
│   ├── h (Planck's constant) = 6.626 × 10⁻³⁴ J·s
│   ├── ℏ = h/(2π) = 1.055 × 10⁻³⁴ J·s   ← "h-bar", used in ALL QM equations
│   ├── mₑ (electron mass) = 9.109 × 10⁻³¹ kg
│   ├── e (elementary charge) = 1.602 × 10⁻¹⁹ C
│   ├── kₑ (Coulomb constant) = 8.988 × 10⁹ N·m²/C²
│   ├── kB (Boltzmann) = 1.381 × 10⁻²³ J/K
│   ├── c (speed of light) = 3.0 × 10⁸ m/s
│   └── 1 eV = 1.602 × 10⁻¹⁹ J (electron-volt, convenient energy unit)
│
├── Quantum chemistry units:
│   ├── 1 Hartree (Ha) = 27.211 eV = 627.5 kcal/mol
│   ├── "Chemical accuracy" = 1.6 mHa ≈ 1 kcal/mol
│   │   If VQE energy error < 1.6 mHa → result is chemically useful
│   ├── Bond lengths: measured in Ångström (Å), 1 Å = 10⁻¹⁰ m
│   └── H₂ bond length = 0.74 Å = 74 pm
│
├── Dimensional analysis trick:
│   ℏ²/(2mₑ) has units of J·m² → this × (d²/dx²) gives J (energy) ✓
│   Always check: does your answer have correct units?
│
└── Exit check: Compute KE of electron at 10⁶ m/s in both Joules and eV.
    KE = ½(9.1×10⁻³¹)(10⁶)² = 4.55×10⁻¹⁹ J = 2.84 eV.

Ph1.1.7  Harmonic Oscillator — The Most Important Classical System
├── Setup: mass m attached to spring with spring constant k
│   ├── Restoring force: F = -kx (Hooke's law)
│   ├── Equation of motion: m(d²x/dt²) = -kx → d²x/dt² = -(k/m)x
│   ├── Define ω = √(k/m) → d²x/dt² = -ω²x
│   └── Solution: x(t) = A·cos(ωt + φ)
│       A = amplitude, φ = initial phase, ω = angular frequency
│
├── Energy of oscillator:
│   ├── KE = ½mv² = ½mω²A²sin²(ωt+φ)
│   ├── PE = ½kx² = ½kA²cos²(ωt+φ) = ½mω²A²cos²(ωt+φ)
│   ├── Total E = ½mω²A²  (constant! KE↔PE exchange)
│   └── Frequency: f = ω/(2π), Period: T = 2π/ω = 2π√(m/k)
│
├── Hamiltonian form:
│   H = p²/(2m) + ½mω²x²
│   Note: we wrote ½kx² as ½mω²x² using k=mω²
│   This is the form used in quantum mechanics
│
├── Quantum version (preview — will solve in Ph2.1):
│   ├── Quantized energies: Eₙ = ℏω(n + ½), n=0,1,2,...
│   ├── Key difference from classical:
│   │   Classical: ANY energy allowed (continuous)
│   │   Quantum: ONLY discrete energies allowed (quantized!)
│   │   E₀=½ℏω → molecule vibrates even at absolute zero temperature
│   ├── Spacing between levels: ΔE = ℏω (constant for harmonic oscillator)
│   └── This is how IR spectroscopy identifies molecule bonds
│
├── Where harmonic oscillator appears in nature:
│   ├── ANY system near a minimum of V(x) ≈ V₀ + ½k(x-x₀)² (Taylor expand!)
│   ├── Molecular bonds (C-H, C-C, O-H vibrations)
│   ├── Crystal lattice vibrations (phonons → heat capacity)
│   ├── Electromagnetic field modes (photons = quantum harmonic oscillator!)
│   └── Protein backbone low-frequency vibrations
│
├── BIO link:
│   ├── IR spectroscopy: C-H stretch ≈ 3000 cm⁻¹, O-H ≈ 3500 cm⁻¹
│   │   Each bond vibration = harmonic oscillator with specific ω
│   │   DNA bases identified by their vibrational fingerprint
│   ├── Molecular dynamics simulations (AMBER, CHARMM):
│   │   V_bond = ½k(r-r₀)² → harmonic oscillator for each chemical bond
│   └── Drug-receptor binding: molecule vibration affects binding affinity
│
├── Worked example with real numbers:
│   H₂ molecule: k = 575 N/m, m_reduced = 8.37×10⁻²⁸ kg
│   ω = √(k/m) = √(575/8.37×10⁻²⁸) = 8.29×10¹⁴ rad/s
│   E₀ = ½ℏω = ½(1.055×10⁻³⁴)(8.29×10¹⁴) = 4.37×10⁻²⁰ J = 0.273 eV
│   E₁ = (3/2)ℏω = 3 × E₀ = 0.819 eV
│   ΔE = E₁-E₀ = ℏω = 0.546 eV (this is the IR absorption energy)
│
└── Code:
    import numpy as np
    k = 575           # spring constant (N/m)
    m = 8.37e-28       # reduced mass of H₂ (kg)
    hbar = 1.055e-34   # ℏ (J·s)
    omega = np.sqrt(k / m)
    for n in range(5):
        E = hbar * omega * (n + 0.5)
        print(f"E_{n} = {E:.3e} J = {E/1.602e-19:.3f} eV")
    # Output: E₀=0.273eV, E₁=0.819eV, E₂=1.365eV, ...

═══════════════════════════════════════════
 GATE TO Ph1.2 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Can write F = -dV/dx and derive force from any V(x)
 □ Can switch between KE = ½mv² and KE = p²/2m instantly
 □ Know what a Hamiltonian H = T + V is (just total energy in p-form)
 □ Can write H for spring system: H = p²/2m + ½mω²x²
 □ Know all 8 constants from Ph1.1.6 (ℏ, mₑ, e, k_e, k_B, h, c, eV↔J)
 □ Computed H₂ energy levels E₀, E₁ with real numbers
 □ Understand: classical → any energy, quantum → only discrete Eₙ
═══════════════════════════════════════════
```

---

## Module Ph1.2: Hamiltonian Mechanics

> **PREREQUISITES: Everything in Ph1.1 gate checklist must be checked off.**
> Specifically: H = T + V, p²/2m form, F = -dV/dx, harmonic oscillator H.
> Also need: partial derivatives from Math M1.2.3 (∂f/∂x, ∂f/∂p).

```
Ph1.2.1  Why This Module Exists — From Newton to Hamilton
├── Newton's way: F = ma (one equation per particle per direction)
│   For H₂ molecule with 2 electrons:
│   6 coupled differential equations (3D × 2 particles)
│   For hemoglobin (4000+ atoms): 12000+ equations. Unmanageable.
│
├── Lagrange's way: L = T - V, then apply Euler-Lagrange equation
│   d/dt(∂L/∂q̇) = ∂L/∂q → automatically gives equations of motion
│   Advantage: works with ANY coordinates (not just x,y,z)
│   Used in: classical mechanics, field theory
│
├── Hamilton's way: H = T + V (in terms of p, not v)
│   Uses POSITION (q) and MOMENTUM (p) as independent variables
│   Gives SYMMETRIC equations: ∂H/∂p and ∂H/∂q
│   Advantage: DIRECTLY translates to quantum mechanics!
│   Schrödinger equation: iℏ∂ψ/∂t = Ĥψ uses the Hamiltonian
│
├── Class 11 level understanding:
│   You already know: H = total energy = T + V
│   New step: write T using p (not v) → H = p²/2m + V
│   That's it! Hamilton's equations are just: ask "how does x change?"
│   and "how does p change?" → each answer comes from ∂H/∂(the other).
│
└── You don't need to be an expert here. The KEY takeaway:
    "The Hamiltonian drives everything in quantum mechanics."

Ph1.2.2  Lagrangian — What It Is (Simplified for Class 11)
├── Definition: L = T - V  (kinetic MINUS potential)
│   Compare: H = T + V (kinetic PLUS potential)
│
├── For a ball falling under gravity:
│   T = ½mv²,  V = mgy (y = height)
│   L = ½mv² - mgy
│   L = ½m(dy/dt)² - mgy   (writing v = dy/dt explicitly)
│
├── Euler-Lagrange equation (don't memorize derivation — just know it exists):
│   d/dt(∂L/∂ẏ) = ∂L/∂y
│   ∂L/∂ẏ = m·ẏ = mv = p   ← this is conjugate momentum!
│   d/dt(mv) = ∂L/∂y = -mg
│   → m(dv/dt) = -mg → a = -g → F = ma recovered ✓
│
├── Conjugate momentum (IMPORTANT definition):
│   pᵢ = ∂L/∂q̇ᵢ  (take derivative of L with respect to velocity)
│   For normal particle: pᵢ = ∂(½mv²)/∂v = mv  ← regular momentum ✓
│   For weird coordinates (angles, etc): p can mean angular momentum etc.
│
├── Worked example — Spring:
│   L = ½mẋ² - ½kx²
│   Conjugate momentum: p = ∂L/∂ẋ = mẋ = mv  ✓
│   Euler-Lagrange: d/dt(mẋ) = ∂L/∂x = -kx → mẍ = -kx → SHM ✓
│
└── You need to know WHAT the Lagrangian is and that p = ∂L/∂q̇.
    Deep theorem proofs NOT needed for quantum track.

Ph1.2.3  From Lagrangian to Hamiltonian (The Legendre Transform)
├── Recipe to build H from L:
│   Step 1: Write L = T - V in terms of q and q̇
│   Step 2: Find conjugate momentum p = ∂L/∂q̇
│   Step 3: Solve for q̇ in terms of p (invert: q̇ = p/m for simple case)
│   Step 4: H = pq̇ - L  (Legendre transform formula)
│
├── Worked example (free particle):
│   L = ½mq̇², V=0
│   Step 2: p = ∂L/∂q̇ = mq̇ → q̇ = p/m
│   Step 4: H = p·(p/m) - ½m(p/m)² = p²/m - p²/(2m) = p²/(2m)  ✓
│   This is just KE = p²/2m.  Nothing new — but now in "Hamiltonian language."
│
├── Worked example (spring):
│   L = ½mq̇² - ½kq²
│   p = mq̇ → q̇ = p/m
│   H = p(p/m) - [½m(p/m)² - ½kq²]
│     = p²/m  - p²/(2m) + ½kq²
│     = p²/(2m) + ½kq²   ← H = T + V  ✓
│
└── Key insight:
    For ALL systems you'll encounter: H = T + V.
    The Legendre transform is the formal justification for this.
    In quantum: just replace H with Ĥ (make it an operator).

Ph1.2.4  Hamilton's Equations of Motion — The Core
├── Two equations for every particle:
│   dq/dt = +∂H/∂p     ("velocity from Hamiltonian")
│   dp/dt = -∂H/∂q     ("force from Hamiltonian")
│
├── Physical meaning:
│   First eq: "How does position change?" → related to momentum (via H)
│   Second eq: "How does momentum change?" → related to force (via H)
│   Together they contain ALL dynamics. Equivalent to Newton's F=ma.
│
├── WORKED EXAMPLE — Spring (DO THIS YOURSELF):
│   Given: H = p²/(2m) + ½kx²
│
│   Equation 1: dx/dt = ∂H/∂p = ∂[p²/(2m)]/∂p = 2p/(2m) = p/m
│               → dx/dt = p/m   → v = p/m  ✓ (makes sense!)
│
│   Equation 2: dp/dt = -∂H/∂x = -∂[½kx²]/∂x = -kx
│               → dp/dt = -kx   → F = -kx  ✓ (Hooke's law!)
│
│   Combine: dp/dt = -kx and dx/dt = p/m
│   → m(d²x/dt²) = -kx → x(t) = A·cos(ωt + φ) where ω=√(k/m)
│
├── WORKED EXAMPLE — Gravitational field:
│   Given: H = p²/(2m) + mgx
│
│   dx/dt = ∂H/∂p = p/m   → velocity ✓
│   dp/dt = -∂H/∂x = -mg  → gravitational force ✓
│   → m(d²x/dt²) = -mg → x(t) = x₀ + v₀t - ½gt²  ← free fall!
│
├── PRACTICE PROBLEM (DO BY HAND):
│   Given: H = p²/(2m) + αx⁴  (anharmonic potential, α=constant)
│   Find: dx/dt and dp/dt using Hamilton's equations.
│   Answer: dx/dt = p/m,  dp/dt = -4αx³
│
└── Why Hamilton's form is better for quantum:
    Classical: q and p are numbers you plug in
    Quantum: q̂ and p̂ become OPERATORS
    Hamilton's equations → Heisenberg equation of motion
    Schrödinger equation is the state-based version of the same physics

Ph1.2.5  Poisson Brackets → Commutators (The Classical-Quantum Bridge)
├── Poisson bracket (classical):
│   {A, B} = Σᵢ (∂A/∂qᵢ · ∂B/∂pᵢ - ∂A/∂pᵢ · ∂B/∂qᵢ)
│
├── You don't need to compute Poisson brackets by hand.
│   Just know these KEY results:
│   {q, p} = 1    (fundamental!)
│   {q, q} = 0    {p, p} = 0
│   {A, H} = dA/dt  (time evolution is a Poisson bracket with H!)
│
├── The quantum translation rule (Dirac's prescription):
│   Every Poisson bracket → commutator divided by iℏ:
│   {A, B} → [Â, B̂]/(iℏ)
│
│   So {q, p} = 1 becomes:
│   [q̂, p̂] = iℏ   ← the CANONICAL COMMUTATION RELATION
│   [x̂, p̂] = iℏ   (using x instead of q)
│   This is the foundation of all quantum mechanics!
│
├── Physical consequence:
│   [x̂, p̂] ≠ 0 means x and p cannot BOTH be measured precisely simultaneously
│   → Heisenberg uncertainty: Δx·Δp ≥ ℏ/2
│   THIS is why energy formulation (not force) is needed in QM
│   You can't track exact trajectory → work with wavefunctions instead
│
└── For now: just memorize {q,p}=1 → [q̂,p̂]=iℏ. Details in Ph2.2.

Ph1.2.6  The Hamiltonian in Biology — Force Fields (Where This Applies)
├── Molecular Dynamics (MD) simulations use CLASSICAL Hamiltonian mechanics:
│
├── MD Hamiltonian: H = KE + V_total
│   V_total is broken into pieces:
│
│   1. V_bond = Σ ½k_bond(r - r₀)²
│      → harmonic oscillator for each chemical bond (C-C, C-H, N-H, etc.)
│      Each bond has its own k and r₀ (tabulated in force field)
│
│   2. V_angle = Σ ½k_angle(θ - θ₀)²
│      → harmonic oscillator for bond angles (C-C-C ≈ 109.5°, etc.)
│
│   3. V_torsion = Σ V_n[1 + cos(nφ - δ)]
│      → energy cost of rotating around bonds (controls protein shape)
│
│   4. V_LJ = Σ 4ε[(σ/r)¹² - (σ/r)⁶]
│      → van der Waals attraction/repulsion between non-bonded atoms
│      σ/r¹² = repulsion (atoms overlap), σ/r⁶ = attraction (London force)
│
│   5. V_elec = Σ kₑ·qᵢqⱼ/rᵢⱼ
│      → Coulomb electrostatic interaction between charged atoms
│
├── AMBER, CHARMM, OPLS, GROMOS = different parameter sets for these terms
│   These simulate: protein folding, drug binding, DNA dynamics, membrane transport
│
├── Why VQE is needed:
│   All of V_bond, V_angle above are APPROXIMATIONS (classical harmonic)
│   Real molecular interactions are QUANTUM (electrons are quantum objects)
│   VQE solves the true quantum Hamiltonian without these approximations
│   Especially critical for: metalloprotein active sites, reaction barriers,
│   excited states (photochemistry → DNA damage, vision, photosynthesis)
│
└── Key takeaway:
    Classical Hamiltonian mechanics = how drug companies simulate proteins TODAY
    Quantum Hamiltonian mechanics (VQE) = the NEXT generation of this

═══════════════════════════════════════════
 GATE TO Ph1.3 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know: Lagrangian L = T - V vs Hamiltonian H = T + V
 □ Can derive conjugate momentum p = ∂L/∂q̇ for simple examples
 □ Can write AND USE Hamilton's equations: dq/dt=∂H/∂p, dp/dt=-∂H/∂q
 □ Solved spring example by hand: got dx/dt=p/m, dp/dt=-kx
 □ Know: {q,p}=1 (classical) → [x̂,p̂]=iℏ (quantum)
 □ Can explain in one sentence: "why does quantum use Hamiltonian not Newton?"
   (Because [x̂,p̂]≠0 means you can't track trajectories → need wavefunctions)
 □ Know what molecular force fields (AMBER/CHARMM) do and their limitations
═══════════════════════════════════════════
```

---

## Module Ph1.3: Wave Mechanics

> **PREREQUISITES: Ph1.2 gate must be passed.**
> Also need from Class 11 ISC Physics: Chapter on Waves (transverse/longitudinal,
> amplitude, wavelength, frequency, velocity). Revise if forgotten.
> Math needed: sin/cos functions, derivatives (for wave equation).

```
Ph1.3.1  Wave Parameters — Class 11 Recap + Extensions
├── What is a wave:
│   ├── A disturbance that transfers ENERGY (not matter) through a medium
│   ├── Transverse: displacement ⊥ to direction of travel (string, EM waves)
│   ├── Longitudinal: displacement ∥ to direction (sound, spring compression)
│   └── Quantum waves: probability amplitude waves (neither transverse nor longitudinal — abstract!)
│
├── All wave parameters (MUST know every one):
│   ├── Amplitude A: maximum displacement from equilibrium (height of wave)
│   │   In quantum: |ψ|² gives probability density, amplitude ψ is complex-valued
│   ├── Wavelength λ: distance between two identical points (peak-to-peak)
│   │   Units: meters (m). Atomic scales: nanometers (nm) to picometers (pm)
│   ├── Frequency f: number of complete oscillations per second
│   │   Units: Hertz (Hz) = 1/s. Light: ~10¹⁴-10¹⁵ Hz
│   ├── Period T: time for one complete oscillation. T = 1/f
│   ├── Angular frequency ω: ω = 2πf (radians per second)
│   │   This is ALWAYS used in quantum equations instead of f
│   ├── Wavenumber k: k = 2π/λ (radians per meter)
│   │   Analogous to ω but for space instead of time
│   └── Wave speed v: v = λf = ω/k
│       Speed of light: c = 3×10⁸ m/s. Speed of sound: ~343 m/s
│
├── Worked example (DO BY HAND):
│   Red light: λ = 700 nm = 700×10⁻⁹ m
│   f = c/λ = 3×10⁸ / 700×10⁻⁹ = 4.29×10¹⁴ Hz
│   ω = 2πf = 2π × 4.29×10¹⁴ = 2.69×10¹⁵ rad/s
│   k = 2π/λ = 2π / 700×10⁻⁹ = 8.98×10⁶ rad/m
│   Energy per photon: E = hf = ℏω = 6.63×10⁻³⁴ × 4.29×10¹⁴ = 2.84×10⁻¹⁹ J = 1.77 eV
│
├── Mathematical description of a traveling wave:
│   y(x,t) = A·sin(kx - ωt + φ)
│   ├── kx: where you are in space (phase from position)
│   ├── ωt: where you are in time (phase from time)
│   ├── φ: initial phase (starting angle at t=0, x=0)
│   ├── (kx - ωt): wave moves in +x direction (positive velocity)
│   └── (kx + ωt): wave moves in -x direction (negative velocity)
│
└── Code (plot a traveling wave):
    import numpy as np, matplotlib.pyplot as plt
    x = np.linspace(0, 10, 500)
    k, omega = 2*np.pi, 2*np.pi  # λ=1, T=1
    for t in [0, 0.25, 0.5]:
        y = np.sin(k*x - omega*t)
        plt.plot(x, y, label=f't={t}')
    plt.xlabel('x'); plt.ylabel('y(x,t)')
    plt.title('Traveling wave at different times')
    plt.legend(); plt.grid(True); plt.show()

Ph1.3.2  The Wave Equation — Why It Matters
├── Classical wave equation (1D):
│   ∂²y/∂t² = v² · ∂²y/∂x²
│
│   Left side: acceleration in time (how fast displacement changes)
│   Right side: curvature in space (how "bent" the wave shape is) × v²
│   Meaning: curvier shape → faster acceleration → higher frequency
│
├── How to verify a solution satisfies the wave equation:
│   Given: y = A·sin(kx - ωt)
│   ∂y/∂t = -Aω·cos(kx-ωt) → ∂²y/∂t² = -Aω²·sin(kx-ωt)
│   ∂y/∂x = Ak·cos(kx-ωt)  → ∂²y/∂x² = -Ak²·sin(kx-ωt)
│   Check: ∂²y/∂t² = v²·∂²y/∂x²
│   -Aω²sin = v²·(-Ak²sin) → ω² = v²k² → v = ω/k ✓
│
├── Connection to Schrödinger:
│   Classical wave eq: ∂²y/∂t² = v²·∂²y/∂x²  (2nd order in time)
│   Schrödinger eq:    iℏ·∂ψ/∂t = Ĥψ          (1st order in time!)
│   Schrödinger is NOT the classical wave equation applied to electrons
│   It was GUESSED by Schrödinger using insight from de Broglie
│   But the wave CONCEPTS (superposition, interference, nodes) carry over
│
└── Self-check: Take y(x,t) = cos(3x - 6t). What is k? ω? v? λ? f?
    k=3, ω=6, v=ω/k=2, λ=2π/3, f=ω/2π=3/π ≈ 0.955 Hz

Ph1.3.3  Superposition — The Principle That Makes Quantum Work
├── Superposition principle:
│   If y₁(x,t) and y₂(x,t) are BOTH solutions of wave equation,
│   then y₁ + y₂ is ALSO a solution. (Linearity!)
│   You can ADD waves and the result is still a valid wave.
│
├── Constructive interference (φ difference = 0 or 2nπ):
│   Waves "in step" → amplitudes ADD → bigger combined wave
│   y₁ = A·sin(kx) + y₂ = A·sin(kx) → y_total = 2A·sin(kx)
│
├── Destructive interference (φ difference = π or (2n+1)π):
│   Waves "out of step" → amplitudes CANCEL → zero!
│   y₁ = A·sin(kx) + y₂ = A·sin(kx+π) = A·sin(kx) - A·sin(kx) = 0
│
├── General case (two waves, phase difference φ):
│   y₁ = A·sin(kx - ωt),  y₂ = A·sin(kx - ωt + φ)
│   y_total = 2A·cos(φ/2)·sin(kx - ωt + φ/2)
│   Envelope = 2A·cos(φ/2): controls how much they add/cancel
│   φ=0 → envelope=2A (max), φ=π → envelope=0 (cancel)
│
├── Double slit experiment analogy (Young's experiment, Class 12 topic):
│   Light through 2 slits → bright and dark bands on screen
│   Bright = constructive (waves from both slits arrive in phase)
│   Dark = destructive (waves arrive out of phase by π)
│   SAME thing happens with ELECTRONS → proves wave nature of matter!
│
├── QUANTUM SUPERPOSITION = same math, deeper meaning:
│   |ψ⟩ = α|0⟩ + β|1⟩ (qubit in superposition)
│   α and β are COMPLEX AMPLITUDES (not just real sin waves)
│   |α|² + |β|² = 1 (normalization = total probability = 1)
│   Measurement: |α|² chance of |0⟩, |β|² chance of |1⟩
│
├── Why interference matters for quantum computing:
│   ├── Grover's algorithm: oracle marks correct answer with -1 phase
│   │   Diffuser creates constructive interference at correct answer
│   │   Destructive interference at wrong answers → they shrink
│   │   After O(√N) iterations → correct answer dominates ≈100%
│   │
│   ├── VQE: ansatz circuit creates interference patterns
│   │   Different θ values → different interference patterns
│   │   Optimizer finds θ where energy expectation is minimized
│   │   THIS is the quantum advantage: exploring interference landscape
│   │
│   └── Key insight: classical computers CANNOT efficiently simulate interference
│       of 2ⁿ amplitudes simultaneously. Quantum computers CAN.
│
└── Code (constructive vs destructive):
    x = np.linspace(0, 10, 500)
    y1 = np.sin(2*np.pi*x)
    fig, axes = plt.subplots(1, 3, figsize=(15,4))
    for i, phi in enumerate([0, np.pi/2, np.pi]):
        y2 = np.sin(2*np.pi*x + phi)
        axes[i].plot(x, y1+y2, 'b-', lw=2)
        axes[i].plot(x, y1, 'r--', alpha=0.5)
        axes[i].plot(x, y2, 'g--', alpha=0.5)
        axes[i].set_title(f'φ={phi:.2f} rad')
    plt.tight_layout(); plt.show()

Ph1.3.4  Standing Waves → Quantization (THE Key Concept)
├── When a wave is trapped between walls (fixed boundaries):
│   Only waves that "fit" are allowed (boundary condition)
│   y(0) = 0 and y(L) = 0 (fixed at both ends)
│
├── Mathematics:
│   y(x,t) = A·sin(kx)·cos(ωt)  (standing wave = product of space × time)
│   Boundary: sin(kL) = 0 → kL = nπ → k = nπ/L, n=1,2,3,...
│   Allowed wavelengths: λₙ = 2L/n (only these fit!)
│   Allowed frequencies: fₙ = nv/(2L) = n·f₁
│
├── Nodes and antinodes:
│   Nodes: points that NEVER move → y=0 always → sin(kx)=0
│   n=1: 0 interior nodes (fundamental)
│   n=2: 1 interior node (first overtone)
│   n=3: 2 interior nodes (second overtone)
│   Antinodes: points of maximum vibration (between nodes)
│
├── THIS IS QUANTIZATION:
│   │   Not every wavelength is allowed — only λₙ = 2L/n
│   │   Not every frequency is allowed — only fₙ = nf₁
│   │   Not every energy is allowed — only Eₙ ∝ n²
│   │   The integer n = QUANTUM NUMBER (first appearance!)
│   └── Classical: any wavelength/frequency is fine
│       Quantum: only discrete values allowed
│       The "box" forces discreteness. Atom = natural "box" for electrons.
│
├── From ISC guitar string to quantum:
│   Guitar string: L=0.65m, standing waves at f₁=330Hz (E4 note)
│   Electron in atom: L≈1Å=10⁻¹⁰m, standing waves at E₁≈13.6eV
│   Same physics, different scale.
│
└── Exit check:
    A string of length L=1m, v=100 m/s.
    What are f₁, f₂, f₃? (Answer: 50, 100, 150 Hz)
    How many nodes does n=4 have? (Answer: 3 interior nodes)

Ph1.3.5  de Broglie Hypothesis — Matter Waves
├── de Broglie (1924): ALL matter has wave-like behavior
│   λ = h/p = h/(mv)  (matter wavelength)
│
├── When does wave nature matter?
│   Baseball (m=0.15kg, v=40m/s): λ=h/p = 6.63e-34/(0.15×40) = 1.1×10⁻³⁴m → WAY too small
│   Electron (m=9.1e-31, v=10⁶): λ = 6.63e-34/(9.1e-31×10⁶) = 0.73nm → comparable to atom!
│
│   Rule: if λ ≈ size of system → quantum effects dominate
│   Atoms are ~0.1-0.5 nm → electron λ ≈ 0.7 nm → QUANTUM EFFECTS!
│   This is why electron behavior in molecules must be solved quantum mechanically
│
├── Experimental proof:
│   Davisson-Germer experiment (1927): electrons diffract off crystal
│   Electrons show interference pattern → WAVE behavior confirmed!
│   If electrons were particles: no interference pattern (just two bands)
│
├── Energy-wavelength relation:
│   E = hf = hv/λ,  and p = h/λ,  and E = p²/2m (kinetic energy)
│   So: E = (h/λ)²/(2m) = h²/(2mλ²)
│   Shorter λ → higher energy (makes sense: faster electron = shorter wavelength)
│
├── BIO link:
│   ├── Electron microscopy: electron λ ≈ 0.005 nm at 50kV
│   │   Much shorter than light (400-700nm) → much higher resolution
│   │   Can image individual atoms in proteins, DNA, viruses
│   ├── X-ray crystallography: X-ray λ ≈ 0.1 nm (comparable to atom spacing)
│   │   Bragg diffraction → determine protein 3D structure
│   │   Rosalind Franklin's Photo 51 → discovered DNA double helix structure!
│   └── Neutron diffraction: finds hydrogen atoms in proteins (X-rays can't)
│
└── Code:
    h = 6.626e-34
    objects = [
        ("Baseball", 0.15, 40),
        ("Bullet", 0.01, 700),
        ("Electron", 9.1e-31, 1e6),
        ("Proton", 1.67e-27, 1e4),
    ]
    for name, m, v in objects:
        lam = h / (m * v)
        print(f"{name:12s}: λ = {lam:.3e} m")
    # Electron: 7.3e-10 m ≈ 0.73 nm → quantum regime!

Ph1.3.6  Photon Energy & Planck's Relation (Bridge to Quantum)
├── Planck (1900): energy of light comes in packets (quanta):
│   E = hf = ℏω  (energy of one photon)
│
├── Photoelectric effect (Einstein 1905):
│   Light hits metal → electrons ejected
│   Classical prediction (WRONG): brighter light = more energy per electron
│   Quantum reality: f > f₀ needed (frequency threshold, not brightness)
│   E_electron = hf - W (work function W = minimum energy to eject)
│
├── Connection to de Broglie:
│   Photon: E = hf, p = h/λ = E/c (massless)
│   Electron: E = p²/2m, p = h/λ (massive)
│   BOTH have wavelength! Both exhibit wave-particle duality.
│
├── Line spectra → discrete energy levels (atoms):
│   Hydrogen emission: only specific wavelengths (colors) emitted
│   656nm (red), 486nm (cyan), 434nm (blue), 410nm (violet)
│   Explanation: electron transitions between discrete energy levels
│   E_photon = E_upper - E_lower = hf → specific f → specific λ
│   This PROVES electrons in atoms have quantized energy levels!
│
├── 🧬 BIO link:
│   ├── DNA absorbs UV at 260nm: π→π* electronic transition
│   │   E = hc/λ = 6.63e-34 × 3e8 / 260e-9 = 7.65e-19 J = 4.77 eV
│   │   This energy breaks molecular bonds → MUTATIONS
│   ├── GFP (Green Fluorescent Protein): absorbs 395nm/475nm, emits 509nm
│   │   Stokes shift: emission λ > absorption λ (some energy → heat)
│   └── Photosynthesis: chlorophyll absorbs 430nm(blue) + 680nm(red)
│       Excited electron drives ATP synthesis → all life depends on E=hf!
│
└── Exit check:
    UV-C germicidal light: λ=254nm.
    Compute photon energy in eV. (Answer: E=hc/λ=4.88eV)
    Is this enough to break a C-C bond (~3.6eV)? YES → kills bacteria by DNA damage.

═══════════════════════════════════════════
 GATE TO Ph2.1 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know all wave parameters: A, λ, f, T, ω, k, v and can compute each
 □ Can write y(x,t) = A·sin(kx-ωt) and explain every term
 □ Understand constructive/destructive interference with worked examples
 □ Know: standing waves → only discrete λₙ = 2L/n → quantization
 □ Can compute de Broglie wavelength for electron: λ = h/(mv)
 □ Know: E = hf = ℏω (photon energy), and when λ matters (λ ≈ system size)
 □ Can explain in one sentence: "why are electron energy levels quantized?"
   (Electron = standing wave in atom, only certain λ fit → only certain E allowed)
 □ Computed DNA UV absorption energy (260nm → 4.77eV) correctly
═══════════════════════════════════════════
```

---

## Module Ph2.1: Schrödinger Equation ⛔ BLOCKER

> **PREREQUISITES: ALL previous gates passed (Ph1.1, Ph1.2, Ph1.3).**
> Specifically MUST know:
> - Complex numbers: e^(iθ) = cosθ + i·sinθ (from Math M1.1)
> - Derivatives and partial derivatives (Math M1.2)
> - Hamiltonian: H = p²/2m + V(x) (from Ph1.2)
> - Wave equation: ∂²y/∂t² = v²·∂²y/∂x² (from Ph1.3)
> - Standing waves + quantization (from Ph1.3.4)
> - de Broglie: λ = h/p → matter has wavelength (Ph1.3.5)
> - E = hf = ℏω (Ph1.3.6)

```
Ph2.1.1  Building the Schrödinger Equation from Things You Already Know
├── This section constructs the equation piece by piece.
│   No memorization of a mysterious formula — we BUILD it.
│
├── Step 1: Start with de Broglie wave for free particle
│   ψ(x,t) = A·e^(i(kx-ωt))   (complex traveling wave)
│   where k = p/ℏ (from λ=h/p → k=2π/λ=p/ℏ)
│   and   ω = E/ℏ (from E=ℏω)
│
├── Step 2: What does ∂ψ/∂t give?
│   ∂ψ/∂t = -iω · ψ = -i(E/ℏ) · ψ
│   → iℏ · ∂ψ/∂t = E · ψ          ← multiply both sides by iℏ
│
├── Step 3: What does ∂²ψ/∂x² give?
│   ∂ψ/∂x = ik · ψ
│   ∂²ψ/∂x² = (ik)² · ψ = -k² · ψ = -(p/ℏ)² · ψ = -p²/ℏ² · ψ
│   → -ℏ²/(2m) · ∂²ψ/∂x² = p²/(2m) · ψ = KE · ψ
│
├── Step 4: Add potential energy V(x)
│   Total energy E = KE + V = p²/(2m) + V(x)
│   From Step 2: iℏ · ∂ψ/∂t = E · ψ
│   From Step 3: E · ψ = [-ℏ²/(2m) · ∂²/∂x² + V(x)] ψ
│   Therefore:
│
│   ┌─────────────────────────────────────────────┐
│   │  iℏ ∂ψ/∂t = [-ℏ²/(2m)(∂²/∂x²) + V(x)] ψ  │
│   │                                              │
│   │  iℏ ∂ψ/∂t = Ĥψ     (compact form)          │
│   └─────────────────────────────────────────────┘
│   THIS IS THE TIME-DEPENDENT SCHRÖDINGER EQUATION (TDSE)
│
├── Every symbol explained (must be able to explain each one):
│   i   = √(-1), imaginary unit from Math M1.1
│   ℏ   = 1.055×10⁻³⁴ J·s, reduced Planck constant from Ph1.1.6
│   ∂/∂t = partial derivative with respect to time (Math M1.2.3)
│   ψ   = wavefunction = complex-valued function of x and t
│   Ĥ   = Hamiltonian OPERATOR (not just a number anymore!)
│       = -ℏ²/(2m)·∂²/∂x² + V(x)
│         ↑ kinetic energy op   ↑ potential energy (stays a function)
│
├── Why the "i" on the left side?
│   Without i: solutions would grow/decay exponentially (unphysical)
│   With i: solutions are OSCILLATING (e^(-iEt/ℏ) = rotation in complex plane)
│   This ensures |ψ|² stays constant → probability is conserved!
│
├── Classical comparison:
│   Newton:      F = ma → d²x/dt² = -dV/dx   (tracks position x)
│   Schrödinger: iℏ∂ψ/∂t = Ĥψ              (tracks wavefunction ψ)
│   Newton → tells you WHERE particle IS
│   Schrödinger → tells you PROBABILITY of where particle could be
│
└── Self-check: Can you write the Schrödinger equation from memory?
    Fill in: iℏ ∂|ψ⟩/∂t = ___|ψ⟩
    Answer: Ĥ = -ℏ²/(2m)·∂²/∂x² + V(x)

Ph2.1.2  What Does ψ (the Wavefunction) Actually Mean?
├── ψ(x,t) is a COMPLEX-VALUED function
│   You CANNOT directly observe ψ. It's not a physical wave you can see.
│   It's a mathematical tool that encodes ALL information about the particle.
│
├── Born's interpretation — THE KEY:
│   |ψ(x,t)|² = probability DENSITY at position x, time t
│   P(finding particle between x and x+dx) = |ψ(x,t)|² dx
│   P(finding particle between a and b) = ∫[a to b] |ψ(x,t)|² dx
│
├── Normalization requirement:
│   ∫[-∞ to +∞] |ψ(x,t)|² dx = 1 (particle MUST be somewhere!)
│   If your calculated ψ doesn't satisfy this → multiply by constant to fix
│   This gives normalization constant A in ψ = A·f(x)
│
├── What |ψ|² looks like for different states:
│   Ground state (n=1): one bump centered in box, max in middle
│   First excited (n=2): two bumps with zero crossing (node) in middle
│   n=3: three bumps, two nodes
│   Higher n → more wiggly → higher energy (more curvature = more KE)
│
├── Probability vs probability density (IMPORTANT distinction):
│   |ψ(x)|² has units of 1/meter (density = per unit length)
│   To get actual probability, multiply by dx:  P = |ψ|² · dx
│   Analogy: mass density ρ(x) [kg/m] → mass = ∫ρ(x)dx
│
├── Complex phase matters!
│   ψ₁ = (1/√2)|0⟩ + (1/√2)|1⟩  → P(|0⟩) = 0.5
│   ψ₂ = (1/√2)|0⟩ - (1/√2)|1⟩  → P(|0⟩) = 0.5 (SAME probability!)
│   But ψ₁ ≠ ψ₂ (different relative phase → different interference → different physics)
│   Measuring in Z-basis: same result
│   Measuring in X-basis: DIFFERENT result (this is where phase matters!)
│
└── Code (visualize |ψ|² for particle in box):
    import numpy as np, matplotlib.pyplot as plt
    L = 1.0
    x = np.linspace(0, L, 500)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for n, ax in zip([1, 2, 3, 4], axes.flat):
        psi = np.sqrt(2/L) * np.sin(n * np.pi * x / L)
        ax.plot(x, psi, 'b-', label=f'ψ_{n}')
        ax.plot(x, psi**2, 'r-', label=f'|ψ_{n}|²')
        ax.set_title(f'n={n}, nodes={n-1}')
        ax.legend(); ax.grid(True)
    plt.tight_layout(); plt.show()

Ph2.1.3  Time-Independent Schrödinger Equation (TISE) — The VQE Target
├── When potential V doesn't change with time (most molecular problems!):
│   We can SEPARATE variables: ψ(x,t) = φ(x) · T(t)
│
├── Separation of variables (FULL derivation, step by step):
│   Start: iℏ ∂[φ(x)T(t)]/∂t = [-ℏ²/(2m)·∂²/∂x² + V(x)] φ(x)T(t)
│
│   Left side:  iℏ · φ(x) · dT/dt        (only T depends on t)
│   Right side: T(t) · [-ℏ²/(2m)·d²φ/dx² + V(x)·φ]  (only φ depends on x)
│
│   Divide both sides by φ(x)·T(t):
│   iℏ(1/T)(dT/dt) = (1/φ)[-ℏ²/(2m)·d²φ/dx² + V·φ]
│
│   Left = function of t ONLY. Right = function of x ONLY.
│   Both sides must equal a CONSTANT → call it E (separation constant = energy!)
│
├── Two separate equations:
│   TIME equation: iℏ(dT/dt) = E·T
│   Solution: T(t) = e^(-iEt/ℏ)  ← just a phase rotation!
│   |T(t)|² = |e^(-iEt/ℏ)|² = 1 → probability doesn't change with time
│   → "stationary state" (probability distribution is time-independent)
│
│   SPACE equation: [-ℏ²/(2m)·d²φ/dx² + V(x)·φ] = E·φ
│
│   ┌──────────────────────────────────┐
│   │  Ĥφ(x) = Eφ(x)    [THE TISE]    │
│   │                                    │
│   │  This is an EIGENVALUE equation!  │
│   └──────────────────────────────────┘
│   Ĥ = operator (analogous to matrix A)
│   φ = eigenfunction/eigenvector (analogous to v in Av = λv)
│   E = eigenvalue (analogous to λ)
│
├── VQE connection (CRITICAL):
│   Molecular Hamiltonian H → big matrix (2ⁿ × 2ⁿ for n qubits)
│   Finding ground state energy = finding SMALLEST eigenvalue of H
│   TISE says: Ĥφ₀ = E₀φ₀ where E₀ = minimum eigenvalue
│   VQE does: min_θ ⟨ψ(θ)|H|ψ(θ)⟩ ≥ E₀ (variational principle)
│   If ansatz is good enough: VQE result ≈ E₀ to chemical accuracy
│
└── Self-check: What does "stationary state" mean?
    |ψ(x,t)|² = |φ(x)|²·|T(t)|² = |φ(x)|² (time-independent!)
    The probability distribution doesn't change with time.
    BUT the complex phase DOES rotate: ψ(x,t) = φ(x)·e^(-iEt/ℏ)

Ph2.1.4  Particle in a Box — Complete Solution (The First Real QM Problem)
├── Setup:
│   Box of length L, infinite potential walls:
│   V(x) = 0     for 0 < x < L  (particle free inside box)
│   V(x) = ∞     for x ≤ 0 or x ≥ L  (cannot escape)
│   Boundary conditions: ψ(0) = 0, ψ(L) = 0 (walls are impenetrable)
│
├── Step 1: Write TISE inside the box (V=0):
│   -ℏ²/(2m) · d²ψ/dx² = Eψ
│   Rearrange: d²ψ/dx² = -(2mE/ℏ²)ψ
│   Define k² = 2mE/ℏ² → d²ψ/dx² = -k²ψ
│   This is the SAME equation as simple harmonic oscillator!
│   General solution: ψ(x) = A·sin(kx) + B·cos(kx)
│
├── Step 2: Apply boundary condition ψ(0) = 0:
│   ψ(0) = A·sin(0) + B·cos(0) = B = 0
│   → B = 0, so ψ(x) = A·sin(kx)
│
├── Step 3: Apply boundary condition ψ(L) = 0:
│   ψ(L) = A·sin(kL) = 0
│   A ≠ 0 (otherwise ψ=0 everywhere = no particle)
│   → sin(kL) = 0 → kL = nπ, n = 1, 2, 3, ...
│   → k = nπ/L
│   (n=0 gives ψ=0 everywhere, not physical, so n starts from 1)
│
├── Step 4: Normalize (find A):
│   ∫₀ᴸ |ψ|² dx = 1
│   ∫₀ᴸ A²·sin²(nπx/L) dx = 1
│   Using ∫₀ᴸ sin²(nπx/L) dx = L/2:     ← standard integral
│   A² · L/2 = 1 → A = √(2/L)
│
├── Result:
│   ψₙ(x) = √(2/L) · sin(nπx/L),   n = 1, 2, 3, ...
│   Eₙ = ℏ²k²/(2m) = n²π²ℏ²/(2mL²)
│
├── Energy level properties:
│   ├── E₁ = π²ℏ²/(2mL²) = ground state (lowest allowed energy, NOT zero)
│   ├── Eₙ = n² · E₁  (energies grow as squares of n)
│   ├── E₁:E₂:E₃ = 1:4:9  (quadratic spacing)
│   ├── Gap: E₂-E₁ = 3E₁ (first excitation energy)
│   ├── Smaller box (smaller L) → larger E₁ → more widely spaced levels
│   │   Atom (L≈10⁻¹⁰m): huge spacing → eV scale → visible/UV light
│   │   Room (L=10m): E₁ ≈ 10⁻⁶⁴ J → practically continuous → classical
│   └── Zero-point energy: E₁ > 0 ALWAYS → particle can never be at rest
│       (Heisenberg: confining to box → Δx=L → Δp ≥ ℏ/(2L) → KE ≥ ℏ²/(8mL²))
│
├── Nodal structure:
│   ψ₁: 0 interior nodes (just one bump)
│   ψ₂: 1 interior node (two bumps, zero crossing in middle)
│   ψₙ: (n-1) interior nodes
│   MORE nodes → MORE curvature → HIGHER KE → HIGHER E
│   This pattern holds for ALL quantum systems (atoms, molecules, etc.)
│
├── Worked numerical example:
│   Electron in box of L = 1 nm = 10⁻⁹ m
│   E₁ = π²(1.055×10⁻³⁴)²/(2×9.109×10⁻³¹×(10⁻⁹)²)
│      = π²×1.113×10⁻⁶⁸/(1.822×10⁻³⁹×10⁻¹⁸)
│      = 6.024×10⁻²⁰ J = 0.376 eV
│   E₂ = 4 × 0.376 = 1.504 eV
│   E₃ = 9 × 0.376 = 3.386 eV
│   Transition E₂→E₁: ΔE = 1.128 eV → λ = hc/ΔE = 1099 nm (infrared)
│   Transition E₃→E₁: ΔE = 3.010 eV → λ = 412 nm (visible violet!)
│
├── Code (complete particle-in-box solver):
│   import numpy as np, matplotlib.pyplot as plt
│   L, hbar, m = 1e-9, 1.055e-34, 9.109e-31  # 1nm box, electron
│   x = np.linspace(0, L, 500)
│   for n in range(1, 5):
│       E_n = (n*np.pi*hbar)**2 / (2*m*L**2)
│       E_eV = E_n / 1.602e-19
│       psi = np.sqrt(2/L) * np.sin(n*np.pi*x/L)
│       print(f"n={n}: E = {E_eV:.3f} eV, nodes = {n-1}")
│       plt.subplot(2,2,n)
│       plt.fill_between(x*1e9, psi**2/1e9, alpha=0.3, color='red')
│       plt.plot(x*1e9, psi/1e4.5, 'b-')
│       plt.title(f"n={n}, E={E_eV:.2f} eV")
│   plt.tight_layout(); plt.show()
│
├── BIO link:
│   ├── π-electrons in benzene ring ≈ particle on a ring (2D box)
│   │   Predicts UV absorption of aromatic molecules (DNA bases!)
│   ├── Quantum dots (nanocrystals): fluorescence color depends on box size L
│   │   Smaller dot = bigger L⁻² = bluer light → used in bio-imaging!
│   ├── HOMO-LUMO gap in molecules:
│   │   HOMO = highest occupied molecular orbital = ground state electron
│   │   LUMO = lowest unoccupied = first excited state
│   │   Gap = E_LUMO - E_HOMO → determines:
│   │     Absorption wavelength (UV-Vis spectroscopy)
│   │     Chemical reactivity (nucleophilic attack at LUMO)
│   │     Conductivity (band gap analogy)
│   └── Drug design: match ligand HOMO with receptor LUMO for optimal binding
│
└── Exit check:
    1. Derive ψₙ and Eₙ for particle in box (write all 4 steps from scratch)
    2. Compute E₁ for electron in L=0.5nm box. Answer: 1.504 eV
    3. Plot |ψ₃|² and mark the 2 interior nodes
    4. Explain why E₁ > 0 (zero-point energy from uncertainty principle)

Ph2.1.5  Hydrogen Atom — Preview (Solved in 3D)
├── Same idea as particle in box, but:
│   V(r) = -kₑe²/r (Coulomb potential, not infinite walls)
│   3D spherical coordinates: ψ(r,θ,φ) = R(r)·Y(θ,φ)
│   Radial part R(r): determines energy levels
│   Angular part Y(θ,φ): determines orbital shapes (s, p, d, f)
│
├── Energy levels (Bohr formula):
│   Eₙ = -13.6 eV / n²,  n = 1,2,3,...
│   E₁ = -13.6 eV (ground state, most tightly bound)
│   E₂ = -3.4 eV  (first excited)
│   E∞ = 0 eV     (ionization = electron free)
│   Ionization energy = 13.6 eV for hydrogen
│
├── Quantum numbers (3D → 3 quantum numbers):
│   n = 1,2,3,...     (principal: determines energy)
│   l = 0,1,...,n-1   (angular momentum: determines orbital shape)
│   mₗ = -l,...,0,...,l  (magnetic: determines orientation)
│   l=0: s orbital (sphere), l=1: p orbital (dumbbell), l=2: d orbital
│
├── Why you need this:
│   VQE for H₂ molecule builds on hydrogen atom orbitals
│   Molecular orbitals = linear combinations of atomic orbitals (LCAO)
│   STO-3G basis set = 3 Gaussians to approximate each hydrogen orbital
│
└── You will NOT solve the hydrogen atom from scratch.
    Just know: Eₙ = -13.6/n², quantum numbers n/l/m, orbital shapes s/p/d/f.
    The 1D particle-in-box teaches you the METHOD. Hydrogen applies it in 3D.

═══════════════════════════════════════════
 GATE TO Ph2.2 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Can write TDSE from memory: iℏ∂ψ/∂t = Ĥψ
 □ Can explain each symbol: i, ℏ, ∂/∂t, ψ, Ĥ
 □ Know: |ψ(x)|² = probability density, ∫|ψ|²dx = 1
 □ Can derive TISE from TDSE via separation of variables
 □ Solved particle-in-box: ψₙ=√(2/L)sin(nπx/L), Eₙ=n²π²ℏ²/(2mL²)
 □ Computed E₁ for electron in 1nm box (≈0.376 eV)
 □ Know: more nodes → higher energy → more curvature of ψ
 □ Know: TISE = eigenvalue equation Ĥψ=Eψ → VQE finds min eigenvalue
 □ Know hydrogen: Eₙ=-13.6/n² eV, quantum numbers n,l,mₗ
═══════════════════════════════════════════
```

---

## Module Ph2.2: Quantum Postulates

> **PREREQUISITES: Ph2.1 gate passed.**
> Must know: Schrödinger equation, ψ meaning, eigenvalue equation Ĥψ=Eψ.
> Must know: Hermitian matrices have real eigenvalues (Math M2.3).
> Must know: inner products ⟨u,v⟩ (Math M3.2).

```
Ph2.2.1  The Six Postulates — Each One Explained With Examples
├── These 6 rules ARE quantum mechanics. Everything else follows from them.
│   Think of them like Newton's 3 laws but for quantum world.
│
├── P1 STATE:
│   "A quantum system is completely described by |ψ⟩ in Hilbert space ℋ."
│   │
│   What this means in simple terms:
│   Classical: state = (position, velocity) → 2 numbers
│   Quantum: state = |ψ⟩ = column vector in ℂⁿ
│   1 qubit: |ψ⟩ = α|0⟩ + β|1⟩ → 2 complex numbers (4 real numbers)
│   2 qubits: |ψ⟩ = α|00⟩ + β|01⟩ + γ|10⟩ + δ|11⟩ → 4 complex numbers
│   n qubits: 2ⁿ complex numbers → exponentially large!
│   Constraint: |α|²+|β|²+... = 1 (normalization)
│
├── P2 OBSERVABLE:
│   "Every measurable quantity = Hermitian operator Â."
│   │
│   WHY Hermitian? Because Hermitian matrices have REAL eigenvalues.
│   Measurements give REAL numbers (you measure 3.5 eV, not 3.5+2i eV).
│   │
│   Examples you already know:
│   Energy → Ĥ (Hamiltonian)
│   Position → x̂ (multiply by x)
│   Momentum → p̂ = -iℏ d/dx
│   Spin-Z → Z = [[1,0],[0,-1]] (eigenvalues +1, -1)
│   │
│   Worked verification: Is Z Hermitian?
│   Z† = Z* transpose = [[1,0],[0,-1]]† = [[1,0],[0,-1]] = Z ✓
│
├── P3 MEASUREMENT OUTCOME:
│   "The ONLY possible results of measuring Â are eigenvalues of Â."
│   │
│   Example: measuring Z on a qubit
│   Z eigenvalues: +1 (for |0⟩) and -1 (for |1⟩)
│   You can ONLY get +1 or -1. NEVER 0.3 or 2.7. ONLY eigenvalues!
│   │
│   For energy measurement of particle in box:
│   Possible results: E₁, E₂, E₃, ... = n²π²ℏ²/(2mL²)
│   You CANNOT measure E = 2.5 × E₁. Only integer-squared multiples.
│
├── P4 STATE COLLAPSE:
│   "After measuring eigenvalue aₙ, state IMMEDIATELY becomes |aₙ⟩."
│   │
│   Before measurement: |ψ⟩ = (3/5)|0⟩ + (4/5)|1⟩ (superposition)
│   Measure Z, get +1 → state is now |0⟩. The |1⟩ component is GONE.
│   Measure Z, get -1 → state is now |1⟩. The |0⟩ component is GONE.
│   │
│   This is IRREVERSIBLE. You cannot reconstruct the original state.
│   This is why quantum computing is tricky: measurement destroys info.
│   VQE measures MANY times to reconstruct ⟨H⟩ statistically.
│
├── P5 BORN RULE (most important for calculations):
│   "P(getting aₙ) = |⟨aₙ|ψ⟩|²"
│   │
│   WORKED EXAMPLE (DO THIS):
│   |ψ⟩ = (3/5)|0⟩ + (4/5)|1⟩
│   P(|0⟩) = |⟨0|ψ⟩|² = |3/5|² = 9/25 = 0.36 = 36%
│   P(|1⟩) = |⟨1|ψ⟩|² = |4/5|² = 16/25 = 0.64 = 64%
│   Check: 0.36 + 0.64 = 1.00 ✓
│   │
│   WORKED EXAMPLE with complex amplitudes:
│   |ψ⟩ = (1/√3)|0⟩ + (i√2/√3)|1⟩
│   P(|0⟩) = |1/√3|² = 1/3 ≈ 33.3%
│   P(|1⟩) = |i√2/√3|² = |i|²·2/3 = 1·2/3 = 2/3 ≈ 66.7%
│   (Note: |i|² = 1, not i²=-1. Modulus squared, not plain squared!)
│
└── P6 TIME EVOLUTION:
    "|ψ(t)⟩ evolves by iℏ d|ψ⟩/dt = Ĥ|ψ⟩"
    Solution: |ψ(t)⟩ = e^(-iĤt/ℏ)|ψ(0)⟩
    This is a UNITARY transformation (preserves normalization)
    Between measurements: evolution is smooth, deterministic, reversible
    AT measurement: collapse is sudden, random, irreversible (P4)
    This duality is the central puzzle of quantum mechanics!

Ph2.2.2  Commutators — What They Tell You About Measurements
├── Definition: [Â,B̂] = ÂB̂ - B̂Â
│   If [Â,B̂] = 0: Â and B̂ COMMUTE → can measure both simultaneously
│   If [Â,B̂] ≠ 0: DON'T COMMUTE → measuring one disturbs the other
│
├── WORKED EXAMPLE — compute [X,Z] step by step:
│   X = [[0,1],[1,0]],  Z = [[1,0],[0,-1]]
│
│   XZ = [[0,1],[1,0]]·[[1,0],[0,-1]] = [[0,-1],[1,0]]
│   ZX = [[1,0],[0,-1]]·[[0,1],[1,0]] = [[0,1],[-1,0]]
│
│   [X,Z] = XZ - ZX = [[0,-1],[1,0]] - [[0,1],[-1,0]]
│         = [[0,-2],[2,0]] = -2·[[0,i],[-i,0]]·i ... wait, let's check:
│         = [[0,-2],[2,0]]
│   Y = [[0,-i],[i,0]]
│   -2iY = -2i·[[0,-i],[i,0]] = [[0,-2],[2,0]] ✓
│   So [X,Z] = -2iY ✓
│
├── All Pauli commutators (MEMORIZE — used daily in QC):
│   [X,Y] = 2iZ    [Y,X] = -2iZ
│   [Y,Z] = 2iX    [Z,Y] = -2iX
│   [Z,X] = 2iY    [X,Z] = -2iY
│   Pattern: cyclic (XYZ → 2i × next), anticyclic → -2i × next
│
├── Code verification:
│   import numpy as np
│   X = np.array([[0,1],[1,0]], dtype=complex)
│   Y = np.array([[0,-1j],[1j,0]])
│   Z = np.array([[1,0],[0,-1]], dtype=complex)
│   comm_XZ = X@Z - Z@X
│   print(comm_XZ)     # [[0, -2], [2, 0]]
│   print(-2j * Y)      # [[0, -2], [2, 0]]  ← same!
│
├── Anti-commutators (bonus — also useful):
│   {Â,B̂} = ÂB̂ + B̂Â
│   {X,Y} = 0, {X,Z} = 0, {Y,Z} = 0  (Paulis anti-commute!)
│   {X,X} = 2I  (XX = I → {X,X} = I+I = 2I)
│
└── VQE measurement link:
    H₂ Hamiltonian = c₁·ZZ + c₂·XX + c₃·YY + c₄·ZI + c₅·IZ + c₆·II
    To measure ⟨H⟩: need to measure each Pauli string separately
    BUT: commuting terms can be measured together in one circuit
    Group 1: {ZZ, ZI, IZ, II} → all pairwise commute → one circuit
    Group 2: {XX} → needs basis rotation → another circuit
    Group 3: {YY} → another basis rotation → another circuit
    Fewer groups = fewer circuit executions = faster VQE!

Ph2.2.3  Heisenberg Uncertainty Principle — Derived from Commutators
├── General uncertainty relation:
│   ΔA · ΔB ≥ ½|⟨[Â,B̂]⟩|
│   where ΔA = √(⟨A²⟩ - ⟨A⟩²) is the standard deviation
│
├── Position-momentum uncertainty:
│   [x̂,p̂] = iℏ (from Ph1.2.5)
│   Δx · Δp ≥ ½|⟨iℏ⟩| = ℏ/2
│   You CANNOT know both x and p precisely simultaneously!
│
├── Energy-time uncertainty:
│   ΔE · Δt ≥ ℏ/2
│   Short-lived states (small Δt) → uncertain energy (large ΔE)
│   This causes spectral line broadening in atoms/molecules
│
├── Physical examples:
│   ├── Electron in atom: Δx ≈ 0.1nm = 10⁻¹⁰m
│   │   Δp ≥ ℏ/(2Δx) = 1.055e-34/(2×10⁻¹⁰) = 5.3×10⁻²⁵ kg·m/s
│   │   Δv = Δp/mₑ = 5.3e-25/9.1e-31 = 5.8×10⁵ m/s
│   │   Electron velocity uncertain by ~10⁶ m/s → can't track trajectory!
│   │
│   └── Baseball: Δx = 1μm = 10⁻⁶m
│       Δp ≥ ℏ/(2×10⁻⁶) = 5.3×10⁻²⁹ kg·m/s
│       Δv = 5.3e-29/0.15 = 3.5×10⁻²⁸ m/s → totally negligible! Classical OK.
│
└── Quantum computing consequence:
    Cannot simultaneously sharply measure non-commuting observables
    [X,Z] ≠ 0 → measuring X on a |0⟩ state (Z-eigenstate) → random result
    This is why VQE needs separate circuits for different Pauli terms

Ph2.2.4  Expectation Values — The VQE Cost Function
├── Classical average: ⟨X⟩ = Σᵢ xᵢ · P(xᵢ)  (weighted average)
│
├── Quantum expectation (SAME idea, quantum notation):
│   ⟨Â⟩ = ⟨ψ|Â|ψ⟩ = Σₙ aₙ · P(aₙ) = Σₙ aₙ · |⟨aₙ|ψ⟩|²
│
├── WORKED EXAMPLE — compute ⟨Z⟩ for |ψ⟩ = (3/5)|0⟩ + (4i/5)|1⟩:
│   Z|0⟩ = +1·|0⟩,  Z|1⟩ = -1·|1⟩
│
│   Method 1 (eigenvalue weighted average):
│   ⟨Z⟩ = (+1)·P(|0⟩) + (-1)·P(|1⟩)
│       = (+1)·(9/25) + (-1)·(16/25)
│       = 9/25 - 16/25 = -7/25 = -0.28
│
│   Method 2 (matrix sandwich — verify):
│   ⟨ψ|Z|ψ⟩ = [3/5, -4i/5]·[[1,0],[0,-1]]·[[3/5],[4i/5]]
│           = [3/5, -4i/5]·[[3/5],[-4i/5]]
│           = (3/5)(3/5) + (-4i/5)(-4i/5)
│           = 9/25 + 16i²/25 = 9/25 - 16/25 = -7/25 ✓
│
│   Both methods give -0.28 ← check: ⟨Z⟩ is REAL (Hermitian operator) ✓
│
├── Code verification:
│   psi = np.array([3/5, 4j/5])
│   Z = np.array([[1,0],[0,-1]])
│   expval = (psi.conj() @ Z @ psi).real
│   print(f"⟨Z⟩ = {expval:.4f}")  # -0.2800
│
├── VQE cost function:
│   E(θ) = ⟨ψ(θ)|H|ψ(θ)⟩   ← this IS the VQE cost function
│   H = Σₖ cₖ Pₖ  (Hamiltonian = sum of weighted Pauli strings)
│   By linearity: ⟨H⟩ = Σₖ cₖ ⟨Pₖ⟩
│   Each ⟨Pₖ⟩ measured separately by running the circuit many times
│   Total ⟨H⟩ = weighted sum of individual Pauli expectation values
│
├── Shot noise (measurement uncertainty):
│   With N shots: ⟨Pₖ⟩_estimated = (count_+1 - count_-1) / N
│   Statistical error: σ ≈ 1/√N
│   N=1024 → error ≈ 3%, N=8192 → error ≈ 1%
│   More shots = more precise ⟨H⟩ but takes more time
│
└── Exit check:
    State |ψ⟩ = (1/√2)|0⟩ + (1/√2)|1⟩ = |+⟩
    1. Compute ⟨Z⟩ analytically. Answer: 0.
    2. Compute ⟨X⟩ analytically. Answer: +1.
    3. Compute ΔZ = √(⟨Z²⟩ - ⟨Z⟩²). Since Z²=I, ⟨Z²⟩=1, ⟨Z⟩=0 → ΔZ=1.
    4. Verify with NumPy using the matrix sandwich method.

═══════════════════════════════════════════
 GATE TO Ph2.3 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Can state all 6 postulates in your own words
 □ Know: observable = Hermitian operator, outcome = eigenvalue
 □ Can compute Born rule probabilities with complex amplitudes
 □ Can calculate commutator [X,Z] = -2iY by matrix multiplication
 □ Know: [Â,B̂]≠0 → can't measure both precisely → uncertainty
 □ Can compute ⟨Z⟩ for any qubit state by BOTH methods (weighted avg + sandwich)
 □ Know: VQE cost = Σₖcₖ⟨Pₖ⟩, each measured in separate circuit group
 □ Understand shot noise: σ ≈ 1/√N
═══════════════════════════════════════════
```

---

## Module Ph2.3: Dirac Notation ⛔ BLOCKER for All QC

> **PREREQUISITES: Ph2.2 gate passed.**
> Must know: inner products (Math M3.2), matrix multiplication (M2.2),
> conjugate transpose A† (Math M2.2.3), Pauli matrices (Math M2.3.3).
> This is the LANGUAGE of quantum computing. Every equation uses it.

```
Ph2.3.1  Ket |ψ⟩ — The Quantum State Vector
├── A ket |ψ⟩ is just a COLUMN VECTOR with a fancy name
│   |0⟩ = [[1],[0]]   (qubit in state 0)
│   |1⟩ = [[0],[1]]   (qubit in state 1)
│   |+⟩ = (1/√2)[[1],[1]]  (equal superposition)
│   |-⟩ = (1/√2)[[1],[-1]] (equal superposition, negative phase)
│
├── Why the notation?
│   Dirac invented |⟩ bracket to make quantum equations look clean.
│   Instead of writing v₁ = [1, 0]ᵀ → write |0⟩
│   Instead of writing c₁v₁ + c₂v₂ → write α|0⟩ + β|1⟩
│   Prettier, faster to read, universally used in quantum.
│
├── General qubit state:
│   |ψ⟩ = α|0⟩ + β|1⟩ = [[α],[β]]
│   α, β ∈ ℂ (complex numbers), |α|²+|β|² = 1
│
└── Code:
    import numpy as np
    ket0 = np.array([[1],[0]], dtype=complex)
    ket1 = np.array([[0],[1]], dtype=complex)
    ket_plus = (ket0 + ket1) / np.sqrt(2)
    print(ket_plus)  # [[0.707], [0.707]]

Ph2.3.2  Bra ⟨ψ| — The Dual (Conjugate Transpose)
├── Bra = conjugate transpose of ket:
│   ⟨ψ| = |ψ⟩† = (|ψ⟩*)ᵀ
│   ⟨0| = [1, 0]    (row vector)
│   ⟨1| = [0, 1]
│
├── For complex state |ψ⟩ = [[1/√2],[i/√2]]:
│   ⟨ψ| = [(1/√2)*, (i/√2)*] = [1/√2, -i/√2]
│   Notice: i → -i (conjugate!)
│
└── Code:
    ket_psi = np.array([[1/np.sqrt(2)], [1j/np.sqrt(2)]])
    bra_psi = ket_psi.conj().T  # conjugate transpose
    print(bra_psi)  # [[0.707, -0.707j]]

Ph2.3.3  Inner Product ⟨φ|ψ⟩ — Bra Times Ket = Scalar
├── ⟨φ|ψ⟩ = bra(row) × ket(column) = dot product = single complex number
│
├── WORKED EXAMPLES (DO ALL BY HAND):
│   ⟨0|0⟩ = [1,0] · [[1],[0]] = 1 (normalized ✓)
│   ⟨1|1⟩ = [0,1] · [[0],[1]] = 1 (normalized ✓)
│   ⟨0|1⟩ = [1,0] · [[0],[1]] = 0 (orthogonal ✓)
│   ⟨1|0⟩ = [0,1] · [[1],[0]] = 0 (orthogonal ✓)
│   ⟨+|+⟩ = (1/√2)[1,1] · (1/√2)[[1],[1]] = ½(1+1) = 1 ✓
│   ⟨+|-⟩ = (1/√2)[1,1] · (1/√2)[[1],[-1]] = ½(1-1) = 0 ✓
│
├── With complex amplitudes:
│   |ψ⟩ = [[1/√2],[i/√2]], |φ⟩ = [[1],[0]] = |0⟩
│   ⟨φ|ψ⟩ = ⟨0|ψ⟩ = [1,0]·[[1/√2],[i/√2]] = 1/√2
│   |⟨0|ψ⟩|² = |1/√2|² = 1/2 → 50% chance of measuring |0⟩ ← Born rule!
│
├── Key properties:
│   ⟨φ|ψ⟩ = ⟨ψ|φ⟩*  (swap = conjugate)
│   ⟨ψ|ψ⟩ = ||ψ||² = real, non-negative (norm squared)
│   ⟨ψ|ψ⟩ = 1 for normalized states
│
└── Code:
    inner = bra_psi @ ket0  # ⟨ψ|0⟩
    prob = np.abs(inner)**2  # Born rule probability
    print(f"⟨ψ|0⟩ = {inner.item()}, P(|0⟩) = {prob.item():.3f}")

Ph2.3.4  Outer Product |ψ⟩⟨φ| — Ket Times Bra = Matrix!
├── |ψ⟩⟨φ| = column × row = MATRIX (operator)
│   This is how you BUILD operators from states
│
├── WORKED EXAMPLES:
│   |0⟩⟨0| = [[1],[0]] · [1,0] = [[1,0],[0,0]]  (projection onto |0⟩)
│   |1⟩⟨1| = [[0],[1]] · [0,1] = [[0,0],[0,1]]  (projection onto |1⟩)
│   |0⟩⟨1| = [[1],[0]] · [0,1] = [[0,1],[0,0]]  (transition |1⟩ → |0⟩)
│   |1⟩⟨0| = [[0],[1]] · [1,0] = [[0,0],[1,0]]  (transition |0⟩ → |1⟩)
│
├── Completeness relation:
│   |0⟩⟨0| + |1⟩⟨1| = [[1,0],[0,0]] + [[0,0],[0,1]] = [[1,0],[0,1]] = I ✓
│   This says: the basis states "fill" the whole space
│
├── Building Pauli X from outer products:
│   X = |0⟩⟨1| + |1⟩⟨0| = [[0,1],[0,0]] + [[0,0],[1,0]] = [[0,1],[1,0]] ✓
│
└── Code:
    proj0 = ket0 @ ket0.conj().T  # |0⟩⟨0|
    proj1 = ket1 @ ket1.conj().T  # |1⟩⟨1|
    print(proj0 + proj1)  # Identity matrix ✓

Ph2.3.5  The Sandwich ⟨ψ|Â|ψ⟩ — Expectation Value
├── This is the MOST IMPORTANT expression in all of VQE:
│   ⟨ψ|Â|ψ⟩ = bra × operator × ket = scalar
│
├── How to evaluate (two equivalent methods):
│   Method 1: First compute |φ⟩ = Â|ψ⟩, then compute ⟨ψ|φ⟩
│   Method 2: Eigendecomposition: ⟨Â⟩ = Σₙ aₙ |⟨aₙ|ψ⟩|² (from Ph2.2.4)
│
├── WORKED EXAMPLE (Method 1):
│   |ψ⟩ = (1/√2)|0⟩ + (1/√2)|1⟩ = |+⟩
│   Compute Â|ψ⟩ = Z|+⟩:
│   Z|+⟩ = (1/√2)Z|0⟩ + (1/√2)Z|1⟩ = (1/√2)|0⟩ - (1/√2)|1⟩ = |-⟩
│   ⟨ψ|Z|ψ⟩ = ⟨+|-⟩ = 0  (they're orthogonal!)
│   ⟨Z⟩ = 0 for |+⟩ state → equally likely to get +1 or -1
│
├── Spectral decomposition (Dirac form of any Hermitian operator):
│   Â = Σₙ aₙ |aₙ⟩⟨aₙ|  (eigenvalue × projection onto eigenvector)
│   Z = (+1)|0⟩⟨0| + (-1)|1⟩⟨1|
│   H_hamiltonian = E₀|E₀⟩⟨E₀| + E₁|E₁⟩⟨E₁| + ...
│   VQE finds: E₀ = minimum eigenvalue in this decomposition
│
└── Code:
    psi_plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    # Method 1:
    phi = Z @ psi_plus       # Z|+⟩
    result = psi_plus.conj() @ phi  # ⟨+|Z|+⟩
    print(f"⟨Z⟩ = {result.real:.4f}")  # 0.0000

Ph2.3.6  Multi-Qubit States — Tensor Product Notation
├── Two-qubit computational basis:
│   |00⟩ = |0⟩⊗|0⟩ = [1,0,0,0]ᵀ
│   |01⟩ = |0⟩⊗|1⟩ = [0,1,0,0]ᵀ
│   |10⟩ = |1⟩⊗|0⟩ = [0,0,1,0]ᵀ
│   |11⟩ = |1⟩⊗|1⟩ = [0,0,0,1]ᵀ
│
├── Product (separable) state — CAN be factored:
│   |+0⟩ = |+⟩⊗|0⟩ = (1/√2)(|00⟩ + |10⟩) = (1/√2)[1,0,1,0]ᵀ
│   This IS a product: (1/√2)(|0⟩+|1⟩) ⊗ |0⟩
│
├── Entangled state — CANNOT be factored:
│   |Φ+⟩ = (1/√2)(|00⟩ + |11⟩) = (1/√2)[1,0,0,1]ᵀ
│
│   Proof it can't factor:
│   Assume |Φ+⟩ = (α|0⟩+β|1⟩)⊗(γ|0⟩+δ|1⟩)
│   = αγ|00⟩ + αδ|01⟩ + βγ|10⟩ + βδ|11⟩
│   Comparing: αγ = 1/√2, αδ = 0, βγ = 0, βδ = 1/√2
│   αδ=0 → either α=0 or δ=0
│   If α=0 → αγ=0 ≠ 1/√2 → contradiction!
│   If δ=0 → βδ=0 ≠ 1/√2 → contradiction!
│   → Cannot factor. |Φ+⟩ is GENUINELY entangled.
│
├── The 4 Bell states (maximally entangled, MEMORIZE):
│   |Φ+⟩ = (1/√2)(|00⟩ + |11⟩)  Circuits: H(0), CNOT(0→1)
│   |Φ-⟩ = (1/√2)(|00⟩ - |11⟩)  Circuits: X(1), H(0), CNOT(0→1)
│   |Ψ+⟩ = (1/√2)(|01⟩ + |10⟩)  Circuits: X(0), H(0), CNOT(0→1)
│   |Ψ-⟩ = (1/√2)(|01⟩ - |10⟩)  Circuits: X(0), X(1), H(0), CNOT(0→1)
│   All 4: orthonormal, maximally entangled, form a complete 2-qubit basis
│
├── Code (Bell state creation + verification):
│   ket00 = np.array([1,0,0,0], dtype=complex)
│   ket11 = np.array([0,0,0,1], dtype=complex)
│   bell = (ket00 + ket11) / np.sqrt(2)
│   print(bell)  # [0.707, 0, 0, 0.707]
│   print(np.abs(bell)**2)  # [0.5, 0, 0, 0.5] → P(|00⟩)=P(|11⟩)=50%
│
└── Why entanglement matters for VQE:
    H₂ molecule ground state has electron CORRELATION
    Hartree-Fock (product state) misses correlation energy ≈ 0.04 Ha
    VQE ansatz with entangling gates → captures correlation
    Without entanglement → VQE cannot reach chemical accuracy

Ph2.3.7  Reading Quantum Equations — The Ultimate Test
├── Can you read this expression and explain EVERY symbol?
│
│   ⟨ψ(θ)|Ĥ|ψ(θ)⟩
│
│   |ψ(θ)⟩: ket, parameterized quantum state (column vector)
│           depends on angles θ = (θ₁,...,θₙ) in the ansatz circuit
│   ⟨ψ(θ)|: bra, conjugate transpose of |ψ(θ)⟩ (row vector)
│   Ĥ: Hamiltonian operator (Hermitian matrix)
│      for molecules: Ĥ = Σₖ cₖ Pₖ (sum of Pauli strings)
│   ⟨ψ(θ)|Ĥ|ψ(θ)⟩: sandwich = expectation value = real number
│                     = "average energy measured for this state"
│
├── The VQE algorithm in one line:
│   E_ground ≈ min_θ ⟨ψ(θ)|Ĥ|ψ(θ)⟩
│   "Find θ that minimizes the energy expectation value."
│
└── FINAL EXIT CHECK (THE DIRAC FLUENCY EXAM):
    Given |ψ⟩ = (√3/2)|0⟩ + (1/2)e^(iπ/4)|1⟩

    1. Write ⟨ψ| (conjugate transpose). Pay attention to e^(iπ/4)→e^(-iπ/4).
    2. Compute ⟨ψ|ψ⟩. Is it 1? → Check normalization.
       |√3/2|²+|1/2|² = 3/4+1/4 = 1 ✓
    3. Compute ⟨Z⟩ = ⟨ψ|Z|ψ⟩.
       = (3/4)(+1) + (1/4)(-1) = 3/4 - 1/4 = 1/2
    4. Compute P(|0⟩) = |⟨0|ψ⟩|² = 3/4 = 75%.
    5. Compute ⟨X⟩ = ⟨ψ|X|ψ⟩. (Hint: X|0⟩=|1⟩, X|1⟩=|0⟩)
       = (√3/2)·(1/2)e^(iπ/4) + (1/2)e^(-iπ/4)·(√3/2)
       = (√3/2)(1/2)[e^(iπ/4)+e^(-iπ/4)] = (√3/2)(1/2)·2cos(π/4)
       = (√3/2)(1/2)(√2) = √6/4 ≈ 0.612
    6. Verify ALL of the above in NumPy. Must match to 10 decimal places.

═══════════════════════════════════════════
 GATE TO PHASE 2 (QC THEORY) — MASTER PHYSICS GATE:
═══════════════════════════════════════════
 □ Can write and read any Dirac expression: ket, bra, inner, outer, sandwich
 □ Know |0⟩,|1⟩,|+⟩,|-⟩ as vectors and can convert between them
 □ Can compute ⟨ψ|Â|ψ⟩ for any 1-qubit state and Pauli operator
 □ Know inner product ⟨φ|ψ⟩ = overlap, |⟨φ|ψ⟩|² = probability
 □ Can build operators from outer products: X = |0⟩⟨1|+|1⟩⟨0|
 □ Know completeness: |0⟩⟨0|+|1⟩⟨1| = I
 □ Can prove |Φ+⟩ is entangled (cannot factor)
 □ Know all 4 Bell states by name and formula
 □ Can read ⟨ψ(θ)|Ĥ|ψ(θ)⟩ and explain every symbol
 □ Passed the FINAL EXIT CHECK (computed ⟨Z⟩, ⟨X⟩, P(|0⟩) for complex state)
═══════════════════════════════════════════
```


---

# TO-DO LIST — PART 2 (PHYSICS PHASE)
> Har topic complete karne ke baad check karo. Ek bhi miss mat karna.

## Phase 1 — Classical to Quantum Foundation

### Ph1.1 Classical Mechanics & Energy
- [ ] Ph1.1.1  Newton's laws in own words (F=ma, inertia, action-reaction)
- [ ] Ph1.1.2  Work = F·d·cosθ; compute for F=10N, d=5m, θ=60°
- [ ] Ph1.1.3  KE = ½mv²; compute for m=2kg, v=3m/s
- [ ] Ph1.1.4  PE = mgh; compute for m=1kg, h=10m
- [ ] Ph1.1.5  Conservation: KE+PE = constant (no friction)
- [ ] Ph1.1.6  Coulomb potential: V(r)=-e²/r (electron in H atom)
- [ ] Ph1.1.7  Total energy: E = KE + PE = p²/2m - e²/r
- [ ] Ph1.1.8  H atom energy: Eₙ=-13.6/n² eV; E₁, E₂, E₃ computed
- [ ] Ph1.1.9  H₂ energy is the sum that VQE must find (explain conceptually)
- [ ] Ph1.1 GATE — passed ✓

### Ph1.2 Hamiltonian Mechanics
- [ ] Ph1.2.1  Lagrangian: L=KE-PE=T-V
- [ ] Ph1.2.2  Euler-Lagrange: d/dt(∂L/∂q̇)-∂L/∂q=0
- [ ] Ph1.2.3  Derive F=ma from EL for L=½mẋ²-V(x)
- [ ] Ph1.2.4  Conjugate momentum: p=∂L/∂q̇ (general)
- [ ] Ph1.2.5  Legendre transform: H=pq̇-L (Hamiltonian)
- [ ] Ph1.2.6  Hamilton's equations: q̇=∂H/∂p, ṗ=-∂H/∂q
- [ ] Ph1.2.7  Solve harmonic oscillator with Hamilton's equations
- [ ] Ph1.2.8  Poisson bracket: {f,g}=Σ(∂f/∂qₖ·∂g/∂pₖ - ∂f/∂pₖ·∂g/∂qₖ)
- [ ] Ph1.2.9  {x,p}=1 (canonical commutation relation, classical)
- [ ] Ph1.2.10 Transition: {·,·} → [·,·]/iℏ (quantum commutator)
- [ ] Ph1.2 GATE — passed ✓

### Ph1.3 Wave Mechanics
- [ ] Ph1.3.1  Wave parameters: wavelength λ, frequency f, speed v=fλ, period T
- [ ] Ph1.3.2  Wave function: y(x,t)=A·sin(kx-ωt); verify is wave equation solution
- [ ] Ph1.3.3  Superposition: add two waves; constructive/destructive interference
- [ ] Ph1.3.4  Standing waves: ψₙ(x)=A·sin(nπx/L); nodes at boundaries
- [ ] Ph1.3.5  de Broglie: λ=h/p (matter waves); compute for electron at 1eV
- [ ] Ph1.3.6  Planck: E=hf (quanta of energy)
- [ ] Ph1.3.7  Photoelectric effect: Eₖ=hf-φ; explains it needs min frequency
- [ ] Ph1.3.8  Atom spectra: Bohr model; ΔE=(13.6)(1/n₁²-1/n₂²) eV
- [ ] Ph1.3.9  UV light causes DNA damage: λ<320nm photons break bonds
- [ ] Ph1.3 GATE — passed ✓

## Phase 2 — Quantum Mechanics Core

### Ph2.1 Schrödinger Equation
- [ ] Ph2.1.1  TDSE: iℏ∂ψ/∂t = Ĥψ — write it and name every symbol
- [ ] Ph2.1.2  Ĥ = -ℏ²/2m · ∂²/∂x² + V(x) (Hamiltonian operator)
- [ ] Ph2.1.3  Kinetic energy operator: p̂→-iℏ∂/∂x; KE=p̂²/2m
- [ ] Ph2.1.4  TISE: Ĥψ=Eψ (eigenvalue equation — time-independent)
- [ ] Ph2.1.5  Born rule: |ψ(x)|² = probability DENSITY
- [ ] Ph2.1.6  ∫|ψ|²dx=1 (normalization) — what this means physically
- [ ] Ph2.1.7  Separation of variables: ψ(x,t)=φ(x)·e^(-iEt/ℏ) — derive TISE
- [ ] Ph2.1.8  Particle in box: ψₙ=√(2/L)sin(nπx/L), Eₙ=n²π²ℏ²/2mL²
- [ ] Ph2.1.9  Compute E₁ for electron in 1nm box
- [ ] Ph2.1.10 Normalization of ψ₁: ∫₀ᴸ|ψ₁|²dx=1 verified analytically
- [ ] Ph2.1.11 P(x∈[0,L/2]) for ground state = 1/2 (by symmetry)
- [ ] Ph2.1 GATE — passed ✓

### Ph2.2 Quantum Postulates
- [ ] Ph2.2.1  Postulate 1: state = ket |ψ⟩ in Hilbert space
- [ ] Ph2.2.2  Postulate 2: observables = Hermitian operators
- [ ] Ph2.2.3  Postulate 3: measurement outcomes = eigenvalues
- [ ] Ph2.2.4  Postulate 4: Born rule P(aₙ)=|⟨aₙ|ψ⟩|²
- [ ] Ph2.2.5  Postulate 5: state collapse after measurement
- [ ] Ph2.2.6  Postulate 6: time evolution = iℏd|ψ⟩/dt=Ĥ|ψ⟩
- [ ] Ph2.2.7  Compute [X̂,Ẑ]=? (commutator XZ-ZX) — write out matrix product
- [ ] Ph2.2.8  [X̂,Ẑ]=2iŶ — verify in NumPy
- [ ] Ph2.2.9  [x̂,p̂]=iℏ (canonical commutation relation, quantum)
- [ ] Ph2.2.10 Uncertainty: ΔA·ΔB ≥ |⟨[Â,B̂]⟩|/2
- [ ] Ph2.2.11 Δx·Δp ≥ ℏ/2 (Heisenberg) — compute for electron vs baseball
- [ ] Ph2.2.12 Expectation ⟨A⟩=⟨ψ|A|ψ⟩ — compute ⟨Z⟩ for |+⟩=0
- [ ] Ph2.2.13 Shot noise: σ≈1/√N; N=1024 shots → 3% error
- [ ] Ph2.2 GATE — passed ✓

### Ph2.3 Dirac Notation
- [ ] Ph2.3.1  Ket |ψ⟩: column vector, state of system
- [ ] Ph2.3.2  Bra ⟨ψ|: row vector, dagger of ket
- [ ] Ph2.3.3  Inner product ⟨φ|ψ⟩ = overlap (complex number)
- [ ] Ph2.3.4  |⟨φ|ψ⟩|² = probability of measuring |φ⟩ in state |ψ⟩
- [ ] Ph2.3.5  ⟨0|ψ⟩=α and ⟨1|ψ⟩=β for |ψ⟩=α|0⟩+β|1⟩
- [ ] Ph2.3.6  Outer product |ψ⟩⟨φ| = matrix (operator)
- [ ] Ph2.3.7  |0⟩⟨0| as projector: compute 2×2 matrix
- [ ] Ph2.3.8  X = |0⟩⟨1|+|1⟩⟨0| verification in NumPy
- [ ] Ph2.3.9  Expectation sandwich: ⟨ψ|Â|ψ⟩ — compute ⟨Z⟩ for |+⟩
- [ ] Ph2.3.10 Spectral: Â=Σaₙ|aₙ⟩⟨aₙ|; Z=(+1)|0⟩⟨0|+(-1)|1⟩⟨1|
- [ ] Ph2.3.11 2-qubit: |00⟩,|01⟩,|10⟩,|11⟩ as 4-component vectors
- [ ] Ph2.3.12 Prove |Φ+⟩=(1/√2)(|00⟩+|11⟩) is entangled (cannot factor)
- [ ] Ph2.3.13 All 4 Bell states from memory (|Φ+⟩,|Φ-⟩,|Ψ+⟩,|Ψ-⟩)
- [ ] Ph2.3.14 Read ⟨ψ(θ)|Ĥ|ψ(θ)⟩ — explain EVERY symbol
- [ ] Ph2.3.15 Final Dirac exam: |ψ⟩=(√3/2)|0⟩+(1/2)e^(iπ/4)|1⟩
  - [ ]   Write ⟨ψ|
  - [ ]   Compute ⟨ψ|ψ⟩ = 1 (normalization check)
  - [ ]   Compute ⟨Z⟩ = 1/2
  - [ ]   Compute P(|0⟩) = 3/4
  - [ ]   Compute ⟨X⟩ = √6/4 ≈ 0.612
  - [ ]   Verify all in NumPy to 10 decimal places
- [ ] Ph2.3 GATE — passed ✓

---

## ⭐ MASTER SIGN-OFF — PART 2

- [ ] All 6 physics module gates passed (Ph1.1→Ph2.3)
- [ ] Can write TDSE and TISE from memory
- [ ] All 6 quantum postulates listed from memory
- [ ] Can prove |Φ+⟩ is entangled (cannot factor proof)
- [ ] All 4 Bell states written from memory
- [ ] Final Dirac exam: ⟨Z⟩=1/2, ⟨X⟩=√6/4, P(|0⟩)=3/4 verified
- [ ] **READY FOR PHASE 2 — QUANTUM COMPUTING THEORY 🚀**
