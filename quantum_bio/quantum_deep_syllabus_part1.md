# ⚛️ Quantum Bioinformatics — Deep Chapter-wise Syllabus
> **Level:** 10th Grade Math/Physics → VQE + Quantum ML for Genetic Engineering
> **Format:** Every Node → Sub-chapters → What/How/Code/Exit Criteria

---

# PHASE 0: BRIDGE SETUP (Weeks 1-2)

---

## Module P0.1: Complex Numbers — 10th to Quantum Bridge

```
P0.1.1  Real → Imaginary → Complex Number System
├── What to learn:
│   ├── Real number line review
│   ├── The problem: √(-1) has no real solution
│   ├── Imaginary unit: i²= -1, so √(-1) = i
│   ├── Complex number: z = a + bi (a = real part, b = imaginary part)
│   └── The complex plane (Argand diagram): x-axis = real, y-axis = imaginary
│
├── Key insight:
│   │  Real numbers = 1D number line
│   │  Complex numbers = 2D plane
│   │  Quantum amplitudes live on this 2D plane
│
└── How to learn:
    ├── Khan Academy: "Intro to complex numbers" (30 min)
    └── Draw: Plot 4 complex numbers on Argand diagram by hand

P0.1.2  Complex Arithmetic
├── What to learn:
│   ├── Addition: (a+bi) + (c+di) = (a+c) + (b+d)i
│   ├── Subtraction: (a+bi) - (c+di) = (a-c) + (b-d)i
│   ├── Multiplication: (a+bi)(c+di) = (ac-bd) + (ad+bc)i
│   └── Division: multiply numerator + denominator by conjugate
│
├── Code:
│   z1 = 3 + 4j         # Python uses j not i
│   z2 = 2 - 1j
│   print(z1 * z2)      # (-10+5j)? verify by hand
│
└── Exit check: Compute (3+4i)(2-i) by hand AND in Python. Must match.

P0.1.3  Conjugate, Modulus, Normalization
├── What to learn:
│   ├── Conjugate: z* = a - bi (flip sign of imaginary part)
│   ├── Modulus: |z| = √(a²+b²)  ← this is the "length" of z
│   ├── |z|² = z * z* = a² + b²
│   └── Normalized: |z|² = 1 (quantum states require this!)
│
├── Quantum link:
│   Qubit state |ψ⟩ = α|0⟩ + β|1⟩
│   Normalization condition: |α|² + |β|² = 1
│   α and β are complex numbers whose moduli² sum to 1
│
└── Exit check:
    Represent |ψ⟩ = (1/√2)|0⟩ + (i/√2)|1⟩ as column vector.
    Verify: |1/√2|² + |i/√2|² = 1. Show work.
```

---

# PHASE 1: MATHEMATICS (Weeks 3-10)

---

## Module M1.1: Complex Numbers — Full Mastery

```
M1.1.1  Polar Form
├── What to learn:
│   ├── Every complex number has two representations:
│   │   Rectangular: z = a + bi
│   │   Polar:       z = r·(cosθ + i·sinθ) = r∠θ
│   ├── Radius r = |z| = √(a²+b²)
│   ├── Angle θ = arctan(b/a)  (use atan2 in code)
│   └── Converting: a = r·cosθ, b = r·sinθ
│
├── Key insight:
│   │  SAME number, different clothing.
│   │  Rectangular → easy to add
│   │  Polar → easy to multiply and rotate
│
├── Code:
│   import numpy as np
│   z = -1 + 1j
│   r = np.abs(z)            # modulus
│   theta = np.angle(z)      # angle in radians
│   z_back = r * np.exp(1j * theta)  # reconstruct
│   print(np.isclose(z, z_back))     # True
│
└── Exit check: Convert (-1+i) to polar form by hand. Verify with NumPy.

M1.1.2  Euler's Formula — Most Important Equation in Quantum
├── What to learn:
│   ├── e^(iθ) = cosθ + i·sinθ  ← Euler's formula
│   ├── Proof: compare Taylor series of e^x, cos(x), sin(x)
│   │   e^x  = 1 + x + x²/2! + x³/3! + ...
│   │   cos  = 1 - x²/2! + x⁴/4! - ...
│   │   sin  = x - x³/3! + x⁵/5! - ...
│   │   So e^(ix) = cos(x) + i·sin(x) ✓
│   ├── Famous cases:
│   │   e^(iπ) = -1       (Euler's identity)
│   │   e^(i·0) = 1
│   │   e^(iπ/2) = i      (rotate 90°)
│   └── Compact polar form: z = r·e^(iθ)
│
├── Quantum link:
│   Every quantum gate uses e^(iθ) internally:
│   Ry(θ) = [[cos(θ/2), -sin(θ/2)], [sin(θ/2), cos(θ/2)]]
│   Written compactly using e^(iθ) in rotation groups
│
└── Exit check: Derive e^(iπ/4) in rectangular form. Show Taylor series sketch.

M1.1.3  Multiplication = Rotation
├── What to learn:
│   ├── Multiplying z₁·z₂ in polar:
│   │   r₁e^(iθ₁) · r₂e^(iθ₂) = r₁r₂ · e^(i(θ₁+θ₂))
│   │   → Magnitudes multiply, angles ADD
│   ├── Multiplying by e^(iφ) rotates z by angle φ
│   ├── De Moivre: (e^(iθ))^n = e^(inθ)
│   └── Roots of unity: z^n = 1 → n solutions on unit circle
│
├── Code:
│   def rotate_complex(z, theta):
│       return z * np.exp(1j * theta)
│   # Test:
│   print(rotate_complex(1+0j, np.pi/2))  # should ≈ 0+1j
│
└── Exit check: Implement rotate_complex(). Verify 4 rotations of π/2 return to start.
```

---

## Module M1.2: Calculus — Derivatives & Partial Derivatives

```
M1.2.1  What is a Derivative (Intuition First)
├── What to learn:
│   ├── Function f(x): machine that converts input x to output f(x)
│   ├── Average rate of change: [f(x+h) - f(x)] / h
│   ├── Derivative = instantaneous rate → limit as h→0
│   ├── Notation: f'(x), df/dx, d/dx[f(x)]
│   └── Geometric meaning: slope of tangent line at point x
│
├── Key insight:
│   │  "How fast is f changing at this exact point?"
│   │  Positive f' → function going up
│   │  Negative f' → function going down
│   │  f' = 0      → peak or valley (minimum or maximum)
│
└── How to learn:
    ├── 3Blue1Brown: "Essence of Calculus" Ep 1-2 (30 min)
    └── Intuition video BEFORE any formula!

M1.2.2  Differentiation Rules
├── What to learn:
│   ├── Power rule:    d/dx[xⁿ] = n·xⁿ⁻¹
│   ├── Constant:      d/dx[c] = 0
│   ├── Sum rule:      d/dx[f+g] = f' + g'
│   ├── Product rule:  d/dx[f·g] = f'g + fg'
│   ├── Chain rule:    d/dx[f(g(x))] = f'(g(x)) · g'(x)
│   └── Exponential:   d/dx[e^x] = e^x, d/dx[e^(ax)] = a·e^(ax)
│
├── Code (numerical verification):
│   def numerical_deriv(f, x, h=1e-7):
│       return (f(x+h) - f(x)) / h
│   # Example: f(x) = x³
│   f = lambda x: x**3
│   # Analytical: f'(2) = 3·4 = 12
│   print(numerical_deriv(f, 2))  # should ≈ 12.0
│
└── Exit check: Differentiate f(x) = x³ + 2x² - 5x. Answer: 3x² + 4x - 5.

M1.2.3  Partial Derivatives — Multi-variable Functions
├── What to learn:
│   ├── f(x,y): function of TWO variables
│   ├── ∂f/∂x = "derivative of f treating y as constant"
│   ├── ∂f/∂y = "derivative of f treating x as constant"
│   ├── Example: f(x,y) = x²y + 3y²
│   │   ∂f/∂x = 2xy      (treat y as constant)
│   │   ∂f/∂y = x² + 6y  (treat x as constant)
│   └── Gradient ∇f = [∂f/∂x, ∂f/∂y] = vector of all partials
│
├── Quantum/ML link:
│   VQE energy E(θ₁, θ₂, ..., θₙ) is function of n parameters
│   To optimize: compute ∂E/∂θₖ for every k
│   This IS gradient descent from your ML track!
│   θₖ_new = θₖ - α · ∂E/∂θₖ
│
├── Code:
│   # Numerical partial derivative
│   def partial_x(f, x, y, h=1e-7):
│       return (f(x+h, y) - f(x, y)) / h
│   f = lambda x, y: x**2 * y + 3*y**2
│   print(partial_x(f, 2, 3))  # should ≈ 2*2*3 = 12
│
└── Exit check:
    L(w) = (wx - y)². Compute ∂L/∂w analytically.
    Verify numerically: (L(w+0.001) - L(w)) / 0.001. Must match.

M1.2.4  Integration (Quantum Probability Context)
├── What to learn:
│   ├── Integration = "area under the curve"
│   ├── Anti-derivative: if F'(x) = f(x), then ∫f(x)dx = F(x) + C
│   ├── Definite integral: ∫[a to b] f(x)dx = F(b) - F(a)
│   ├── Key integrals:
│   │   ∫xⁿdx = xⁿ⁺¹/(n+1) + C
│   │   ∫e^x dx = e^x + C
│   │   ∫sin(x)dx = -cos(x) + C
│   └── u-substitution (chain rule reversed)
│
├── Quantum link:
│   P(particle in region a→b) = ∫[a to b] |ψ(x)|² dx
│   Normalization: ∫[-∞ to ∞] |ψ(x)|² dx = 1
│   Expectation value: ⟨E⟩ = ∫ ψ*(x) H ψ(x) dx
│
└── Exit check:
    Compute ∫(3x²+2)dx from 0 to 1. Answer: [x³+2x] from 0→1 = 3.
```

---

## Module M1.3: Ordinary Differential Equations (ODE)

```
M1.3.1  First Order ODE — Separation of Variables
├── What to learn:
│   ├── ODE: equation involving function AND its derivative
│   ├── dy/dt = ky  (simplest ODE)
│   ├── Separate: dy/y = k·dt
│   ├── Integrate both sides: ln|y| = kt + C
│   └── Solution: y(t) = y(0)·e^(kt)
│
├── Examples:
│   k > 0: exponential growth (bacteria, viral replication)
│   k < 0: exponential decay (radioactive decay)
│   k = -iE/ℏ: QUANTUM TIME EVOLUTION ← this is the key case
│
└── Exit check: Solve dy/dt = -3y, y(0)=2. Answer: y(t) = 2e^(-3t).

M1.3.2  Complex ODE → Schrödinger Connection
├── What to learn:
│   ├── What if k is imaginary? k = -iE/ℏ
│   ├── Solution: ψ(t) = ψ(0)·e^(-iEt/ℏ)
│   ├── Magnitude: |e^(-iEt/ℏ)| = 1 always
│   │   → Total probability conserved! (physical requirement)
│   ├── This IS the time-dependent Schrödinger equation:
│   │   iℏ·dψ/dt = E·ψ
│   └── TISE: Ĥψ = Eψ (eigenvalue equation)
│
├── Quantum link:
│   Every quantum state evolves as e^(-iEt/ℏ)
│   Superposition: Ψ(t) = c₁ψ₁e^(-iE₁t/ℏ) + c₂ψ₂e^(-iE₂t/ℏ)
│
└── Exit check:
    Verify by substitution: ψ(t) = ψ(0)·e^(-iEt/ℏ) satisfies iℏ(dψ/dt) = Eψ.
    Show step by step.
```

---

## Module M2.1: Vectors & Vector Spaces

```
M2.1.1  What is a Vector
├── What to learn:
│   ├── Vector = ordered list of numbers: v = [v₁, v₂, ..., vₙ]
│   ├── Geometric interpretation: arrow from origin to point
│   ├── Column vector (used in quantum): [1, 0]ᵀ
│   ├── Operations:
│   │   Addition: u+v = [u₁+v₁, u₂+v₂]
│   │   Scalar mul: αv = [αv₁, αv₂]
│   └── Zero vector: 0 = [0,0,...,0]
│
├── Quantum link:
│   |0⟩ = [1, 0]ᵀ   (qubit basis state)
│   |1⟩ = [0, 1]ᵀ   (qubit basis state)
│   |ψ⟩ = α[1,0]ᵀ + β[0,1]ᵀ = [α, β]ᵀ  ← superposition!
│
└── Code:
    import numpy as np
    ket0 = np.array([1, 0], dtype=complex)
    ket1 = np.array([0, 1], dtype=complex)
    psi = (1/np.sqrt(2))*ket0 + (1/np.sqrt(2))*ket1  # |+⟩ state

M2.1.2  Linear Combinations, Span, Linear Independence
├── What to learn:
│   ├── Linear combination: αu + βv + γw + ...
│   ├── Span: ALL possible linear combinations of a set
│   │   Span{[1,0],[0,1]} = all of ℝ² (any 2D vector)
│   ├── Linear independence: no vector = combo of others
│   │   {[1,0],[0,1]}: independent ✓
│   │   {[1,0],[2,0]}: dependent ✗ (second = 2 × first)
│   └── Basis: linearly independent set that spans the space
│
├── Quantum link:
│   {|0⟩, |1⟩} = computational basis of qubit space
│   They are linearly independent and span the 2D qubit space
│   Any qubit state = α|0⟩ + β|1⟩ (linear combination of basis)
│
└── Exit check:
    Are {[1,2,3], [4,5,6], [7,8,9]} linearly independent?
    Test: does row reduction give a zero row? Yes → dependent.

M2.1.3  Complex Vector Spaces
├── What to learn:
│   ├── Vectors where entries can be complex numbers
│   ├── ℂⁿ: n-dimensional complex vector space
│   ├── Qubit lives in ℂ²: vectors [α, β]ᵀ with α,β ∈ ℂ
│   └── All linear algebra extends to ℂⁿ
│
├── Dimension — why it matters:
│   1 qubit → ℂ²: 2D space
│   2 qubits → ℂ⁴: 4D space
│   n qubits → ℂ^(2ⁿ): grows EXPONENTIALLY
│   50 qubits → 2⁵⁰ ≈ 10¹⁵ dimensions (classical can't store this!)
│
└── Exit check:
    Represent |+⟩ = (1/√2)|0⟩+(1/√2)|1⟩ and |-⟩ = (1/√2)|0⟩-(1/√2)|1⟩.
    Verify they are orthogonal (inner product = 0) in NumPy.
```

---

## Module M2.2: Matrices & Operations

```
M2.2.1  Matrix as Linear Transformation
├── What to learn:
│   ├── Matrix M is a function: takes vector v → new vector Mv
│   ├── Rows × Columns notation: M is (m×n) if m rows, n cols
│   ├── Matrix-vector product (the MOST important operation):
│   │   [Mv]ᵢ = Σⱼ Mᵢⱼ vⱼ
│   └── Geometric meaning: rotation, scaling, reflection, shear
│
├── Quantum link:
│   Quantum gate = matrix. Applying gate to qubit = Mv multiplication.
│   X|0⟩: [[0,1],[1,0]] · [1,0]ᵀ = [0,1]ᵀ = |1⟩  ← NOT gate!
│
└── Code:
    X = np.array([[0,1],[1,0]])  # Pauli X gate
    ket0 = np.array([1,0])
    print(X @ ket0)              # [0, 1] = |1⟩

M2.2.2  Matrix Multiplication
├── What to learn:
│   ├── (A·B)ᵢⱼ = Σₖ Aᵢₖ Bₖⱼ
│   │   → row i of A dotted with column j of B
│   ├── NOT commutative: AB ≠ BA in general
│   ├── Associative: (AB)C = A(BC)
│   ├── Identity: AI = IA = A
│   └── Dimensions must match: (m×k)·(k×n) = (m×n)
│
├── Drill (no calculator):
│   [[1,2],[3,4]] × [[5,6],[7,8]] = ?
│   Row 1 · Col 1: 1·5 + 2·7 = 19
│   Row 1 · Col 2: 1·6 + 2·8 = 22
│   Row 2 · Col 1: 3·5 + 4·7 = 43
│   Row 2 · Col 2: 3·6 + 4·8 = 50
│   Result: [[19,22],[43,50]]
│
└── Exit check: Do the above by hand in 90 seconds. Then verify with @.

M2.2.3  Transpose, Conjugate, Dagger
├── What to learn:
│   ├── Transpose Aᵀ: rows become columns, Aᵀᵢⱼ = Aⱼᵢ
│   ├── Complex conjugate A*: replace every entry aᵢⱼ with aᵢⱼ*
│   ├── Dagger (conjugate transpose): A† = (Aᵀ)* = (A*)ᵀ
│   │   → MOST IMPORTANT operation in quantum mechanics
│   └── For real matrices: A† = Aᵀ
│
├── Code:
│   A = np.array([[1+1j, 2], [3, 4-2j]])
│   A_dag = A.conj().T    # conjugate transpose = dagger
│   print(A_dag)
│
└── Exit check: Compute A† for A=[[1+i, 2],[3, 4-2i]] by hand.

M2.2.4  Determinant, Inverse, Trace
├── What to learn:
│   ├── Determinant det(A):
│   │   2×2: det([[a,b],[c,d]]) = ad - bc
│   │   Geometric meaning: scaling factor of transformation
│   │   det=0 → matrix is singular (no inverse)
│   ├── Inverse A⁻¹: A·A⁻¹ = I
│   │   2×2: A⁻¹ = (1/det)·[[d,-b],[-c,a]]
│   └── Trace Tr(A): sum of diagonal elements Σᵢ Aᵢᵢ
│
└── Exit check:
    det([[2,1],[5,3]]) = 2·3-1·5 = 1.
    Verify A·A⁻¹=I for this matrix in NumPy.
```

---

## Module M2.3: Special Matrices — Hermitian & Unitary

```
M2.3.1  Hermitian Matrices
├── Definition: A = A†  (self-adjoint)
├── Properties:
│   ├── All eigenvalues are REAL (critical for quantum)
│   ├── Eigenvectors for different eigenvalues are orthogonal
│   └── Can always be diagonalized (spectral theorem)
│
├── Why quantum needs Hermitian:
│   ANY physically measurable quantity (energy, position, momentum)
│   is represented by a Hermitian operator.
│   WHY? Because measurements must give REAL numbers,
│   and Hermitian matrices have real eigenvalues.
│
├── Test: is A = [[3, 1+i],[1-i, 2]] Hermitian?
│   A† = [[3, 1+i],[1-i, 2]]ᵀ* = [[3, 1-i],[1+i, 2]]* = [[3, 1+i],[1-i, 2]]
│   A† = A ✓ → YES, Hermitian
│
└── Exit check: is_hermitian(A) → A† = A? Implement and test on 3 matrices.

M2.3.2  Unitary Matrices
├── Definition: U†U = UU† = I  (preserves inner product)
├── Properties:
│   ├── |det(U)| = 1
│   ├── Columns are orthonormal
│   ├── Maps unit vectors to unit vectors (preserves length/probability)
│   └── Rows are also orthonormal
│
├── Why quantum needs Unitary:
│   Quantum gates MUST be unitary:
│   1. Quantum evolution preserves total probability = 1
│   2. Unitary maps preserve norm: ||Uv|| = ||v||
│   3. Quantum must be reversible: U⁻¹ = U† (another gate)
│
└── Exit check: is_unitary(U) → U†U ≈ I? Implement. Test on Hadamard.

M2.3.3  Pauli Matrices — Memorize Cold
├── The four Paulis (foundation of all qubit gates):
│
│   I = [[1, 0],    X = [[0, 1],    Y = [[0, -i],    Z = [[1,  0],
│        [0, 1]]         [1, 0]]         [i,  0]]         [0, -1]]
│
├── Properties (ALL four):
│   ├── Hermitian: σᵢ = σᵢ†
│   ├── Unitary: σᵢ†σᵢ = I
│   ├── Self-inverse: σᵢ² = I
│   ├── Traceless (except I): Tr(X)=Tr(Y)=Tr(Z)=0
│   └── Algebra: XY=iZ, YZ=iX, ZX=iY, XYZ=iI
│
├── Physical meaning:
│   X: bit flip (NOT gate)   |0⟩↔|1⟩
│   Z: phase flip            |0⟩→|0⟩, |1⟩→-|1⟩
│   Y: both flip             combines X and Z
│   Z eigenvalues: +1 (for |0⟩), -1 (for |1⟩)
│
├── Code:
│   I = np.eye(2, dtype=complex)
│   X = np.array([[0,1],[1,0]], dtype=complex)
│   Y = np.array([[0,-1j],[1j,0]])
│   Z = np.array([[1,0],[0,-1]], dtype=complex)
│
└── Exit check: Write all 4 from memory in 60 sec. Verify XYZ=iI in NumPy.

M2.3.4  Hadamard & Rotation Gates
├── Hadamard H = (1/√2)[[1,1],[1,-1]]:
│   ├── H|0⟩ = |+⟩ = (1/√2)(|0⟩+|1⟩)  ← creates superposition
│   ├── H|1⟩ = |-⟩ = (1/√2)(|0⟩-|1⟩)
│   ├── H² = I  (self-inverse)
│   └── Mnemonic: Hadamard puts qubit "on the equator" of Bloch sphere
│
├── Rotation gates (key for VQE ansatz):
│   Rx(θ) = [[cos(θ/2), -i·sin(θ/2)],[-i·sin(θ/2), cos(θ/2)]]
│   Ry(θ) = [[cos(θ/2), -sin(θ/2)],  [sin(θ/2),   cos(θ/2)]]
│   Rz(θ) = [[e^(-iθ/2), 0],          [0, e^(iθ/2)]]
│
└── Exit check:
    Decompose [[3,1+i],[1-i,2]] = c₀I + c₁X + c₂Y + c₃Z.
    Find c₀,c₁,c₂,c₃. (Hint: cₖ = Tr(σₖM)/2)
    This IS how Qiskit's SparsePauliOp works.
```

---

## Module M3.1: Eigenvalues & Eigenvectors

```
M3.1.1  Definition and Physical Meaning
├── What to learn:
│   ├── Eigenvector equation: Av = λv
│   │   A: matrix (operator), v: eigenvector, λ: eigenvalue (scalar)
│   ├── Av just SCALES v by λ (doesn't rotate it!)
│   ├── Every matrix has special "preferred directions" = eigenvectors
│   └── Each preferred direction has a "stretch factor" = eigenvalue
│
├── Quantum postulate (CRITICAL):
│   When you MEASURE a quantum observable Â,
│   the ONLY possible outcomes are the eigenvalues of Â.
│   After measurement, the state collapses to the eigenvector.
│   VQE goal = find minimum eigenvalue of molecular Hamiltonian.
│
└── How to learn:
    3Blue1Brown EoLA Ep 13-14 [MANDATORY BEFORE ANYTHING ELSE]

M3.1.2  Computing Eigenvalues — Characteristic Equation
├── Method:
│   ├── det(A - λI) = 0   ← characteristic equation
│   ├── Solve for λ (get polynomial in λ)
│   ├── Roots of polynomial = eigenvalues
│   └── For each λ, solve (A-λI)v = 0 to find eigenvector
│
├── Worked example (2×2):
│   A = [[4,1],[2,3]]
│   det([[4-λ,1],[2,3-λ]]) = (4-λ)(3-λ) - 2·1 = 0
│   λ² - 7λ + 10 = 0
│   (λ-5)(λ-2) = 0 → λ = 5 or λ = 2
│   For λ=5: (A-5I)v=0 → [[-1,1],[2,-2]]v=0 → v=[1,1]ᵀ
│   For λ=2: (A-2I)v=0 → [[2,1],[2,1]]v=0  → v=[1,-2]ᵀ
│
└── Code:
    A = np.array([[4,1],[2,3]])
    eigenvalues, eigenvectors = np.linalg.eig(A)
    print(eigenvalues)     # [5. 2.]
    print(eigenvectors)    # columns are eigenvectors

M3.1.3  Hermitian Matrices — Special Properties
├── What to learn:
│   ├── Hermitian A → all eigenvalues are REAL
│   ├── Hermitian A → eigenvectors of different eigenvalues are ORTHOGONAL
│   ├── Spectral theorem: A = Σᵢ λᵢ |vᵢ⟩⟨vᵢ|
│   └── Diagonalization: A = PDP† where D=diag(λ₁,λ₂,...), P=eigenvectors
│
├── VQE connection:
│   Molecular Hamiltonian H is Hermitian
│   → All energy eigenvalues E₀, E₁, E₂,... are real numbers ✓
│   → Ground state energy E₀ = MINIMUM eigenvalue
│   → VQE's job: find E₀ = min ⟨ψ(θ)|H|ψ(θ)⟩
│
└── Exit check:
    H = 0.5·Z + 0.5·X (Pauli matrices).
    Compute eigenvalues by hand via characteristic equation.
    Verify with np.linalg.eig(). Which is the ground state energy?
```

---

## Module M3.2: Inner Products & Orthogonality

```
M3.2.1  Inner Product in ℝⁿ and ℂⁿ
├── Real inner product: ⟨u,v⟩ = uᵀv = Σᵢ uᵢvᵢ
├── Complex inner product: ⟨u,v⟩ = u†v = Σᵢ uᵢ*vᵢ
│   IMPORTANT: conjugate the FIRST vector, not second!
│   (Dirac notation: ⟨φ|ψ⟩ = conjugate of |φ⟩ dotted with |ψ⟩)
├── Norm (length): ||v|| = √⟨v,v⟩
├── Orthogonal: ⟨u,v⟩ = 0
└── Orthonormal: pairwise orthogonal + each normalized (||v||=1)

M3.2.2  Born Rule = Inner Product Squared
├── Quantum measurement probability:
│   P(measuring |φ⟩ when in state |ψ⟩) = |⟨φ|ψ⟩|²
│   This IS the Born rule — derived from inner product!
│
├── Example:
│   State |ψ⟩ = (√3/2)|0⟩ + (1/2)|1⟩
│   P(|0⟩) = |⟨0|ψ⟩|² = |√3/2|² = 3/4 = 75%
│   P(|1⟩) = |⟨1|ψ⟩|² = |1/2|² = 1/4 = 25%
│   Sum = 1 ✓ (normalization preserved)
│
└── Exit check: For |ψ⟩=(√3/2)|0⟩+(1/2)|1⟩, compute P(|0⟩) and P(|1⟩). Verify sum=1.

M3.2.3  Gram-Schmidt Orthogonalization
├── Problem: given any basis, make it orthonormal
├── Algorithm (for 2 vectors u,v):
│   e₁ = u / ||u||
│   v₂ = v - ⟨e₁,v⟩·e₁    ← remove component along e₁
│   e₂ = v₂ / ||v₂||
│
└── Why care: Measurement bases must be orthonormal.
    Gram-Schmidt constructs valid measurement bases.
```

---

## Module M3.3: Tensor Products

```
M3.3.1  Kronecker Product for Matrices
├── A⊗B: replace each entry aᵢⱼ with the block aᵢⱼ·B
│
├── Example (2×2 ⊗ 2×2 = 4×4):
│   [[a,b],[c,d]] ⊗ [[e,f],[g,h]]
│   = [[a·[[e,f],[g,h]], b·[[e,f],[g,h]]],
│      [c·[[e,f],[g,h]], d·[[e,f],[g,h]]]]
│   = [[ae,af,be,bf],[ag,ah,bg,bh],[ce,cf,de,df],[cg,ch,dg,dh]]
│
└── Code: np.kron(A, B)

M3.3.2  Tensor Product of Quantum States
├── |ψ⟩ ⊗ |φ⟩:
│   |0⟩⊗|0⟩ = |00⟩ = [1,0,0,0]ᵀ
│   |0⟩⊗|1⟩ = |01⟩ = [0,1,0,0]ᵀ
│   |1⟩⊗|0⟩ = |10⟩ = [0,0,1,0]ᵀ
│   |1⟩⊗|1⟩ = |11⟩ = [0,0,0,1]ᵀ
│
├── Gate on multi-qubit: apply gate to qubit k with I on others
│   H on qubit 0 of 2-qubit system: H⊗I (4×4 matrix)
│
├── Why exponential:
│   1 qubit → ℂ² (2 numbers)
│   2 qubits → ℂ² ⊗ ℂ² = ℂ⁴ (4 numbers)
│   n qubits → ℂ^(2ⁿ) (2ⁿ numbers)
│   50 electrons in molecule → 2⁵⁰ ≈ 10¹⁵ numbers
│   Classical computers CANNOT store this → quantum advantage!
│
└── Exit check:
    Compute |+⟩⊗|0⟩ manually. Result: (1/√2)[1,0,1,0]ᵀ.
    Then compute (H⊗I)|00⟩ using np.kron and @ operator. Verify.
```

---

## Module M4.1-M4.3: Hilbert Spaces, Probability, ODEs

```
M4.1  Hilbert Space — The Container of Quantum States
├── Abstract definition: complete inner product space
├── Practical (finite QC): ℂⁿ with standard inner product
├── Orthonormal basis: {|e₁⟩,...,|eₙ⟩} where ⟨eᵢ|eⱼ⟩ = δᵢⱼ
├── Basis expansion: any |ψ⟩ = Σᵢ cᵢ|eᵢ⟩ where cᵢ = ⟨eᵢ|ψ⟩
└── Completeness: Σᵢ|eᵢ⟩⟨eᵢ| = I  (resolution of identity)

M4.1.1  Completeness Relation in Practice
├── |0⟩⟨0| + |1⟩⟨1| = I₂    (qubit completeness)
├── Inserting I: ⟨φ|ψ⟩ = ⟨φ|I|ψ⟩ = Σₙ ⟨φ|n⟩⟨n|ψ⟩
└── Exit check: Prove |+⟩⟨+| + |-⟩⟨-| = I₂ in NumPy.

M4.2  Probability & Statistics (Quantum Context)
├── 4.2.1 Expectation values:
│   Classical: E[X] = Σ xᵢ P(xᵢ)
│   Quantum:   ⟨A⟩ = ⟨ψ|A|ψ⟩  ← this is VQE's cost function
│
├── 4.2.2 Shot noise:
│   Running circuit N times → estimate ⟨A⟩ from average of outcomes
│   Statistical error ≈ 1/√N  (need large N for precision)
│   VQE typically uses N = 1024 to 8192 shots per evaluation
│
└── Exit check: Compute ⟨Z⟩ for |+⟩ analytically. Verify with 10000-shot Qiskit sim.

M4.3  ODE Review — Schrödinger Equation Form
├── TISE: Ĥψ = Eψ (eigenvalue equation — matrix form!)
├── TDSE: iℏ∂|ψ⟩/∂t = H|ψ⟩ (how state evolves in time)
├── Solution: |ψ(t)⟩ = e^(-iHt/ℏ)|ψ(0)⟩
└── Particle in box: ψₙ(x) = √(2/L)·sin(nπx/L), Eₙ = n²π²ℏ²/2mL²
```
