# ⚛️ Quantum Bioinformatics — Deep Chapter-wise Syllabus
> **Level:** ########## → VQE + Quantum ML for Genetic Engineering
> **Format:** Every Node → Sub-chapters → What/How/Code/Exit Criteria

---

# PHASE 0: BRIDGE SETUP (Weeks 1-2)

---

## Module P0.1: Complex Numbers — 10th to Quantum Bridge

> **PREREQUISITES from Class 9/10 you MUST know:**
> Number system: Natural, Integer, Rational, Irrational, Real numbers.
> Basic algebra: solving x² = 4 → x = ±2.
> Coordinate geometry: x-y plane, plotting points (x, y).
> Pythagoras theorem: a² + b² = c².
> If any of these are shaky → revise NCERT Class 9 Ch1 and 10 Ch1 first.

```
P0.1.1  The Problem with Real Numbers — Why We Need Complex Numbers
├── What you already know (Class 9/10):
│   ├── Natural: 1, 2, 3, ...
│   ├── Integer: ..., -2, -1, 0, 1, 2, ...
│   ├── Rational: p/q form (1/2, 3/4, -7/3)
│   ├── Irrational: √2, π, e (cannot be written as p/q)
│   └── Real: ALL of the above combined → the number LINE
│
├── The problem: what is √(-1)?
│   √4 = 2  because 2² = 4  ✓
│   √9 = 3  because 3² = 9  ✓
│   √(-1) = ?   because (??)² = -1  → NO real number works!
│   Squaring any real number gives a POSITIVE result.
│   So: √(-1) has NO solution in real numbers.
│
├── The solution: INVENT a new number
│   Define: i = √(-1), which means i² = -1
│   i is called the IMAGINARY unit.
│   This isn't "fake" — it's just extending our number system,
│   exactly like how we "invented" negative numbers to solve x+5=3.
│
├── Complex number z = a + bi:
│   a = real part (Re(z))     ← normal real number
│   b = imaginary part (Im(z)) ← coefficient of i
│   Examples:
│   ├── z = 3 + 4i  → a=3, b=4
│   ├── z = 2 - i   → a=2, b=-1
│   ├── z = 5       → a=5, b=0  (pure real)
│   └── z = 3i      → a=0, b=3  (pure imaginary)
│
├── Argand diagram (Complex PLANE):
│   Instead of 1D number LINE → we have 2D PLANE
│   x-axis = real part, y-axis = imaginary part
│   Plot z = 3 + 4i → point at (3, 4) on this plane
│   Plot z = -1 + 2i → point at (-1, 2)
│   Plot z = 4 → point at (4, 0) (on x-axis = real axis)
│   Plot z = 3i → point at (0, 3) (on y-axis = imaginary axis)
│
├── WHY this matters for quantum:
│   Quantum amplitudes (α, β in |ψ⟩=α|0⟩+β|1⟩) are COMPLEX numbers
│   They live in 2D complex plane, not 1D real line
│   The "phase" of a complex number = angle on this plane
│   Phase is what creates quantum interference effects!
│
├── Code:
│   import numpy as np
│   z1 = 3 + 4j              # Python: j not i
│   z2 = complex(-1, 2)      # alternative: complex(real, imag)
│   print(z1.real, z1.imag)  # 3.0 and 4.0
│   # Plot:
│   import matplotlib.pyplot as plt
│   numbers = [3+4j, -1+2j, 2-3j, -2-1j]
│   plt.axhline(0, color='k'); plt.axvline(0, color='k')
│   for z in numbers:
│       plt.plot(z.real, z.imag, 'ro', ms=8)
│       plt.annotate(f'{z}', (z.real, z.imag))
│   plt.xlabel('Real'); plt.ylabel('Imaginary')
│   plt.title('Argand Diagram'); plt.grid(True); plt.show()
│
└── Self-check: Plot these by hand → (2+3i), (-4+i), (3-2i), -5i.
    What axis does a pure real number lie on? (Real/x-axis)
    What axis does a pure imaginary number lie on? (Imaginary/y-axis)

P0.1.2  Complex Arithmetic — Step by Step
├── ADDITION & SUBTRACTION (combine like terms):
│   (a+bi) + (c+di) = (a+c) + (b+d)i
│   Treat real parts and imaginary parts SEPARATELY
│
│   WORKED EXAMPLE:
│   (3+4i) + (2-i) = (3+2) + (4+(-1))i = 5 + 3i ✓
│   (3+4i) - (2-i) = (3-2) + (4-(-1))i = 1 + 5i ✓
│
├── MULTIPLICATION (use FOIL then i²=-1):
│   (a+bi)(c+di) = ac + adi + bci + bdi²
│                = ac + adi + bci + bd(-1)     [since i²=-1]
│                = (ac - bd) + (ad + bc)i
│
│   WORKED EXAMPLE step by step:
│   (3+4i)(2-i)
│   = 3·2 + 3·(-i) + 4i·2 + 4i·(-i)
│   = 6 - 3i + 8i - 4i²
│   = 6 + 5i - 4(-1)              [since i²=-1]
│   = 6 + 5i + 4
│   = 10 + 5i  ✓
│
│   WORKED EXAMPLE: multiply two imaginary numbers:
│   (2i)(3i) = 6i² = 6(-1) = -6  (real number!)
│
├── DIVISION (multiply by conjugate):
│   Conjugate of (c+di) = (c-di) [flip sign of imaginary part]
│   Trick: (c+di)(c-di) = c² + d² (always real!)
│
│   WORKED EXAMPLE: (3+4i) ÷ (2-i)
│   = (3+4i)/(2-i) × (2+i)/(2+i)     [multiply by conjugate/conjugate]
│   = (3+4i)(2+i) / [(2-i)(2+i)]
│   = (6+3i+8i+4i²) / (4+1)
│   = (6+11i-4) / 5
│   = (2+11i) / 5
│   = 0.4 + 2.2i  ✓
│
├── Code:
│   z1 = 3 + 4j
│   z2 = 2 - 1j
│   print(f"Add: {z1+z2}")   # (5+3j)
│   print(f"Sub: {z1-z2}")   # (1+5j)
│   print(f"Mul: {z1*z2}")   # (10+5j)
│   print(f"Div: {z1/z2}")   # (0.4+2.2j)
│
└── Exit check (DO BY HAND then verify in Python):
    1. (1+2i)(3-4i) = ?  Answer: (11+2i)
    2. (2+3i)² = ?       Answer: (-5+12i)
    3. (1+i)/(1-i) = ?   Answer: i  (pure imaginary!)

P0.1.3  Conjugate, Modulus, and Normalization
├── CONJUGATE z* = a - bi:
│   Just flip the sign of imaginary part.
│   z = 3 + 4i → z* = 3 - 4i
│   z = 2 - 5i → z* = 2 + 5i
│   z = 7 (real) → z* = 7 (real stays same)
│   z = 4i → z* = -4i
│
│   Key property: z·z* = (a+bi)(a-bi) = a² + b² (ALWAYS real, positive!)
│   Proof: (a+bi)(a-bi) = a² - abi + abi - b²i² = a² - b²(-1) = a² + b²  ✓
│
├── MODULUS |z| = √(a² + b²):
│   This is the DISTANCE from origin to z on the complex plane
│   Same as Pythagoras theorem! (right triangle with legs a and b)
│   |z|² = a² + b² = z·z* (quick formula)
│
│   WORKED EXAMPLES:
│   |3+4i| = √(3²+4²) = √(9+16) = √25 = 5
│   |1+i|  = √(1²+1²) = √2 ≈ 1.414
│   |5|    = √(25+0)  = 5  (real numbers: modulus = absolute value)
│   |3i|   = √(0+9)   = 3  (pure imaginary: modulus = coefficient)
│
│   Code verification:
│   z = 3 + 4j
│   print(abs(z))           # 5.0  (Python's abs() works on complex)
│   print(np.abs(z))        # 5.0
│   print(np.sqrt(z.real**2 + z.imag**2))  # 5.0 (manual)
│
├── NORMALIZATION (z on the unit circle):
│   If |z| = 1 → z lies on unit circle (circle of radius 1)
│   z = e^(iθ) is always on the unit circle (|e^(iθ)| = 1)
│   To normalize z: divide by |z| → gives z/|z| with unit modulus
│
├── WHY THIS MATTERS — Quantum connection:
│   Qubit state: |ψ⟩ = α|0⟩ + β|1⟩
│   α = amplitude for |0⟩, β = amplitude for |1⟩
│   Born rule: P(measuring |0⟩) = |α|², P(measuring |1⟩) = |β|²
│   Total probability = 1: |α|² + |β|² = 1 (normalization condition!)
│
│   WORKED EXAMPLE:
│   |ψ⟩ = (1/√2)|0⟩ + (i/√2)|1⟩
│   α = 1/√2, β = i/√2
│   |α|² = |1/√2|² = (1/√2)·(1/√2) = 1/2
│   |β|² = |i/√2|² = |i|²·|1/√2|² = 1·(1/2) = 1/2
│           ↑ |i|=1 because i has modulus 1!
│   |α|²+|β|² = 1/2 + 1/2 = 1 ✓ (properly normalized)
│
│   P(|0⟩) = 1/2 = 50%, P(|1⟩) = 1/2 = 50%
│   Even though β = i/√2 has a complex value, probability is still real!
│
└── Code:
    alpha = 1/np.sqrt(2)
    beta  = 1j/np.sqrt(2)
    print(f"|α|² = {abs(alpha)**2:.3f}")   # 0.500
    print(f"|β|² = {abs(beta)**2:.3f}")    # 0.500
    print(f"Sum  = {abs(alpha)**2 + abs(beta)**2:.3f}")  # 1.000 ✓

═══════════════════════════════════════════
 GATE TO M1.1 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Can define i and explain i² = -1
 □ Can plot any complex number z = a+bi on Argand diagram
 □ Can add, subtract, multiply complex numbers by hand (no calculator)
 □ Can compute conjugate z* and modulus |z| = √(a²+b²) for any z
 □ Know: |z|² = z·z* = a²+b² (ALWAYS real)
 □ Verified: |1+2i| = √5 ≈ 2.236 by hand AND Python
 □ Verified: |i/√2|² = 1/2 (complex amplitude normalization)
 □ Can explain: "why are quantum amplitudes complex and not just real?"
   (Phase: complex numbers have an angle θ that determines interference.
    Real-only amplitudes cannot create destructive interference.)
═══════════════════════════════════════════
```

---

# PHASE 1: MATHEMATICS (Weeks 3-10)

---

## Module M1.1: Complex Numbers — Full Mastery

> **PREREQUISITES: P0.1 gate passed.**
> Also need from Class 10: sin, cos, tan definitions (SOHCAHTOA),
> and basic graph shapes of sin/cos. If forgotten → revise NCERT Class 10 Ch8.
> Also: Taylor series idea (not full proof needed — just the idea that
> a function can be written as an infinite polynomial).

```
M1.1.1  Polar Form — Two Ways to Describe the Same Number
├── The two representations:
│   Rectangular: z = a + bi   (good for: addition, subtraction)
│   Polar:       z = r·e^(iθ) (good for: multiplication, rotation)
│   They describe the SAME complex number, just different clothing.
│
├── Converting Rectangular → Polar:
│   Given: z = a + bi
│   r = |z| = √(a²+b²)    ← modulus (from P0.1.3)
│   θ = angle of z from positive real axis
│     = arctan(b/a)  [careful about quadrant! use atan2 in code]
│
│   WORKED EXAMPLE:
│   z = 1 + i    → a=1, b=1
│   r = √(1²+1²) = √2 ≈ 1.414
│   θ = arctan(1/1) = arctan(1) = π/4 = 45°
│   Polar: z = √2 · e^(iπ/4)
│
│   z = -1 + i   → a=-1, b=1
│   r = √(1+1) = √2
│   θ = arctan(1/-1) = arctan(-1) but in 2nd quadrant → θ = 3π/4 = 135°
│   Polar: z = √2 · e^(i·3π/4)
│
│   z = -1        → a=-1, b=0
│   r = 1, θ = π (pointing left on real axis)
│   Polar: z = 1·e^(iπ) = e^(iπ)
│   NOTE: this is Euler's famous result: e^(iπ) = -1 !
│
├── Converting Polar → Rectangular:
│   Given: z = r·e^(iθ) = r(cosθ + i·sinθ)   ← Euler's formula
│   a = r·cosθ   (real part)
│   b = r·sinθ   (imaginary part)
│
│   WORKED EXAMPLE:
│   z = 2·e^(iπ/3) = 2·(cos60° + i·sin60°)
│                  = 2·(0.5 + i·(√3/2))
│                  = 1 + i√3
│   Verify: r = √(1² + (√3)²) = √(1+3) = 2 ✓, θ = arctan(√3/1) = π/3 ✓
│
├── Code:
│   import numpy as np
│   z = -1 + 1j
│   r = np.abs(z)             # modulus = √2
│   theta = np.angle(z)       # angle in radians = 3π/4
│   print(f"r = {r:.4f}, θ = {theta:.4f} rad = {np.degrees(theta):.1f}°")
│   z_back = r * np.exp(1j * theta)   # reconstruct from polar
│   print(f"Reconstructed: {z_back}")  # (-1+1j)
│   print(np.isclose(z, z_back))       # True ✓
│
└── Exit check (DO BY HAND):
    Convert these to polar (find r and θ):
    1. z = 1 + 0i    → r=1, θ=0
    2. z = 0 + 1i    → r=1, θ=π/2
    3. z = -1 + 0i   → r=1, θ=π
    4. z = 0 - 1i    → r=1, θ=-π/2 (or 3π/2)
    5. z = 3 + 4i    → r=5, θ=arctan(4/3)≈53.13°
    Verify all 5 with numpy.

M1.1.2  Euler's Formula — The Most Important Equation in All of Quantum
├── The formula: e^(iθ) = cosθ + i·sinθ
│   (Euler pronounced "Oiler", Swiss mathematician 1707-1783)
│   This single equation connects:
│   e (base of natural log), i (imaginary unit), cos, sin
│   It's so beautiful that when θ=π: e^(iπ) + 1 = 0
│   (combines 5 most fundamental constants: e, i, π, 1, 0)
│
├── WHERE DOES IT COME FROM? (Taylor series proof):
│   Any well-behaved function can be written as infinite polynomial:
│   e^x = 1 + x + x²/2! + x³/3! + x⁴/4! + x⁵/5! + ...
│   cos(x) = 1 - x²/2! + x⁴/4! - x⁶/6! + ...     (even powers)
│   sin(x) = x - x³/3! + x⁵/5! - x⁷/7! + ...      (odd powers)
│
│   Now plug in x = iθ into e^x:
│   e^(iθ) = 1 + iθ + (iθ)²/2! + (iθ)³/3! + (iθ)⁴/4! + ...
│           = 1 + iθ + i²θ²/2! + i³θ³/3! + i⁴θ⁴/4! + ...
│   Use i²=-1, i³=-i, i⁴=+1, i⁵=i, ... (pattern repeats every 4):
│           = 1 + iθ - θ²/2! - iθ³/3! + θ⁴/4! + iθ⁵/5! - ...
│   Group real and imaginary parts:
│           = (1 - θ²/2! + θ⁴/4! - ...) + i(θ - θ³/3! + θ⁵/5! - ...)
│           =    cos(θ)                  + i·sin(θ)   ✓
│
├── KEY VALUES (must memorize all 8):
│   ┌────────────┬──────────────────────────┬────────────────┐
│   │   θ        │  e^(iθ) = cosθ + i·sinθ │  Location      │
│   ├────────────┼──────────────────────────┼────────────────┤
│   │   0        │  1 + 0i = 1              │  Right (+x)    │
│   │   π/6 (30°)│  (√3/2) + (1/2)i        │  1st quadrant  │
│   │   π/4 (45°)│  (1/√2) + (1/√2)i       │  1st quadrant  │
│   │   π/3 (60°)│  (1/2) + (√3/2)i        │  1st quadrant  │
│   │   π/2 (90°)│  0 + i                   │  Top (+y)      │
│   │   π  (180°)│  -1 + 0i = -1            │  Left (-x)     │
│   │   3π/2(270°)│ 0 - i                   │  Bottom (-y)   │
│   │   2π (360°)│  1 + 0i = 1              │  Back to start │
│   └────────────┴──────────────────────────┴────────────────┘
│
├── Why DOES this appear in quantum physics?
│   Schrödinger: ψ(t) = ψ(0) · e^(-iEt/ℏ)
│   This is PHASE ROTATION: |e^(-iEt/ℏ)|=1 always (probability conserved!)
│   The state "rotates" in the complex plane as time passes.
│   Different energy levels rotate at different speeds (different E).
│   Two levels rotating at different speeds → their PHASE DIFFERENCE changes
│   → INTERFERENCE pattern → this is what VQE's cost function captures!
│
├── Code (verify all 8 values):
│   import numpy as np
│   angles = [0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]
│   names  = ['0', 'π/6', 'π/4', 'π/3', 'π/2', 'π', '3π/2', '2π']
│   for a, n in zip(angles, names):
│       z = np.exp(1j * a)
│       print(f"e^(i·{n}) = {z.real:.4f} + {z.imag:.4f}i, |z|={abs(z):.4f}")
│   # Notice: all have |z| = 1 exactly! ← they're on the unit circle
│
└── Exit check:
    1. Compute e^(iπ/4) in rectangular form:
       cos(π/4) + i·sin(π/4) = (1/√2) + i/√2 = 0.7071 + 0.7071i
    2. Verify |e^(iπ/4)| = 1.  (0.7071²+0.7071² = 0.5+0.5 = 1 ✓)
    3. What is e^(iπ/4) × e^(iπ/4)?
       = e^(i·2π/4) = e^(iπ/2) = i  (rotation of 90° twice = 90°+90°=180°... wait no, 45°+45°=90°)

M1.1.3  Multiplication = Pure Rotation
├── In polar form: z₁·z₂ = (r₁e^(iθ₁))·(r₂e^(iθ₂)) = r₁r₂·e^(i(θ₁+θ₂))
│   Magnitudes MULTIPLY, angles ADD.
│
├── The beautiful implication:
│   Multiplying any z by e^(iφ) (a unit-modulus complex number):
│   |z| stays the same (r×1 = r unchanged)
│   θ increases by φ (rotation by angle φ)
│   → "Multiplying by e^(iφ)" = PURE ROTATION by angle φ!
│
├── WORKED EXAMPLES:
│   Start with z = 1 (at angle 0), multiply by e^(iπ/2) = i:
│   i × 1 = i  (now at angle π/2 = 90°)  → rotated 90° ✓
│
│   Start with z = 1+i (at 45°), rotate by π/4:
│   (1+i) × e^(iπ/4) = (√2 · e^(iπ/4)) × e^(iπ/4)
│                     = √2 · e^(iπ/2)
│                     = √2 · i  (now at 90°, radius still √2)  ✓
│
├── De Moivre's theorem:
│   (e^(iθ))^n = e^(inθ)   ← rotating n times = multiplying angle by n
│   (cosθ + i·sinθ)^n = cos(nθ) + i·sin(nθ)
│
│   Application: What is i⁴?
│   i = e^(iπ/2), so i⁴ = e^(i·4π/2) = e^(i·2π) = 1 ✓
│   i¹=i, i²=-1, i³=-i, i⁴=1  (cycle of 4!)
│
├── Quantum gate connection:
│   Rotation gates Rx(θ), Ry(θ), Rz(θ) in quantum computing
│   literally ROTATE the qubit state on the Bloch sphere by angle θ
│   They use e^(iθ/2) internally:
│   Rz(θ) = [[e^(-iθ/2), 0], [0, e^(iθ/2)]]  ← rotation in complex plane!
│   VQE ansatz = sequence of rotation gates with adjustable θ values
│
├── Roots of unity (bonus — beautiful math):
│   z^n = 1 has EXACTLY n solutions → nth roots of unity
│   e^(2πik/n) for k=0,1,...,n-1 (equally spaced on unit circle!)
│   Square roots of 1: e^(0)=1, e^(iπ)=-1 → {1, -1} ✓
│   4th roots of 1: k=0,1,2,3 → {1, i, -1, -i} ✓ (north/south/east/west)
│   Used in: quantum Fourier transform, Grover algorithm phase kickback
│
├── Code:
│   def rotate_complex(z, phi):
│       return z * np.exp(1j * phi)
│   z = 1 + 0j  # start at (1+0i)
│   print(f"Start: {z}")
│   for k in range(1, 5):
│       z = rotate_complex(z, np.pi/2)
│       print(f"After {k}×90°: {np.round(z, 4)}")
│   # Should give: i, -1, -i, 1 (returns to start!)
│
└── Exit check:
    1. z = -1+i, multiply by e^(-iπ/4). What angle does result have?
       z has θ=3π/4, rotated by -π/4 → new θ=3π/4-π/4 = π/2 → result is i√2
    2. Compute (1+i)^8 using De Moivre. r=√2, θ=π/4.
       (√2·e^(iπ/4))^8 = (√2)^8·e^(i·2π) = 16·1 = 16
    3. Verify in Python: print((1+1j)**8) → (16+0j) ✓

═══════════════════════════════════════════
 GATE TO M1.2 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Can convert any z=a+bi to polar form r·e^(iθ) by hand
 □ Can convert polar r·e^(iθ) back to a+bi using cosθ, sinθ
 □ Know all 8 key values of e^(iθ) (table in M1.1.2)
 □ Understand Euler's formula e^(iθ)=cosθ+i·sinθ (can sketch Taylor proof)
 □ Know: multiplying by e^(iφ) = rotating by angle φ
 □ Verified: 4 rotations of 90° return to start
 □ Can compute (1+i)^8 = 16 using De Moivre
 □ Know: quantum rotation gates Rz(θ) use e^(±iθ/2) → polar form is essential
═══════════════════════════════════════════
```

---

## Module M1.2: Calculus — Derivatives & Partial Derivatives

> **PREREQUISITES: M1.1 gate passed.**
> Also need from Class 11 PCB Math (if you had it) or Class 11 Physics:
> - Basic idea of rate of change (speed = distance/time)
> - Graph reading: slope of a line = rise/run
> - Functions: what f(x) means (machine that takes x, gives output)
> If completely new to calculus: watch 3Blue1Brown "Essence of Calculus" Ep 1-3 FIRST.

```
M1.2.1  What is a Derivative — The Key Idea
├── Start with something you KNOW: speed
│   If you travel 100 km in 2 hours → average speed = 50 km/h
│   Average rate of change = (change in distance) / (change in time) = Δd/Δt
│
├── Instantaneous speed (derivative):
│   Suppose d(t) = t² (distance at time t, in meters)
│   Average speed from t=2 to t=2+h:
│   = [d(2+h) - d(2)] / h
│   = [(2+h)² - 4] / h
│   = [4 + 4h + h² - 4] / h
│   = [4h + h²] / h
│   = 4 + h
│   As h→0 (smaller and smaller time interval): speed → 4 m/s
│   This is the DERIVATIVE at t=2: d'(t)=2t → d'(2)=4 ✓
│
├── Formal definition:
│   f'(x) = lim[h→0] [f(x+h) - f(x)] / h
│   = "instantaneous rate of change of f at x"
│   = "slope of tangent line to graph of f at point x"
│
├── Notation (all mean the same thing):
│   f'(x) = df/dx = d/dx[f(x)] = Ḟ (physics dot notation)
│
├── Visual meaning:
│   f'(x) > 0: function going UP (slope positive)
│   f'(x) < 0: function going DOWN (slope negative)
│   f'(x) = 0: FLAT (peak, valley, or inflection point)
│   f''(x) > 0 at f'=0: MINIMUM (valley) ← VQE optimization target!
│
├── Biology analogy (since you have PCB):
│   Population P(t) of bacteria at time t
│   P'(t) = growth rate (cells per hour)
│   P'(t) > 0: population growing
│   P'(t) = 0: population at steady state
│   P'(t) < 0: population dying
│   Drug effect = changing P'(t) → calculus tells you drug efficacy!
│
└── Code (visualize derivative as slope):
    import numpy as np, matplotlib.pyplot as plt
    x = np.linspace(0, 3, 100)
    f = x**2
    f_prime = 2*x     # analytical derivative
    plt.plot(x, f, 'b-', label='f(x)=x²')
    plt.plot(x, f_prime, 'r-', label="f'(x)=2x")
    plt.axhline(0, color='k', lw=0.5)
    plt.legend(); plt.grid(True); plt.title('Function and its derivative')
    plt.show()

M1.2.2  Differentiation Rules — Complete Table with Worked Examples
├── RULE 1: Power Rule  d/dx[xⁿ] = n·xⁿ⁻¹
│   The most used rule. "Bring exponent down, reduce exponent by 1"
│   d/dx[x²]  = 2x¹ = 2x
│   d/dx[x³]  = 3x²
│   d/dx[x]   = 1·x⁰ = 1
│   d/dx[1]   = 0    (constant → derivative = 0)
│   d/dx[x⁻¹] = -x⁻² = -1/x² (negative exponents work too!)
│   d/dx[√x]  = d/dx[x^(1/2)] = (1/2)x^(-1/2) = 1/(2√x)
│
├── RULE 2: Sum Rule  d/dx[f+g] = f' + g'
│   Derivatives distribute over addition.
│   d/dx[x³ + 2x² - 5x + 7]
│   = 3x² + 2(2x) - 5(1) + 0
│   = 3x² + 4x - 5  ✓
│
├── RULE 3: Product Rule  d/dx[f·g] = f'g + fg'
│   d/dx[x² · e^x] = 2x·e^x + x²·e^x = e^x(2x + x²) = x(x+2)e^x
│
├── RULE 4: Chain Rule  d/dx[f(g(x))] = f'(g(x)) · g'(x)
│   "Derivative of outer (leaving inner unchanged) × derivative of inner"
│   d/dx[e^(3x)] = e^(3x) · 3 = 3e^(3x)   [outer=e^u, inner=3x]
│   d/dx[(x²+1)⁵] = 5(x²+1)⁴ · 2x = 10x(x²+1)⁴
│   d/dx[sin(2x)] = cos(2x) · 2 = 2cos(2x)
│
├── RULE 5: Exponential and Trig
│   d/dx[e^x]    = e^x  (unique: e^x is its own derivative!)
│   d/dx[e^(ax)] = a·e^(ax)  (chain rule)
│   d/dx[sin(x)] = cos(x)
│   d/dx[cos(x)] = -sin(x)  (note: negative!)
│   d/dx[ln(x)]  = 1/x
│
├── WORKED EXAMPLE — differentiate completely:
│   f(x) = x³·sin(x) + e^(2x)
│   Term 1: d/dx[x³sin(x)] = 3x²·sin(x) + x³·cos(x)   [product rule]
│   Term 2: d/dx[e^(2x)] = 2e^(2x)                      [chain rule]
│   f'(x) = 3x²sin(x) + x³cos(x) + 2e^(2x)
│
├── Quantum relevance:
│   E(θ) VQE cost function has terms like e^(iθ) and trigonometric
│   "Parameter shift rule" (quantum gradient): ∂E/∂θ = [E(θ+π/2) - E(θ-π/2)]/2
│   This is a quantum analog of the chain rule!
│
├── Code (verify numerically):
│   def numerical_deriv(f, x, h=1e-7):
│       return (f(x+h) - f(x)) / h
│   # Test: f(x) = x³, analytical: f'(2) = 3×4 = 12
│   f = lambda x: x**3
│   print(f"Numerical: {numerical_deriv(f, 2):.6f}")  # 12.000000
│   print(f"Analytical: {3*4}")                        # 12
│
└── Exit check (DO BY HAND):
    Differentiate each:
    1. f(x) = 4x⁵ - 3x² + 7         → Answer: 20x⁴ - 6x
    2. g(x) = e^(5x)                  → Answer: 5e^(5x)
    3. h(x) = sin(x²)                 → Answer: 2x·cos(x²)  [chain rule]
    4. k(x) = x²·e^x                  → Answer: e^x(x²+2x)  [product+chain]
    Verify all 4 numerically with numerical_deriv.

M1.2.3  Partial Derivatives — Functions of Multiple Variables
├── Why needed:
│   VQE energy function: E(θ₁, θ₂, ..., θₙ) — depends on n angles
│   To minimize E: need gradient → need ∂E/∂θₖ for each k
│   This is EXACTLY gradient descent (from your ML track!)
│
├── What is a partial derivative:
│   f(x,y) = function of TWO variables (graph = surface in 3D)
│   ∂f/∂x: "how does f change as x changes, keeping y FIXED"
│   Same process as normal derivative, just treat other variables as constants
│
├── WORKED EXAMPLES — full step by step:
│   f(x,y) = x²y + 3y² + 5x
│
│   ∂f/∂x = treat y as a constant:
│   ∂/∂x[x²y] = y·∂/∂x[x²] = y·2x = 2xy   [y is constant, so outside]
│   ∂/∂x[3y²] = 0                             [no x terms → zero]
│   ∂/∂x[5x]  = 5                             [power rule]
│   → ∂f/∂x = 2xy + 5
│
│   ∂f/∂y = treat x as a constant:
│   ∂/∂y[x²y] = x²·∂/∂y[y] = x²·1 = x²       [x² is constant]
│   ∂/∂y[3y²] = 6y                              [power rule]
│   ∂/∂y[5x]  = 0                               [no y term]
│   → ∂f/∂y = x² + 6y
│
│   Verify at point (x=2, y=3):
│   ∂f/∂x|(2,3) = 2(2)(3) + 5 = 12 + 5 = 17
│   ∂f/∂y|(2,3) = (2)² + 6(3) = 4 + 18 = 22
│
├── Gradient vector:
│   ∇f = [∂f/∂x, ∂f/∂y] = [17, 22] at point (2,3)
│   Points in direction of STEEPEST ASCENT of f
│   Gradient descent: move in -∇f direction (go downhill)
│
├── VQE gradient connection:
│   E(θ) measured by quantum hardware (circuit + measurement)
│   ∂E/∂θₖ computed by "parameter shift rule": run circuit twice
│   θₖ_new = θₖ - α · ∂E/∂θₖ   ← gradient descent update!
│   Repeat until E(θ) converges to minimum (= ground state energy)
│
├── Code:
│   def partial_deriv(f, var_idx, point, h=1e-7):
│       point_plus = list(point)
│       point_plus[var_idx] += h
│       return (f(*point_plus) - f(*point)) / h
│   f = lambda x, y: x**2 * y + 3*y**2 + 5*x
│   print(partial_deriv(f, 0, [2, 3]))  # ∂f/∂x at (2,3) ≈ 17
│   print(partial_deriv(f, 1, [2, 3]))  # ∂f/∂y at (2,3) ≈ 22
│
└── Exit check:
    f(x,y,z) = x²y + yz² + 3xyz
    Compute: ∂f/∂x, ∂f/∂y, ∂f/∂z at point (1,2,-1)
    Answers:
    ∂f/∂x = 2xy + 3yz = 2(1)(2) + 3(2)(-1) = 4 - 6 = -2
    ∂f/∂y = x² + z² + 3xz = 1 + 1 + 3(1)(-1) = 2-3 = -1
    ∂f/∂z = 2yz + 3xy = 2(2)(-1) + 3(1)(2) = -4+6 = 2

M1.2.4  Integration — Area Under the Curve
├── Integration = REVERSE of differentiation (anti-derivative)
│   If F'(x) = f(x), then ∫f(x)dx = F(x) + C
│
├── Visual meaning:
│   Definite integral ∫[a to b]f(x)dx = area under curve from a to b
│   Think: sum of infinitely thin rectangles (height=f(x), width=dx)
│
├── KEY INTEGRALS (memorize for quantum work):
│   ∫xⁿdx = xⁿ⁺¹/(n+1) + C  (power rule reversed)
│   ∫e^x dx = e^x + C
│   ∫e^(ax)dx = e^(ax)/a + C
│   ∫sin(x)dx = -cos(x) + C
│   ∫cos(x)dx = sin(x) + C
│   ∫1/x dx = ln|x| + C
│
├── WORKED EXAMPLE 1 (polynomial):
│   ∫(3x²+2)dx from 0 to 1
│   Step 1: find antiderivative: F(x) = x³ + 2x
│   Step 2: F(1) - F(0) = (1+2) - (0+0) = 3  ✓
│
├── WORKED EXAMPLE 2 (quantum normalization):
│   Verify ψ₁(x) = √2·sin(πx) is normalized on [0,1]:
│   ∫₀¹ |ψ₁|² dx = ∫₀¹ 2·sin²(πx)dx
│   Use identity: sin²(u) = (1-cos(2u))/2
│   = ∫₀¹ 2·(1-cos(2πx))/2 dx
│   = ∫₀¹ (1 - cos(2πx)) dx
│   = [x - sin(2πx)/(2π)]₀¹
│   = (1 - sin(2π)/(2π)) - (0 - 0)
│   = 1 - 0 = 1 ✓  (properly normalized)
│
├── Quantum uses of integration:
│   ├── Probability: P(a<x<b) = ∫[a to b] |ψ(x)|² dx
│   ├── Normalization: ∫[-∞ to ∞] |ψ(x)|² dx = 1 (MUST be satisfied)
│   ├── Expectation: ⟨E⟩ = ∫ψ*(x)·Ĥ·ψ(x) dx (in position basis)
│   └── Overlap: ⟨φ|ψ⟩ = ∫φ*(x)ψ(x)dx = inner product!
│
├── Code:
│   from scipy import integrate
│   import numpy as np
│   # Normalize ψ₁ = √2 sin(πx) on [0,1]:
│   psi = lambda x: np.sqrt(2) * np.sin(np.pi * x)
│   prob_density = lambda x: abs(psi(x))**2
│   norm, _ = integrate.quad(prob_density, 0, 1)
│   print(f"Norm = {norm:.6f}")  # 1.000000 ✓ (normalized!)
│
└── Exit check:
    1. ∫x²dx from 1 to 3  → Answer: [x³/3]₁³ = 27/3 - 1/3 = 26/3 ≈ 8.667
    2. ∫e^(2x)dx from 0 to 1 → Answer: [e^(2x)/2]₀¹ = e²/2 - 1/2 ≈ 3.195
    3. Verify ψ₂ = √2·sin(2πx) is normalized on [0,1] in code.

═══════════════════════════════════════════
 GATE TO M1.3 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know derivative = "rate of change" / "slope of tangent"
 □ Can apply power, chain, product rules without errors
 □ Differentiated x³sin(x) + e^(2x) correctly
 □ Can compute ∂f/∂x and ∂f/∂y for any 2-variable polynomial
 □ Computed gradient [∂f/∂x, ∂f/∂y] at a specific point
 □ Know: gradient descent → θ_new = θ_old - α·∂E/∂θ (VQE update rule)
 □ Can integrate basic functions (power, exponential, trig)
 □ Verified that ψ₁=√2·sin(πx) normalizes to 1 over [0,1]
═══════════════════════════════════════════
```

---

## Module M1.3: Ordinary Differential Equations (ODE)

> **PREREQUISITES: M1.2 gate passed.**
> Need: derivatives and antiderivatives.
> Real-world intuition: bacterial growth, radioactive decay, population dynamics.

```
M1.3.1  First Order ODEs — Separation of Variables
├── What is an ODE?
│   An equation that involves BOTH a function AND its derivative.
│   Not: "find x when f(x) = 5" (algebraic equation)
│   But: "find f(t) when df/dt = 3f" (differential equation)
│
├── WHY this matters for quantum:
│   Schrödinger equation IS an ODE (or PDE in multiple variables):
│   iℏ dψ/dt = Eψ → solution: ψ(t) = ψ(0)·e^(-iEt/ℏ)
│   Every quantum state evolves according to an ODE!
│
├── Simplest ODE: dy/dt = ky  (exponential growth/decay)
│   This says: "rate of change of y is proportional to y itself"
│   Biology: bacteria double → dP/dt = kP (k>0 growth)
│   Physics: radioactive decay → dN/dt = -λN (λ>0, k=-λ)
│
│   SOLVING by separation of variables:
│   Step 1: Write dy/y = k·dt  (separate y's and t's on different sides)
│   Step 2: Integrate both sides:
│           ∫dy/y = ∫k·dt
│           ln|y| = kt + C
│   Step 3: Exponentiate:
│           |y| = e^(kt+C) = e^C · e^(kt)
│           Rename e^C as A (just a constant):
│           y(t) = A·e^(kt)
│   Step 4: Apply initial condition y(0) = y₀:
│           y(0) = A·e^0 = A = y₀
│           → y(t) = y₀·e^(kt)
│
├── WORKED EXAMPLE 1 — Bacteria:
│   dP/dt = 0.5P  (growth rate 50% per hour)
│   P(0) = 1000 bacteria
│   Solution: P(t) = 1000·e^(0.5t)
│   At t=3h: P(3) = 1000·e^1.5 = 1000×4.48 ≈ 4482 bacteria
│
├── WORKED EXAMPLE 2 — Radioactive decay:
│   dN/dt = -0.1N  (10% decay per unit time)
│   N(0) = 10000 atoms
│   Solution: N(t) = 10000·e^(-0.1t)
│   Half-life: when N = 5000 → 10000e^(-0.1t) = 5000 → t = ln(2)/0.1 ≈ 6.93
│
├── Code:
│   import numpy as np, matplotlib.pyplot as plt
│   t = np.linspace(0, 10, 100)
│   # Bacteria: y(t) = 1000 e^(0.5t)
│   P = 1000 * np.exp(0.5 * t)
│   # Radioactive: N(t) = 10000 e^(-0.1t)
│   N = 10000 * np.exp(-0.1 * t)
│   fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,4))
│   ax1.plot(t, P); ax1.set_title('Bacterial Growth'); ax1.set_xlabel('t (hours)')
│   ax2.plot(t, N); ax2.set_title('Radioactive Decay'); ax2.set_xlabel('t')
│   plt.show()
│
└── Exit check:
    Solve: dy/dt = -3y, y(0) = 2.
    Answer: y(t) = 2e^(-3t).
    At t=1: y(1) = 2e^(-3) ≈ 0.0996. Verify with scipy.

M1.3.2  Complex ODE — The Quantum Case (Schrödinger Time Evolution)
├── What happens if k is IMAGINARY? k = -iE/ℏ
│   (i = imaginary unit, E = energy, ℏ = Planck's constant)
│
├── The Schrödinger time equation:
│   iℏ dψ/dt = Eψ
│   Rearrange: dψ/dt = (-iE/ℏ)ψ  ← same form as dy/dt = ky with k=-iE/ℏ
│
├── Solving by the same method:
│   Solution: ψ(t) = ψ(0) · e^(-iEt/ℏ)
│
├── What does e^(-iEt/ℏ) look like?
│   Using Euler's formula: e^(-iEt/ℏ) = cos(Et/ℏ) - i·sin(Et/ℏ)
│   This oscillates! Like sin/cos, not growing or shrinking.
│   |e^(-iEt/ℏ)|² = cos²(Et/ℏ) + sin²(Et/ℏ) = 1 ALWAYS
│   → Total probability = |ψ(t)|² = |ψ(0)|² · 1 = CONSERVED ✓
│   This is why quantum state evolution preserves probability!
│
├── WORKED VERIFICATION (substitute back):
│   Claim: ψ(t) = ψ(0)·e^(-iEt/ℏ) satisfies iℏ dψ/dt = Eψ
│   LHS: iℏ · dψ/dt = iℏ · ψ(0) · (-iE/ℏ) · e^(-iEt/ℏ)
│                    = iℏ · (-iE/ℏ) · ψ(t)
│                    = -i²E · ψ(t)  [since -i²=+1, i²=-1]
│                    = E · ψ(t)
│   RHS: Eψ(t)
│   LHS = RHS ✓  (substitution works)
│
├── Superposition (two energy levels):
│   If ψ₁ and ψ₂ are solutions, so is their sum:
│   Ψ(t) = c₁·ψ₁(0)·e^(-iE₁t/ℏ) + c₂·ψ₂(0)·e^(-iE₂t/ℏ)
│   Each component rotates at its own frequency ω=E/ℏ
│   Beat frequency = (E₂-E₁)/ℏ = oscillation in ⟨ψ|X|ψ⟩ etc.
│   This is quantum INTERFERENCE in time!
│
└── TISE connection:
    When looking for steady states: ψ(x,t) = φ(x)·e^(-iEt/ℏ)
    Plug into TDSE: [iℏ∂/∂t][φ(x)e^(-iEt/ℏ)] = Ĥ[φ(x)e^(-iEt/ℏ)]
    LHS: iℏ·φ(x)·(-iE/ℏ)e^(-iEt/ℏ) = E·φ(x)e^(-iEt/ℏ)
    RHS: e^(-iEt/ℏ)·Ĥφ(x)
    Cancel e^(-iEt/ℏ): Ĥφ(x) = Eφ(x)    ← THIS IS TISE!
    Eigenvalue equation: finding φ(x) and E is VQE's entire job.

═══════════════════════════════════════════
 GATE TO M2.1 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Can solve dy/dt = ky → y(t) = y₀·e^(kt) by separation of variables
 □ Solved bacteria example: P(t) = 1000e^(0.5t), P(3) ≈ 4482
 □ Know: Schrödinger equation dψ/dt = (-iE/ℏ)ψ is SAME FORM as dy/dt=ky
 □ Know solution: ψ(t) = ψ(0)·e^(-iEt/ℏ) (phase rotation, not growth!)
 □ Verified by substitution that ψ(t)=ψ(0)e^(-iEt/ℏ) satisfies iℏdψ/dt=Eψ
 □ Know: |e^(-iEt/ℏ)|=1 always → probability conserved  
 □ Can derive TISE (Ĥφ=Eφ) by separation of variables from TDSE
═══════════════════════════════════════════
```

---

## Module M2.1: Vectors & Vector Spaces

> **PREREQUISITES: M1.3 gate passed.**
> Also from Class 10 Coordinate Geometry: plotting points (x,y), distance formula.
> Physics Class 11: vector vs scalar, adding vectors (parallelogram law).
> If shaky → revise NCERT Class 11 Physics Ch4 (vectors overview).

```
M2.1.1  What is a Vector — From Arrows to Lists of Numbers
├── You already know (Class 10/11):
│   Position in 2D: point (x, y) → arrow from (0,0) to (x,y)
│   This arrow IS a vector: v = [x, y]ᵀ
│   Length = √(x²+y²) (Pythagoras!)
│
├── But vectors can have MORE than 2 entries:
│   3D: v = [x, y, z]ᵀ  (arrow in 3D space)
│   4D: v = [a, b, c, d]ᵀ (can't visualize, but math works the same!)
│   n-D: v = [v₁, v₂, ..., vₙ]ᵀ
│
├── COLUMN vector (used in quantum):
│   v = [[v₁],   ← this is a 3×1 matrix (3 rows, 1 column)
│        [v₂],
│        [v₃]]
│   Python writes it as: np.array([[1],[0]]) or np.array([1, 0])
│
├── Vector OPERATIONS (worked examples):
│   Addition: [3,4]ᵀ + [1,-2]ᵀ = [3+1, 4+(-2)]ᵀ = [4, 2]ᵀ
│   Scalar multiply: 3·[2, -1]ᵀ = [6, -3]ᵀ
│   Scale-then-add: 2·[1,0]ᵀ + 3·[0,1]ᵀ = [2,0]ᵀ + [0,3]ᵀ = [2,3]ᵀ
│
├── Quantum connection:
│   |0⟩ = [1, 0]ᵀ   (qubit "pointing up")
│   |1⟩ = [0, 1]ᵀ   (qubit "pointing down")
│   Superposition: |ψ⟩ = α|0⟩ + β|1⟩ = [α, β]ᵀ
│   Example: |+⟩ = (1/√2)|0⟩ + (1/√2)|1⟩ = [1/√2, 1/√2]ᵀ
│
├── DOT PRODUCT (inner product in ℝ²):
│   u·v = u₁v₁ + u₂v₂ + ... (multiply entry by entry, sum)
│   [3,4]·[1,2] = 3·1 + 4·2 = 3+8 = 11
│   Geometric meaning: measures how ALIGNED two vectors are
│   u·v = |u||v|cos(θ) → θ=0 (same dir): maximum, θ=90°: ZERO
│   Orthogonal (perpendicular): u·v = 0
│
├── Code:
│   import numpy as np
│   ket0 = np.array([1, 0], dtype=complex)
│   ket1 = np.array([0, 1], dtype=complex)
│   # Superposition |+⟩ state:
│   ket_plus = (1/np.sqrt(2))*ket0 + (1/np.sqrt(2))*ket1
│   print(ket_plus)         # [0.707, 0.707]
│   print(np.dot(ket0, ket1))  # 0  (orthogonal!)
│   print(np.dot(ket0, ket0))  # 1  (normalized!)
│
└── Exit check:
    1. Add [2+3i, 1-i]ᵀ + [1, 2i]ᵀ by hand.
    2. Compute |ψ⟩ = (0.6)|0⟩ + (0.8)|1⟩. Verify |0.6|²+|0.8|² = 1.

M2.1.2  Linear Combinations, Span, Independence
├── Linear combination: αu + βv (stretch u by α, stretch v by β, add)
│   Example: 2·[1,0]ᵀ + 3·[0,1]ᵀ = [2,3]ᵀ
│           -1·[1,0]ᵀ + 4·[0,1]ᵀ = [-1,4]ᵀ
│   → By choosing different α,β we can reach EVERY point in 2D!
│
├── Span: the SET of all possible linear combinations
│   Span{[1,0], [0,1]} = ALL of ℝ² (reach every 2D point)
│   Span{[1,0], [2,0]} = just the x-axis (only reach points [t,0])
│
├── Linear INDEPENDENCE: none of the vectors is redundant
│   {[1,0], [0,1]}: INDEPENDENT ✓ (can't make [1,0] from [0,1])
│   {[1,0], [2,0]}: DEPENDENT ✗  ([2,0] = 2×[1,0], redundant!)
│   Test: put as rows, row-reduce. Any ZERO row → dependent.
│
├── Basis: a linearly independent set that spans the whole space
│   {[1,0], [0,1]}: basis of ℝ² (standard basis)
│   {|0⟩, |1⟩}: basis of qubit space ℂ²  (computational basis)
│   {|+⟩, |-⟩}: also a basis of ℂ² (Hadamard basis)
│   → Any basis works for expansion, different bases useful for different tasks
│
├── Why this matters for VQE:
│   Measuring in different BASES:
│   Z-basis {|0⟩,|1⟩}: measure Pauli Z
│   X-basis {|+⟩,|-⟩}: measure Pauli X (apply H gate first, then measure Z)
│   Y-basis: measure Pauli Y (apply specific rotation first)
│
├── Code:
│   v1 = np.array([1, 2, 3])
│   v2 = np.array([4, 5, 6])
│   v3 = np.array([7, 8, 9])
│   M = np.array([v1, v2, v3])  # as rows
│   rank = np.linalg.matrix_rank(M)
│   print(f"Rank = {rank}")  # 2 (not 3) → linearly DEPENDENT
│
└── Exit check:
    Are {[2,1], [1,3]} linearly independent?
    Method: det([[2,1],[1,3]]) = 6-1 = 5 ≠ 0 → INDEPENDENT ✓

M2.1.3  Complex Vector Spaces ℂⁿ
├── Extension: allow complex entries in vectors
│   ℝ² = pairs of real numbers (ordinary 2D plane)
│   ℂ² = pairs of COMPLEX numbers (qubit lives here)
│   Example: v = [(1+i)/√2, (1-i)/√2]ᵀ ∈ ℂ²
│
├── Inner product in ℂⁿ (IMPORTANT difference from real):
│   Real: ⟨u,v⟩ = Σᵢ uᵢvᵢ   (multiply entries, sum)
│   Complex: ⟨u,v⟩ = Σᵢ uᵢ*vᵢ  (conjugate FIRST vector!)
│   Why conjugate? So ⟨v,v⟩ = Σ|vᵢ|² is always REAL and ≥ 0.
│
│   WORKED EXAMPLE:
│   u = [1, i]ᵀ, v = [1, 1]ᵀ
│   ⟨u,v⟩ = (1)*·1 + (i)*·1 = 1·1 + (-i)·1 = 1 - i
│   ⟨v,u⟩ = (1)*·1 + (1)*·i = 1·1 + 1·i = 1 + i = (1-i)* ✓ (conjugate swap)
│
├── Why exponential dimension:
│   1 qubit → ℂ²: 2 complex numbers (4 real parameters)
│   2 qubits → ℂ⁴: 4 complex numbers
│   n qubits → ℂ^(2ⁿ): 2ⁿ complex numbers
│   50 qubits → 2⁵⁰ ≈ 10¹⁵ complex numbers (!!)
│   H₂ molecule: 4 spin-orbitals → 2⁴=16 basis states → 16×16 Hamiltonian
│
├── Code:
│   u = np.array([1, 1j])          # u = [1, i]ᵀ
│   v = np.array([1, 1])           # v = [1, 1]ᵀ
│   inner = np.dot(u.conj(), v)    # ⟨u|v⟩ = u†v
│   print(f"⟨u|v⟩ = {inner}")      # (1-1j)
│
└── Exit check:
    Compute ⟨+|−⟩ using complex inner product.
    |+⟩ = [1/√2, 1/√2]ᵀ, |−⟩ = [1/√2, -1/√2]ᵀ
    ⟨+|−⟩ = (1/√2)*(1/√2) + (1/√2)*(-1/√2) = 1/2 - 1/2 = 0 (orthogonal ✓)

═══════════════════════════════════════════
 GATE TO M2.2 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know: vector = ordered list of numbers, column vector notation
 □ Can add and scalar-multiply vectors by hand
 □ Know: dot product u·v = Σ uᵢvᵢ, orthogonal when u·v=0
 □ Know: span, linear independence, basis (with examples)
 □ Know: |0⟩=[1,0]ᵀ, |1⟩=[0,1]ᵀ, |+⟩=[1/√2, 1/√2]ᵀ, |−⟩=[1/√2, -1/√2]ᵀ
 □ Know: complex inner product conjugates the FIRST vector
 □ Verified ⟨+|−⟩ = 0 by hand AND NumPy
 □ Know: n qubits → 2ⁿ-dimensional space (exponential!)
═══════════════════════════════════════════
```

---

## Module M2.2: Matrices & Operations

> **PREREQUISITES: M2.1 gate passed.**
> Think of a matrix as a "function machine" that transforms vectors.
> From Class 10: coordinate geometry, grid transformations (rotation).

```
M2.2.1  Matrix as Linear Transformation — Gate = Matrix
├── Matrix M: a rectangular grid of numbers (m rows × n columns)
│   Applied to a vector: new_vec = M · old_vec
│   Every row of M "mixes" the entries of the input vector
│
├── Matrix-vector product — STEP BY STEP:
│   M = [[a, b], [c, d]],  v = [x, y]ᵀ
│   Mv = [a·x + b·y,   ← (row 1)·v
│          c·x + d·y]  ← (row 2)·v
│
│   WORKED EXAMPLE (Pauli X gate):
│   X = [[0,1],[1,0]],  |0⟩ = [1,0]ᵀ
│   X|0⟩ = [0·1 + 1·0,  = [0, 1]ᵀ = |1⟩
│           1·1 + 0·0]
│   → X gate FLIPS |0⟩ to |1⟩ (NOT gate) ✓
│
│   WORKED EXAMPLE (Hadamard gate):
│   H = (1/√2)[[1,1],[1,-1]],  |0⟩ = [1,0]ᵀ
│   H|0⟩ = (1/√2)[1·1+1·0, 1·1+(-1)·0]ᵀ = (1/√2)[1, 1]ᵀ = |+⟩ ✓
│   H|1⟩ = (1/√2)[1·0+1·1, 1·0+(-1)·1]ᵀ = (1/√2)[1,-1]ᵀ = |−⟩ ✓
│
├── Code:
│   X = np.array([[0,1],[1,0]], dtype=complex)
│   H = np.array([[1,1],[1,-1]], dtype=complex) / np.sqrt(2)
│   ket0 = np.array([1,0], dtype=complex)
│   ket1 = np.array([0,1], dtype=complex)
│   print("X|0⟩ =", X @ ket0)   # [0, 1] = |1⟩
│   print("H|0⟩ =", H @ ket0)   # [0.707, 0.707] = |+⟩
│   print("H|1⟩ =", H @ ket1)   # [0.707, -0.707] = |−⟩
│
└── Exit check:
    Compute Z|1⟩ by hand.  Z = [[1,0],[0,-1]], |1⟩ = [0,1]ᵀ
    Answer: [0,-1]ᵀ = -|1⟩  (phase flip: |1⟩ gets a -1 phase)

M2.2.2  Matrix Multiplication — Applying Two Gates in Sequence
├── Two gates: first apply A, then B → total gate = B·A (NOT A·B!)
│   (A·B)ᵢⱼ = row i of A · column j of B = Σₖ Aᵢₖ Bₖⱼ
│
├── NOT commutative: AB ≠ BA in general!
│   X·Z ≠ Z·X (this is the commutator [X,Z] ≠ 0 from Ph2.2!)
│
├── WORKED EXAMPLE — X·Z vs Z·X:
│   XZ = [[0,1],[1,0]] · [[1,0],[0,-1]]
│      = [[0·1+1·0, 0·0+1·(-1)],[1·1+0·0, 1·0+0·(-1)]]
│      = [[0,-1],[1,0]]
│
│   ZX = [[1,0],[0,-1]] · [[0,1],[1,0]]
│      = [[0,1],[-1,0]]
│
│   XZ ≠ ZX ✓ (they DON'T commute — explains quantum uncertainty!)
│
├── Dimension rule: (m×k)·(k×n) = (m×n) — middle dims must match!
│   (2×2)·(2×1) = (2×1)  ← gate applied to qubit state ✓
│   (2×2)·(2×2) = (2×2)  ← two gates composed ✓
│   (2×2)·(3×3) = ERROR  ← dimensions don't match!
│
├── Code:
│   X = np.array([[0,1],[1,0]], dtype=complex)
│   Z = np.array([[1,0],[0,-1]], dtype=complex)
│   print("XZ =", X @ Z)
│   print("ZX =", Z @ X)
│   print("XZ == ZX?", np.allclose(X@Z, Z@X))  # False ✓
│
└── Exit check (DO BY HAND in 90 seconds):
    HXH = ?   (H = (1/√2)[[1,1],[1,-1]], X = [[0,1],[1,0]])
    Step 1: XH = [[0,1],[1,0]]·(1/√2)[[1,1],[1,-1]] = (1/√2)[[1,-1],[1,1]]
    Step 2: H·(XH) = (1/√2)[[1,1],[1,-1]]·(1/√2)[[1,-1],[1,1]] = [[1,0],[0,-1]] = Z
    Result: HXH = Z (H transforms X into Z basis!) ✓

M2.2.3  Transpose, Conjugate, Dagger (†)
├── Transpose Aᵀ: reflect across main diagonal (rows ↔ columns)
│   A  = [[1, 2],     Aᵀ = [[1, 3],
│         [3, 4]]           [2, 4]]
│
├── Complex conjugate A*: negate all imaginary parts
│   A  = [[1+i, 2],    A* = [[1-i, 2],
│         [3, 4-2i]]        [3, 4+2i]]
│
├── Dagger A† = (A*)ᵀ = conjugate transpose (MOST USED in quantum):
│   Step 1: conjugate (flip signs of imaginary parts)
│   Step 2: transpose (swap rows and columns)
│
│   WORKED EXAMPLE:
│   A = [[1+i, 2 ],
│        [3,   4-2i]]
│   Step 1 conjugate: [[1-i, 2], [3, 4+2i]]
│   Step 2 transpose: A† = [[1-i, 3 ],
│                            [2,   4+2i]]
│
│   Key property: (AB)† = B†A† (order REVERSES)
│   For bra: ⟨ψ| = |ψ⟩†  (bra = dagger of ket)
│
├── Code:
│   A = np.array([[1+1j, 2], [3, 4-2j]])
│   A_dag = A.conj().T    # dagger = conjugate transpose
│   print("A ="); print(A)
│   print("A† ="); print(A_dag)
│
└── Exit check:
    V = [[1+i, 0], [0, 1-i]] — compute V†V. Is it I? (Is V unitary?)
    V† = [[1-i, 0],[0, 1+i]]
    V†V = [[(1-i)(1+i), 0],[0, (1+i)(1-i)]] = [[2, 0],[0, 2]] ≠ I
    → NOT unitary! (|1+i| = √2 ≠ 1 → would need to divide by √2)

M2.2.4  Determinant, Trace, Inverse — Key Properties
├── DETERMINANT det(A) — 2×2:
│   det([[a,b],[c,d]]) = ad - bc
│
│   Geometric meaning: scaling factor (area scaling)
│   det = 0: matrix is SINGULAR (no inverse, not reversible)
│   For quantum gates: |det(U)| = 1 (unitary gates preserve volume)
│
│   WORKED EXAMPLES:
│   det(X) = det([[0,1],[1,0]]) = 0·0 - 1·1 = -1  (|det| = 1 ✓ unitary)
│   det(Z) = det([[1,0],[0,-1]]) = (1)(-1) - 0 = -1 ✓
│   det(H) = det((1/√2)[[1,1],[1,-1]]) = (1/2)(1(-1)-1·1) = (1/2)(-2) = -1 ✓
│
├── TRACE Tr(A) — sum of diagonal elements:
│   Tr([[a,b],[c,d]]) = a + d
│   Key property: Tr(AB) = Tr(BA) (cyclic!)
│   Tr(X) = 0 + 0 = 0,  Tr(Z) = 1 + (-1) = 0
│   Tr(I) = 1 + 1 = 2
│   Quantum use: ⟨A⟩ = Tr(ρA) for mixed states
│   VQE: cₖ = Tr(σₖH)/2  (extract Pauli coefficient from Hamiltonian)
│
├── INVERSE A⁻¹ — for 2×2:
│   A⁻¹ = (1/det)·[[d,-b],[-c,a]]
│   For quantum gates: U⁻¹ = U† (dagger IS the inverse!)
│   This is why quantum gates are reversible by applying U† after U.
│
├── Code:
│   X = np.array([[0,1],[1,0]], dtype=complex)
│   print(f"det(X) = {np.linalg.det(X):.4f}")  # -1.0000
│   print(f"Tr(X)  = {np.trace(X):.4f}")        # 0.0000
│   print("X†X =", X.conj().T @ X)              # Identity ✓
│
└── Exit check:
    A = [[3,1],[2,2]]. Compute det(A), Tr(A), A⁻¹.
    det = 3·2 - 1·2 = 4
    Tr  = 3 + 2 = 5
    A⁻¹ = (1/4)[[2,-1],[-2,3]]
    Verify A·A⁻¹ = I in NumPy.

═══════════════════════════════════════════
 GATE TO M2.3 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Can compute matrix-vector product Mv by hand for 2×2 case
 □ Verified: X|0⟩=|1⟩, H|0⟩=|+⟩, Z|1⟩=-|1⟩ by matrix multiplication
 □ Can multiply 2×2 matrices by hand; know AB≠BA in general
 □ Computed XZ and ZX and confirmed they differ
 □ Can compute A† (dagger) for any 2×2 complex matrix
 □ Know: det=0 → singular; |det(U)|=1 for quantum gates; U⁻¹=U†
 □ Know Tr(A) = sum of diagonal; Tr(X)=Tr(Y)=Tr(Z)=0, Tr(I)=2
 □ Verified HXH = Z by matrix multiplication
═══════════════════════════════════════════
```

---

## Module M2.3: Special Matrices — Hermitian & Unitary

> **PREREQUISITES: M2.2 gate passed.**
> Need: dagger A†, matrix multiplication, trace, determinant.

```
M2.3.1  Hermitian Matrices — Real Measurements
├── Definition: A is Hermitian if A = A†  (equals its own dagger)
│
├── Why Hermitian matrices give REAL eigenvalues:
│   Suppose Av = λv (eigenvalue equation)
│   Take dagger of both sides: v†A† = λ*v†
│   Since A†=A: v†A = λ*v†
│   Now multiply original (Av=λv) on left by v†:
│   v†Av = λ(v†v) = λ||v||²
│   Multiply transformed on right by v: v†Av = λ*||v||²
│   So: λ||v||² = λ*||v||² → λ = λ* → λ is REAL ✓
│
├── WORKED CHECK — is Z Hermitian?
│   Z = [[1,0],[0,-1]]   (all real, diagonal)
│   Z† = (Z*)ᵀ = Z*ᵀ = Zᵀ = Z  (transpose of real diagonal = itself)
│   Z† = Z ✓ → Hermitian
│
├── WORKED CHECK — is A = [[3, 1+i],[1-i, 2]] Hermitian?
│   A† = step1: A* = [[3, 1-i],[1+i, 2]]
│        step2: transpose = [[3, 1+i],[1-i, 2]] = A ✓
│   → Hermitian. Eigenvalues will be real.
│
├── Code:
│   def is_hermitian(A, tol=1e-10):
│       return np.allclose(A, A.conj().T, atol=tol)
│   I = np.eye(2, dtype=complex)
│   X = np.array([[0,1],[1,0]], dtype=complex)
│   Y = np.array([[0,-1j],[1j,0]])
│   Z = np.array([[1,0],[0,-1]], dtype=complex)
│   for name, M in [('I',I),('X',X),('Y',Y),('Z',Z)]:
│       print(f"{name} Hermitian: {is_hermitian(M)}")  # all True
│
└── Exit check: A = [[2, 1+2i],[1-2i, 3]]. Is it Hermitian? Find eigenvalues.
    Answer: Yes Hermitian. Eigenvalues: λ = 2.5 ± √(0.25+5) = 2.5 ± 2.29...
    Verify with np.linalg.eigvalsh (returns REAL for Hermitian).

M2.3.2  Unitary Matrices — Probability Preserving Gates
├── Definition: U†U = UU† = I  (dagger = inverse)
├── Consequences:
│   1) ||Uv|| = ||v||  (length preserved → probability preserved!)
│   2) Columns are orthonormal (inner products preserved)
│   3) U⁻¹ = U† → every quantum gate is reversible!
│
├── WORKED CHECK — is H Hadamard unitary?
│   H = (1/√2)[[1,1],[1,-1]]
│   H† = (1/√2)[[1,1],[1,-1]] (H is real → H†=Hᵀ = H itself!)
│   HH† = (1/2)[[1,1],[1,-1]]·[[1,1],[1,-1]]
│        = (1/2)[[1+1, 1-1],[1-1, 1+1]] = (1/2)[[2,0],[0,2]] = I ✓
│
├── Code:
│   def is_unitary(U, tol=1e-10):
│       return np.allclose(U.conj().T @ U, np.eye(len(U)), atol=tol)
│   H = np.array([[1,1],[1,-1]], dtype=complex) / np.sqrt(2)
│   print(f"H unitary: {is_unitary(H)}")  # True
│   for M in [X, Y, Z, H]:
│       print(is_unitary(M))  # all True for Paulis and H
│
└── Exit check: Rx(π) = [[0, -i],[-i, 0]]. Verify this is unitary by hand.
    Rx(π)† = [[0, i],[i, 0]]. Product: [[0,-i],[-i,0]]·[[0,i],[i,0]] = I ✓

M2.3.3  Pauli Matrices — Memorize All Four
├── The four Paulis (MEMORIZE — used every day in QC):
│   I = [[1, 0],    X = [[0, 1],    Y = [[0,-i],    Z = [[1, 0],
│        [0, 1]]         [1, 0]]         [i, 0]]         [0,-1]]
│
├── Key properties of ALL four:
│   Hermitian (σ=σ†) AND Unitary (σ†σ=I)  → both observable AND gate!
│   Self-inverse: X²=Y²=Z²=I  (applying twice = doing nothing)
│   Traceless: Tr(X)=Tr(Y)=Tr(Z)=0, Tr(I)=2
│
├── Algebra (memorize these products):
│   XY = iZ,  YZ = iX,  ZX = iY
│   YX = -iZ, ZY = -iX, XZ = -iY
│   XX = YY = ZZ = I  (self-inverse)
│
├── Physical meaning:
│   X: bit flip: X|0⟩=|1⟩, X|1⟩=|0⟩  (quantum NOT gate)
│   Z: phase flip: Z|0⟩=|0⟩, Z|1⟩=-|1⟩
│   Y: both: Y|0⟩=i|1⟩, Y|1⟩=-i|0⟩
│   Z eigenvalues: +1 for |0⟩, -1 for |1⟩  ← measuring Z = reading qubit!
│
├── Code:
│   I = np.eye(2, dtype=complex)
│   X = np.array([[0,1],[1,0]], dtype=complex)
│   Y = np.array([[0,-1j],[1j,0]])
│   Z = np.array([[1,0],[0,-1]], dtype=complex)
│   print(np.allclose(X@Y, 1j*Z))  # True: XY=iZ
│   print(np.allclose(X@Y@Z, 1j*I))  # True: XYZ=iI
│
└── Exit check: Write all 4 from memory in 60 sec. Verify XYZ=iI in NumPy.

M2.3.4  Hadamard & Rotation Gates
├── Hadamard H = (1/√2)[[1,1],[1,-1]]:
│   H|0⟩ = |+⟩ = (1/√2)(|0⟩+|1⟩)  ← creates equal superposition
│   H|1⟩ = |-⟩ = (1/√2)(|0⟩-|1⟩)
│   H² = I  (applying H twice = identity)
│   H changes basis: Z-basis ↔ X-basis
│
├── Rotation gates (key for VQE ansatz — these are VQE's parameters):
│   Ry(θ): [[cos(θ/2), -sin(θ/2)], [sin(θ/2), cos(θ/2)]]
│   Rz(θ): [[e^(-iθ/2), 0], [0, e^(iθ/2)]]
│   VQE: choose θ₁,θ₂,...,θₙ → minimize ⟨H⟩
│
├── Pauli decomposition (HOW Hamiltonians are built):
│   Any 2×2 Hermitian M = c₀I + c₁X + c₂Y + c₃Z
│   Coefficients: cₖ = Tr(σₖM)/2
│
│   WORKED EXAMPLE: decompose Z:
│   c₀ = Tr(I·Z)/2 = Tr([[1,0],[0,-1]])/2 = 0/2 = 0
│   c₁ = Tr(X·Z)/2 = Tr([[0,-1],[1,0]])/2 = 0/2 = 0
│   c₃ = Tr(Z·Z)/2 = Tr(I)/2 = 2/2 = 1
│   → Z = 1·Z ✓ (only one term, coefficient 1)
│
└── Exit check:
    Decompose M = [[3,1+i],[1-i,2]] = c₀I + c₁X + c₂Y + c₃Z.
    c₀=Tr(IM)/2=5/2=2.5, c₁=Tr(XM)/2=1, c₂=Tr(YM)/2=1, c₃=Tr(ZM)/2=0.5
    Verify: 2.5I + 1X + 1Y + 0.5Z = M using NumPy.

═══════════════════════════════════════════
 GATE TO M3.1 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know definition: Hermitian A=A†, Unitary U†U=I
 □ Verified: I,X,Y,Z are ALL Hermitian AND Unitary
 □ Know WHY Hermitian → real eigenvalues (see proof above)
 □ Know WHY Unitary → probability conserved (||Uv||=||v||)
 □ Can write all 4 Pauli matrices from memory
 □ Verified XY=iZ and XYZ=iI in NumPy
 □ Know: H maps |0⟩↔|+⟩, |1⟩↔|-⟩ (Hadamard switches bases)
 □ Know Pauli decomposition formula: cₖ = Tr(σₖM)/2
═══════════════════════════════════════════
```

---

## Module M3.1: Eigenvalues & Eigenvectors

> **PREREQUISITES: M2.3 gate passed.**
> Need: matrix multiplication, determinant, solving linear equations.

```
M3.1.1  What is an Eigenvector — The Special Directions
├── For most vectors v: Av points in a DIFFERENT direction from v
│   But some special vectors v stay in same direction after A: Av = λv
│   These are EIGENVECTORS, λ is the EIGENVALUE (stretch factor)
│
├── The central insight:
│   Z|0⟩ = +1·|0⟩  → |0⟩ is eigenvector of Z with eigenvalue +1
│   Z|1⟩ = -1·|1⟩  → |1⟩ is eigenvector of Z with eigenvalue -1
│   Measuring Z on |0⟩ always gives +1. On |1⟩ always gives -1.
│   Measuring Z on |+⟩ = gives +1 or -1 with equal probability (50/50)
│
├── VQE goal restated:
│   H|ground⟩ = E₀|ground⟩  (TISE — eigenvalue equation!)
│   E₀ = MINIMUM eigenvalue of molecular Hamiltonian H
│   VQE: approximate |ground⟩ without solving eigenvalue equation directly
│
└── Watch FIRST: 3Blue1Brown EoLA Ep 13-14 (eigenvectors/values visual)

M3.1.2  Computing Eigenvalues — Characteristic Equation
├── Method: det(A - λI) = 0
│   This gives a polynomial in λ whose roots = eigenvalues
│   For each eigenvalue λₖ: solve (A-λₖI)v=0 to find eigenvector vₖ
│
├── WORKED EXAMPLE — eigenvalues of Z:
│   Z - λI = [[1-λ, 0],[0, -1-λ]]
│   det = (1-λ)(-1-λ) - 0 = -(1-λ)(1+λ) = -(1-λ²) = λ²-1 = 0
│   λ² = 1 → λ = ±1
│   λ=+1: (Z-I)v = [[0,0],[0,-2]]v = 0 → v = [1,0]ᵀ = |0⟩ ✓
│   λ=-1: (Z+I)v = [[2,0],[0,0]]v = 0 → v = [0,1]ᵀ = |1⟩ ✓
│
├── WORKED EXAMPLE — eigenvalues of H_molecule_toy = 0.5Z + 0.5X:
│   H = [[0.5, 0.5],[0.5, -0.5]]
│   det(H-λI) = (0.5-λ)(-0.5-λ) - 0.25 = -(0.25-λ²) - 0.25 = λ²-0.5 = 0
│   λ = ±1/√2 ≈ ±0.707
│   Ground state energy E₀ = -1/√2 ≈ -0.707 (minimum eigenvalue)
│
├── Code:
│   Z = np.array([[1,0],[0,-1]], dtype=complex)
│   X = np.array([[0,1],[1,0]], dtype=complex)
│   H = 0.5*Z + 0.5*X
│   vals, vecs = np.linalg.eigh(H)  # eigh for Hermitian (real eigenvalues)
│   print(f"Eigenvalues: {vals}")   # [-0.707, +0.707]
│   print(f"Ground state E₀ = {min(vals):.4f}")  # -0.7071
│
└── Exit check:
    A = [[4,1],[2,3]]. Compute eigenvalues by hand.
    det([[4-λ,1],[2,3-λ]]) = (4-λ)(3-λ)-2 = λ²-7λ+10 = 0
    → (λ-5)(λ-2) = 0 → λ=5, λ=2

M3.1.3  Spectral Theorem — The Most Important Result
├── For Hermitian A with eigenvalues λₙ and orthonormal eigenvectors |vₙ⟩:
│   A = Σₙ λₙ |vₙ⟩⟨vₙ|  (spectral decomposition)
│   → Any Hermitian = sum of projections scaled by eigenvalues
│
├── WORKED EXAMPLE — spectral decomposition of Z:
│   Z = (+1)|0⟩⟨0| + (-1)|1⟩⟨1|
│      = [[1,0],[0,0]] - [[0,0],[0,1]] = [[1,0],[0,-1]] ✓
│
├── Link to expectation values:
│   ⟨A⟩ = ⟨ψ|A|ψ⟩ = ⟨ψ|(Σₙ λₙ|vₙ⟩⟨vₙ|)|ψ⟩ = Σₙ λₙ|⟨vₙ|ψ⟩|²
│   = weighted average of eigenvalues with Born rule probabilities ✓
│
└── VQE variational principle:
    For ANY |ψ⟩: ⟨ψ|H|ψ⟩ ≥ E₀ (ground state energy)
    Proof: ⟨ψ|H|ψ⟩ = Σₙ Eₙ|cₙ|² ≥ E₀·Σₙ|cₙ|² = E₀·1 = E₀
    → Minimizing ⟨ψ|H|ψ⟩ over all |ψ⟩ gives E₀ exactly!

═══════════════════════════════════════════
 GATE TO M3.2 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Know: Av=λv (eigenvector = special direction, eigenvalue = stretch)
 □ Know: Z eigenstates are |0⟩(λ=+1) and |1⟩(λ=-1)
 □ Can compute eigenvalues using det(A-λI)=0 for 2×2
 □ Computed eigenvalues of 0.5Z+0.5X → ±1/√2 ≈ ±0.707
 □ Know spectral theorem: A = Σλₙ|vₙ⟩⟨vₙ|
 □ Verified Z = |0⟩⟨0| - |1⟩⟨1| in NumPy
 □ Know variational principle: ⟨ψ|H|ψ⟩ ≥ E₀ for all |ψ⟩
═══════════════════════════════════════════
```

---

## Module M3.2: Inner Products & Orthogonality

> **PREREQUISITES: M3.1 gate passed + M2.1.3 (complex inner product).**

```
M3.2.1  Inner Product — The Overlap Between States
├── Real: ⟨u,v⟩ = u·v = Σᵢ uᵢvᵢ
├── Complex (quantum): ⟨u,v⟩ = u†v = Σᵢ uᵢ*vᵢ  (conjugate FIRST!)
├── Norm: ||v|| = √⟨v,v⟩ = √Σ|vᵢ|²
├── Orthogonal: ⟨u,v⟩ = 0  (states with no overlap)
├── Orthonormal: ⟨eᵢ,eⱼ⟩ = δᵢⱼ (1 if i=j, 0 otherwise)
│
├── WORKED EXAMPLES:
│   ⟨0|0⟩ = [1,0]·[1,0] = 1 ✓
│   ⟨0|1⟩ = [1,0]·[0,1] = 0 ✓ (orthogonal!)
│   ⟨+|+⟩ = (1/√2)[1,1]·(1/√2)[1,1] = ½(1+1) = 1 ✓
│   ⟨+|-⟩ = (1/√2)[1,1]·(1/√2)[1,-1] = ½(1-1) = 0 ✓
│
└── Code:
    ket0 = np.array([1,0], dtype=complex)
    ket1 = np.array([0,1], dtype=complex)
    ket_p = np.array([1,1], dtype=complex)/np.sqrt(2)
    ket_m = np.array([1,-1], dtype=complex)/np.sqrt(2)
    print(np.dot(ket0.conj(), ket1))  # 0
    print(np.dot(ket_p.conj(), ket_m))  # 0

M3.2.2  Born Rule = Inner Product Squared
├── P(measuring |φ⟩ when in state |ψ⟩) = |⟨φ|ψ⟩|²
│
├── WORKED EXAMPLES:
│   |ψ⟩ = (√3/2)|0⟩ + (1/2)|1⟩
│   P(|0⟩) = |⟨0|ψ⟩|² = |√3/2|² = 3/4 = 75%
│   P(|1⟩) = |⟨1|ψ⟩|² = |1/2|² = 1/4 = 25%
│   Sum: 75%+25% = 100% ✓
│
│   |ψ⟩ = (1/√2)|0⟩ + (i/√2)|1⟩  (complex amplitude!)
│   P(|0⟩) = |1/√2|² = 1/2
│   P(|1⟩) = |i/√2|² = |i|²/2 = 1/2 ✓ (equal despite complex amplitude)
│
│   P(|+⟩) for |ψ⟩=(√3/2)|0⟩+(1/2)|1⟩:
│   ⟨+|ψ⟩ = (1/√2)[1,1]·[√3/2, 1/2] = (1/√2)(√3/2+1/2)
│   P(|+⟩) = |(√3+1)/(2√2)|² = (√3+1)²/8 = (4+2√3)/8 = (2+√3)/4 ≈ 0.933
│
└── Exit check:
    State |ψ⟩=(3/5)|0⟩+(4/5)|1⟩. Compute P(|0⟩), P(|1⟩), P(|+⟩).
    P(|0⟩)=9/25, P(|1⟩)=16/25. P(|+⟩)=|(3/5+4/5)/(√2)|² = (49/25)/2 = 49/50.

M3.2.3  Gram-Schmidt — Building Orthonormal Bases
├── Problem: given any two linearly independent vectors, make them orthonormal
├── Algorithm:
│   Step 1: normalize u → e₁ = u/||u||
│   Step 2: subtract projection: v₂ = v - ⟨e₁,v⟩·e₁
│   Step 3: normalize v₂ → e₂ = v₂/||v₂||
│   Result: {e₁, e₂} — orthonormal!
│
├── WORKED EXAMPLE: u=[1,1], v=[1,0]
│   e₁ = [1,1]/√2 = [1/√2, 1/√2]
│   ⟨e₁,v⟩ = [1/√2, 1/√2]·[1,0] = 1/√2
│   v₂ = [1,0] - (1/√2)·[1/√2, 1/√2] = [1,0] - [1/2, 1/2] = [1/2, -1/2]
│   e₂ = [1/2,-1/2]/||[1/2,-1/2]|| = [1/2,-1/2]/(1/√2) = [1/√2, -1/√2]
│   → {e₁, e₂} = {|+⟩, |-⟩}  ← this IS the Hadamard basis!
│
└── Code:
    def gram_schmidt(u, v):
        e1 = u / np.linalg.norm(u)
        v2 = v - np.dot(e1.conj(), v) * e1
        e2 = v2 / np.linalg.norm(v2)
        return e1, e2

═══════════════════════════════════════════
 GATE TO M3.3 — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Can compute complex inner product ⟨φ|ψ⟩ = Σ φᵢ*ψᵢ by hand
 □ Verified ⟨0|1⟩=0, ⟨+|-⟩=0 by hand AND NumPy
 □ Can apply Born rule: P(|φ⟩) = |⟨φ|ψ⟩|² with complex amplitudes
 □ Computed P(|0⟩)=9/25, P(|1⟩)=16/25 for |ψ⟩=(3/5)|0⟩+(4/5)|1⟩
 □ Know: |⟨φ|ψ⟩|² = 1 means states are equal; = 0 means orthogonal
 □ Can apply Gram-Schmidt to two vectors to get orthonormal basis
═══════════════════════════════════════════
```

---

## Module M3.3: Tensor Products

> **PREREQUISITES: M3.2 gate passed.**
> Core idea: multiple qubits → combine their spaces using ⊗.

```
M3.3.1  Kronecker Product — Combining Spaces
├── A⊗B: each element aᵢⱼ of A gets replaced by aᵢⱼ·B
│   (m×n)⊗(p×q) = (mp × nq) — dimensions multiply!
│
├── WORKED EXAMPLE — X⊗Z (2-qubit gate):
│   X⊗Z = [[0,1],[1,0]] ⊗ [[1,0],[0,-1]]
│        = [[0·Z, 1·Z],[1·Z, 0·Z]]
│        = [[0, 0, 1, 0],
│           [0, 0, 0,-1],
│           [1, 0, 0, 0],
│           [0,-1, 0, 0]]  (4×4 matrix)
│
├── Code:
│   X = np.array([[0,1],[1,0]], dtype=complex)
│   Z = np.array([[1,0],[0,-1]], dtype=complex)
│   XZ_tensor = np.kron(X, Z)
│   print(XZ_tensor.shape)  # (4, 4)
│
└── Why this matters: H₂ Hamiltonian = c₁(Z⊗Z) + c₂(X⊗X) + c₃(Y⊗Y) + ...
    Each term is a Kronecker product of single-qubit Paulis!

M3.3.2  Tensor Product of States — Multi-qubit
├── |ψ⟩⊗|φ⟩ = column vector of size (dim_ψ × dim_φ):
│   |0⟩⊗|0⟩ = [[1],[0]] ⊗ [[1],[0]] = [[1],[0],[0],[0]] = |00⟩
│   |0⟩⊗|1⟩ = [[1],[0]] ⊗ [[0],[1]] = [[0],[1],[0],[0]] = |01⟩
│   |1⟩⊗|0⟩ = [[0],[1]] ⊗ [[1],[0]] = [[0],[0],[1],[0]] = |10⟩
│   |1⟩⊗|1⟩ = [[0],[1]] ⊗ [[0],[1]] = [[0],[0],[0],[1]] = |11⟩
│
├── WORKED: compute |+⟩⊗|0⟩ manually:
│   |+⟩⊗|0⟩ = (1/√2)(|0⟩+|1⟩)⊗|0⟩ = (1/√2)(|00⟩+|10⟩)
│            = (1/√2)[1,0,1,0]ᵀ
│
├── Gate on qubit k: use I on all others
│   H on qubit 0, of 2-qubit system: H⊗I (4×4)
│   H on qubit 1:  I⊗H (4×4)
│   (H⊗I)|00⟩ = (1/√2)[1,0,1,0]ᵀ = |+0⟩ (only qubit 0 goes to |+⟩) ✓
│
├── n-qubit exponential scaling:
│   n=1: 2 amplitudes  (ℂ²)
│   n=2: 4 amplitudes  (ℂ⁴)
│   n=4: 16 amplitudes (H₂ molecule!)
│   n=50: 2⁵⁰ ≈ 10¹⁵ → classical computers cannot store this
│   This exponential is WHY quantum computers can solve problems classical can't!
│
├── Code:
│   ket0 = np.array([1,0], dtype=complex)
│   ket1 = np.array([0,1], dtype=complex)
│   ket_plus = (ket0+ket1)/np.sqrt(2)
│   print(np.kron(ket_plus, ket0))  # [0.707, 0, 0.707, 0] = |+0⟩
│   H = np.array([[1,1],[1,-1]], dtype=complex)/np.sqrt(2)
│   I = np.eye(2, dtype=complex)
│   HI = np.kron(H, I)   # H on qubit 0
│   ket00 = np.kron(ket0, ket0)
│   print(HI @ ket00)    # [0.707, 0, 0.707, 0] = |+0⟩ ✓
│
└── Exit check:
    Compute (H⊗H)|00⟩. Expected result: (1/2)[1,1,1,1]ᵀ = (1/2)(|00⟩+|01⟩+|10⟩+|11⟩).
    Verify in NumPy using np.kron(H,H) @ ket00.

═══════════════════════════════════════════
 GATE TO M4.x — Do NOT proceed until ALL boxes checked:
═══════════════════════════════════════════
 □ Can compute Kronecker product A⊗B for 2×2 matrices by hand
 □ Know |00⟩=[1,0,0,0]ᵀ, |01⟩=[0,1,0,0]ᵀ, |10⟩=[0,0,1,0]ᵀ, |11⟩=[0,0,0,1]ᵀ
 □ Can compute |+⟩⊗|0⟩ = (1/√2)[1,0,1,0]ᵀ by hand
 □ Know: applying gate to qubit k = I⊗...⊗G⊗...⊗I (I on all others)
 □ Verified (H⊗I)|00⟩ = |+0⟩ in NumPy
 □ Know: n qubits → 2ⁿ amplitudes → exponential quantum advantage
═══════════════════════════════════════════
```

---

## Module M4.1-M4.3: Hilbert Space, Probability & Summary

> **PREREQUISITES: All previous gates passed.**
> This module ties everything together — how it all connects to VQE.

```
M4.1  Hilbert Space — The Container of Quantum States
├── Formal definition: a COMPLETE inner product space
│   Practical (for finite QC): ℂⁿ with standard complex inner product
│   Qubit: ℂ² (2-dimensional Hilbert space)
│   2-qubit: ℂ⁴ (4-dimensional Hilbert space via tensor product)
│
├── Orthonormal basis {|e₁⟩,...,|eₙ⟩}: basis vectors where:
│   ⟨eᵢ|eⱼ⟩ = δᵢⱼ = 1 if i=j, else 0
│   ANY state |ψ⟩ = Σᵢ cᵢ|eᵢ⟩  where cᵢ = ⟨eᵢ|ψ⟩
│
├── Completeness relation (VERY useful in derivations):
│   Σᵢ |eᵢ⟩⟨eᵢ| = I  (resolution of identity)
│   Qubit: |0⟩⟨0| + |1⟩⟨1| = I  (verify: [[1,0],[0,0]]+[[0,0],[0,1]]=I ✓)
│   Also: |+⟩⟨+| + |-⟩⟨-| = I  (ANY orthonormal basis works!)
│
├── Uses of completeness (insert "1=I" anywhere in a calculation):
│   ⟨φ|ψ⟩ = ⟨φ|I|ψ⟩ = Σₙ⟨φ|n⟩⟨n|ψ⟩  (expand in basis)
│   ⟨ψ|A|ψ⟩ = Σₙ λₙ|⟨n|ψ⟩|²           (spectral theorem + completeness)
│
└── Code:
    ket0 = np.array([[1],[0]], dtype=complex)
    ket1 = np.array([[0],[1]], dtype=complex)
    P0 = ket0 @ ket0.conj().T  # |0⟩⟨0|
    P1 = ket1 @ ket1.conj().T  # |1⟩⟨1|
    print(np.allclose(P0 + P1, np.eye(2)))  # True ✓

M4.2  Probability & Measurement (Quantum Context)
├── Classical expectation: E[X] = Σ xᵢ P(xᵢ)
├── Quantum expectation value: ⟨A⟩ = ⟨ψ|A|ψ⟩
│   = Σₙ aₙ|⟨aₙ|ψ⟩|² (eigenvalue × Born probability)
│   = classical average of eigenvalues weighted by Born probabilities
│
├── Shot noise when measuring quantum hardware:
│   Each circuit run (shot) gives ONE of the eigenvalues (+1 or -1 for Paulis)
│   Estimate ⟨A⟩ from N shots: mean of outcomes
│   Statistical error σ = std/√N ≈ 1/√N
│   N=1024 → error ≈ 3%; N=8192 → error ≈ 1%
│
├── VQE measurement strategy:
│   H = Σₖ cₖ Pₖ  (Hamiltonian as Pauli sum)
│   ⟨H⟩ = Σₖ cₖ⟨Pₖ⟩  (linearity of expectation)
│   Measure ⟨Z⊗Z⟩: run circuit N shots, record ±1, average → ⟨Z⊗Z⟩
│   Measure ⟨X⊗X⟩: apply H⊗H first (basis rotation), then measure Z⊗Z
│   Total ⟨H⟩ = c₁⟨Z⊗Z⟩ + c₂⟨X⊗X⟩ + ... (linear combination)
│
└── Exit check:
    |ψ⟩ = |+⟩ = [1/√2, 1/√2]ᵀ.
    Compute ⟨Z⟩ analytically: (+1)P(|0⟩)+(-1)P(|1⟩) = (+1)(1/2)+(-1)(1/2) = 0.
    Compute ⟨X⟩: (+1)P(|+⟩)+(-1)P(|-⟩) = (+1)(1)+(-1)(0) = 1.
    Verify both with NumPy matrix sandwich.

M4.3  The Complete Math→Quantum Bridge
├── Everything you learned maps directly to VQE:
│   P0.1/M1.1  Complex numbers, e^(iθ)  → Quantum amplitudes, rotation gates
│   M1.2/M1.3  Derivatives, ODEs        → Schrödinger equation, optimization
│   M2.1       Vectors                  → Quantum states |ψ⟩
│   M2.2       Matrices                 → Quantum gates and operators
│   M2.3       Hermitian, Unitary       → Observables and gates
│   M3.1       Eigenvalues              → Energy levels, measurement outcomes
│   M3.2       Inner products           → Born rule, probabilities
│   M3.3       Tensor products          → Multi-qubit states
│   M4.1/M4.2  Hilbert space, stats     → Complete quantum framework
│
└── VQE in one line (you can now understand EVERY symbol):
    E₀ ≈ min_θ ⟨ψ(θ)|Ĥ|ψ(θ)⟩
    where Ĥ = Σₖ cₖ Pₖ (Paulis ← M2.3)
    and |ψ(θ)⟩ = Ry(θ₁)⊗Ry(θ₂)...|00...0⟩ (tensor product ← M3.3)
    and min by gradient descent: θₖ → θₖ - α·∂E/∂θₖ (partial deriv ← M1.2)

═══════════════════════════════════════════
 ⭐ MASTER MATH GATE — YOU ARE READY FOR PHYSICS PHASE ⭐
═══════════════════════════════════════════
 □ P0.1: Complex numbers, modulus, normalization
 □ M1.1: Polar form, Euler's formula, 8 key e^(iθ) values
 □ M1.2: Derivatives (5 rules), partial derivatives, integrals
 □ M1.3: ODE separation of variables, ψ(t)=ψ(0)e^(-iEt/ℏ)
 □ M2.1: Vectors, dot product, independence, complex inner product
 □ M2.2: Matrix-vector, matrix-matrix, dagger, det, trace
 □ M2.3: Hermitian (real eigenvalues), Unitary (preserves norm), Paulis
 □ M3.1: Eigenvalue equation Av=λv, characteristic equation, spectral theorem
 □ M3.2: Inner product, Born rule, Gram-Schmidt
 □ M3.3: Tensor/Kronecker product, multi-qubit states
 □ M4.x: Hilbert space, completeness, measurement strategy
 □ Can read: E₀ ≈ min_θ ⟨ψ(θ)|Ĥ|ψ(θ)⟩ and explain EVERY symbol
═══════════════════════════════════════════ 
```


---

# TO-DO LIST — PART 1 (MATH PHASE)
> Har topic complete karne ke baad check karo. Ek bhi miss mat karna.

## Phase 0 — Bridge

- [ ] P0.1.1  Number types: N⊂Z⊂Q⊂R⊂C — explain each
- [ ] P0.1.2  i=√-1; compute i², i³, i⁴
- [ ] P0.1.3  Add: (2+3i)+(1-4i) by hand
- [ ] P0.1.4  Multiply: (2+3i)(1-4i) by hand
- [ ] P0.1.5  Conjugate: z* = a-bi
- [ ] P0.1.6  Modulus: |z|=√(a²+b²) — derive from Pythagoras
- [ ] P0.1.7  Prove |z|² = z·z* by hand
- [ ] P0.1.8  Argand diagram: plot and label z=3+4i
- [ ] P0.1.9  Quantum normalization: verify |0.6|²+|0.8i|²=1
- [ ] P0.1.10 Phase: e^(iθ) rotates z, |z| unchanged
- [ ] P0.1 GATE — passed ✓

## Phase 1 — Calculus & ODEs

### M1.1 Complex Numbers Full Mastery
- [ ] M1.1.1  z=1+i → polar (r=√2, θ=π/4) by hand
- [ ] M1.1.2  z=-1+i → polar (θ=3π/4, 2nd quadrant) by hand
- [ ] M1.1.3  z=3+4i → polar (r=5) by hand
- [ ] M1.1.4  Convert 2·e^(iπ/3) back to rectangular (1+i√3)
- [ ] M1.1.5  Verify NumPy: np.isclose(z, r*np.exp(1j*theta))
- [ ] M1.1.6  Euler: e^(iθ)=cosθ+i·sinθ — write from memory
- [ ] M1.1.7  Taylor proof: plug ix into e^x, group even/odd powers
- [ ] M1.1.8  All 8 values: e^(i·0)=1, e^(iπ/6), e^(iπ/4), e^(iπ/3), e^(iπ/2)=i, e^(iπ)=-1, e^(i·3π/2)=-i, e^(i·2π)=1
- [ ] M1.1.9  Verify |e^(iθ)|=1 for all 8 values in NumPy
- [ ] M1.1.10 Multiplication = rotation: angles ADD, moduli MULTIPLY
- [ ] M1.1.11 4 rotations of π/2 return to start — code verify
- [ ] M1.1.12 De Moivre: (1+i)^8=16 by hand, verify Python
- [ ] M1.1.13 Rz(θ)=[[e^(-iθ/2),0],[0,e^(iθ/2)]] — polar form connection
- [ ] M1.1 GATE — passed ✓

### M1.2 Calculus
- [ ] M1.2.1  Derivative = instantaneous rate; f'(x)=lim[f(x+h)-f(x)]/h
- [ ] M1.2.2  Derive d(t²)/dt=2t from first principles
- [ ] M1.2.3  Power rule: d/dx[xⁿ]=n·xⁿ⁻¹ (x⁵, x⁻¹, √x)
- [ ] M1.2.4  Sum rule: differentiate x³+2x²-5x+7
- [ ] M1.2.5  Product rule: d/dx[x²·eˣ]=eˣ(x²+2x)
- [ ] M1.2.6  Chain rule: d/dx[e^(3x)], d/dx[sin(x²)], d/dx[(x²+1)⁵]
- [ ] M1.2.7  Full: differentiate x³sin(x)+e^(2x)
- [ ] M1.2.8  Partial ∂f/∂x for f=x²y+3y²+5x → 2xy+5
- [ ] M1.2.9  Partial ∂f/∂y same → x²+6y
- [ ] M1.2.10 Gradient ∇f=[17,22] at (2,3); gradient descent update rule
- [ ] M1.2.11 Antiderivatives: ∫xⁿdx, ∫eˣdx, ∫sin(x)dx, ∫cos(x)dx
- [ ] M1.2.12 Definite integral: ∫(3x²+2)dx from 0→1 = 3
- [ ] M1.2.13 Verify ∫₀¹ 2sin²(πx)dx=1 (normalization check, use identity)
- [ ] M1.2 GATE — passed ✓

### M1.3 ODEs
- [ ] M1.3.1  ODE definition: equation involving f and f' together
- [ ] M1.3.2  Solve dy/dt=ky by separation of variables → y=y₀e^(kt)
- [ ] M1.3.3  Bacteria: P(t)=1000e^(0.5t), compute P(3)≈4482
- [ ] M1.3.4  Radioactive: N(t)=10000e^(-0.1t), half-life≈6.93
- [ ] M1.3.5  Solve dy/dt=-3y, y(0)=2 → 2e^(-3t); verify scipy
- [ ] M1.3.6  Quantum ODE: dψ/dt=(-iE/ℏ)ψ — imaginary k case
- [ ] M1.3.7  Solution: ψ(t)=ψ(0)·e^(-iEt/ℏ) — verify by substitution
- [ ] M1.3.8  |e^(-iEt/ℏ)|=1 always → probability conserved
- [ ] M1.3.9  Derive TISE (Ĥφ=Eφ) by separation of variables from TDSE
- [ ] M1.3 GATE — passed ✓

## Phase 2 — Linear Algebra

### M2.1 Vectors
- [ ] M2.1.1  Vector = ordered list; column notation [v₁,v₂]ᵀ
- [ ] M2.1.2  Add [3,4]+[1,-2]=[4,2]; scale 3·[2,-1]=[6,-3]
- [ ] M2.1.3  |0⟩=[1,0]ᵀ, |1⟩=[0,1]ᵀ, |+⟩=[1/√2,1/√2]ᵀ, |-⟩=[1/√2,-1/√2]ᵀ
- [ ] M2.1.4  Dot product: [3,4]·[1,2]=11; orthogonal when dot=0
- [ ] M2.1.5  Verify ⟨0|1⟩=0 by hand, AND NumPy
- [ ] M2.1.6  Span, linear independence, basis — 2D examples
- [ ] M2.1.7  rank([1,2,3],[4,5,6],[7,8,9])=2 → dependent (code verify)
- [ ] M2.1.8  Complex inner product: ⟨u,v⟩=Σuᵢ*vᵢ (conjugate FIRST!)
- [ ] M2.1.9  Compute ⟨u|v⟩ for u=[1,i]ᵀ, v=[1,1]ᵀ → 1-i
- [ ] M2.1.10 Verify ⟨+|-⟩=0 by hand and NumPy
- [ ] M2.1.11 n qubits → 2ⁿ dimensions; H₂ has 4 qubits → 16-dim
- [ ] M2.1 GATE — passed ✓

### M2.2 Matrices
- [ ] M2.2.1  Mv product: each row of M dotted with v
- [ ] M2.2.2  X|0⟩=[0,1]ᵀ (NOT gate) by hand
- [ ] M2.2.3  H|0⟩=|+⟩, H|1⟩=|-⟩ by hand
- [ ] M2.2.4  Z|1⟩=-|1⟩ (phase flip) by hand
- [ ] M2.2.5  AB≠BA: compute XZ and ZX, confirm different
- [ ] M2.2.6  HXH=Z by hand (two matrix multiplications)
- [ ] M2.2.7  Transpose Aᵀ, conjugate A*, dagger A†=(A*)ᵀ
- [ ] M2.2.8  Compute A† for [[1+i,2],[3,4-2i]] step by step
- [ ] M2.2.9  (AB)†=B†A†; bra ⟨ψ|=|ψ⟩†
- [ ] M2.2.10 det([[a,b],[c,d]])=ad-bc; det(X)=-1, det(Z)=-1
- [ ] M2.2.11 Tr diagonal sum; Tr(X)=Tr(Y)=Tr(Z)=0, Tr(I)=2
- [ ] M2.2.12 U⁻¹=U† for quantum gates (reversible via dagger)
- [ ] M2.2 GATE — passed ✓

### M2.3 Special Matrices
- [ ] M2.3.1  Hermitian: A=A†; proof → all eigenvalues real
- [ ] M2.3.2  Verify Z, X, Y, I all Hermitian in Python
- [ ] M2.3.3  Unitary: U†U=I; ||Uv||=||v||; verify HH†=I by hand
- [ ] M2.3.4  All 4 Paulis from memory in 60 sec
- [ ] M2.3.5  Verify XY=iZ, XYZ=iI in NumPy
- [ ] M2.3.6  Pauli decomp: cₖ=Tr(σₖM)/2 — apply to [[3,1+i],[1-i,2]]
- [ ] M2.3 GATE — passed ✓

## Phase 3 — Advanced

### M3.1 Eigenvalues & Eigenvectors
- [ ] M3.1.1  Av=λv; eigenvector = special direction (no rotation)
- [ ] M3.1.2  Z|0⟩=+1·|0⟩, Z|1⟩=-1·|1⟩ — verify
- [ ] M3.1.3  det(A-λI)=0; compute eigenvalues of Z → ±1
- [ ] M3.1.4  Compute eigenvalues of 0.5Z+0.5X by hand → ±1/√2
- [ ] M3.1.5  Spectral: Z=(+1)|0⟩⟨0|+(-1)|1⟩⟨1| — verify NumPy
- [ ] M3.1.6  Variational: ⟨ψ|H|ψ⟩≥E₀ for all |ψ⟩ (prove briefly)
- [ ] M3.1 GATE — passed ✓

### M3.2 Inner Products
- [ ] M3.2.1  ⟨u,v⟩=Σuᵢ*vᵢ; norm ||v||; orthogonal when 0
- [ ] M3.2.2  ⟨0|0⟩=1, ⟨0|1⟩=0, ⟨+|+⟩=1, ⟨+|-⟩=0 by hand
- [ ] M3.2.3  Born rule P(|φ⟩)=|⟨φ|ψ⟩|²
- [ ] M3.2.4  P(|0⟩)=3/4, P(|1⟩)=1/4 for |ψ⟩=(√3/2)|0⟩+(1/2)|1⟩
- [ ] M3.2.5  Complex amplitude: P(|0⟩)=P(|1⟩)=1/2 for (1/√2)|0⟩+(i/√2)|1⟩
- [ ] M3.2.6  Gram-Schmidt: {[1,1],[1,0]} → {|+⟩,|-⟩} Hadamard basis
- [ ] M3.2 GATE — passed ✓

### M3.3 Tensor Products
- [ ] M3.3.1  A⊗B: replace aᵢⱼ with aᵢⱼ·B block
- [ ] M3.3.2  Compute X⊗Z by hand (4×4 matrix)
- [ ] M3.3.3  |00⟩, |01⟩, |10⟩, |11⟩ as 4-element column vectors
- [ ] M3.3.4  |+⟩⊗|0⟩=(1/√2)[1,0,1,0]ᵀ by hand
- [ ] M3.3.5  Gate on qubit 0: H⊗I; gate on qubit 1: I⊗H
- [ ] M3.3.6  (H⊗H)|00⟩=(1/2)[1,1,1,1]ᵀ — NumPy verify
- [ ] M3.3.7  H₂ Hamiltonian = Σcₖ(Pauli⊗Pauli) strings
- [ ] M3.3 GATE — passed ✓

## Phase 4 — Framework

### M4.x Hilbert Space & VQE Bridge
- [ ] M4.1.1  Hilbert space = ℂⁿ with inner product
- [ ] M4.1.2  Completeness: |0⟩⟨0|+|1⟩⟨1|=I; also |+⟩⟨+|+|-⟩⟨-|=I
- [ ] M4.1.3  Verify both in NumPy
- [ ] M4.2.1  ⟨A⟩=⟨ψ|A|ψ⟩ — VQE cost function
- [ ] M4.2.2  ⟨Z⟩=0, ⟨X⟩=1 for |+⟩ — computed analytically
- [ ] M4.2.3  Shot noise: σ≈1/√N; 1024 shots → ~3% error
- [ ] M4.2.4  Measure ⟨X⊗X⟩: rotate basis with H⊗H first, then Z-measure
- [ ] M4.3.1  Read E₀≈min_θ⟨ψ(θ)|Ĥ|ψ(θ)⟩ — explain every symbol
- [ ] M4.3.2  Map each math module to its quantum concept

---

## ⭐ MASTER SIGN-OFF — PART 1

- [ ] All 10 module gates passed
- [ ] All 4 Paulis written from memory in 60 sec
- [ ] Eigenvalues of 0.5Z+0.5X = ±0.707 by hand
- [ ] VQE equation explained symbol by symbol
- [ ] **READY FOR PART 2 — PHYSICS PHASE 🚀**

