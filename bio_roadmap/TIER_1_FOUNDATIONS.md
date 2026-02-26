# 🔴 TIER 1: ABSOLUTE FOUNDATIONS
## DNA → Central Dogma → Amino Acids & Biochemistry → Gene Regulation → Cell Biology

> **PRIORITY:** CRITICAL — Learn FIRST, everything else depends on this
> **TIME:** Month 1-3
> **PREREQUISITE:** Basic Python, school-level biology (cells, DNA exists)

---

## 1.1 DNA — THE CODE OF LIFE
**Why:** Every bioinformatics tool processes DNA. You must understand the molecule itself.

```
A. STRUCTURE:
├── Double helix — Watson & Crick 1953
├── Nucleotides = Phosphate + Deoxyribose sugar + Nitrogenous Base
│   ├── Purine bases (2 rings): Adenine (A), Guanine (G)
│   └── Pyrimidine bases (1 ring): Thymine (T), Cytosine (C)
├── Base pairing rules: A-T (2 H-bonds), G-C (3 H-bonds)
│   └── G-C pairs are STRONGER → high-GC regions are more stable (harder to denature)
├── Antiparallel strands: one runs 5'→3', other runs 3'→5'
│   └── DNA polymerase ONLY reads 3'→5', ONLY synthesizes 5'→3'
├── Phosphodiester bonds: link nucleotides (sugar-phosphate backbone)
├── Major groove & Minor groove → proteins bind here to read DNA
│   └── TFs, Cas9 use major groove for sequence recognition
├── Supercoiling: DNA wound around histones → nucleosome → chromatin → chromosome
│   ├── Topoisomerases: enzymes that relieve supercoiling tension
│   └── Euchromatin (loose, active) vs Heterochromatin (tight, silent)
└── Chromosomes: 46 in humans (23 pairs = 22 autosomes + XX/XY)

B. DNA IN CONTEXT:
├── Nuclear DNA: ~3.2 billion base pairs (haploid genome = 1 copy)
├── Mitochondrial DNA: 16,569 bp, circular, maternally inherited, 37 genes
├── Only ~1.5% is protein-coding! (~20,000 genes)
├── Non-coding DNA (the other 98.5%):
│   ├── Regulatory sequences (enhancers, silencers, promoters)
│   ├── Introns (removed during mRNA processing)
│   ├── Non-coding RNAs (miRNA, lncRNA, piRNA)
│   ├── Repetitive elements (SINEs, LINEs, ALU elements)
│   ├── Transposons ("jumping genes" — activated in aging!)
│   └── Telomeric DNA (TTAGGG repeats at chromosome ends)
├── Genome size ≠ complexity (onion genome > human genome)
└── CpG dinucleotides: where methylation happens → aging clocks!

C. WHAT YOU DON'T NEED (skip):
├── Detailed X-ray crystallography history
└── Chargaff's derivation — just know A=T%, G=C%

🔗 BIOINFORMATICS CONNECTION:
├── FASTA/FASTQ files = text files of A/T/G/C strings (you code this!)
├── Complement: A↔T, G↔C → reverse complement matters in alignment
├── GC content: high GC = stable, low GC = easy to denature (PCR design!)
├── Restriction sites: palindromic sequences → CRISPR PAM sites are similar
└── k-mer counting: break DNA into short words → genome assembly algorithms
```

---

## 1.2 THE CENTRAL DOGMA — DNA → RNA → Protein
**Why:** Gene expression is what aging disrupts. You need this to understand what goes wrong.

```
THE FLOW:
DNA → [Transcription] → pre-mRNA → [Processing] → mRNA → [Translation] → Protein

═══════════════════════════════════════════
A. TRANSCRIPTION (DNA → mRNA):
═══════════════════════════════════════════

Machinery:
├── RNA Polymerase II: enzyme that reads DNA template, builds mRNA
├── Promoter: DNA region where Pol II attaches (contains TATA box ~25bp upstream)
├── Template strand: read 3'→5' by Pol II
├── Coding strand: same sequence as mRNA (but T instead of U)
└── RNA uses: A, U (Uracil replaces Thymine), G, C

Steps:
├── 1. Initiation: TFs + Pol II assemble at promoter → transcription bubble forms
├── 2. Elongation: Pol II moves along template, adds ribonucleotides (5'→3')
│   └── Speed: ~40 nucleotides/second in eukaryotes
├── 3. Termination: Pol II reaches termination signal → releases mRNA
└── Result: pre-mRNA (still has introns)

mRNA Processing (eukaryotes only — CRITICAL):
├── 5' cap: 7-methylguanosine added → protects from degradation, ribosome recognition
├── 3' Poly-A tail: ~200 adenines added → stability, nuclear export
├── Splicing: INTRONS removed, EXONS joined by spliceosome
│   ├── Introns: non-coding sequences WITHIN gene (can be huge, >100kb!)
│   ├── Exons: coding sequences that remain in mature mRNA
│   ├── Alternative splicing: same gene → different exon combinations → different proteins!
│   │   └── Human: ~20,000 genes but >100,000 protein variants via alt splicing!
│   └── 🔗 WHY cDNA exists: cDNA = mRNA reverse-transcribed back to DNA (no introns!)
│       RNA-seq libraries use cDNA, not genomic DNA
└── Mature mRNA: 5'cap + 5'UTR + CDS (coding) + 3'UTR + polyA tail

🔗 BIOINFORMATICS:
├── RNA-seq = measuring transcription levels genome-wide
├── STAR aligner: maps RNA-seq reads to genome (handles intron-spanning reads!)
├── BWA: maps DNA reads (does NOT handle splicing → wrong tool for RNA-seq)
└── GTF/GFF files: contain exon/intron boundaries for every gene

═══════════════════════════════════════════
B. TRANSLATION (mRNA → Protein):
═══════════════════════════════════════════

Machinery:
├── Ribosome: 80S in eukaryotes (60S large + 40S small subunits)
│   └── Made of rRNA + proteins → assembled in nucleolus
├── mRNA: the message (codons)
├── tRNA: adapter molecules, carry amino acids
│   ├── Anticodon loop: base-pairs with mRNA codon
│   └── 3' end: amino acid attachment site (CCA sequence)
└── Aminoacyl-tRNA synthetases: 20 enzymes, each loads correct AA onto correct tRNA

The Genetic Code:
├── 3 nucleotides = 1 codon = 1 amino acid
├── 64 possible codons → 20 amino acids + 3 stop signals
├── AUG = Methionine = START codon (every protein starts with Met)
├── UAA, UAG, UGA = STOP codons (no amino acid loaded)
├── REDUNDANCY: multiple codons → same AA (e.g., Leu has 6 codons!)
│   └── "Wobble position": 3rd base of codon is flexible
├── Near-universal: same code in almost all life (evidence of common origin)
└── Exceptions: mitochondrial code slightly different (UGA = Trp, not Stop)

Steps:
├── 1. Initiation: 40S binds 5'cap, scans for AUG → 60S joins → 80S ribosome
├── 2. Elongation: tRNA brings AA → peptide bond forms → ribosome moves 1 codon
│   └── Amino acids linked: N-terminus (first) → C-terminus (last)
├── 3. Termination: Stop codon → release factor enters → protein released
└── Post-translational: protein folds, gets modified, goes to correct location

🔗 BIOINFORMATICS:
├── ORF finding: scan for AUG...STOP in all 6 reading frames (3 forward + 3 reverse)
├── Codon usage bias: organisms prefer certain codons → codon optimization for synthetic genes
├── Rosalind problems: "Translating RNA into Protein" — you already coded this (BS7_PROT!)
└── Ribosome profiling (Ribo-seq): measure which mRNAs are actively being translated
```

---

## 1.3 AMINO ACIDS — THE BUILDING BLOCKS OF PROTEINS
**Why:** Proteins do everything in your body. Every drug, every aging mechanism targets specific amino acids.

```
⚠️ GAP FIX: This was completely missing from the original roadmap.

THE 20 STANDARD AMINO ACIDS:

Classified by side chain (R-group) properties:

NONPOLAR (Hydrophobic) — hide INSIDE protein core:
├── Glycine     (Gly, G) — Smallest AA, only H as side chain, very flexible
├── Alanine     (Ala, A) — Small methyl group (-CH₃)
├── Valine      (Val, V) — Branched chain ← essential AA (must eat!)
├── Leucine     (Leu, L) — Branched chain ← essential AA
├── Isoleucine  (Ile, I) — Branched chain ← essential AA
├── Proline     (Pro, P) — RING structure → forces KINK in protein chain!
│   └── Found in collagen (age-related decline → wrinkles!)
├── Phenylalanine (Phe, F) — Aromatic ring
├── Tryptophan  (Trp, W) — Largest AA, double ring, rare in proteins
├── Methionine  (Met, M) — Contains sulfur, START codon (AUG)
│   └── Methionine restriction extends lifespan in model organisms!
└── WHY MATTERS: Hydrophobic collapse = main driving force of protein folding!

POLAR UNCHARGED (Hydrophilic) — on protein SURFACE:
├── Serine      (Ser, S) — -OH group → PHOSPHORYLATION target!
├── Threonine   (Thr, T) — -OH group → PHOSPHORYLATION target! Essential AA
├── Tyrosine    (Tyr, Y) — -OH + aromatic → PHOSPHORYLATION target!
│   └── Kinases add phosphate to S/T/Y → signal transduction!
├── Asparagine  (Asn, N) — Amide group, glycosylation site
├── Glutamine   (Gln, Q) — Amide group, important in metabolism
└── Cysteine    (Cys, C) — -SH thiol group → forms DISULFIDE BONDS (S-S)!
    └── 2 Cys residues → S-S bridge → stabilizes protein 3D structure

POSITIVELY CHARGED (+) at pH 7:
├── Lysine      (Lys, K) — -NH₃⁺ → ACETYLATION, METHYLATION, UBIQUITINATION target!
│   └── Histone Lys modifications = epigenetics! (H3K4me3, H3K27ac, etc.)
├── Arginine    (Arg, R) — Guanidinium group (+), strong positive charge
└── Histidine   (His, H) — Imidazole ring, pKa ~6 → charge changes near pH 7!
    └── Acts as pH sensor in proteins, found in enzyme active sites

NEGATIVELY CHARGED (-) at pH 7:
├── Aspartate   (Asp, D) — -COO⁻, short side chain
└── Glutamate   (Glu, E) — -COO⁻, longer side chain
    └── Glutamate = major neurotransmitter precursor!

MEMORIZE THE 1-LETTER CODES (you use them in every Rosalind problem):
G A V L I P F W M — nonpolar
S T Y N Q C       — polar uncharged
K R H             — positive
D E               — negative

POST-TRANSLATIONAL MODIFICATIONS (PTMs):
├── Phosphorylation: kinase adds PO₄³⁻ to Ser/Thr/Tyr → ON/OFF switch for signaling
│   └── ~30% of all proteins are phosphorylated at some point!
├── Ubiquitination: ubiquitin (76 AA protein) attached to Lys → targets for degradation
│   └── Proteasome degrades ubiquitinated proteins → proteostasis!
├── Acetylation: acetyl group on Lys → neutralizes positive charge
│   └── Histone acetylation → gene activation (H3K27ac)
│   └── Sirtuins = deacetylases → aging pathway!
├── Methylation: methyl groups on Lys/Arg → epigenetic marks
├── Glycosylation: sugar chains attached to Asn/Ser/Thr → cell surface markers
└── SUMOylation: small ubiquitin-like modifier → nuclear transport, DNA repair

🔗 BIOINFORMATICS CONNECTION:
├── BLOSUM62, PAM250 matrices: amino acid substitution scores for alignment
│   └── Based on how often one AA replaces another in evolution
├── Protein language models (ESM-2) embed amino acid sequences
├── PTM prediction: ML models predict phosphorylation sites from sequence
├── PDB files: 3D coordinates of ATOMS in each amino acid
└── Drug design: drugs bind to specific AA residues in protein active sites
```

---

## 1.4 BIOCHEMISTRY ESSENTIALS — THE CHEMISTRY OF LIFE
**Why:** Drug design, metabolic aging, enzyme targets — all require biochemistry.

```
⚠️ GAP FIX: This was completely missing from the original roadmap.

═══════════════════════════════════════════
A. WATER, pH, BUFFERS — THE ENVIRONMENT:
═══════════════════════════════════════════
├── Water: polar solvent, H-bonds → dissolves ions & polar molecules
├── pH = -log₁₀[H⁺]: measures acidity
│   ├── pH 7 = neutral (blood pH 7.35-7.45)
│   ├── pH < 7 = acidic (stomach = pH 2)
│   └── pH > 7 = basic (intestine = pH 8)
├── pKa: pH at which an amino acid side chain is 50% protonated
│   └── Histidine pKa ~6 → charge changes near physiological pH!
├── Buffers: resist pH change (phosphate buffer, bicarbonate buffer in blood)
└── WHY: Enzymes ONLY work in narrow pH range → lab conditions must be precise!
    CRISPR reactions done in specific buffers!

═══════════════════════════════════════════
B. THERMODYNAMICS — WHEN REACTIONS HAPPEN:
═══════════════════════════════════════════
├── Gibbs Free Energy: ΔG = ΔH - TΔS
│   ├── ΔG < 0: spontaneous (exergonic) — reaction proceeds forward
│   ├── ΔG > 0: non-spontaneous (endergonic) — needs energy input
│   └── ΔG = 0: equilibrium
├── ATP hydrolysis: ATP → ADP + Pᵢ, ΔG = -30.5 kJ/mol
│   └── ATP = universal energy currency → drives non-spontaneous reactions
├── Coupled reactions: unfavorable reaction + ATP hydrolysis = favorable overall
└── WHY: Protein folding is driven by ΔG (hydrophobic collapse lowers free energy)
    Drug binding to target protein = governed by ΔG of binding

═══════════════════════════════════════════
C. ENZYME KINETICS — HOW ENZYMES WORK:
═══════════════════════════════════════════
├── Enzymes: biological catalysts → speed up reactions without being consumed
├── Active site: 3D pocket where substrate binds → specificity!
│   └── Lock-and-key model vs Induced fit model
├── Michaelis-Menten equation: v = Vmax·[S] / (Km + [S])
│   ├── v: reaction velocity (rate)
│   ├── Vmax: maximum rate (when all enzyme is saturated)
│   ├── Km: substrate concentration at half-Vmax
│   │   └── Low Km = high affinity (enzyme grabs substrate easily)
│   │   └── High Km = low affinity
│   └── Kcat (turnover number): how many substrates per enzyme per second
│       └── Kcat/Km = "catalytic efficiency" → perfect enzyme ~10⁸ M⁻¹s⁻¹
│
├── Enzyme Inhibition:
│   ├── Competitive: inhibitor binds active site (competes with substrate)
│   │   └── Increases apparent Km, Vmax unchanged
│   ├── Non-competitive (Allosteric): binds elsewhere → changes enzyme shape
│   │   └── Km unchanged, decreases Vmax
│   └── Uncompetitive: binds enzyme-substrate complex
│
├── AGING-RELEVANT ENZYMES:
│   ├── Telomerase: TERT (protein) + TERC (RNA) → extends telomeres
│   ├── Caspases: executioner proteases → apoptosis
│   ├── Kinases (CDK, mTOR, AMPK): phosphorylate targets → signaling
│   ├── Sirtuins (SIRT1-7): NAD⁺-dependent deacetylases
│   └── Proteasome: degrades ubiquitinated proteins → declines 50% with age!
│
└── WHY: Drug design = designing enzyme INHIBITORS
    IC50 = concentration at which enzyme is 50% inhibited → key drug metric
    Rapamycin is an ALLOSTERIC inhibitor of mTOR → extends lifespan!

═══════════════════════════════════════════
D. METABOLIC PATHWAYS — ENERGY OF LIFE:
═══════════════════════════════════════════
├── Glycolysis: Glucose → 2 Pyruvate + 2 ATP + 2 NADH (cytoplasm)
│   └── 10 enzymatic steps, ancient pathway (all life has it)
├── TCA Cycle (Krebs): Acetyl-CoA → CO₂ + 3 NADH + FADH₂ + GTP (mitochondria matrix)
│   └── 8 enzymatic steps, circular → complete oxidation of carbon
├── Oxidative Phosphorylation (ETC): NADH/FADH₂ → ~34 ATP (inner mito membrane)
│   ├── Complex I: NADH → ubiquinone (4 H⁺ pumped)
│   ├── Complex II: FADH₂ → ubiquinone (0 H⁺ pumped, enters at lower energy)
│   ├── Complex III: ubiquinone → cytochrome c (4 H⁺ pumped)
│   ├── Complex IV: cytochrome c → O₂ → H₂O (2 H⁺ pumped)
│   └── ATP Synthase (Complex V): H⁺ gradient drives ATP synthesis (rotary motor!)
│       └── ~3 H⁺ per ATP → total ~34 ATP per glucose
│
├── ROS (Reactive Oxygen Species):
│   ├── Leaked electrons at Complex I and III → partial reduction of O₂
│   ├── Superoxide (O₂⁻), H₂O₂, hydroxyl radical (·OH)
│   ├── Damage: DNA mutations, protein oxidation, lipid peroxidation
│   ├── Antioxidant defense: SOD, catalase, glutathione, NRF2 pathway
│   └── AGING: ROS production increases with age → mitochondrial damage → more ROS → vicious cycle!
│
├── NAD⁺/NADH ratio:
│   ├── NAD⁺ = oxidized form (electron acceptor)
│   ├── NADH = reduced form (electron carrier → gives to ETC)
│   ├── NAD⁺ declines 50% from age 20→80!
│   ├── Sirtuins NEED NAD⁺ to function → less NAD⁺ = less sirtuin activity = aging!
│   └── INTERVENTION: NMN, NR (NAD⁺ precursors) → restore levels
│
├── Other key metabolites:
│   ├── Acetyl-CoA: fuel for TCA, building block for fatty acids, acetylation donor!
│   ├── α-ketoglutarate: TCA intermediate → also cofactor for TET enzymes (demethylation!)
│   │   └── Supplementation extends lifespan in worms!
│   ├── S-adenosylmethionine (SAM): methyl group donor for DNA methylation
│   └── Glutathione: major antioxidant, declines with age
│
└── 🔗 CONNECTION:
    ├── Metabolomics (Module 13 DEEP SYLLABUS): measure all metabolites
    ├── NAD⁺/NADH = central to aging interventions
    ├── Quantum computing: simulating ETC complex I exactly = key future application
    └── Caloric restriction → less glucose → less mTOR → more autophagy → longer life!
```

---

## 1.5 GENE EXPRESSION & REGULATION — THE MASTER SWITCH
**Why:** Aging = gene expression gone wrong. Reprogramming = resetting gene expression.

```
HOW GENES GET TURNED ON/OFF:

═══════════════════════════════════════════
A. PROMOTERS & ENHANCERS:
═══════════════════════════════════════════
├── Promoter: region upstream of gene where RNA Pol II binds
│   ├── Core promoter: TATA box (~-25bp), INR (+1), DPE
│   └── Proximal promoter: -100 to -250bp, TF binding sites
├── Enhancer: can be 10-1000 kb away, loops via chromatin to contact promoter
│   └── How? Mediator complex + cohesin facilitate enhancer-promoter loop
├── Silencer: opposite of enhancer, recruits repressors
├── Insulator: prevents enhancer from activating wrong gene
│   └── CTCF protein = main insulator binding factor
└── CpG islands: CG-rich regions near promoters
    ├── Normally UNmethylated → gene active
    └── Methylated → gene silenced → changes with aging!

═══════════════════════════════════════════
B. TRANSCRIPTION FACTORS (TFs) — THE COMMANDERS:
═══════════════════════════════════════════
├── Proteins that bind specific DNA sequences to activate/repress genes
├── DNA-binding domains:
│   ├── Zinc finger: Cys₂His₂ or Cys₄ → each finger reads 3bp
│   ├── Helix-turn-helix (HTH): one helix fits major groove
│   ├── bHLH (basic helix-loop-helix): dimerize → bind E-box (CANNTG)
│   ├── bZIP (leucine zipper): dimerize → bind DNA
│   └── Homeodomain: 60aa domain, development genes
│
├── CRITICAL TFs FOR AGING:
│   ├── NF-κB → inflammation → inflammaging (Hallmark #11)
│   ├── p53 → DNA damage response → senescence (Hallmark #8)
│   ├── FOXO3 → longevity! (highest in centenarians)
│   │   └── Activated by AMPK, inhibited by AKT/insulin
│   ├── NRF2 → antioxidant response → protects mitochondria
│   ├── TFEB → autophagy activation → lysosomal biogenesis
│   ├── HIF-1α → hypoxia response → angiogenesis
│   └── HSF1 → heat shock response → chaperone expression
│
├── YAMANAKA FACTORS (Oct4, Sox2, Klf4, c-Myc = OSKM):
│   ├── These are TFs that reprogram differentiated cells → iPSCs
│   ├── Partial expression → age reversal without losing cell identity
│   ├── This is THE mechanism Altos Labs, Retro Biosciences are using!
│   └── Retro Biosciences + OpenAI → 50x reprogramming efficiency with AI!
│
└── Coactivators & Corepressors:
    ├── CBP/p300 → acetyltransferases recruited by TFs → open chromatin
    ├── NCoR/SMRT → recruited by repressors → close chromatin
    └── Mediator complex: bridge between TFs and RNA Pol II

═══════════════════════════════════════════
C. NON-CODING RNAs (ncRNAs):
═══════════════════════════════════════════
├── microRNA (miRNA): ~22nt
│   ├── Mature miRNA loaded into RISC complex (Argonaute protein)
│   ├── Binds 3'UTR of target mRNA → degrades it or blocks translation
│   ├── One miRNA can target ~200 different mRNAs!
│   ├── miR-21: oncomiR (promotes cancer)
│   ├── miR-34a: p53 target, tumor suppressor
│   └── miR-29: increases with age → targets collagen genes → skin aging!
│
├── lncRNA (long non-coding RNA): >200nt
│   ├── XIST: X-chromosome inactivation
│   ├── HOTAIR: represses HOX genes, misregulated in cancer
│   ├── TERRA: telomeric RNA, regulates telomere length
│   └── Many aging-associated lncRNAs discovered in last 5 years
│
├── siRNA (small interfering RNA): ~21nt, double-stranded
│   ├── Used in RNAi gene silencing therapeutics!
│   ├── ONPATTRO (patisiran): first FDA-approved siRNA drug (2018)
│   └── Mechanism: Dicer cleaves dsRNA → RISC loads → mRNA degradation
│
├── piRNA: germline-specific
│   └── Protect genome from transposons (jumping genes)
│
└── circRNA (circular RNA): newly discovered regulatory RNAs
    └── miRNA sponges — sequester miRNAs away from targets

🔗 BIOINFORMATICS CONNECTION:
├── ChIP-seq: where TFs bind across genome → TF binding analysis
├── ATAC-seq: which DNA regions are open/accessible
├── miRNA-seq: measure microRNA expression
├── Motif analysis: find TF binding site sequences (MEME, HOMER tools)
└── Gene regulatory network inference: build TF→target gene maps
```

---

## 1.6 CELL BIOLOGY — WHAT YOU ACTUALLY NEED TO KNOW
**Why:** You're engineering cells. Know your machine.

```
═══════════════════════════════════════════
A. CELL MEMBRANE & TRANSPORT:
═══════════════════════════════════════════
⚠️ GAP FIX: Membrane was missing, needed for gene delivery understanding.

├── Phospholipid bilayer:
│   ├── Hydrophilic heads face water (outside & cytoplasm)
│   ├── Hydrophobic tails face each other (interior of membrane)
│   ├── Fluid mosaic model: membrane proteins float in lipid sea
│   └── Cholesterol: stiffens membrane, controls fluidity
│
├── Membrane proteins:
│   ├── Integral (transmembrane): span entire bilayer (channels, receptors)
│   ├── Peripheral: attached to one surface (signaling scaffolds)
│   └── GPI-anchored: lipid-linked to outer surface
│
├── Transport types:
│   ├── Passive (no energy):
│   │   ├── Simple diffusion: O₂, CO₂ cross freely (small, nonpolar)
│   │   ├── Facilitated diffusion: channels/carriers (glucose via GLUT)
│   │   └── Osmosis: water through aquaporins
│   ├── Active (needs ATP):
│   │   ├── Na⁺/K⁺ ATPase: 3 Na⁺ out, 2 K⁺ in → membrane potential
│   │   └── ABC transporters: drug efflux (cancer drug resistance!)
│   └── Vesicular transport:
│       ├── Endocytosis: cell takes in material (phagocytosis, pinocytosis)
│       │   └── Receptor-mediated endocytosis → THIS IS HOW LNPs ENTER CELLS!
│       └── Exocytosis: cell releases material (secretion, SASP factors!)
│
└── WHY: Gene delivery vectors must cross membranes!
    ├── LNPs: ionizable lipid fuses with endosomal membrane → cargo escapes
    ├── AAVs: enter via receptor-mediated endocytosis → escape endosome
    └── Electroporation: electric pulse makes temporary pores in membrane

═══════════════════════════════════════════
B. ORGANELLES — AGING RELEVANT:
═══════════════════════════════════════════

Mitochondria (🔥 MOST IMPORTANT):
├── Double membrane: outer (permeable), inner (cristae → ETC lives here)
├── Matrix: TCA cycle enzymes, mtDNA, mitochondrial ribosomes
├── Own genome: 16,569bp circular DNA, 37 genes, 13 ETC subunits
├── Fission (Drp1) & Fusion (Mfn1/2, OPA1): quality control balance
├── Mitophagy: PINK1 accumulates on damaged mito → recruits Parkin → autophagy
└── AGING: dysfunction = Hallmark #7 (ROS↑, ATP↓, mtDNA mutations↑)

Nucleus:
├── Double membrane with nuclear pores (~2000 pores per nucleus!)
├── Nuclear pore complex: controls import/export (NLS/NES signals on proteins)
├── Nucleolus: rRNA transcription + ribosome assembly
├── Lamins (A/B): nuclear skeleton → LAMINOPATHIES = premature aging!
│   └── Progeria (Hutchinson-Gilford): mutation in Lamin A → child ages rapidly
└── Chromatin: DNA + histones → levels of compaction

Lysosomes:
├── pH ~4.5, filled with >60 hydrolytic enzymes (acid hydrolases)
├── Autophagosome fuses with lysosome → degrades cellular waste
├── Lipofuscin: undegradable waste accumulates in lysosomes with age!
└── Lysosomal storage diseases: Gaucher, Tay-Sachs → enzyme deficiency

Endoplasmic Reticulum:
├── Rough ER: ribosomes → secretory/membrane protein synthesis + folding
├── Smooth ER: lipid synthesis, Ca²⁺ storage, detoxification
├── ER stress → Unfolded Protein Response (UPR):
│   ├── IRE1, PERK, ATF6 sensors detect misfolded proteins
│   └── If stress too high → apoptosis (death signal from ER!)
└── UPR activation increases with age → chronic ER stress!

Golgi Apparatus:
├── cis (receiving) → medial → trans (shipping) faces
├── Sorts, modifies (glycosylation), packages proteins for secretion
└── SASP factors from senescent cells processed through Golgi!

Proteasome:
├── 26S proteasome: barrel-shaped protein degradation machine
├── Recognizes ubiquitin tags → unfolds protein → cleaves into peptides
├── Activity DECLINES 50% with age → proteotoxic stress!
└── 🔗 AGING: Loss of proteostasis = Hallmark #4

═══════════════════════════════════════════
C. CELL CYCLE — AGING CRITICAL:
═══════════════════════════════════════════
├── Phases: G1 (growth) → S (DNA synthesis) → G2 (prep) → M (mitosis)
├── G0: quiescence (reversible exit from cycle)
│
├── Accelerators (Cyclin-CDK):
│   ├── Cyclin D / CDK4,6 → G1 progression
│   ├── Cyclin E / CDK2 → G1/S transition
│   ├── Cyclin A / CDK2 → S phase
│   └── Cyclin B / CDK1 → G2/M transition
│
├── Brakes (Tumor suppressors):
│   ├── p53: "Guardian of genome" → DNA damage → arrest or apoptosis
│   ├── p21 (CDKN1A): CDK inhibitor, activated by p53
│   ├── p16INK4a: inhibits CDK4/6 → blocks G1 → SENESCENCE marker!
│   └── Rb: sequesters E2F → no S phase genes transcribed
│
├── Checkpoints: G1/S, intra-S, G2/M, spindle assembly
│   └── DNA damage activates ATM/ATR kinases → p53 → p21 → arrest
│
└── Hayflick Limit: normal cells divide ~50 times then senesce
    └── Caused by telomere shortening at each division!

═══════════════════════════════════════════
D. CELL DEATH vs SENESCENCE:
═══════════════════════════════════════════

Apoptosis (Programmed cell death — CLEAN, GOOD):
├── Intrinsic: DNA damage → p53 → BAX/BAK → mitochondria releases cytochrome c
│   → Apoptosome forms → Caspase-9 → Caspase-3/7 → DEATH
├── Extrinsic: Death ligand (FasL) → Fas receptor → FADD → Caspase-8 → Caspase-3/7
├── Anti-apoptotic: Bcl-2, Bcl-xL (block cytochrome c release)
└── Result: cell shrinks, DNA fragments, packaged into bodies, phagocytes eat

Necrosis (Uncontrolled death — MESSY, BAD):
├── Cell swells → membrane ruptures → contents leak → INFLAMMATION
└── Causes: severe injury, ischemia, infection

Senescence (ZOMBIE CELLS — KEY AGING TARGET):
├── TRIGGERS: telomere shortening, DNA damage, ROS, oncogene activation
├── Cell stops dividing BUT remains metabolically active
├── SASP (Senescence-Associated Secretory Phenotype):
│   └── Secretes: IL-6, IL-8, TNF-α, IL-1β, MMP-3, VEGF → inflammatory soup!
├── Markers: p16INK4a↑, p21↑, SA-β-gal↑, γH2AX↑, LMNB1↓
├── Few = useful (wound healing, embryo development)
├── Many = TOXIC (tissue damage, chronic inflammation, cancer promotion)
└── INTERVENTIONS:
    ├── Senolytics: kill senescent cells (dasatinib+quercetin, navitoclax)
    ├── Senomorphics: suppress SASP without killing (rapamycin, metformin)
    └── CAR-T senolytics: T cells engineered to target senescent cells (Amor 2024)

Ferroptosis (Iron-dependent death — newly discovered 2012):
├── Accumulation of lipid peroxides → membrane damage
├── GPX4 enzyme prevents it (uses glutathione)
└── Emerging role in neurodegeneration and aging
```

---

# ✅ TIER 1 COMPLETION CHECKLIST

> **Instructions:** Do NOT move to Tier 2 until ALL boxes are checked.
> For each item, you must be able to EXPLAIN it out loud or WRITE it from memory.

```
═══════════════════════════════════════════════════════════════════════
 SECTION 1.1 — DNA STRUCTURE
═══════════════════════════════════════════════════════════════════════
 □ Can draw DNA double helix and label: sugar-phosphate backbone,
   bases, H-bonds, major groove, minor groove, 5' and 3' ends
 □ Know the 4 bases and pairing: A-T (2 bonds), G-C (3 bonds)
 □ Know: antiparallel means strands run in opposite directions
 □ Understand: euchromatin (open, active) vs heterochromatin (closed, silent)
 □ Know: only ~1.5% of genome codes for proteins
 □ Know: mitochondrial DNA is 16,569bp, circular, maternally inherited
 □ Can explain what CpG dinucleotides are and why they matter for aging
 □ Can write Python code to compute GC content of a DNA sequence
 □ Can compute reverse complement of any short DNA string

═══════════════════════════════════════════════════════════════════════
 SECTION 1.2 — CENTRAL DOGMA
═══════════════════════════════════════════════════════════════════════
 □ Can trace: DNA → transcription → pre-mRNA → processing → mature mRNA → translation → protein
 □ Know what RNA Polymerase II does and where it binds (promoter, TATA box)
 □ Can explain: 5' cap, poly-A tail, splicing, introns vs exons
 □ Know WHY alternative splicing matters (20K genes → 100K+ proteins)
 □ Know the difference between cDNA and genomic DNA
 □ Can explain why STAR aligner (not BWA) is needed for RNA-seq
 □ Know: AUG = start (Met), UAA/UAG/UGA = stop codons
 □ Can translate a short mRNA sequence into amino acids using codon table
 □ Know: 64 codons → 20 amino acids + 3 stops (redundancy via wobble)
 □ Can explain codon optimization and why it matters for synthetic genes

═══════════════════════════════════════════════════════════════════════
 SECTION 1.3 — AMINO ACIDS
═══════════════════════════════════════════════════════════════════════
 □ Can name all 20 amino acids with 1-letter AND 3-letter codes
 □ Can classify: nonpolar, polar uncharged, positive, negative
 □ Know which AAs are phosphorylation targets (S, T, Y) and why it matters
 □ Know Lysine = target for acetylation, methylation, ubiquitination → epigenetics!
 □ Can explain disulfide bonds (2 Cys → S-S bridge → protein stability)
 □ Know Proline causes kinks, Glycine is most flexible
 □ Can explain what BLOSUM62 matrix is and why amino acid properties matter for alignment
 □ Know at least 3 post-translational modifications and their biological role

═══════════════════════════════════════════════════════════════════════
 SECTION 1.4 — BIOCHEMISTRY
═══════════════════════════════════════════════════════════════════════
 □ Understand pH scale, can explain pKa, know blood pH is 7.35-7.45
 □ Know ΔG < 0 = spontaneous, ΔG > 0 = needs energy
 □ Can explain: ATP hydrolysis provides energy for non-spontaneous reactions
 □ Can draw/explain Michaelis-Menten curve (v vs [S], identify Km and Vmax)
 □ Know the 3 types of enzyme inhibition and how they affect Km/Vmax
 □ Know: rapamycin = allosteric inhibitor of mTOR → extends lifespan
 □ Can trace glucose through: glycolysis → TCA → ETC → ATP
 □ Know where ROS are produced (Complex I and III of ETC)
 □ Can explain why NAD⁺ declines with age and why NMN/NR might help
 □ Know: α-ketoglutarate = TCA intermediate AND cofactor for TET enzymes (demethylation!)

═══════════════════════════════════════════════════════════════════════
 SECTION 1.5 — GENE REGULATION
═══════════════════════════════════════════════════════════════════════
 □ Can explain: promoter, enhancer, silencer, insulator, CpG island
 □ Know what CTCF protein does (insulator function)
 □ Know at least 5 aging-relevant TFs (NF-κB, p53, FOXO3, NRF2, TFEB)
 □ Can explain Yamanaka factors (OSKM) and how partial reprogramming reverses aging
 □ Know at least 3 types of ncRNAs (miRNA, lncRNA, siRNA) and their function
 □ Can explain RNAi mechanism: dsRNA → Dicer → RISC → mRNA degradation
 □ Know: ChIP-seq maps TF binding, ATAC-seq maps chromatin accessibility

═══════════════════════════════════════════════════════════════════════
 SECTION 1.6 — CELL BIOLOGY
═══════════════════════════════════════════════════════════════════════
 □ Can describe phospholipid bilayer and fluid mosaic model
 □ Know: receptor-mediated endocytosis is how LNPs and AAVs enter cells
 □ Can explain passive vs active transport with examples
 □ Can list aging-relevant organelles and their role:
   mitochondria, lysosomes, ER, Golgi, proteasome, nucleus
 □ Know: Lamin A mutations → Progeria (premature aging disease)
 □ Can trace: ATP synthesis from ETC (Complex I→IV + ATP synthase)
 □ Know mitophagy pathway: PINK1 → Parkin → selective autophagy
 □ Can explain all cell cycle phases (G1→S→G2→M) and G0
 □ Know: Cyclin D/CDK4,6 drives G1; p16INK4a inhibits CDK4,6 → senescence
 □ Can explain difference between apoptosis (clean) and necrosis (messy)
 □ Can list 5 SASP factors and explain why senescent cells are dangerous
 □ Know senolytics (dasatinib+quercetin) and senomorphics (rapamycin) difference
 □ Know: Baker 2016 proved senescent cell clearance extends lifespan in mice

═══════════════════════════════════════════════════════════════════════
 TOTAL: 55 items. All checked = TIER 1 COMPLETE → proceed to TIER 2
═══════════════════════════════════════════════════════════════════════
```
