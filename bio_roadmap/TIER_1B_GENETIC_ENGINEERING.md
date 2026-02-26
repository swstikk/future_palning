# 🔴 TIER 1B: GENETIC ENGINEERING — THE TOOLS
## CRISPR → Gene Delivery → Recombinant DNA → Synthetic Biology

> **PRIORITY:** CRITICAL — Learn in parallel from Month 4-9
> **PREREQUISITE:** Tier 1 sections 1.1-1.4 (DNA, Central Dogma, Amino Acids, Biochemistry)
> **NOTE:** This runs PARALLEL with Tier 2, not after it

---

## GE.1 CRISPR-Cas SYSTEMS — IN DEPTH
**Why:** This is your primary weapon. Know it at implementation level.

```
═══════════════════════════════════════════
A. HOW Cas9 WORKS (step-by-step molecular mechanism):
═══════════════════════════════════════════

The Players:
├── Cas9 protein: 1,368 amino acids, 2 nuclease domains
│   ├── HNH domain: cuts the guide-complementary (target) strand
│   └── RuvC domain: cuts the non-target strand
├── sgRNA (single guide RNA): ~100nt chimeric molecule
│   ├── crRNA portion (~20nt): complementary to target DNA → SPACER
│   └── tracrRNA portion (~80nt): scaffold that Cas9 recognizes
└── PAM (Protospacer Adjacent Motif): short DNA sequence REQUIRED for recognition
    ├── SpCas9 PAM: 5'-NGG-3' (N = any base, G = guanine)
    ├── PAM must be in the TARGET DNA, not in the guide RNA
    ├── On the NON-target strand
    └── PAM is what distinguishes self (CRISPR array, no PAM) from foreign (phage, has PAM)

The Mechanism:
├── Step 1: sgRNA loads into Cas9 → conformational change → "armed" complex
├── Step 2: Complex scans DNA via 3D diffusion
│   └── Searches for PAM first (PAM-interacting domain reads NGG on DNA)
├── Step 3: PAM found → Cas9 opens 10-12bp of adjacent DNA → R-loop formation
│   └── RNA-DNA hybrid forms IF complementary
├── Step 4: If full 20nt RNA-DNA match → HNH + RuvC domains activated
│   └── HNH cuts target strand (3bp upstream of PAM)
│   └── RuvC cuts non-target strand → BLUNT-ENDED DSB!
├── Step 5: DSB triggers cellular repair:
│   ├── NHEJ (most common): Ku70/80 bind → error-prone ligation → INDELs
│   │   └── USE FOR: Gene KNOCKOUTS (frameshift → no functional protein)
│   └── HDR (homology-directed repair): donor template → precise insertion
│       ├── USE FOR: Gene knock-ins, point mutations, tag insertion
│       ├── Requires: donor template with homology arms (~500bp-1kb each side)
│       └── Efficiency: ~5-50% (lower than NHEJ, S/G2 phase only)
└── Result: gene disruption (NHEJ) or precise edit (HDR)

Off-Target Effects:
├── Cas9 tolerates 1-3 mismatches → can cut unintended sites!
├── Mismatches tolerated better at PAM-distal end (positions 1-8 of spacer)
├── Bulk mismatches at PAM-proximal end (positions 10-20) → usually no cutting
├── COMPUTATIONAL PREDICTION:
│   ├── CRISPOR (crispor.tefor.net): score on-target + off-target
│   ├── CRISPick (Broad Institute): ML-based guide selection
│   ├── Cas-OFFinder: exhaustive off-target search
│   └── GUIDE-seq, CIRCLE-seq: experimental off-target detection
├── High-fidelity variants (reduced off-targets):
│   ├── eSpCas9 (1.1): 3 mutations → less off-target, slightly less on-target
│   ├── SpCas9-HF1: 4 mutations → specificity improved
│   └── HiFi Cas9 (IDT): single R691A mutation → commercial standard
└── Most off-target concern: clinical gene therapy (can't afford even 1 off-target cut!)

═══════════════════════════════════════════
B. CRISPR VARIANTS — YOUR EXPANDING TOOLKIT:
═══════════════════════════════════════════

dCas9 (dead Cas9): D10A + H840A mutations → NO cutting, but still binds DNA
├── CRISPRi (interference): dCas9 + KRAB repressor domain
│   └── Silences gene WITHOUT cutting DNA → reversible, no permanent change!
├── CRISPRa (activation): dCas9 + VP64, p65, Rta (VPR) activator domains
│   └── Activates gene WITHOUT editing → turn ON endogenous genes!
├── Epigenome editing:
│   ├── dCas9-DNMT3A → methylate specific CpGs → silence specific gene
│   ├── dCas9-TET1 → demethylate specific CpGs → activate specific gene  
│   ├── dCas9-p300 → add H3K27ac → activate enhancer
│   └── dCas9-KRAB-MeCP2 → combo repression → strongest silencing
└── AGING USE: dCas9-TET1 → demethylate age-methylated CpGs → reset clock!

Base Editors (NO double-strand break!):
├── CBE (Cytosine Base Editor):
│   ├── nCas9 (D10A, nicks 1 strand) + APOBEC1 cytidine deaminase
│   ├── Converts C → U (reads as T) within editing window (positions 4-8)
│   ├── UGI (uracil glycosylase inhibitor) prevents cell from reverting U→C
│   └── Net result: C•G → T•A base pair conversion
├── ABE (Adenine Base Editor):
│   ├── nCas9 + TadA* (evolved adenine deaminase)
│   ├── Converts A → Inosine (reads as G) within editing window
│   └── Net result: A•T → G•C base pair conversion
├── Coverage: CBE + ABE together can make 4 of 12 possible transitions
│   └── C→T, G→A (CBE), A→G, T→C (ABE) = the 4 transition mutations
├── ABE8e: latest high-efficiency variant from David Liu lab
└── DdCBE: base editor for MITOCHONDRIAL DNA!
    ├── Uses TALE proteins (not Cas9) because RNA can't enter mitochondria
    ├── First tool ever to edit mtDNA (Mok et al., Nature 2020)
    └── AGING: fix mtDNA mutations that accumulate with age!

Prime Editing (THE most precise):
├── pegRNA: extended sgRNA with reverse transcript template + primer binding site
├── PE2: nCas9-H840A + M-MLV reverse transcriptase
├── PE3: PE2 + additional nicking guide for opposite strand
├── Mechanism:
│   ├── Nick target strand → RT uses pegRNA template → writes new sequence
│   ├── Cell incorporates new sequence during repair
│   └── Net result: any edit you want at that location!
├── Can do: all 12 substitutions + small insertions up to ~40bp + small deletions
├── Advantage: no DSB, no donor template DNA needed
└── Disadvantage: lower efficiency (~10-30%), complex guide design

Cas Variants Beyond Cas9:
├── SaCas9 (S. aureus): 1,053 aa → fits in single AAV! PAM: NNGRRT
├── Cas12a (Cpf1): PAM = TTTN (upstream!), creates staggered cuts (5' overhangs)
│   ├── Self-processes crRNA array → multiplex editing from one transcript!
│   └── Better for HDR (staggered cuts → sticky ends → easier ligation)
├── Cas12b (C2c1): smaller, thermostable variants exist
├── Cas13 (RNA targeting): cuts RNA not DNA → reversible knockdown!
│   ├── No genome modification → safer for therapeutic applications
│   └── REPAIR: RNA editing → change specific nucleotides in transcript
├── CasΦ (Cas-Phi): very small, found in bacteriophages
├── Cas14 (CasX): ultra-compact, <1000 aa
└── IscB/IsrB: ancestors of Cas9, even smaller → future engineering potential

═══════════════════════════════════════════
C. GUIDE RNA DESIGN (computational):
═══════════════════════════════════════════
├── Spacer: 20nt sequence complementary to target
│   └── Must have adjacent PAM (NGG for SpCas9) in target DNA
├── Design criteria:
│   ├── GC content: 40-70% optimal (too low = weak binding, too high = off-target)
│   ├── Avoid polyT runs (>4 T's = Pol III termination signal!)
│   ├── Avoid secondary structure in spacer (hairpins block target recognition)
│   ├── Off-target score: minimize sites with <3 mismatches genome-wide
│   └── Position within gene: for KO → target early exons (constitutive ones)
│
├── Computational Tools:
│   ├── CRISPOR: gold standard web tool → on-target + off-target scores
│   ├── CRISPick (Broad): ML-optimized guide selection
│   ├── Benchling: integrated design + ordering platform
│   └── Python: crispy, guidescan libraries
│
└── 🔗 ML CONNECTION: 
    ├── DeepCRISPR: deep learning model for guide efficiency prediction
    ├── CRISPR-ML: GNN-based off-target prediction
    └── THIS is a research area YOU can contribute to right now!

═══════════════════════════════════════════
D. CRISPR FOR AGING (applications you'll build):
═══════════════════════════════════════════
├── Telomere elongation:
│   └── CRISPRa for TERT → activate telomerase in somatic cells
├── Epigenetic reprogramming:
│   ├── dCas9-TET1 → demethylate age-associated CpGs
│   └── dCas9-p300 → add activating marks to silenced longevity genes
├── Senolytic gene circuits:
│   ├── CRISPRi p16/p21 → push zombie cells toward apoptosis
│   └── Combine with p16-promoter for specificity (Tier GE.4)
├── Mitochondrial repair:
│   └── DdCBE → fix pathogenic mtDNA mutations (3243A>G, 8344A>G)
├── Longevity gene activation:
│   └── CRISPRa for: SIRT6, FOXO3, TERT, TFEB
├── Disease gene correction:
│   └── Base editing for: progeria (LMNA C1824T), sickle cell (HBB), etc.
└── Immune rejuvenation:
    └── Edit CAR-T cells to target senescent cell markers → immuno-senolytics!
```

---

## GE.2 GENE DELIVERY — GET IT INTO THE CELL
**Why:** Designed therapy is worthless if it can't reach target cells.

```
═══════════════════════════════════════════
A. VIRAL VECTORS:
═══════════════════════════════════════════

AAV (Adeno-Associated Virus) — GOLD STANDARD:
├── Non-integrating: stays as episome in nucleus → lower insertional risk
├── Genome: 4.7kb ssDNA → PACKAGING LIMIT ~4.7kb including ITRs
│   └── PROBLEM: SpCas9 cDNA alone is 4.1kb → barely fits with promoter!
│   └── SOLUTION: use SaCas9 (smaller), split-intein Cas9, or use mRNA-LNP
├── Serotypes (capsid variants → different tissue tropism):
│   ├── AAV1: skeletal muscle, heart
│   ├── AAV2: retina, liver, CNS → MOST studied
│   ├── AAV5: lung, CNS
│   ├── AAV6: muscle, lung, hematopoietic
│   ├── AAV8: liver (highest transduction) → metabolic diseases
│   ├── AAV9: CNS (crosses blood-brain barrier!), heart, muscle
│   └── AAVrh10: CNS → emerging for neurodegeneration
├── Production: HEK293T cells + triple transfection:
│   ├── Plasmid 1: ITR-flanked transgene (your gene of interest)
│   ├── Plasmid 2: Rep/Cap (replication + capsid proteins)
│   └── Plasmid 3: Adenovirus helper functions
├── FDA-approved: Luxturna (AAV2, retina), Zolgensma (AAV9, SMA), Hemgenix (AAV5, hemophilia B)
└── Limitations:
    ├── Pre-existing immunity: ~40-60% of humans have anti-AAV antibodies
    ├── Immunogenicity: can trigger immune response at high doses
    └── Re-dosing difficult → one-shot therapy for most current approaches

Lentivirus:
├── Integrates into host genome → permanent modification
│   └── 3rd generation: self-inactivating (SIN) → safer design
├── Packaging: ~8kb insert capacity → larger than AAV!
├── Risk: insertional mutagenesis → can activate oncogenes
│   └── Occurred in early SCID-X1 gene therapy (X-SCID trial, 2002)
├── Use: ex vivo therapy (take cells out → transduce → put back)
│   └── CAR-T cell manufacturing uses lentivirus heavily
└── Pseudotyping: VSV-G envelope protein → broad tropism

Adenovirus:
├── NON-integrating, high transduction, large insert (~36kb!)
├── Strong immune response → inflammation → NOT for chronic therapy
└── Use: oncolytic virotherapy (exploit immune activation against cancer!)

═══════════════════════════════════════════
B. NON-VIRAL DELIVERY:
═══════════════════════════════════════════

LNPs (Lipid Nanoparticles) — THE FUTURE:
├── Components:
│   ├── Ionizable lipid: positive at low pH (endosomal escape), neutral at blood pH
│   ├── PEG-lipid: stealth layer, prevents immune clearance
│   ├── Cholesterol: structural stability
│   └── Phospholipid (DSPC): bilayer structure
├── Payload: mRNA, siRNA, Cas9 mRNA + sgRNA, base editor mRNA
├── Entry mechanism:
│   ├── ApoE coating in blood → LDL receptor-mediated endocytosis → liver
│   ├── Endosome → acidification → ionizable lipid becomes positive → membrane disruption
│   └── Cargo escapes into cytoplasm → translated by ribosomes
├── Advantages: no genome integration, no size limit, transient (safety), scalable
├── Disadvantages: liver tropism (hard to target other organs), transient expression
├── FDA-approved: ONPATTRO (patisiran, siRNA-LNP), COVID vaccines (Pfizer/Moderna)
├── EMERGING: SORT technology → organ-selective LNPs (lung, spleen, etc.)
│   └── Adding charged lipids changes biodistribution!
└── For CRISPR: deliver Cas9 mRNA + sgRNA via LNP → transient editing → safer!

Other non-viral:
├── Electroporation: electric pulse → temporary pores → cargo enters
│   └── For: ex vivo cell therapy (T cells, iPSCs)
├── Polymeric nanoparticles: PEI, PLGA → sustained release
├── Exosomes: natural cell-derived vesicles → low immunogenicity
├── VLPs (Virus-Like Particles): viral capsid without genome → safe delivery
└── Cell-penetrating peptides (CPPs): TAT, penetratin → conjugate to cargo

═══════════════════════════════════════════
C. DELIVERY STRATEGY FOR AGING (tissue-specific):
═══════════════════════════════════════════
├── LIVER: LNP (mRNA) or AAV8 → metabolic interventions, NAD⁺
├── MUSCLE: AAV1/6/9 → mitochondrial diseases, sarcopenia
├── BRAIN: AAV9/rh10 or intrathecal → neurodegeneration, cognitive aging
├── BLOOD CELLS: Lentivirus ex vivo or LNP → immune rejuvenation
├── EYE: AAV2 subretinal → retinal degeneration (Luxturna model)
├── SKIN: Electroporation or topical LNP → senolytic testing, clock validation
└── SYSTEMIC: IV LNP → liver first, then engineer for pan-tissue delivery

🔗 COMPUTATIONAL CONNECTION:
├── Capsid engineering: ML used to design novel AAV capsids (directed evolution + AI)
│   └── Dyno Therapeutics → AI-designed AAV capsids for CNS
├── Codon optimization: optimize gene sequence for expression in target species/cell
├── LNP formulation: ML screening of lipid combinations for optimal delivery
└── Promoter design: computational regulatory element analysis
```

---

## GE.3 RECOMBINANT DNA TECHNOLOGY — BUILDING BLOCKS
**Why:** Still used in every molecular biology lab. Foundation for all genetic engineering.

```
═══════════════════════════════════════════
A. RESTRICTION ENZYMES:
═══════════════════════════════════════════
├── Type II: most useful, recognize palindromic sequences (4-8bp)
├── Examples:
│   ├── EcoRI: 5'...G↓AATTC...3' → 5' overhang (sticky end)
│   ├── HindIII: 5'...A↓AGCTT...3' → 5' overhang
│   ├── BamHI: 5'...G↓GATCC...3' → 5' overhang
│   ├── SmaI: 5'...CCC↓GGG...3' → blunt end
│   └── NotI: 5'...GC↓GGCCGC...3' → 8bp recognition (rare cutter)
├── Sticky ends: complementary overhangs → self-anneal → easier ligation
├── Blunt ends: no overhangs → less efficient ligation
├── Isoschizomers: different enzymes, same recognition site
└── Methylation-sensitive: some cannot cut methylated DNA → biological significance!

═══════════════════════════════════════════
B. CLONING WORKFLOW:
═══════════════════════════════════════════
├── 1. Cut: vector + insert with SAME restriction enzymes
├── 2. Ligate: T4 DNA Ligase → joins vector + insert (phosphodiester bonds)
│   └── Molar ratio: 3:1 insert:vector optimal
├── 3. Transform: heat shock (42°C, 45sec) or electroporation → into E. coli
├── 4. Select: grow on antibiotic plates → only transformed bacteria survive
├── 5. Screen:
│   ├── Colony PCR: PCR with insert-specific primers
│   ├── Blue-white selection: lacZ disruption → white colony = insert present
│   └── Restriction digest: cut plasmid, check band pattern on gel
└── 6. Verify: Sanger sequencing of insert

Modern alternatives:
├── Gibson Assembly: seamless cloning, overlapping fragments, no restriction enzymes!
├── Golden Gate Assembly: Type IIS enzymes, directional, multi-fragment
├── Gateway cloning: recombination-based (attB/attP sites)
└── TOPO cloning: topoisomerase-based, fast, for PCR products

═══════════════════════════════════════════
C. PLASMID VECTORS:
═══════════════════════════════════════════
├── Core components:
│   ├── Origin of replication (ORI): determines copy number
│   │   ├── High-copy: pUC ori (~500 copies/cell)
│   │   └── Low-copy: pSC101 ori (~5 copies/cell)
│   ├── Selectable marker: antibiotic resistance (AmpR, KanR, CmR)
│   ├── Multiple Cloning Site (MCS): cluster of unique restriction sites
│   └── Promoter (for expression vectors): drives insert transcription
├── Expression vector extras:
│   ├── Inducible promoter: T7/lac (IPTG-inducible), Tet-On/Off, CMV
│   ├── Fusion tags: His6 (purification), GFP (visualization), FLAG (detection)
│   ├── Ribosome binding site (RBS) or Kozak sequence
│   └── Terminator: stops transcription
└── Reading vector maps: circular diagram with all elements annotated
    └── GenBank format: standard text format for sequences + annotations

═══════════════════════════════════════════
D. PCR (Polymerase Chain Reaction) — CRITICAL:
═══════════════════════════════════════════
├── Purpose: amplify specific DNA region exponentially
├── Components: template DNA, forward primer, reverse primer, DNA Pol (Taq), dNTPs, Mg²⁺, buffer
├── Steps per cycle:
│   ├── Denaturation: 94-98°C → separate DNA strands
│   ├── Annealing: 50-65°C → primers bind template (Tm-dependent)
│   └── Extension: 72°C → Taq polymerase synthesizes new strand (5'→3')
├── Exponential: n cycles → 2ⁿ copies (30 cycles → ~10⁹ copies!)
│
├── Primer design:
│   ├── Length: 18-25bp
│   ├── Tm: 55-65°C (both primers should have similar Tm)
│   ├── GC content: 40-60%
│   ├── 3' end: end with G or C (GC clamp → stronger binding)
│   └── Tool: Primer3, Primer-BLAST (NCBI)
│
├── PCR Variants:
│   ├── RT-PCR: RNA → cDNA (reverse transcriptase) → PCR
│   ├── qPCR (real-time): quantify DNA/RNA amount in real-time
│   │   └── SYBR Green or TaqMan probes → Ct value = cycle threshold
│   ├── Digital PCR (dPCR): absolute quantification, no standard curve
│   ├── Nested PCR: two rounds with different primers → higher specificity
│   ├── Overlap extension PCR: join two fragments → mutations, fusions
│   └── Emulsion PCR: individual molecules in droplets → NGS library prep
│
└── BIOINFORMATICS:
    ├── in silico PCR: UCSC Genome Browser → predict amplicon from primer pair
    ├── Primer3: automated primer design given target region
    └── Tm calculation: nearest-neighbor thermodynamic model

═══════════════════════════════════════════
E. GEL ELECTROPHORESIS:
═══════════════════════════════════════════
├── Agarose gel (0.5-2%): separate DNA by size (100bp-10kb range)
├── Polyacrylamide gel (PAGE): higher resolution, smaller fragments or proteins
├── SDS-PAGE: proteins denatured + SDS → separate by molecular weight
├── Staining: ethidium bromide (UV, mutagenic!) or SYBR Safe (safer)
├── DNA ladder: size markers (1kb ladder, 100bp ladder)
└── Western blot: transfer protein from gel to membrane → probe with antibody
    └── Used to verify protein expression after gene editing!
```

---

## GE.4 SYNTHETIC BIOLOGY — PROGRAMMING CELLS
**Why:** You want to build a senolytic gene circuit. This section teaches how.

```
═══════════════════════════════════════════
A. BIOLOGICAL PARTS & CIRCUIT COMPONENTS:
═══════════════════════════════════════════
├── Input sensors (promoters):
│   ├── Constitutive: always ON (CMV, EF1α, PGK in mammalian; T7, lac in E. coli)
│   ├── Inducible:
│   │   ├── Tet-On/Tet-Off: doxycycline-controlled (most used in mammalian research)
│   │   ├── Cumate: cumate-inducible → orthogonal to Tet system
│   │   ├── LacI/IPTG: standard in E. coli expression
│   │   ├── AraC/L-arabinose: another E. coli standard
│   │   └── Light-inducible (optogenetics): blue light → CRY2/CIB1 dimerization
│   └── Disease-responsive:
│       ├── p16INK4a promoter → active only in SENESCENT cells!
│       ├── HIF1α-responsive → active in hypoxic (tumor) environments
│       └── NF-κB responsive → active during inflammation
│
├── Processing (regulatory elements):
│   ├── Riboswitches: RNA structure that responds to small molecules
│   ├── IRES: internal ribosome entry site → multiple proteins from one mRNA
│   ├── 2A peptides (T2A, P2A): self-cleaving → multiple proteins from one ORF
│   └── Destabilization domains (DD): protein degraded unless stabilized by drug
│
└── Outputs (effector genes):
    ├── Reporters: GFP, mCherry, luciferase → visualize circuit activity
    ├── Therapeutic: pro-apoptotic (truncated BID, PUMA), cytokines, antibodies
    ├── Kill switches: iCaspase9 → dimerized by AP1903 → apoptosis ON DEMAND
    └── Selection: antibiotic resistance, puromycin, blasticidin

═══════════════════════════════════════════
B. GENETIC LOGIC GATES:
═══════════════════════════════════════════
├── NOT gate: repressor A → when A present, output OFF
│   └── Implementation: TetR protein represses promoter → add doxycycline → TetR inactive → gene ON
├── AND gate: both inputs needed
│   └── Implementation: split T7 RNAP → each half under different promoter → functional only if BOTH present
├── OR gate: either input sufficient
│   └── Implementation: two different promoters driving same output gene
├── NAND gate: NOT AND → output OFF only if BOTH inputs present
├── NOR gate: NOT OR → output ON only if NEITHER input present
├── Toggle switch: bistable, two stable states
│   └── Two repressors mutually inhibit → cell "remembers" which state it's in
│   └── Collins lab (2000): first synthetic toggle switch in E. coli
└── Band-pass filter: output ON only in medium range of input
    └── Use cascading NOT gates with different thresholds

═══════════════════════════════════════════
C. SENOLYTIC GENE CIRCUIT (YOUR TARGET PROJECT):
═══════════════════════════════════════════

Design specification:
├── Goal: Kill senescent cells, leave normal cells ALIVE
├── Input: p16INK4a expression level (HIGH in senescent, LOW in normal)
│
├── VERSION 1 (simple):
│   ├── p16 promoter → Caspase 3 (apoptosis gene)
│   ├── p16 HIGH → Caspase 3 expressed → senescent cell dies ✓
│   ├── p16 LOW → no Caspase 3 → normal cell survives ✓
│   └── Problem: even LOW p16 expression might activate → LEAKY! → kills normal cells!
│
├── VERSION 2 (with safety AND gate):
│   ├── Input 1: p16 promoter → half of split Caspase
│   ├── Input 2: SA-β-gal reporter → other half of split Caspase
│   ├── AND logic: BOTH p16 HIGH + SA-β-gal HIGH → Caspase reconstituted → death
│   └── Much safer: requires TWO senescence markers
│
├── VERSION 3 (with external kill switch):
│   ├── AND gate V2 + iCaspase9 controlled by AP1903 drug
│   ├── Even if circuit misfires → drug must be administered to execute
│   └── Safety layers: senescence markers + drug trigger
│
└── COMPUTATIONAL DESIGN:
    ├── Cello (MIT): compiler that takes logic specification → outputs DNA sequence
    ├── iBioSim: circuit simulation (ODE-based → predict circuit behavior)
    ├── Benchling: sequence design + plasmid visualization
    └── Parts: igem.org/Registry → standard biological parts catalog

═══════════════════════════════════════════
D. CELL REPROGRAMMING — RESET THE CLOCK:
═══════════════════════════════════════════
├── Yamanaka Factors (OSKM): Oct4, Sox2, Klf4, c-Myc
│   ├── Full expression → iPSC (fully dedifferentiated) → cancer risk!
│   └── Partial/pulsatile expression → age reversal WITHOUT dedifferentiation!
│
├── Delivery methods for OSKM:
│   ├── mRNA transfection → transient, safe → expensive, repeated dosing
│   ├── Sendai virus → non-integrating RNA virus → approved for iPSC generation
│   ├── AAV with inducible promoter (Tet-On) → controllable, ONE injection
│   └── LNP-mRNA → transient, scalable → liver-tropic (challenge for other tissues)
│
├── Partial reprogramming results:
│   ├── Ocampo et al. 2016: cyclic OSKM in progeria mice → 30% lifespan extension!
│   ├── Lu et al. 2020: OSK (no c-Myc, safer) → restored vision in aged mice
│   ├── Retro Biosciences + OpenAI 2024: 50× efficiency with ML optimization
│   └── Altos Labs: massive funding ($3B) for partial reprogramming → human trials
│
├── Key safety concern:
│   ├── c-Myc is an ONCOGENE → drives cancer if overexpressed
│   ├── Solution: use OSK only (drop Myc) or use pulsatile delivery
│   └── Oscillator circuit (GE.4-C) → pulsatile Yamanaka delivery → safest approach
│
└── Non-Yamanaka approaches:
    ├── Chemical reprogramming: small molecules mimic OSKM effect
    │   └── 7-compound cocktail (Deng Lab, 2022) → iPSC without TF overexpression!
    ├── Direct transdifferentiation: fibroblast → neuron (Ascl1, Brn2, Myt1l)
    └── In vivo reprogramming: reprogram cells INSIDE living organism (no ex vivo)
```

---

# ✅ TIER 1B COMPLETION CHECKLIST

```
═══════════════════════════════════════════════════════════════════════
 SECTION GE.1 — CRISPR
═══════════════════════════════════════════════════════════════════════
 □ Can trace Cas9 mechanism: sgRNA loading → PAM scan → R-loop → DSB
 □ Know: HNH cuts target strand, RuvC cuts non-target strand
 □ Know: PAM = NGG for SpCas9, must be in target DNA (not guide)
 □ Can explain NHEJ (error-prone, knockouts) vs HDR (precise, needs template)
 □ Know off-target risks and at least 2 computational tools to predict them
 □ Know at least 3 high-fidelity Cas9 variants (eSpCas9, HF1, HiFi)
 □ Can explain dCas9 applications: CRISPRi, CRISPRa, epigenome editing
 □ Can explain base editing: CBE (C→T), ABE (A→G), and DdCBE (mitochondrial!)
 □ Can explain prime editing: pegRNA, reverse transcriptase, all 12 substitutions
 □ Know 4 Cas variants beyond Cas9: SaCas9, Cas12a, Cas13, CasΦ
 □ Can design a guide RNA: check GC%, off-target, polyT, PAM, position
 □ Know at least 3 tools for guide design (CRISPOR, CRISPick, Benchling)
 □ Can list 5 specific CRISPR applications for aging research

═══════════════════════════════════════════════════════════════════════
 SECTION GE.2 — GENE DELIVERY
═══════════════════════════════════════════════════════════════════════
 □ Know AAV: non-integrating, 4.7kb limit, serotypes determine tissue tropism
 □ Can match serotype to tissue: AAV2 (retina), AAV8 (liver), AAV9 (CNS/muscle)
 □ Know: pre-existing immunity (~40-60%) is a challenge for AAV
 □ Can explain lentivirus: integrating, ex vivo, larger capacity, CAR-T use
 □ Know LNP components: ionizable lipid, PEG-lipid, cholesterol, phospholipid
 □ Can explain how LNPs enter cells: ApoE → LDL receptor → endocytosis → escape
 □ Know: COVID vaccines = mRNA-LNP technology (Pfizer/Moderna)
 □ Know 3 FDA-approved gene therapies and which vector they use
 □ Can match delivery strategy to aging target tissue (liver, brain, muscle, blood, skin)
 □ Know: electroporation used for ex vivo cell transfection

═══════════════════════════════════════════════════════════════════════
 SECTION GE.3 — RECOMBINANT DNA
═══════════════════════════════════════════════════════════════════════
 □ Know at least 4 restriction enzymes and their sites (EcoRI, HindIII, BamHI, NotI)
 □ Can explain sticky ends vs blunt ends → which ligates more efficiently
 □ Can trace cloning workflow: cut → ligate → transform → select → screen → verify
 □ Know modern cloning alternatives: Gibson Assembly, Golden Gate, Gateway
 □ Can read a plasmid vector map: ORI, marker, MCS, promoter
 □ Can explain PCR mechanism: denature → anneal → extend (3 temperatures)
 □ Know PCR variants: RT-PCR, qPCR (Ct), digital PCR, nested PCR
 □ Can design primers: 18-25bp, Tm 55-65°C, GC 40-60%, GC clamp
 □ Know gel electrophoresis: agarose (DNA), PAGE (protein), Western blot
 □ Can use Primer3 and NCBI Primer-BLAST for primer design

═══════════════════════════════════════════════════════════════════════
 SECTION GE.4 — SYNTHETIC BIOLOGY
═══════════════════════════════════════════════════════════════════════
 □ Can name 5 inducible promoter systems (Tet-On, cumate, LacI, AraC, optogenetics)
 □ Can explain at least 4 genetic logic gates (NOT, AND, OR, toggle switch)
 □ Can design a senolytic circuit: p16 promoter → apoptosis gene (version 1)
 □ Can explain why AND gate version is safer (requires 2 senescence markers)
 □ Can explain kill switch safety mechanism (iCaspase9 + AP1903)
 □ Know Yamanaka factors (OSKM) and the difference between full vs partial reprogramming
 □ Know: pulsatile OSKM delivery safer than continuous (oscillator circuit!)
 □ Can explain why c-Myc is dangerous and why OSK-only is preferred
 □ Know: chemical reprogramming (7-compound cocktail) as alternative to TF overexpression
 □ Know design tools: Cello (compiler), iBioSim (simulation), Benchling (design)
 □ Can explain 2A peptides and IRES as multi-gene expression strategies

═══════════════════════════════════════════════════════════════════════
 TOTAL: 44 items. All checked = TIER 1B COMPLETE → proceed to TIER 3
═══════════════════════════════════════════════════════════════════════
```
