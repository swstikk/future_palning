# 🟠 TIER 2: ESSENTIAL KNOWLEDGE
## Genetics → Cell Signaling → Epigenetics → Microbiology → Phylogenetics

> **PRIORITY:** HIGH — Learn within Month 3-6
> **PREREQUISITE:** Tier 1 checklist 100% complete
> **NOTE:** Protein Structure moved here from Tier 3 (ordering fix — needed before Epigenetics)

---

## 2.1 PROTEIN STRUCTURE & FUNCTION
**Why:** Epigenetics discusses protein domains (bromodomains, chromodomains). Drug design targets proteins. AlphaFold = your key AI tool.

```
⚠️ ORDERING FIX: Moved from Tier 3 → Tier 2 (needed BEFORE epigenetics)

═══════════════════════════════════════════
A. FOUR LEVELS OF PROTEIN STRUCTURE:
═══════════════════════════════════════════
├── Primary (1°): Linear amino acid sequence
│   └── Written N→C terminal, determined by gene sequence
│   └── Single mutation can change ONE AA → disease (sickle cell: Glu→Val at position 6)
│
├── Secondary (2°): Local folding patterns (H-bonds between backbone atoms)
│   ├── α-helix: right-handed spiral, H-bond every 3.6 residues
│   │   └── Common in transmembrane proteins (hydrophobic AAs in helix)
│   ├── β-sheet: extended strands side-by-side, H-bonds between strands
│   │   ├── Parallel (same direction) or Antiparallel (opposite) β-sheets
│   │   └── Amyloid fibrils = misfolded β-sheets = Alzheimer's, Parkinson's
│   ├── Loops/Turns: connecting regions, often on protein surface
│   └── Intrinsically Disordered Regions (IDRs): no fixed 2° structure!
│       └── ~40% of human proteins have IDRs → involved in signaling, regulation!
│
├── Tertiary (3°): Complete 3D fold of single polypeptide
│   ├── Driven by: hydrophobic collapse (nonpolar AAs hide inside)
│   ├── Stabilized by: H-bonds, ionic bonds, disulfide (S-S), Van der Waals
│   ├── Domains: independently folding units within a protein
│   │   └── Same domain can appear in different proteins (modular design!)
│   └── Misfolding → aggregation → toxicity (prions, amyloid)
│
└── Quaternary (4°): Multiple polypeptide subunits
    ├── Hemoglobin: 4 subunits (2α + 2β) = tetramer
    ├── Cas9: single polypeptide (no quaternary structure)
    └── Ribosomes: dozens of protein + RNA subunits

═══════════════════════════════════════════
B. KEY PROTEIN CLASSES FOR YOU:
═══════════════════════════════════════════
├── Enzymes: catalyze reactions (you learned kinetics in Tier 1.4)
│   ├── Example: Cas9 is an ENDONUCLEASE (cuts DNA)
│   └── Example: Telomerase is a REVERSE TRANSCRIPTASE (extends telomeres)
│
├── Transcription factors: bind DNA → regulate genes (Tier 1.5)
│   ├── DNA-binding domains: zinc finger, HTH, bHLH, bZIP
│   └── Example: p53, FOXO3, NF-κB
│
├── Structural proteins:
│   ├── Collagen: most abundant protein in body, triple helix → ECM
│   │   └── Cross-linking increases with age → tissue stiffness → wrinkles, fibrosis
│   ├── Lamins A/B: nuclear envelope → mutations = progeria
│   └── Actin, Tubulin: cytoskeleton
│
├── Signaling proteins:
│   ├── Kinases: add phosphate (mTOR, AMPK, AKT, ERK, CDKs)
│   ├── Phosphatases: remove phosphate (PTEN, PP2A)
│   ├── GTPases: molecular switches (RAS family — oncogenes!)
│   └── Receptors: transmembrane proteins that sense external signals
│
├── Chaperones (protein quality control):
│   ├── HSP70: prevents aggregation, refolds partially unfolded proteins
│   ├── HSP90: stabilizes signaling proteins (kinases, steroid receptors)
│   ├── HSP60 (chaperonin): barrel-shaped, helps fold proteins inside chamber
│   └── AGING: chaperone expression DECLINES → proteotoxic stress!
│
└── Ubiquitin-Proteasome System:
    ├── E1 → E2 → E3 ligases: tag proteins with ubiquitin chains
    ├── RING domain: E3 ligase catalytic domain
    ├── K48 polyubiquitin: target for proteasome degradation
    ├── K63 polyubiquitin: signaling (not degradation)
    └── AGING: proteasome activity drops 50% with age!

═══════════════════════════════════════════
C. PROTEIN DOMAINS & MOTIFS (used in bioinformatics):
═══════════════════════════════════════════
├── SH2 domain: binds phosphotyrosine → signal transduction
├── SH3 domain: binds proline-rich peptides
├── PH domain: binds membrane phospholipids (PIP2/PIP3) → membrane targeting
├── Bromodomain: reads acetylated lysine → epigenetic reader!
├── Chromodomain: reads methylated lysine → epigenetic reader!
├── RING domain: E3 ubiquitin ligase activity
├── Kinase domain: phosphotransferase catalytic domain
└── Databases: Pfam (domains), InterPro (families), PROSITE (motifs)

🔗 AI/BIOINFORMATICS CONNECTION:
├── AlphaFold2/3 (DeepMind): predicts 3D structure from sequence!
│   └── 200M structures predicted → revolutionized structural biology
├── ESMFold (Meta): protein language model → even faster structure prediction
├── RFdiffusion (Baker Lab): DESIGN new proteins from scratch!
├── Protein-protein interaction: STRING database
└── Molecular dynamics simulation → future Quantum Computing application!
```

---

## 2.2 GENETICS — THE INHERITANCE RULES
**Why:** Gene variants, polymorphisms, GWAS — all of genomics builds on this.

```
═══════════════════════════════════════════
A. MENDELIAN GENETICS (fast track — 1-2 days):
═══════════════════════════════════════════
├── Alleles: different versions of same gene (e.g., eye color variants)
├── Dominant: expressed when 1 copy present (A_)
├── Recessive: expressed only when homozygous (aa)
├── Genotype: genetic makeup (AA, Aa, aa)
├── Phenotype: observable trait
├── Punnett square: predict offspring ratios
├── Homozygous (AA/aa) vs Heterozygous (Aa)
└── Hardy-Weinberg Equilibrium: p² + 2pq + q² = 1 (no evolution scenario)
    └── Used in GWAS to check for population stratification

═══════════════════════════════════════════
B. BEYOND MENDEL (actually important for bioinformatics):
═══════════════════════════════════════════
├── Incomplete dominance: heterozygote intermediate (red + white → pink flower)
├── Codominance: both alleles fully expressed (AB blood type)
├── Epistasis: one gene masks/modifies another gene's effect
│   └── Common in complex disease genetics, detected in GWAS
├── Polygenic traits: hundreds/thousands of genes each contribute small effect
│   ├── Height, BMI, intelligence, LIFESPAN!
│   └── Polygenic Risk Score (PRS): sum of all contributing SNP effects
├── Pleiotropy: one gene → multiple traits
│   └── Antagonistic pleiotropy theory of aging: genes beneficial young → harmful old!
├── Linkage: genes close on same chromosome inherited together
├── Recombination: crossing over during meiosis → breaks linkage
└── Linkage Disequilibrium (LD): statistical correlation between nearby variants
    └── KEY in GWAS: SNPs in LD with causal variant also show association!

═══════════════════════════════════════════
C. MUTATIONS — THE ERRORS:
═══════════════════════════════════════════
├── Point mutations:
│   ├── Silent/Synonymous: codon changes but SAME amino acid (wobble position)
│   ├── Missense: different amino acid → may affect protein function
│   │   └── Conservative: similar AA (Leu→Ile) vs Non-conservative (Glu→Lys)
│   ├── Nonsense: creates STOP codon → truncated protein (usually loss of function)
│   └── Transition vs Transversion:
│       ├── Transition: purine↔purine (A↔G) or pyrimidine↔pyrimidine (C↔T) — MORE common
│       └── Transversion: purine↔pyrimidine — LESS common
│
├── Frameshift: insertion/deletion of non-multiple-of-3 bases → shifts reading frame
│   └── Usually devastating: changes all downstream amino acids
├── Splice site mutations: affect exon/intron boundaries → aberrant mRNA
├── Copy Number Variants (CNVs): large segments duplicated/deleted (>1kb)
├── Structural Variants (SVs): inversions, translocations, large insertions
├── Microsatellite instability: tandem repeat expansion/contraction
│   └── Huntington's disease: CAG repeat expansion in HTT gene
└── Somatic vs Germline:
    ├── Germline: inherited, in every cell (relevant for genetic diseases)
    └── Somatic: acquired, only in affected cell lineage (cancer, aging!)
        └── Somatic mutation burden increases with age → Hallmark #1!

═══════════════════════════════════════════
D. DNA REPAIR MECHANISMS:
═══════════════════════════════════════════
├── Base Excision Repair (BER): small single-base damage
│   └── Glycosylase → AP endonuclease → polymerase → ligase
├── Nucleotide Excision Repair (NER): bulky lesions (UV-induced pyrimidine dimers)
│   └── Xeroderma Pigmentosum: NER-deficient → extreme UV sensitivity → premature skin aging
├── Mismatch Repair (MMR): replication errors (wrong base inserted)
│   └── MLH1, MSH2 deficiency → Lynch syndrome (hereditary cancer)
├── Homologous Recombination (HR): accurate DSB repair using sister chromatid
│   └── BRCA1, BRCA2 → HR-deficient = cancer predisposition
├── Non-Homologous End Joining (NHEJ): error-prone DSB repair (any phase)
│   └── THIS is what repairs CRISPR-Cas9 cuts → causes INDELs for knockouts!
│
└── AGING: DNA repair capacity DECLINES → unrepaired damage accumulates
    ├── Werner syndrome: WRN helicase → premature aging
    ├── Cockayne syndrome: NER deficiency → accelerated aging
    └── Most progeroid syndromes = DNA repair gene mutations!

🔗 BIOINFORMATICS CONNECTION:
├── SNP = Single Nucleotide Polymorphism: most common variant type
├── GWAS: millions of SNPs tested for association with trait (longevity!)
├── Variant Calling: GATK pipeline finds mutations in sequencing data
├── VCF files: store all detected variants with quality scores
├── ClinVar database: clinical significance of variants (pathogenic? benign?)
└── gnomAD: population allele frequencies (is this variant rare or common?)
```

---

## 2.3 CELL SIGNALING PATHWAYS — THE AGING CONTROL SYSTEM
**Why:** Every longevity intervention targets a signaling pathway.

```
THE 6 MUST-KNOW AGING PATHWAYS:

═══════════════════════════════════════════
A. mTOR (Mechanistic Target of Rapamycin):
═══════════════════════════════════════════
├── Master regulator of cell growth & metabolism
├── mTORC1 (complex 1): promotes growth, INHIBITS autophagy
│   ├── Activated by: high amino acids, insulin, growth factors, glucose
│   ├── Downstream: S6K1 → ribosome biogenesis (protein synthesis↑)
│   ├── Inhibits ULK1 → blocks autophagy initiation
│   └── Inhibits TFEB → blocks lysosome biogenesis
├── mTORC2: regulates AKT, cytoskeleton, metabolism
├── AGING: mTOR stays hyperactivated → cells keep building, not recycling
├── INTERVENTION: Rapamycin (allosteric mTOR inhibitor)
│   └── Extends lifespan in yeast, worms, flies, mice → MOST reproducible result!
└── Key crosstalk: AMPK INHIBITS mTORC1 (they're antagonistic!)

═══════════════════════════════════════════
B. AMPK (AMP-activated Protein Kinase):
═══════════════════════════════════════════
├── Energy sensor: activated when AMP:ATP is HIGH (cellular energy LOW)
├── When activated:
│   ├── TURNS ON: autophagy, fatty acid oxidation, glucose uptake, mitophagy
│   └── TURNS OFF: protein synthesis, lipid synthesis, mTORC1
├── AGING: AMPK becomes LESS responsive with age → cannot sense low energy
├── INTERVENTION: Metformin (activates AMPK indirectly, blocks Complex I)
│   └── TAME Trial (Targeting Aging with Metformin) → first real aging clinical trial!
├── Exercise → activates AMPK → #1 reason exercise extends lifespan
└── AMP:ATP → AMPK → phosphorylates TSC2 → inhibits Rheb → inhibits mTORC1

═══════════════════════════════════════════
C. INSULIN/IGF-1 SIGNALING (IIS):
═══════════════════════════════════════════
├── Insulin → receptor → IRS1 → PI3K → PIP2→PIP3 → AKT → mTOR
├── IGF-1 (Insulin-like Growth Factor 1): promotes cell growth
├── AGING: DOWNREGULATING this pathway extends lifespan in ALL organisms tested!
│   ├── C. elegans daf-2 mutant (IGF-1R) lives 2× longer!
│   ├── Small dogs (low IGF-1) live longer than large dogs
│   └── Laron syndrome (GH receptor deficiency): extremely long-lived, cancer-resistant
├── INTERVENTION: Caloric restriction → lowers insulin → lowers IGF-1 → longer life
├── AKT phosphorylates FOXO3 → excludes from nucleus → FOXO3 CAN'T activate targets
│   └── Low insulin → FOXO3 stays in nucleus → activates longevity genes!
└── Paradox: chronic high insulin (metabolic syndrome) → accelerated aging!

═══════════════════════════════════════════
D. SIRTUIN PATHWAY (SIRT1-7):
═══════════════════════════════════════════
├── NAD⁺-dependent deacetylases & ADP-ribosyltransferases
├── Require NAD⁺ as co-factor → NAD⁺ declines 50% age 20→80!
├── SIRT1 (nuclear): metabolism, DNA repair, FOXO activation, NF-κB inhibition
├── SIRT2 (cytoplasmic): cell cycle, myelination
├── SIRT3 (mitochondrial): ETC Complex I regulation, antioxidant (MnSOD)
├── SIRT4 (mitochondrial): lipid metabolism
├── SIRT5 (mitochondrial): succinylation, malonylation regulation
├── SIRT6 (nuclear): telomere maintenance, DNA repair, OVEREXPRESSION → longer life!
│   └── SIRT6 KO mice die at 4 weeks (premature aging!)
├── SIRT7 (nucleolar): rRNA transcription
├── AGING: NAD⁺↓ → sirtuins less active → everything downstream degrades
└── INTERVENTION:
    ├── NMN (nicotinamide mononucleotide): NAD⁺ precursor
    ├── NR (nicotinamide riboside): another NAD⁺ precursor
    └── CD38 inhibitors: CD38 enzyme degrades NAD⁺ → block it = more NAD⁺!

═══════════════════════════════════════════
E. p53/p21/p16 PATHWAY (DNA Damage Response):
═══════════════════════════════════════════
├── DNA damage → ATM/ATR kinases detect → phosphorylate p53
├── p53 stabilizes (normally degraded by MDM2) → enters nucleus
├── p53 transcribes:
│   ├── p21 (CDKN1A): CDK inhibitor → cell cycle arrest → TIME to repair
│   ├── BAX: pro-apoptotic → if damage too severe → cell death
│   ├── GADD45: DNA repair genes
│   └── 14-3-3σ: G2/M arrest
├── If damage repaired → p53 drops → cell resumes cycle
├── If damage chronic → permanent p21 expression → SENESCENCE
├── p16INK4a → CDK4/6 inhibition → Rb stays active → E2F blocked → SENESCENCE
│   └── p16 accumulates with age → THE biomarker for senescent cells
└── AGING: chronic DNA damage → chronic DDR → senescence → SASP → tissue damage

═══════════════════════════════════════════
F. NF-κB (Inflammaging):
═══════════════════════════════════════════
├── Master regulator of inflammatory response
├── Components: p65/RelA + p50 heterodimer (most common form)
├── Normally held in cytoplasm by IκB (inhibitor)
├── ACTIVATION: ROS, TNF-α, IL-1β, LPS, DNA damage → IKK phosphorylates IκB → degraded
│   → NF-κB translocates to nucleus
├── Target genes: IL-6, TNF-α, IL-1β, IL-8, COX-2, iNOS, MMP-9
│   └── THESE are the SASP factors!
├── AGING: NF-κB constitutively active → INFLAMMAGING (Hallmark #11)
│   └── Even without infection → baseline inflammation elevated in elderly
└── INTERVENTION:
    ├── Senolytics: remove SASP-producing cells → NF-κB activity drops
    ├── Anti-inflammatory drugs: aspirin, but side effects for chronic use
    └── SIRT1 directly inhibits NF-κB → NAD⁺ restoration → anti-inflammatory

🔗 BIOINFORMATICS CONNECTION:
├── KEGG Pathway Database: all pathways mapped with genes
├── Gene Set Enrichment Analysis (GSEA): which pathways enriched in your data?
├── STRING database: protein-protein interactions in these pathways
├── Drug-target databases: DrugBank, ChEMBL → map drugs to pathway nodes
├── Network analysis: build signaling networks → find intervention points
└── ODE modeling: simulate pathway dynamics (Systems Biology Module 14)
```

---

## 2.4 EPIGENETICS — THE AGING CLOCK
**Why:** Epigenetics IS the aging clock. Reprogramming = epigenetic reset.

```
═══════════════════════════════════════════
A. DNA METHYLATION:
═══════════════════════════════════════════
├── Addition of CH₃ group to Cytosine (5-methylcytosine, 5mC)
├── Mostly at CpG dinucleotides (C followed by G)
│   └── Human genome: ~28 million CpGs, ~70-80% are methylated
├── Methylated CpG → gene SILENCED (blocks TF access)
├── Writers: DNMT1 (maintenance), DNMT3A/B (de novo)
│   ├── DNMT1: copies methylation pattern during replication (faithful copy)
│   └── DNMT3A/B: add new methylation marks (development, reprogramming)
├── Erasers: TET1/2/3 enzymes → 5mC → 5hmC → 5fC → 5caC → unmethylated
│   └── α-ketoglutarate is TET cofactor → metabolism → epigenetics link!
├── CpG islands: GC-rich regions near gene promoters (~60% of human genes)
│   ├── Normally UNmethylated → gene active
│   └── Aberrant methylation → gene silencing (cancer, aging)
│
├── AGING PATTERN:
│   ├── Global hypomethylation: overall methylation decreases with age
│   ├── Focal hypermethylation: specific CpG islands GET methylated → silence genes
│   ├── Methylation DRIFT: pattern becomes increasingly noisy/random
│   └── These changes are HIGHLY PREDICTABLE → aging clocks!
│
├── AGING CLOCKS (most powerful aging biomarker):
│   ├── Horvath Clock (2013): 353 CpG sites → biological age across ALL tissues
│   │   └── Accuracy: ±3.6 years (measures DNA from blood, saliva, etc.)
│   ├── Hannum Clock (2013): 71 CpGs → blood-specific, more accurate for blood
│   ├── PhenoAge (2018): optimized for mortality prediction
│   ├── GrimAge (2019): predicts time to death, incorporates smoking/proteins
│   ├── DunedinPACE (2022): measures RATE of aging (not absolute age)
│   └── YOUR PROJECT: implement clock using Python + methylation array data
│       → ElasticNet regression on 450K/EPIC array data!
│
└── WHY: Epigenetic reprogramming (OSKM) RESETS methylation → age reversal!
    After Yamanaka reprogramming, Horvath clock shows younger age!

═══════════════════════════════════════════
B. HISTONE MODIFICATIONS:
═══════════════════════════════════════════
├── Histones (H2A, H2B, H3, H4): N-terminal tails protrude → get modified
│
├── KEY MARKS TO KNOW (naming: H3K4me3 = Histone3, Lysine4, trimethylation):
│   ├── H3K4me3: Active promoters → gene ON
│   ├── H3K27me3: Polycomb-repressed → gene OFF (facultative heterochromatin)
│   ├── H3K9me3: constitutive heterochromatin (centromeres, transposon silencing)
│   ├── H3K27ac: active enhancers → gene ON
│   ├── H3K36me3: actively transcribed gene bodies
│   ├── H4K16ac: higher-order chromatin decompaction
│   │   └── Lost during aging → chromatin becomes less organized!
│   └── H3K9ac: active transcription start sites
│
├── Enzymes (remember from Tier 1.3 amino acid PTMs):
│   ├── Writers: HATs (acetyltransferases: CBP/p300), HMTs (methyltransferases: EZH2, MLL)
│   ├── Erasers: HDACs (deacetylases: HDAC1-11 + Sirtuins!), HDMs (LSD1, JMJD)
│   └── Readers: Bromodomain (reads acetyl), Chromodomain (reads methyl), Tudor domain
│
├── Bivalent Domains: H3K4me3 + H3K27me3 on SAME promoter
│   └── "Poised" state: gene ready to activate OR repress → embryonic stem cells!
│   └── Lost during differentiation → aging disrupts these marks!
│
└── AGING:
    ├── H3K9me3 boundaries erode → transposons reactivate → genomic instability!
    ├── H3K27me3 patterns redistribute → wrong genes silenced/activated
    └── Global loss of H4K16ac → chromatin disorganization

═══════════════════════════════════════════
C. CHROMATIN REMODELING:
═══════════════════════════════════════════
├── ATP-dependent remodeler families:
│   ├── SWI/SNF (BAF): tumor suppressor, opens chromatin → gene activation
│   ├── ISWI: nucleosome spacing, maintains regular arrays
│   ├── CHD: NuRD complex, deacetylation + remodeling → gene silencing
│   └── INO80: DNA repair, variant histone exchange
├── Variant histones: H3.3, H2A.Z, macroH2A → replace canonical histones
│   └── H2A.Z at promoters → transcription regulation
└── AGING: remodeling becomes less efficient → gene expression noise increases

═══════════════════════════════════════════
D. 3D GENOME ORGANIZATION:
═══════════════════════════════════════════
├── TADs (Topologically Associating Domains): ~100kb-1Mb loops
│   └── Enhancers only contact promoters within same TAD
├── CTCF + Cohesin: form TAD boundaries → loop extrusion model
├── A/B Compartments:
│   ├── A (active): euchromatin, gene-rich → early replicating
│   └── B (inactive): heterochromatin, gene-poor → late replicating
├── Lamina-Associated Domains (LADs): heterochromatin at nuclear periphery
│   └── Loss of LADs → gene derepression in aging!
└── AGING: TAD boundaries weaken → enhancers contact wrong genes!
    Hi-C studies show topology reorganization in aged cells

🔗 BIOINFORMATICS CONNECTION:
├── WGBS/RRBS: whole genome / reduced representation bisulfite sequencing → methylation
├── Infinium 450K/EPIC Array: measure ~850K CpGs → clock data (cheaper than WGBS)
├── ChIP-seq: histone mark mapping genome-wide
├── ATAC-seq: open chromatin (transposase-based)
├── Hi-C / Micro-C: genome-wide 3D contact maps
├── CUT&RUN / CUT&TAG: lower input alternatives to ChIP-seq
└── All feed into Module 10 (Epigenomics) and Module 13 (Multi-Omics) of DEEP SYLLABUS
```

---

## 2.5 MICROBIOLOGY ESSENTIALS — THE MISSING PIECE
**Why:** All cloning uses bacteria. CRISPR evolved from bacterial immunity. Microbiome = aging.

```
⚠️ GAP FIX: This was completely missing from the original roadmap.

═══════════════════════════════════════════
A. THREE DOMAINS OF LIFE:
═══════════════════════════════════════════
├── Bacteria: prokaryotes, no nucleus, circular genome, peptidoglycan wall
├── Archaea: prokaryotes, extremophiles, different membrane lipids
└── Eukarya: eukaryotes (protists, fungi, plants, animals)
    └── YOU are Eukarya, but you USE Bacteria for genetic engineering!

═══════════════════════════════════════════
B. E. COLI — YOUR LAB WORKHORSE:
═══════════════════════════════════════════
├── Why E. coli? fast growth (20 min doubling), well-characterized, cheap
├── Genome: 4.6 Mb, circular, ~4300 genes
├── Growth curve: lag → exponential (log) → stationary → death phase
├── Transformation: heat shock (42°C, 45 sec) or electroporation → uptake plasmid
├── Selection: antibiotic plates → only transformed bacteria survive
├── KEY STRAINS:
│   ├── DH5α: cloning strain (high transformation efficiency)
│   ├── BL21(DE3): expression strain (T7 promoter system)
│   └── TOP10: general purpose cloning
└── Expression systems hierarchy:
    ├── E. coli: fast, cheap, no PTMs (can't glycosylate!)
    ├── Yeast (S. cerevisiae, Pichia): some PTMs, eukaryotic, slow
    ├── Insect cells (Sf9): better PTMs, medium
    └── Mammalian (HEK293, CHO): best PTMs, expensive, slow → for therapeutics

═══════════════════════════════════════════
C. BACTERIAL GENE REGULATION — THE OPERON:
═══════════════════════════════════════════
├── Operon: cluster of genes under ONE promoter (bacteria only)
│   ├── Structure: Promoter → Operator → [Gene1, Gene2, Gene3...]
│   └── Polycistronic mRNA: one mRNA → multiple proteins (unlike eukaryotes!)
│
├── Lac Operon (THE first genetic circuit ever studied!):
│   ├── lacZ: β-galactosidase (breaks down lactose)
│   ├── lacY: permease (imports lactose)
│   ├── lacA: transacetylase
│   ├── LacI (repressor): binds operator → BLOCKS transcription
│   ├── When lactose present: allolactose binds LacI → LacI falls off → genes ON
│   ├── CAP/cAMP: when glucose LOW → cAMP high → CAP binds → activates transcription
│   └── LOGIC: Gene ON when (lactose present) AND (glucose absent)
│       → This IS a genetic AND gate! (same concept as synthetic biology!)
│
├── Trp Operon (negative feedback):
│   ├── Genes for tryptophan synthesis
│   ├── When Trp HIGH: Trp binds repressor → repressor binds operator → genes OFF
│   └── When Trp LOW: repressor inactive → genes ON → more Trp made
│
└── WHY THIS MATTERS:
    ├── Synthetic biology GE.4 uses SAME LOGIC from operons → gene circuits
    ├── IPTG/LacI system → widely used inducible promoter in research!
    └── Understanding operons = understanding how to design expression systems

═══════════════════════════════════════════
D. CRISPR — A BACTERIAL IMMUNE SYSTEM:
═══════════════════════════════════════════
├── CRISPR evolved as adaptive immunity against PHAGES (viruses)
├── How bacteria use CRISPR:
│   ├── Phage infects → bacteria cuts phage DNA → stores fragment in CRISPR array
│   ├── CRISPR array: spacers (phage DNA) separated by repeats
│   ├── On re-infection: spacer transcribed → crRNA → guides Cas to phage DNA → cuts!
│   └── THIS is what Doudna & Charpentier adapted for gene editing!
│
├── Phage biology (brief):
│   ├── Bacteriophages: viruses that infect bacteria
│   ├── Lytic cycle: phage hijacks cell → makes copies → cell bursts
│   ├── Lysogenic cycle: phage DNA integrates → prophage → dormant
│   └── Phage therapy: using phages to kill antibiotic-resistant bacteria (re-emerging!)
│
└── Horizontal Gene Transfer (HGT):
    ├── Conjugation: direct cell-to-cell transfer via pilus (F-plasmid)
    ├── Transformation: uptake of naked DNA from environment
    ├── Transduction: phage carries bacterial DNA to new host
    └── WHY: antibiotic resistance spreads via HGT → global health crisis!

═══════════════════════════════════════════
E. MICROBIOME & AGING:
═══════════════════════════════════════════
├── Human microbiome: ~38 trillion bacteria (roughly equal to human cells!)
├── Gut microbiome: >1000 species, majority in large intestine
├── Functions: digestion, vitamin synthesis, immune training, metabolism
├── AGING:
│   ├── Microbiome diversity DECREASES with age
│   ├── Inflammatory species INCREASE (Firmicutes↑, Bacteroidetes↓)
│   ├── Gut barrier integrity declines → "leaky gut" → systemic inflammation
│   └── Fecal microbiota transplant (FMT) from young → old mice extends lifespan!
├── 🔗 MODULE 26 in DEEP SYLLABUS: metagenomics analysis
└── 16S rRNA sequencing: identify bacterial species (QIIME2 pipeline)
```

---

## 2.6 PHYLOGENETICS & MOLECULAR EVOLUTION — COMPARE ACROSS SPECIES
**Why:** Finding longevity genes requires comparing genomes across species. GWAS needs evolutionary context.

```
⚠️ GAP FIX: This was completely missing from the original roadmap.

═══════════════════════════════════════════
A. SEQUENCE HOMOLOGY:
═══════════════════════════════════════════
├── Homologs: genes with shared evolutionary ancestor
│   ├── Orthologs: homologs separated by SPECIATION (same function likely)
│   │   └── Human p53 ↔ Mouse Trp53 are orthologs (same function in both)
│   └── Paralogs: homologs separated by GENE DUPLICATION (may diverge in function)
│       └── Human SIRT1 and SIRT2 are paralogs
│
├── Sequence identity vs similarity:
│   ├── Identity: exact same residue at that position
│   └── Similarity: chemically similar residue (conservative substitution)
│       → This is where BLOSUM62 and PAM250 matrices come in!
│
├── Multiple Sequence Alignment (MSA):
│   ├── Aligns 3+ sequences simultaneously to find conserved regions
│   ├── Tools: ClustalOmega, MUSCLE, MAFFT, T-Coffee
│   └── Output: consensus sequence + conservation scores
│       → Highly conserved = functionally important = good CRISPR target!
│
└── Databases:
    ├── OrthoMCL, OrthoDB: ortholog databases
    ├── KEGG Orthology: functional annotation of orthologs
    └── Ensembl Compara: multi-species genome comparison

═══════════════════════════════════════════
B. PHYLOGENETIC TREES:
═══════════════════════════════════════════
├── What: branching diagrams showing evolutionary relationships
├── Methods:
│   ├── Distance-based:
│   │   ├── UPGMA: simplest, assumes molecular clock (constant rate)
│   │   └── Neighbor-Joining (NJ): faster, no clock assumption → commonly used!
│   ├── Character-based:
│   │   ├── Maximum Parsimony (MP): fewest changes → most likely tree
│   │   ├── Maximum Likelihood (ML): statistically most likely tree given data + model
│   │   │   └── RAxML, IQ-TREE: standard tools
│   │   └── Bayesian: posterior probability of tree → MrBayes, BEAST
│   └── Choosing method: ML or Bayesian for publication, NJ for quick exploration
│
├── Key concepts:
│   ├── Bootstrap: resample alignment → how confident is each branch? (>70% = reliable)
│   ├── Outgroup: distant species to root the tree
│   ├── Molecular clock: mutations accumulate at roughly constant rate → date divergence
│   └── Rooted vs Unrooted trees
│
└── Tools: MEGA, FigTree (visualization), iTOL (interactive Tree of Life)

═══════════════════════════════════════════
C. SELECTION & EVOLUTION OF AGING:
═══════════════════════════════════════════
├── dN/dS (ω ratio): non-synonymous / synonymous substitution rate
│   ├── ω < 1: purifying selection (gene is conserved → important function!)
│   ├── ω = 1: neutral evolution (no selection pressure)
│   └── ω > 1: positive selection (gene evolving faster than expected → adaptation!)
│       → Use PAML codeml program to detect
│
├── AGING APPLICATION:
│   ├── Compare longevity genes across species:
│   │   ├── Naked mole rat: 30+ year lifespan (vs 2-3 yrs for similar-sized mice)
│   │   │   → BETTER DNA repair, DIFFERENT p53 responses, HIGH hyaluronic acid
│   │   ├── Bowhead whale: 200+ years → largest mammal, slow life history
│   │   │   → ERCC1 gene variants, unique p53 regulation
│   │   ├── Turritopsis dohrnii: "immortal jellyfish" → reverts to polyp stage!
│   │   └── Hydra vulgaris: continuous stem cell replacement → unlimited lifespan
│   ├── Genes under purifying selection across long-lived species = validated targets
│   └── Comparative genomics approach: what's unique in long-lived species?
│
└── WHY:
    ├── Identifies WHICH longevity genes are most important (conserved across evolution)
    ├── CRISPR targets should ideally be in conserved regions (less off-target risk)
    ├── Protein structure comparison: conserved domains = functional cores
    └── Module 17 (Comparative Genomics) in DEEP SYLLABUS

🔗 BIOINFORMATICS CONNECTION:
├── BLAST: compare your sequence to database → find homologs
├── Ensembl / UCSC Genome Browser: multi-species genome alignment
├── PhyloP, phastCons: per-base conservation scores (genome browser tracks)
├── Synteny analysis: are genes in same order across species?
└── Rosalind has several phylogenetics problems (tree distances, character tables)
```

---

# ✅ TIER 2 COMPLETION CHECKLIST

```
═══════════════════════════════════════════════════════════════════════
 SECTION 2.1 — PROTEIN STRUCTURE
═══════════════════════════════════════════════════════════════════════
 □ Can describe all 4 levels of protein structure with examples
 □ Know α-helix (3.6 residues/turn) and β-sheet (parallel/antiparallel)
 □ Can explain hydrophobic collapse as main folding driving force
 □ Know what protein domains are and can name 3 (SH2, bromodomain, kinase domain)
 □ Know: bromodomain reads acetylation, chromodomain reads methylation
 □ Can list 3 aging-relevant chaperones (HSP70, HSP90, HSP60)
 □ Know the proteasome pathway: ubiquitin → E1/E2/E3 → proteasome → peptides
 □ Know AlphaFold predicts 3D structure from sequence
 □ Know: misfolded proteins → aggregation → amyloid → Alzheimer's

═══════════════════════════════════════════════════════════════════════
 SECTION 2.2 — GENETICS
═══════════════════════════════════════════════════════════════════════
 □ Understand dominant/recessive, genotype/phenotype, Hardy-Weinberg
 □ Can explain epistasis, pleiotropy, polygenic traits, linkage disequilibrium
 □ Know all mutation types: silent, missense, nonsense, frameshift, splice site
 □ Know transition (purine↔purine) vs transversion (purine↔pyrimidine)
 □ Can explain somatic vs germline mutations and relevance to aging
 □ Know all 5 DNA repair pathways (BER, NER, MMR, HR, NHEJ) and when each is used
 □ Know: NHEJ repairs CRISPR cuts → causes INDELs for knockouts
 □ Know: Werner syndrome, Cockayne syndrome = DNA repair → premature aging
 □ Can explain what a GWAS study does and what VCF files contain
 □ Know: Polygenic Risk Score (PRS) = sum of many small SNP effects

═══════════════════════════════════════════════════════════════════════
 SECTION 2.3 — CELL SIGNALING
═══════════════════════════════════════════════════════════════════════
 □ Can trace: mTOR pathway from nutrients → mTORC1 → S6K1 / ULK1 inhibition
 □ Know: rapamycin inhibits mTOR → extends lifespan in ALL organisms
 □ Can explain AMPK: energy sensor, activated by low AMP:ATP, inhibits mTORC1
 □ Know: metformin activates AMPK → TAME trial
 □ Can trace: insulin → receptor → PI3K → AKT → mTOR
 □ Know the Sirtuin family (SIRT1-7), their locations, and NAD⁺ dependency
 □ Can explain p53 → p21 → senescence pathway
 □ Know NF-κB pathway and how chronic activation → inflammaging
 □ Can explain the crosstalk: AMPK ⊣ mTOR, Insulin → AKT ⊣ FOXO3
 □ Know: KEGG, GSEA, STRING are tools for pathway analysis

═══════════════════════════════════════════════════════════════════════
 SECTION 2.4 — EPIGENETICS
═══════════════════════════════════════════════════════════════════════
 □ Can explain DNA methylation: 5mC at CpGs, DNMT writers, TET erasers
 □ Know at least 3 Horvath-style aging clocks and what makes each different
 □ Can name 5 histone marks and classify as activating or repressing
 □ Know: bromodomain = acetyl reader, chromodomain = methyl reader
 □ Can explain bivalent domains (H3K4me3 + H3K27me3) and their role in stem cells
 □ Know: TADs organize genome, CTCF+cohesin form boundaries
 □ Can explain how aging disrupts: methylation drift, heterochromatin erosion, TAD weakening
 □ Know: Yamanaka reprogramming resets Horvath clock → biological age reversal
 □ Know which assays measure what: WGBS (methylation), ChIP-seq (histones), ATAC-seq (accessibility), Hi-C (3D)

═══════════════════════════════════════════════════════════════════════
 SECTION 2.5 — MICROBIOLOGY
═══════════════════════════════════════════════════════════════════════
 □ Know 3 domains of life and key differences (Bacteria, Archaea, Eukarya)
 □ Know WHY E. coli is used for cloning (fast, cheap, well-characterized)
 □ Can explain the Lac operon: LacI repressor, allolactose inducer, CAP/cAMP
 □ Understand: Lac operon = biological AND gate → foundation of synthetic biology
 □ Know: CRISPR literally evolved as bacterial adaptive immunity against phages
 □ Can explain horizontal gene transfer: conjugation, transformation, transduction
 □ Know microbiome aging pattern: diversity↓, inflammation↑, barrier integrity↓
 □ Know: 16S rRNA sequencing identifies bacterial species

═══════════════════════════════════════════════════════════════════════
 SECTION 2.6 — PHYLOGENETICS
═══════════════════════════════════════════════════════════════════════
 □ Can define orthologs vs paralogs with examples
 □ Know MSA tools: ClustalOmega, MUSCLE, MAFFT
 □ Can explain 3 tree-building methods: NJ, ML, Bayesian
 □ Know what bootstrap values mean (>70% = reliable branch)
 □ Can explain dN/dS ratio: <1 (purifying), =1 (neutral), >1 (positive selection)
 □ Know 3 long-lived species and what makes them special
   (naked mole rat, bowhead whale, Hydra)
 □ Know: conserved regions across species = functionally important = good CRISPR targets
 □ Know BLAST as primary tool for finding homologous sequences

═══════════════════════════════════════════════════════════════════════
 TOTAL: 60 items. All checked = TIER 2 COMPLETE → proceed to TIER 1B (GE)
═══════════════════════════════════════════════════════════════════════
```
