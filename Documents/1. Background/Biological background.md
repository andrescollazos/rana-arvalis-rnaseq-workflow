# Thesis Background and Introduction – Outline

## 1. Ecological and evolutionary framework: thermal adaptation along latitudinal gradients

### 1.1 Latitudinal gradients as drivers of phenotypic divergence
- Temperature and season length vary predictably with latitude
- High-latitude environments impose strong time constraints on development
- Ectotherms are especially sensitive to ambient temperature during early life stages

### 1.2 Phenotypic plasticity and local adaptation
- Definitions and distinction:
  - phenotypic plasticity
  - genetic adaptation
- Plasticity as:
  - a short-term buffer to environmental variation
  - a potential facilitator or constraint on adaptive evolution
- Genotype-by-environment interactions at the phenotypic level

### 1.3 Countergradient and cogradient variation
- Conceptual definitions
- Relationship between genetic effects and environmental effects
- Relevance for interpreting latitudinal life-history patterns
- Importance of common garden experiments for disentangling causes

### 1.4 Post-glacial recolonization and historical contingency
- European post-glacial recolonization routes
- Demographic history as a determinant of present-day divergence
- Interaction between colonization history and selection along gradients

### 1.5 The moor frog (*Rana arvalis*) as a study system
- Geographic distribution across Europe
- Replicated latitudinal gradient (Latvia/Kalmar → Norrbotten)
- Life-history traits shaped by thermal environment
- Summary of findings from Papers 1–4:
  - divergence in developmental rate and life-history traits
  - scale dependence of divergence (local vs regional vs historical)
  - alignment or mismatch between plasticity and population divergence

**Purpose of this section:**  
Justifies why thermal adaptation along latitudinal gradients is a biologically meaningful problem in *Rana arvalis*.

### Reading list

**Papers 1–4 (Laurila group, *Rana arvalis*)**  
*Primary, non-negotiable*

**Focus on:**

- exact geographic range of sampled populations (Latvia/Kalmar → Norrbotten)
- how latitude maps onto temperature and season length
- which traits diverge (developmental rate, growth, life-history timing)
- which patterns persist under common garden conditions
- how divergence changes with spatial scale (local vs regional vs historical)
- where plasticity aligns with divergence and where it does not
- explicit statements about colonization history and its role

Do **not** focus on:

- statistical minutiae
- effect size values unless used to make a biological argument

**Mallick 2022 MSc thesis (countergradient variation in *R. arvalis*)**

**Focus on:**

- clean definitions of countergradient vs cogradient variation
- how time constraints at high latitude are framed biologically
- how common garden logic is explained pedagogically
- how recolonization routes are described and linked to phenotypic patterns
- how results are discussed in ecological, not molecular, terms

Do **not** focus on:

- R code
- detailed statistical tests

---

## 2. From phenotype to molecular regulation

### 2.1 Limits of phenotypic data alone
- Phenotypic convergence can mask underlying regulatory divergence
- Similar traits may arise from different molecular mechanisms
- Need to examine regulatory responses to temperature directly

### 2.2 Gene regulation as a target of selection
- Regulation as a primary mode of adaptive change
- Expression patterns as regulatory phenotypes
- Population-associated vs plastic regulatory responses

### 2.3 Genotype-by-environment interactions at the molecular level
- Expression-level G×E
- Plastic responses vs population-associated expression differences
- Evolution of regulatory robustness and plasticity

### 2.4 Cis- and trans-regulatory variation
- Conceptual distinction
- Relative stability and environmental sensitivity
- Implications for adaptive divergence along thermal gradients

### 2.5 Tissue-specific regulatory responses
- Liver as a central metabolic and stress-responsive tissue
- Relevance for thermal adaptation and developmental processes

**Purpose of this section:**  
Establishes gene regulation as the mechanistic link between environment, phenotype, and adaptation.

### Reading list

**Ballinger et al. 2023, PNAS (house mouse)**  
*Comparative anchor*

**Focus on:**

- conceptual separation of cis vs trans regulation
- robustness vs plasticity of expression across temperatures
- how expression divergence is linked to adaptation
- how G×E is interpreted at the regulatory level
- how tissue-specific responses are justified biologically

Read this as:
- a **conceptual model**, not a system to replicate

Do **not** focus on:

- mouse-specific phenotypes
- selection scan technical details

**Xu et al. 2025, Journal of Thermal Biology (frog, liver, thermal stress)**

**Focus on:**

- how temperature treatments are biologically motivated
- why liver is chosen as the focal tissue
- how population-specific transcriptional responses are described
- how plastic vs population-associated expression is distinguished
- how results are framed in terms of thermal resilience and adaptation

Use this paper to:
- calibrate expectations for what expression responses to temperature look like in amphibians

Do **not** focus on:

- software versions
- enrichment plots beyond their biological interpretation

---

## 3. Reference-based RNA-seq as a framework for ecological genomics

### 3.1 Regulatory signals accessible through reference-based analysis
- Population-level comparability across samples
- Integration of:
  - gene expression
  - exon usage
  - transposable element expression

### 3.2 Structure of regulatory contrasts
- Temperature effects
- Population-origin effects
- Temperature × population interactions

### 3.3 Statistical logic of expression analysis
- Count-based modeling
- Variance and dispersion as biological signals
- Interpretation limits of differential expression

### 3.4 Role of differential expression in ecological inference
- Identifying candidate regulatory responses
- Generating hypotheses rather than causal claims

**Purpose of this section:**  
Justifies how regulatory variation is detected and interpreted in a population and environmental context.

## Reading list

**DESeq2 original paper (Love et al. 2014)**

**Focus on:**

- biological interpretation of count-based models
- why variance matters biologically
- how contrasts map to experimental design
- limits of differential expression inference
- normalization and variance modeling (negative binomial logic)  

Do **not** focus on:

- implementation details

**Any review on RNA-seq in ecological or evolutionary genomics**  
(e.g., gene expression in natural populations)

**Focus on:**

- what expression differences can and cannot say about adaptation
- common pitfalls in ecological interpretation
- how population structure affects inference

Stop reading once:
- conceptual limits are clear

---

## 4. Co-expression and regulatory architecture

### 4.1 Polygenic and coordinated nature of thermal adaptation
- Complex traits arise from coordinated gene activity
- Single-gene approaches are insufficient

### 4.2 Co-expression modules as regulatory units
- Conceptual meaning of a module
- Modules as emergent properties of regulation
- Relationship between modules and environmental or population effects

### 4.3 Network-based perspectives on adaptation
- Regulatory architecture rather than individual genes
- Population differences in network structure
- Environmentally responsive modules

### 4.4 Limits and interpretation of co-expression
- Correlation vs causation
- Biological interpretation of network structure

### 4.5 Causal inference layered on co-expression
- Conceptual meaning of causality in regulatory networks
- Why causal frameworks are attractive for adaptive hypotheses
- Positioning of WGCNA-based and alternative network approaches

**Purpose of this section:**  
Explains why network-level analyses are required to understand regulatory responses to temperature.

## Reading list (Network-level regulation)

**WGCNA original paper or tutorial review**

**Focus on:**

- what a co-expression module represents biologically
- why modules are interpreted as regulatory units
- how trait or environment associations are used
- explicit caveats about correlation-based inference

Do **not** focus on:

- adjacency functions
- parameter tuning

**Paper on causal inference layered on WGCNA output**  
(e.g., the NAR Genomics and Bioinformatics paper you were given)

**Focus on:**

- what “causal” means in this context
- how causality is inferred from observational data
- what biological claims are allowed vs not allowed
- why this approach is proposed as an extension, not a replacement

**GWENA documentation or methods paper**

**Focus on:**

- how it conceptually differs from WGCNA
- what biological questions it claims to answer
- when a simpler formulation is advantageous

---

## 5. Transposable elements as components of the regulatory response

### 5.1 Transposable elements in eukaryotic genomes
- Abundance and genomic distribution
- Relationship to gene regulation

### 5.2 Environmental stress and TE expression
- Stress-induced transcriptional activation
- Temperature as a trigger for TE expression changes

### 5.3 TEs as regulatory elements
- Influence on nearby gene expression
- Contribution to regulatory innovation
- Population and environment-dependent TE activity

### 5.4 TE expression as part of thermal adaptation
- TE transcription vs mobilization
- TE activity as a component of regulatory plasticity
- Integration of TE expression with gene-level and network-level analyses

**Purpose of this section:**  
Places TE expression within the same regulatory framework as genes and co-expression networks.

## Reading list (TE biology and stress responses)

**Review on transposable elements and environmental stress**  
(any recent TE regulation review)

**Focus on:**

- stress-induced TE transcription
- distinction between transcription and mobilization
- TE-derived regulatory elements
- population and environment-specific TE activity

Do **not** focus on:

- TE classification schemes
- evolutionary history of TE families unless tied to regulation

**TE expression method papers (TEtranscripts, Telescope)**

**Focus on:**

- what signal is being measured (TE transcription)
- how TE expression is disentangled from gene expression
- biological interpretation of increased TE expression

Do **not** focus on:

- algorithmic internals


---

## 6. Integrated regulatory perspective on thermal adaptation

- Thermal adaptation as a multi-layered regulatory process
- Interaction between:
  - gene expression
  - co-expression architecture
  - transposable element activity
- Relevance for understanding population divergence under climate change
- Positioning of the present thesis within this framework

**Purpose of this section:**  
Leads directly into the thesis aims and research questions.
