# Thermal Plasticity and Adaptation Along a Latitudinal Gradient in Moor Frog (Rana arvalis)

## Description

This repository contains the RNA-seq analysis workflows and scripts developed for the MSc thesis:

**“Thermal Plasticity and Adaptation Along a Latitudinal Gradient in Moor Frog (Rana arvalis)”**

The project investigates how thermal plasticity and population divergence are reflected at the transcriptional regulatory level. Using liver RNA-seq data from multiple populations sampled along a latitudinal gradient and reared under two temperature treatments, the analyses focus on:

- Differential gene expression  
- Gene co-expression network structure  
- Expression of transposable elements  

The aim is to characterize regulatory variation associated with temperature, population origin, and their interaction.

---

## Workflow

The analysis follows a sequential reference-based RNA-seq workflow:

1. **Read alignment**  
   Paired-end FASTQ files are aligned to the *Rana arvalis* reference genome using STAR.

2. **Quantification**  
   Gene-level counts are generated with featureCounts to produce a count matrix.

3. **Differential expression analysis**  
   DESeq2 is used to model:
   - Temperature effects  
   - Population-origin effects  
   - Temperature × population interactions  

4. **Co-expression network analysis**

5. **Transposable element expression analysis**

