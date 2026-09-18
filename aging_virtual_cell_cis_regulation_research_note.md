# Aging, AlphaGenome, cis-Regulation, and Virtual Cell

> Working note for a proposed research direction: **Virtual Cell (ML) + scRNA-seq + Stereo-seq for precision aging and longevity**.  
> Updated: 2026-09-09. URLs are included so this file can be used directly in VS Code / a reference manager workflow.

> **Project-specific decision for the current Repro-CM/HDF study:** use the existing bulk RNA-seq only for exploratory splice-candidate discovery; do not use AlphaGenome in the core paper because no DNA-variant question is present; consider a narrowly defined “virtual fibroblast” only after generating replicated, matched perturbation data and validating unseen predictions prospectively. Stereo-seq becomes relevant only after expansion to a 3D skin, explant, or in vivo niche. See `aging_repro_cm_paper_strategy_review_ko.md` for the full logic audit, experimental design, and figure-by-figure paper plan.

## 1. One-sentence research framing

Build a cell- and tissue-context-aware model that starts with DNA sequence/genotype, molecular cell state, and spatial niche, then predicts how a genetic, chemical, or environmental intervention changes an aging-related cellular state.

\[
P(\text{future cell state} \mid \text{genotype/sequence},\ \text{initial state},\ \text{niche},\ \text{perturbation},\ \text{time})
\]

This is more useful than a static aging atlas only if it can propose and rank interventions that are then experimentally tested.

## 2. Concept map: atlas, AlphaGenome, and Virtual Cell

| Layer | Main question | Useful output | Typical data/model |
|---|---|---|---|
| Genome map / atlas | What is located at this locus and in which cells is it active? | Gene, enhancer, TF, chromatin and cell-type annotation | UCSC, ENCODE, cell atlases |
| AlphaGenome | What molecular **cis-regulatory effect** could a DNA-sequence variant produce? | Predicted change in expression, splicing, chromatin accessibility, TF binding and chromatin contacts | Sequence-to-function deep model |
| Virtual Cell | Given a cell state and an intervention, how will its state change? | Post-perturbation cell-state distribution / response | scRNA-seq, Perturb-seq, foundation and state-transition models |
| Spatial Virtual Cell | How does the same intervention work in its tissue microenvironment? | Niche-conditioned response | Stereo-seq / spatial transcriptomics plus cell-state models |

The layers are complementary rather than substitutes:

`variant -> cis-regulatory delta -> cell-state transition -> aging/function phenotype`.

## 3. cis-regulatory effect

### Definition

A **cis-regulatory effect** is the effect of a DNA sequence feature or variant on regulation of a gene on the same physical DNA molecule/allele. “Cis” does not mean only the immediately adjacent gene: enhancer–promoter contacts can bridge tens of kilobases to megabases in three-dimensional chromatin.

### cis versus trans

| | cis regulation | trans regulation |
|---|---|---|
| Driver | Local promoter, enhancer, insulator or splicing sequence | Diffusible TF, splicing factor, cytokine, pathway state |
| Scope | Usually a local gene or gene set | Potentially many genes on many chromosomes |
| Allele specificity | Often present | Usually affects both alleles similarly |
| Example | A SNP weakens an enhancer motif and lowers one allele's RNA | NF-kB activation induces inflammatory genes genome-wide |

### Major mechanisms

1. **Promoter / enhancer activity**: a sequence change alters transcription-factor binding and the output of a promoter or enhancer.
2. **Chromatin and 3D genome**: effects can be mediated by accessibility, histone modifications, CTCF/insulator sites, and enhancer–promoter contact.
3. **Splicing**: variants in splice donor/acceptor sites, branch points, or exonic/intronic splicing regulatory elements alter exon usage or intron retention. Total gene expression can stay constant while functional isoform balance changes.
4. **Allele-specific regulation**: in a heterozygote, unequal expression/accessibility of the two alleles is especially informative because both alleles share the same trans environment.

### Measurements and causal validation

| Biological question | Common evidence |
|---|---|
| Does genotype alter total RNA? | eQTL, allele-specific expression |
| Does genotype alter isoform use? | sQTL, junction QTL |
| Does genotype alter accessibility? | caQTL, scATAC-seq |
| Does it change TF binding / chromatin marks? | ChIP-seq or CUT&Tag QTL, motif analysis |
| Which gene does an enhancer regulate? | Hi-C/HiChIP, ABC model, CRISPRi |
| Is the exact base causal? | MPRA, base editing, CRISPR perturbation |

Association alone is insufficient: LD may point to a linked, non-causal variant, and enhancer-to-gene assignment can be incorrect without perturbational evidence.

### Why cis effects matter in aging

The DNA variant is fixed, but its observed effect is conditional on cell type, age, TF availability, chromatin, inflammation, and local niche:

\[
\text{observed cis effect}=f(\text{sequence},\ \text{cell type},\ \text{TF abundance},\ \text{chromatin},\ \text{age/niche})
\]

Thus an enhancer variant might be nearly silent in young cells yet become functional in old/inflamed cells after NF-kB, AP-1, C/EBP, STAT, FOXO, or NRF2 programs and chromatin accessibility change. This motivates age-stratified and cell-type-resolved eQTL/sQTL analyses.

## 4. AlphaGenome

AlphaGenome is a sequence-to-function AI model for predicting the molecular consequences of DNA sequence variation. The practical use is to compare predictions from reference and alternate alleles and obtain a **predicted cis-regulatory delta**.

Potential outputs include changes in:

- RNA expression-related tracks
- splice junction / transcript processing signals
- chromatin accessibility and histone marks
- transcription-factor binding
- chromatin contact-related signals

It is a candidate-prioritization model, not proof of causal biology. It should be combined with cell-context data and wet-lab validation.

### Suggested variant-to-phenotype workflow

1. Start from aging/longevity GWAS, rare variants, or regulatory loci of interest.
2. Use AlphaGenome to prioritize variants and possible affected regulatory modalities.
3. Intersect predictions with cell-type-specific scRNA-seq/scATAC-seq and ENCODE cCRE evidence.
4. Connect enhancer candidates to likely target genes using chromatin contact / ABC evidence.
5. Test top candidates with MPRA or base editing, followed by Perturb-seq or functional aging assays.

## 5. What “Virtual Cell” means

An AI Virtual Cell is a model intended to represent and simulate molecular/cellular state across multiple biological contexts and to predict state transitions following intervention.

It is useful to separate three levels:

| Model type | What it does | What it does **not** establish by itself |
|---|---|---|
| Cell foundation/embedding model | Learns a useful representation of cell and gene state | A causal intervention response |
| Perturbation/state-transition model | Predicts transcriptional response to gene/drug/environmental perturbation | Full cellular function, protein activity, or organismal lifespan |
| Spatial model | Adds cell neighbourhood / tissue niche | Dynamic multicellular causal response unless trained on such data |

At present, most models are more accurately called **virtual transcriptomic cells** than full virtual cells. They mainly predict RNA-level states and are not complete simulators of proteins, metabolism, organelles, electrophysiology, or tissue mechanics.

## 6. Representative papers and materials

### Foundational vision and classical whole-cell modelling

1. Karr JR et al. **A Whole-Cell Computational Model Predicts Phenotype from Genotype.** *Cell* (2012).  
   https://doi.org/10.1016/j.cell.2012.05.044  
   - Classical mechanistic whole-cell model of *Mycoplasma genitalium*.
   - Important historical contrast: it explicitly encodes many biological modules, whereas modern AI approaches learn patterns from massive data.

2. Bunne C et al. **How to build the virtual cell with artificial intelligence: Priorities and opportunities.** *Cell* (2024).  
   https://doi.org/10.1016/j.cell.2024.11.015  
   Preprint: https://arxiv.org/abs/2409.11654  
   - Essential conceptual reading for AI Virtual Cell.
   - Proposes multiscale/multimodal representations, virtual instruments, perturbation prediction, interpretability, and experiment–model feedback loops.

### Single-cell foundation / representation models

3. Theodoris CV et al. **Transfer learning enables predictions in network biology.** *Nature* (2023).  
   https://doi.org/10.1038/s41586-023-06139-9  
   - Geneformer; pretrained on roughly 30 million single-cell profiles.
   - Good for gene/cell representation and network-level prediction, but not automatically a complete perturbation simulator.

4. Cui H et al. **scGPT: toward building a foundation model for single-cell multi-omics using generative AI.** *Nature Methods* (2024).  
   https://www.nature.com/articles/s41592-024-02201-0  
   - Generative single-cell foundation model trained on >33 million cells.
   - Supports annotation, batch integration, multiomic integration, GRN inference, and perturbation-related tasks.

5. Rosen Y et al. **Universal cell embedding provides a foundation model for cell biology.** *Nature* (2026).  
   https://www.nature.com/articles/s41586-026-10689-z  
   - UCE: cross-tissue and cross-species universal cell embedding.
   - Strong reference for broad generalizable cell representations, not a direct state-transition simulator.

### Perturbation / state-transition models

6. Lotfollahi M et al. **scGen predicts single-cell perturbation responses.** *Nature Methods* (2019).  
   https://www.nature.com/articles/s41592-019-0494-8  
   - Early variational model for transferring perturbation responses across cell types and settings.

7. Lotfollahi M et al. **Predicting cellular responses to complex perturbations in high-throughput screens.** *Molecular Systems Biology* (2023).  
   https://doi.org/10.15252/msb.202211517  
   - CPA (compositional perturbation autoencoder): predicts responses across doses, cell types, times, species, and combinations.

8. Roohani Y, Huang K, Leskovec J. **Predicting transcriptional outcomes of novel multigene perturbations with GEARS.** *Nature Biotechnology* (2024).  
   https://www.nature.com/articles/s41587-023-01905-6  
   - Uses a gene–gene knowledge graph to predict single and combinatorial genetic perturbations, including previously unseen genes/combinations.
   - Particularly relevant for prioritizing combinations within mTOR, NF-kB, FOXO, autophagy, DNA repair, and senescence networks.

9. Adduri AK et al. **Predicting cellular responses to perturbation across diverse contexts with State.** *Cell* (online 31 August 2026).  
   https://doi.org/10.1016/j.cell.2026.07.052  
   Preprint: https://www.biorxiv.org/content/10.1101/2025.06.26.661135v2  
   Code: https://github.com/ArcInstitute/state  
   - Maps initial cell state plus perturbation to a predicted post-perturbation transcriptomic state/distribution.
   - One of the clearest current implementations of a practical Virtual Cell direction.

### Spatial context and virtual cells

10. **Nicheformer: a foundation model for single-cell and spatial omics.** *Nature Methods* (2025).  
    https://doi.org/10.1038/s41592-025-02814-z  
    - Incorporates cellular microenvironment / neighbourhood information into representation learning.
    - Most directly relevant to adding Stereo-seq information to a cell-state model.

11. Wang Z et al. **NicheTrans: spatial-aware cross-omics translation.** *Nature Methods* (2026).  
    https://doi.org/10.1038/s41592-026-03153-3  
    - Integrates spatial microenvironment information with multimodal translation.

12. Wang J, Huang Y, Winther O. **SpatialFormer: universal spatial representation learning from subcellular molecular to multicellular landscapes.** *Nature Computational Science* (2026).  
    https://doi.org/10.1038/s43588-026-01016-7  
    - Relevant if the project expands from cell-level annotation toward multiscale tissue organisation.

### Aging-specific perturbation resource

13. Zhu S et al. **Multiomic single-cell perturbation screens reveal critical lncRNA regulators of senescence.** *Nature Aging* (2026).  
    https://doi.org/10.1038/s43587-026-01100-7  
    - Direct example of using single-cell multiomic perturbation screens to identify causal regulators of senescence.
    - Useful template for an aging-specific Virtual Cell training/validation dataset.

### Benchmarking and limitations

14. **Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines.** *Nature Methods* (2025).  
    https://www.nature.com/articles/s41592-025-02772-6  
    - Required critical reading: apparent model gains can disappear under robust baselines and difficult out-of-distribution evaluation.

15. Viñas Torné R et al. **Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation.** *Nature Biotechnology* (2026).  
    https://doi.org/10.1038/s41587-025-02777-8  
    - Focuses on evaluation that better isolates perturbation-specific effects.

16. Arc Virtual Cell Initiative Team, Goodarzi H. **Virtual Cell Challenge 2026: Benchmarking zero-shot generalization across cellular contexts.** *Cell* (2026).  
    https://doi.org/10.1016/j.cell.2026.08.004  
    Challenge overview: https://arcinstitute.org/news/virtual-cell-challenge-2026  
    - A useful reality check: unseen-cell-context, zero-shot perturbation prediction remains an open challenge.

### cis-regulatory references and resources

17. de Boer CG, Taipale J. **Hold out the genome: a roadmap to solving the cis-regulatory code.** *Nature* (2024).  
    https://www.nature.com/articles/s41586-023-06661-w  
    - Conceptual and evaluation-oriented reading on learning the cis-regulatory code.

18. ENCODE SCREEN / Registry of candidate cis-regulatory elements (cCREs).  
    https://screen.wenglab.org/about  
    - Search and visualise candidate promoters, enhancers, TF binding, chromatin state, and 3D genome relationships across human/mouse biosamples.

19. Avsec Z et al. **Advancing regulatory variant effect prediction with AlphaGenome.** *Nature* (2026).  
    Paper: https://www.nature.com/articles/s41586-025-10014-0  
    DeepMind overview: https://deepmind.google/blog/alphagenome-ai-for-better-understanding-the-genome/  
    Code and research materials: https://github.com/google-deepmind/alphagenome_research  
    - Sequence-to-function model for prioritising non-coding regulatory variant effects.

## 7. Concrete aging Virtual Cell study design

### Minimal viable research aim

> Predict which genetic or pharmacological perturbations move specific aged cell populations toward a pre-defined healthier molecular and functional state, conditional on cell type and spatial niche.

### Data layers

| Layer | Recommended role |
|---|---|
| scRNA-seq | Define cell identity, state, trajectory, and aging programs |
| scATAC-seq or multiome | Establish accessible regulatory elements and TF programs |
| Stereo-seq | Add spatial neighbourhood, inflammation/fibrosis niche, and regional context |
| WGS/genotyping | Supply candidate variants and enable eQTL/sQTL/ASE analyses |
| Perturb-seq / CRISPRi/a | Supply causal state-transition training labels |
| Time course | Separate early response, adaptation, and stable trajectory changes |
| Functional assays | Validate RNA predictions against senescence, SASP, proliferation, mitochondrial, DNA damage, and tissue-function endpoints |

### Stepwise analysis plan

1. Create a cell-type-resolved young–middle–old atlas and control donor/batch/sex effects.
2. Derive aging state scores and trajectories; do not assume all cells move along one linear aging axis.
3. Use Stereo-seq to define niche states (for example, inflammatory, fibrotic, vascular, immune-rich).
4. Map genetic variants to likely cis-regulatory mechanisms with AlphaGenome plus cCRE/scATAC/3D evidence.
5. Train or adapt a perturbation model with age-relevant Perturb-seq and time-course data.
6. Rank interventions by predicted movement toward healthy cell state **and** preservation of cell identity/function.
7. Validate the top single and combinatorial interventions experimentally.

### Candidate quality criteria

A high-priority aging marker/target should ideally satisfy several criteria:

- reproducible across donors and cohorts;
- cell-type and spatial-niche specificity are explicit;
- supported by RNA and, where appropriate, ATAC/protein/splicing evidence;
- associated with a functional phenotype, not only chronological age classification;
- changes predictably after perturbation;
- has orthogonal causal support (CRISPR, base editing, MPRA, organoid, or in vivo assay).

## 8. Important cautions

- Young-versus-old observational scRNA-seq can train an age classifier; it cannot, by itself, identify a rejuvenation intervention.
- A transcriptomic “reversal” can be non-beneficial if it compromises cell identity, proliferation control, immune surveillance, or cancer protection.
- Standard 3-prime spatial transcriptomics is generally weak for detailed isoform/splicing analysis. If splicing is a central endpoint, add junction-depth RNA-seq, targeted isoform assays, or long-read RNA-seq.
- Sequence-based prediction does not automatically capture donor environment, chromatin state, protein activity, metabolism, or tissue-scale feedback.
- Evaluate against simple linear/additive baselines and test on held-out donors, perturbations, and cellular contexts.

## 9. Recommended reading order

1. Bunne et al. 2024 — what the AI Virtual Cell should be.
2. GEARS 2024 — how genetic perturbation prediction is operationalised.
3. AlphaGenome paper — sequence-to-cis-regulatory prediction.
4. Nicheformer 2025 — spatial context.
5. Zhu et al. 2026 — aging/senescence perturbation example.
6. Nature Methods 2025 perturbation benchmark — limitations before committing to a model.

## 10. Search terms for follow-up

```text
(aging OR senescence OR longevity) AND (Perturb-seq OR CRISPRi OR single-cell perturbation)
age-dependent eQTL sQTL single-cell
spatial transcriptomics aging cell niche Stereo-seq
cis-regulatory variant enhancer promoter splicing aging
virtual cell foundation model perturbation prediction
AlphaGenome regulatory variant aging GWAS
```
