# Bioinformatics notes

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com)

Bioinformatics learning and data analysis tips and tricks. Please, [contribute and get in touch](CONTRIBUTING.md)! 

# Table of content

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Awesome](#awesome)
  - [Notes by Ming Tang](#notes-by-ming-tang)
- [Pipelines](#pipelines)
  - [liftOver](#liftover)
  - [k-mers](#k-mers)
- [Courses](#courses)
- [Videos](#videos)
- [Bioinformatics core organization](#bioinformatics-core-organization)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# Genomics notes

- [RNA-seq notes](https://github.com/mdozmorov/RNA-seq)
- [scRNA-seq notes](https://github.com/mdozmorov/scRNA-seq_notes)
- [ChIP-seq notes](https://github.com/mdozmorov/ChIP-seq_notes)
- [Methylation notes](https://github.com/mdozmorov/Methylation_notes)
- [scATAC-seq notes](https://github.com/mdozmorov/scATAC-seq_notes)
- [SNP notes](https://github.com/mdozmorov/SNP_notes)
- [Hi-C tools](https://github.com/mdozmorov/HiC_tools)
- [Hi-C data](https://github.com/mdozmorov/HiC_data)
- [scHi-C notes](https://github.com/mdozmorov/scHiC_notes)
- [Cancer notes](https://github.com/mdozmorov/Cancer_notes)
- [Immuno notes](https://github.com/mdozmorov/Immuno_notes)

See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes.

## Awesome

- Historical perspective on genome sequencing technologies. From the landmark 1953 Watson and Crick publication through sequencing of nucleic acids, shotgun sequencing (Messing, Sanger), Human Genome Project (HGP) & Celera Genomics, milestones in genome assembly, Next and Third generation sequencing, single molecule sequencing (SMRT, ONT, PacBio), long-read genome assemblers (FALCON, Canu). Box1 - companies developing sequencing technologies, Box 2 - more details on HGP, Box 3 - bioinformatics tools for genome assembly. <details>
    <summary>Paper</summary>
    Giani, Alice Maria, Guido Roberto Gallo, Luca Gianfranceschi, and Giulio Formenti. “Long Walk to Genomics: History and Current Approaches to Genome Sequencing and Assembly.” Computational and Structural Biotechnology Journal, November 2019, S2001037019303277. https://doi.org/10.1016/j.csbj.2019.11.002.
</details>

- The complete assembly of human genome (haploid CHM13 cell line). 3.055 billion base pairs, no gaps for all 22 chromosomes plus ChrX, new genes. Resolving ribosomal rDNA sequences. PacBio, Oxford Nanopore, other technologies. [The Telomere-to-Telomere (T2T) consortium](https://sites.google.com/ucsc.edu/t2tworkinggroup), [UCSC](http://genome.ucsc.edu/cgi-bin/hgTracks?genome=t2t-chm13-v1.0&hubUrl=http://t2t.gi.ucsc.edu/chm13/hub/hub.txt), [NCBI PRJNA559484](https://www.ncbi.nlm.nih.gov/bioproject/559484), [GitHub with download links to FASTA, gff3, liftover chains](https://github.com/marbl/CHM13). [Assembly issues tracker](https://github.com/marbl/CHM13-issues) <details>
    <summary>Paper</summary>
    Nurk, Sergey, Sergey Koren, Arang Rhie, Mikko Rautiainen, Andrey V. Bzikadze, Alla Mikheenko, Mitchell R. Vollger et al. "The complete sequence of a human genome." bioRxiv (2021). https://doi.org/10.1101/2021.05.26.445798
</details>

- [T2T-Y](https://github.com/marbl/CHM13) - Complete sequencing of human Y chromosome, as a part of the T2T project. Lots of information about corrected errors, complete structures of sex-determining genes, additional protein-coding genes, alternating satellite patterns, centromere structure, repeats. Illumina 33X, PacBio HiFi 42X, ONT 125X sequencing. Table 1 - stats comparison with GRCh38-Y. Improves alignment, variant calling, helps removing human contamination from microbial studies. [Code used in the paper](https://github.com/arangrhie/T2T-HG002Y). <details>
    <summary>Paper</summary>
    Rhie, Arang, Sergey Nurk, Monika Cechova, Savannah J. Hoyt, Dylan J. Taylor, Nicolas Altemose, Paul W. Hook, et al. “The Complete Sequence of a Human Y Chromosome.” Preprint. Genomics, December 1, 2022. https://doi.org/10.1101/2022.12.01.518724.
</details>

- [Milestones in Genomic Sequencing](https://www.nature.com/immersive/d42859-020-00099-0/index.html) by Nature, 2000-2021 period, interactive infographics

- [Genome Coordinate Cheat Sheet](http://alternateallele.blogspot.com/2012/03/genome-coordinate-cheat-sheet.html)

- [Awesome-Bioinformatics](https://github.com/danielecook/Awesome-Bioinformatics) - A curated list of awesome Bioinformatics libraries and software, by Daniel Cook and community-contributed

- [awosome-bioinformatics](https://github.com/openbiox/awosome-bioinformatics) - A curated list of resources for learning bioinformatics

- [awesome-bioinformatics-education](https://github.com/lskatz/awesome-bioinformatics-education) - A curated list of resources specific to learning bioinformatics. Courses, hands-on training, tools.

- [awesome-cancer-variant-databases](https://github.com/seandavi/awesome-cancer-variant-databases) - A community-maintained repository of cancer clinical knowledge bases and databases focused on cancer variants, by Sean Davis and community-maintained

- [awesome-single-cell](https://github.com/seandavi/awesome-single-cell) - List of software packages for single-cell data analysis, including RNA-seq, ATAC-seq, etc. By Sean Davis and community-maintained

- [awesome-alternative-splicing](https://github.com/HussainAther/awesome-alternative-splicing) - Alternative splicing resources

- [awesome-10x-genomics](https://github.com/johandahlberg/awesome-10x-genomics) - List of tools and resources related to the 10x genomics GEMCode/Chromium system, by Johan Dahlberg and community-contributed

- [awesome-multi-omics](https://github.com/mikelove/awesome-multi-omics) - List of software packages for multi-omics analysis, by Mike Love and community-maintained

- [awesome-bioinformatics-benchmarks](https://github.com/j-andrews7/awesome-bioinformatics-benchmarks) - A curated list of bioinformatics bench-marking papers and resources

- [awesome_genome_browsers](https://github.com/davemcg/awesome_genome_browsers) - genome browsers and genomic visualization tools. By David McGaughey

- [awesome-expression-browser](https://github.com/federicomarini/awesome-expression-browser) - A curated list of software and resources for exploring and visualizing (browsing) expression data, and more. By Federico Marini and community-maintained

- [awesome-genome-visualization](https://github.com/cmdcolin/awesome-genome-visualization) - A list of interesting genome browser or genome-browser-like implementations. By Colin Diesh

- [awesome-microbes](https://github.com/stevetsa/awesome-microbes) - List of software packages (and the people developing these methods) for microbiome (16S), metagenomics (WGS, Shot-gun sequencing), and pathogen identification/detection/characterization. By Steve Tsang and community-contributed

- [aws-for-bioinformatics](https://github.com/lynnlangit/aws-for-bioinformatics) - Amazon Web Services for Bioinformatics Researchers

- [algorithmsInBioinformatics](https://github.com/joachimwolff/algorithmsInBioinformatics) - Bioinformatics algorithms: Needleman-Wunsch, Feng-Doolittle, Gotoh and Nussinov implemented in Python. By Joachim Wolff, with lecture notes. Also, [rklib](https://github.com/agalitsyna/rklib) - Rabin-Karp implementation of sequence substring search for DNA/RNA, and [lv89](https://github.com/lh3/lv89) - C implementation of the Landau-Vishkin algorithm to compute the edit distance between two strings.

- [BAM-tricks](https://github.com/IARCbioinfo/BAM-tricks) - Tip and tricks for BAM files

- [Bedtools Cheatsheet](https://gist.github.com/ilevantis/6d6ecf8718a5803acff736c2dffc933e#subtract)

- [bioicons](https://bioicons.com/) - A library of free open source icons for science illustrations in biology and chemistry. SVG downloads. [GitHub](https://github.com/duerrsimon/bioicons)

- [bioconvert](https://github.com/bioconvert/bioconvert) - converter between various data formats. Python, command line.

- [biotools](https://github.com/jdidion/biotools) - A massive collection of references on the topics of bioinformatics, sequencing technologies, programming, machine learning, and more. By John Didion

- [bio.tools](https://bio.tools/) - a community-driven comprehensive and consistent registry of bioinformatics resources. Supported by ELIXIR - the European infrastructure for biological information. <details>
    <summary>Paper</summary>
    Ison, Jon, Kristoffer Rapacki, Hervé Ménager, Matúš Kalaš, Emil Rydza, Piotr Chmura, Christian Anthon, et al. “Tools and Data Services Registry: A Community Effort to Document Bioinformatics Resources.” Nucleic Acids Research 44, no. D1 (January 4, 2016): D38–47. https://doi.org/10.1093/nar/gkv1116.
</details>

- [BioNumPy](https://github.com/bionumpy/bionumpy) - Numpy-like analysis of biological data (FASTQ, BED, BAM, etc.). Indexing, vectorized functions, high efficiency on BED intersection, kmer counting, and VCF operations. [Documentation](https://bionumpy.github.io/bionumpy/) with examples of filtering FASTQ reads, working with BAM files, computing GC content, general handling of sequences, intervals, PWMs. <details>
    <summary>Paper</summary>
    Rand, Knut Dagestad, Ivar Grytten, Milena Pavlovic, Chakravarthi Kanduri, and Geir Kjetil Sandve. “BioNumPy: Fast and Easy Analysis of Biological Data with Python.” Preprint. Bioinformatics, December 22, 2022. https://doi.org/10.1101/2022.12.21.521373.
</details>

- [Circa](https://omgenomics.com/circa/) - online tool for Circos plot creation

- [Jbrowse 2](https://jbrowse.org/jb2/) - Java-based genome browser that includes linear and advanced visualization for synteny, dotplots, breakpoints, gene fusions, whole genome overview, Circos, Spreadsheet view, SV inspector. Many file formats, including Hi-C (hic), BAM/CRAM (various sort pileup options), plain text (VFC, BED, etc.), some specific to specific views. Bookmarks (can be imported from BED), sessions (shareable). Web, desktop versions, can be run from [Jupyter](https://gmod.github.io/jbrowse-jupyter/docs/html/index.html) notebooks, R ([JBrovseR](https://gmod.github.io/JBrowseR/)), [command line](https://jbrowse.org/jb2/docs/cli/), extended with [plugins](https://jbrowse.org/jb2/plugin_store/). [Documentation](https://jbrowse.org/jb2/docs/), [Demos/tutorials](https://jbrowse.org/jb2/demos/). <details>
    <summary>Paper</summary>
    Diesh, Colin, Garrett J Stevens, Peter Xie, Teresa De Jesus Martinez, Elliot A. Hershberg, Angel Leung, Emma Guo, et al. “JBrowse 2: A Modular Genome Browser with Views of Synteny and Structural Variation.” Preprint. Bioinformatics, July 31, 2022. https://doi.org/10.1101/2022.07.28.501447.
</details>

- [genome_assembly_tools](https://github.com/nadegeguiglielmoni/genome_assembly_tools) - List of genome assembly tools. Links to software, publications

- [For all your seq... DNA & RNA](https://www.illumina.com/content/dam/illumina-marketing/documents/applications/ngs-library-prep/for-all-you-seek-dna.pdf) - Illumina flyer with infographics of all sequencing-by-synthesis technologies. [RNA](https://www.illumina.com/content/dam/illumina-marketing/documents/applications/ngs-library-prep/for-all-you-seek-rna.pdf) and [DNA](https://www.illumina.com/content/dam/illumina-marketing/documents/applications/ngs-library-prep/for-all-you-seek-dna.pdf) versions

- [learngenomics.dev](https://github.com/stjude/learngenomics.dev) - A guided, intuitive introduction to genomics for software engineers. Curated by the community. [Website](https://learngenomics.dev/)

- [Learning Bioinformatics At Home](https://github.com/harvardinformatics/learning-bioinformatics-at-home)

- [List of software/websites/databases/other stuff for genome engineering](https://github.com/davidliwei/awesome-CRISPR)

- [multimodal-scRNA-seq](https://github.com/arnavm/multimodal-scRNA-seq) - Figure depicting the breadth of multimodal scRNA-seq technologies. References to technology-specific papers

- [ref-gen](https://github.com/lh3/ref-gen) - Human reference genome analysis sets, by Heng Li. Recommended human genome sequences for various genome asseblies, justification.

- [SequencEnG](http://education.knoweng.org/sequenceng/) - Hierarchical summary of 66 sequencing technologies, computational algorithms, references to papers. 
    - Zhang, Y., Manjunath, M., Kim, Y., Heintz, J., and Song, J.S. (2019). SequencEnG: an interactive knowledge base of sequencing techniques. Bioinformatics 35, 1438–1440.

- [SRPlot](https://bioinformatics.com.cn/srplot) - free online visualization platform of many bioinformatics-specific plots. Tang D, Chen M, Huang X, Zhang G, Zeng L, Zhang G, Wu S, Wang Y. [SRplot: A free online platform for data visualization and graphing.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0294236) PLoS One. 2023 Nov 9;18(11):e0294236. doi: 10.1371/journal.pone.0294236. PMID: 37943830.

- [Wasm](https://webassembly.org/) - WebAssembly, a framework for deploying tools and code on the web.
  - [webR/WebAssembly Demo](https://jperkel.github.io) - running R in web browser
  - [biowasm](https://biowasm.com/) - running C/C++ genomics tools in the browser
  - [sandbox.bio](https://sandbox.bio/tutorials/?id=bedtools-intro) - an interactive bedtools tutorial developed by the Quinlan Lab
  - Perkel, Jeffrey M. "[No installation required: how WebAssembly is changing scientific computing](https://doi.org/10.1038/d41586-024-00725-1)." Nature 627, no. 8003 (2024): 455-456.

- [The Leek group guide to genomics papers](https://github.com/jtleek/genomicspapers)- Jeff Leek recommended list of genomics papers

### Notes by Ming Tang

- [Unix, R and python tools for genomics and data science](https://github.com/crazyhottommy/getting-started-with-genomics-tools-and-resources) - links and references to many computational biology resources, by Ming Tang 

- [bioinformatics-one-liners](https://github.com/crazyhottommy/bioinformatics-one-liners) - collection of bioinformatics-genomics bash one liners, using awk, sed etc., by Ming Tang

- [cloud_computing_resources](https://github.com/crazyhottommy/cloud_computing_resources) - genomics-oriented cloud resources.

- [RNA-seq-analysis](https://github.com/crazyhottommy/RNA-seq-analysis) - RNAseq analysis notes, by Ming Tang

- [ChIP-seq-analysis](https://github.com/crazyhottommy/ChIP-seq-analysis) - ChIP-seq analysis notes, by Ming Tang

- [DNA-seq-analysis](https://github.com/crazyhottommy/DNA-seq-analysis) - Notes on whole exome and whole genome sequencing analysis, by Ming Tang

- [DNA-methylation-analysis](https://github.com/crazyhottommy/DNA-methylation-analysis) - DNA methylation analysis notes from Ming Tang

- [scATAC-seq-analysis-notes](https://github.com/crazyhottommy/scATAC-seq-analysis-notes) - single-cell ATAC-seq notes, by Ming Tang


## Pipelines

- [List of bioinformatics tools developed by IHEC Int'l Human Epigenome Consortium researchers](http://ihec-epigenomes.org/research/tools/) - tools for all types of genomic analyses

- [Learning Nextflow in 2020](https://www.nextflow.io/blog/2020/learning-nextflow-in-2020.html) - Materials to learn NextFlow
- [Nextflow Tower](https://tower.nf/) - management of Nextflow data pipelines. [Video overview](https://youtu.be/gtwS3A7qLj0)

- [nf-core](https://www.nextflow.io/) - a framework for Nextflow-based pipeline creation, community-driven. Integrated with Conda, Docker, Biocontainers. Scalable to a cloud level. Pipeline assemblers: [Flowcraft](https://github.com/assemblerflow/flowcraft), [Pipeliner](https://github.com/montilab/pipeliner). [nf-core GitHub](https://github.com/nf-core/). All pipelines are on the [nf-core hub](https://nf-co.re/). 
    - Ewels, Philip, Alexander Peltzer, Sven Fillinger, Johannes Alneberg, Harshil Patel, Andreas Wilm, Maxime Garcia, Paolo Di Tommaso, and Sven Nahnsen. “[Nf-Core: Community Curated Bioinformatics Pipelines](https://doi.org/10.1101/610741).” Preprint. Bioinformatics, April 16, 2019. 

- [bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen) - Validated, scalable, community developed variant calling and RNA-seq analysis. [Documentation](https://bcbio-nextgen.readthedocs.org)

- [Appyters Catalog](https://appyters.maayanlab.cloud) of over 40 Appyters. Appyters extend Jupyter Notebooks into standalone web-based bioinformatics applications. Adds Jupyter "magics" to construct output Python code with jinja, the UI is built with svelte and Bootstrap and served with nginx, while data-persistent capabilities such as page hits are deferred to a PostgreSQL database that can be accessed via a PostgREST microservice API. Many bioinformatics workflows, including machine learning, omics data analysis. Entry form for data upload and settings, the Appyter executes the corresponding notebook in the cloud or locally. Comparison with Google Colab, Binder, Papermill and similar tools. Ability to convert Appyters to Docker images. Examples of [bulk RNA-seq Appyter](https://appyters.maayanlab.cloud/#/Bulk_RNA_seq) (preprocessing, differential expression using limma, edgeR, DESeq2, etc., enrichment using GO, KEGG, Reactome etc., correlation with LINCS L1000 databases), [scRNA-seq Appyter]( https://appyters.maayanlab.cloud/#/scRNA_seq) (filtering by mitochondria, counts, imputation, clustering, differential analysis, pseudotime, cell type prediction), [the Harmonizome-ML Appyter](https://appyters.maayanlab.cloud/#/harmonizome_ml) (organized gene-centric data from multiple databases, used to benchmark ML algorithms to predict gene-gene functional associations), [the Drugmonizome-ML Appyter](https://appyters.maayanlab.cloud/#/Drugmonizome_ML) (ML interface to the [Drugmonizome](https://maayanlab.cloud/drugmonizome/) database for predicting drug and small-molecule attributes), [the patient-cohort RNA-seq data viewer](https://appyters.maayanlab.cloud/#/Patient_Cohorts_RNASeq_Viewer) for exploring patient's characteristics, e.g., from TCGA, [data processing Appyters](https://appyters.maayanlab.cloud/#/?tags=Harmonizome&q=ETL) to extract, transform, and load data for Harmonizome, [the Enrichr visualization Appyter](https://appyters.maayanlab.cloud/#/?tags=Enrichr&q=skylar), [the set comparison Appyter](https://appyters.maayanlab.cloud/#/CompareSets). [Table 1](https://www.cell.com/action/showFullTableHTML?isHtml=true&tableId=tbl1&pii=S2666-3899%2821%2900023-4) - overview of 27 Appyters. <details>
    <summary>Paper</summary>
    Clarke, Daniel J.B., Minji Jeon, Daniel J. Stein, Nicole Moiseyev, Eryk Kropiwnicki, Charles Dai, Zhuorui Xie, et al. “Appyters: Turning Jupyter Notebooks into Data-Driven Web Apps.” Patterns, March 2021, 100213. https://doi.org/10.1016/j.patter.2021.100213.
</details>

### liftOver

- [Benchmarking of six liftover tools](https://github.com/phuoc362/Lifted) (UCSC liftOver, rtracklayer::liftOver, CrossMap, NCBI Remap, flo and segment liftover) on converting CpG sites and WGBS samples. Considering integrity-preserving and non integrity-preserving. Explanation of chain file format (Figure 1). segment liftover generates more reliable results than USCS liftOver. <details>
    <summary>Paper</summary>
    Luu, Phuc-Loi, Phuc-Thinh Ong, Thanh-Phuoc Dinh, and Susan J Clark. “Benchmark Study Comparing Liftover Tools for Genome Conversion of Epigenome Sequencing Data.” NAR Genomics and Bioinformatics 2, no. 3 (September 1, 2020): lqaa054. https://doi.org/10.1093/nargab/lqaa054.
</details>

- [AirLift](https://github.com/CMU-SAFARI/AirLift) - command-line tool (C) for remapping reads between genome assemblies. Very fast, highly accurate in identifying ground truth SNP/InDels in lifted data. Describe limitations of current tools. Methods, lookup tables, eight steps, each read remapped using one of four independent cases. Supports multiple file formats, from BAM to BED. Compared with liftOver, CrossMap. Supplementary - description of other tools (liftOver, CrossMap, MCBI Genome Remapping Service, Segment_liftover, PyLiftover, nf-LO, LevioSAM, Liftoff). <details>
    <summary>Paper</summary>
    Kim, Jeremie S., Can Firtina, Meryem Banu Cavlak, Damla Senol Cali, Mohammed Alser, Nastaran Hajinazar, Can Alkan, and Onur Mutlu. “AirLift: A Fast and Comprehensive Technique for Remapping Alignments between Reference Genomes.” arXiv, November 21, 2022. http://arxiv.org/abs/1912.08735.
</details>

- [CrossMap](https://github.com/liguowang/CrossMap) - genome coordinates conversion between different assemblies (such as hg18 (NCBI36) <=> hg19 (GRCh37)). It supports commonly used file formats including BAM, CRAM, SAM, Wiggle, BigWig, BED, GFF, GTF, MAF VCF, and gVCF. [Documentation](https://crossmap.readthedocs.io/en/latest/)

- [HiCLift](https://github.com/XiaoTaoWang/HiCLift) - coordinate liftover for region pairs. Converts chromatin contacts in [pairs](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) format, [allValidPairs](https://nservant.github.io/HiC-Pro/RESULTS.html) outputted by HiC-Pro, cool format and hic (highest resolution) format. Uses UCSC chain files, represented as IntervalTrees, and the T2T chain files. Benchmarked using HiCRep on the original and converted matrices. <details>
    <summary>Paper</summary>
    Wang, Xiaotao, and Feng Yue. “HiCLift: A Fast and Efficient Tool for Converting Chromatin Interaction Data between Genome Assemblies.” Edited by Tobias Marschall. Bioinformatics 39, no. 6 (June 1, 2023): btad389. https://doi.org/10.1093/bioinformatics/btad389.
</details>

- [Liftoff](https://github.com/agshumate/Liftoff) - An accurate GFF3/GTF lift over pipeline

- [liftover](https://github.com/jeremymcrae/liftover) - liftover for python, made fast with cython

- [liftOverBedpe](https://github.com/dphansti/liftOverBedpe) - A liftOver wrapper to accomodate BEDPE files, by Doug Phanstiel, requires Python 2.7

- [pairLiftOver](https://pypi.org/project/pairLiftOver/) - Python package that converts the two-dimensional genomic coordinates of chromatin contact pairs between assemblies, by Xiaotao Wang. Input - 4D-DCIC [pairs](https://github.com/4dn-dcic/pairix) format or allValidPairs defined by HiC-Pro. Based on [pyliftover](https://github.com/konstantint/pyliftover)

### k-mers

- [KAT](https://github.com/TGAC/KAT) - The K-mer Analysis Toolkit, compares k-mer spectra. Analyzes jellyfish hashes or FASTQ/FASTA files.

- [kcmbt_mt](https://github.com/abdullah009/kcmbt_mt) - A very fast k-mer counter

- [kmdiff](https://github.com/tlemane/kmdiff) - Differential k-mer analysis

## Courses

- [computational_biology](https://github.com/fayjustin/computational_biology) - computational biology lectures (pdf) and labs (ipynb). Alignment algorithms (BLAST, BW transform, dynamic programming), molecular evolution (Markov models, likelihood), phylogenetics (trees, pruning, MCMC), annotations (HMM), motifs (PWMs), machine learning (clustering, classification), bioinformatics (pipelines), biophysics (protein structure, prediction).

- [training-collection](https://github.com/sib-swiss/training-collection) - Bioinformatics training materials: Scripting and languages (UNIX shell, Python, R, git), Sequence data analysis (General, Miscellaneous, RNA-seq, ChIP-seq, Single cell, Variant analysis, Genome assembly, Metagenomics), Computational methods and pipelines (Containers, Nextflow, Snakemake, Galaxy, CWL, High performance computing), Statistics and machine learning, Reproducibility and data management. By [SIB Swiss Institute of Bioinformatics](https://github.com/sib-swiss)

- [2020-GGG298](https://github.com/ngs-docs/2020-GGG298) - Course materials for GGG298 - Tools to support data-intensive research, by Titus Brown. Unix, Conda, Snakemake, project organization, Git, Slurm, R/Rmarkdown

- [BI-BE-CS-183-2023](https://github.com/pachterlab/BI-BE-CS-183-2023) - Introduction to Computational Biology and Bioinformatics Course at Caltech, 2023.

- [biocomp-book](https://github.com/camilogarciabotero/biocomp-book) - a course on Fundamentals of Computational Biology. Sequencing analysis and genome assembly/analysis. [Website](https://camilogarciabotero.github.io/biocomp-book/)

- [bioinformatics](https://github.com/ossu/bioinformatics) - Path to a free self-taught education in Bioinformatics! Systematically organized links to courses.

- [biomedicalresearch2021](https://github.com/schatzlab/biomedicalresearch2021) - Course Materials for EN.601.452 / AS.020.415 Computational Biomedical Research & Advanced Biomedical Research, by Michael Schatz. Links to other courses, papers.

- [bggn213_W23](https://github.com/bioboot/bggn213_W23) - A hands-on introduction to the computer-based analysis of genomic and biomolecular data from the Division of Biological Sciences, UCSD. By prof. Barry Grant. Lecture slides, lab, homework material, videos. [Website](https://bioboot.github.io/bggn213_W23/)

- [RNA-seq theory and analysis](https://www.sudosight.com/RNA-seek/RNA-seq/Theory/) - theory and practical guidelines for RNA-seq data analysis, by Skyler Kuhn

- [Introduction to Bioinformatics and Computational Biology](https://liulab-dfci.github.io/bioinfo-combio/) - the course material for STAT115/215 BIO/BST282 at Harvard University. [YouTube](https://www.youtube.com/playlist?list=PLeB-Dlq-v6taAXK6ZCGfqImrNWJzFt3p3)

- [Introduction to Bioimage Analysis](https://bioimagebook.github.io/) - This book tries explain the main ideas of image analysis in a practical and engaging way. By Pete Bankhead. [GitHub](https://github.com/bioimagebook/bioimagebook.github.io)

- [GRanges tutorial](https://seandavi.github.io/ITR/RangesAndSignal.html) - Working with Genomic Ranges, by Sean Davis

- [Data science for economists](https://github.com/uo-ec607/lectures) - from R, tidyverse to GitHub, web scraping, Docker, Google Cloud and more. By [Grant McDermott](http://grantmcdermott.com/)

- [Reproducible research and data analysis with Linux containers and Nextflow pipelines](https://biocorecrg.github.io/CoursesCRG_Containers_Nextflow_May_2021/index.html), [GitHub](https://github.com/biocorecrg/CoursesCRG_Containers_Nextflow_May_2021)

- [Bioinformatics Coffee Hour](https://github.com/harvardinformatics/bioinformatics-coffee-hour) - Short lessons from FAS Informatics coffee hour,data science, command line, R basics, Snakemake

- [The Biostar Handbook: 2nd Edition](https://www.biostarhandbook.com/index.html) - Biostar Handbook - bioinformatics survival guide. A practical overview for the data analysis methods of bioinformatics. From Unix/command line to each type of sequencing data ana analysis

- [Beginner's Handbook to Next Generation Sequencing by GenoHub](https://genohub.com/next-generation-sequencing-handbook/) - omics technologies, experimental descriptions

- [Computational Genomics with R](http://compgenomr.github.io/book/) book by Altuna Akalin. From R basics to different types of bioinformatics analyses. [GitHub](https://github.com/compgenomr/book)

- The [Bioinformatics algorithms](http://bioinformaticsalgorithms.com/index.htm) web site. [Videos](http://bioinformaticsalgorithms.com/videos.htm) and [Texts](http://bioinformaticsalgorithms.com/faqs.htm) covering topics from biology to genome assembly algorithms, motif finding, sequence comparisons and more bioinformatics tasks and solutions. By Phillip Compeau and Pavel Pevzner

- [Algorithms for DNA Sequencing](http://www.langmead-lab.org/teaching-materials/) - Ben Langmead's course material. [Slides on GitHub](https://github.com/BenLangmead/ads1-slides.git), Python/Colab code examples, Youtube videos 

- [Applied Computational Genomics Course at UU](https://github.com/quinlan-lab/applied-computational-genomics) by Aaron Quinlan. [All slides on Google drive](https://drive.google.com/drive/folders/0B5Jmsvw39gJkT0dZOExrR3EwVjg), Youtube videos for each lecture

- [Introduction to Computational Biology](https://biodatascience.github.io/compbio/) course by Mike Love. [GitHub](https://github.com/biodatascience/compbio_src)

- [Data Analysis in Genome Biology](http://girke.bioinformatics.ucr.edu/GEN242/) course by Thomas Girke. Bioinformatics of NGS data analysis [GitHub](https://github.com/tgirke/GEN242)

- [Bioconductor for Genomic Data Science](https://kasperdanielhansen.github.io/genbioconductor/) course by Kasper Hansen. Includes videos, code examples and lecture material. [GitHub](https://github.com/kasperdanielhansen/genbioconductor)

- [BPA-CSIRO Workshops](https://bpa-csiro-workshops.github.io/) - Cancer Genomics, Introduction to Next Generation Sequencing Hands-on Workshop. Links to topic-oriented GitHub repositories, [PDF handouts](http://bpa-csiro-workshops.github.io/btp-workshop-ngs/). [GitHub](https://github.com/BPA-CSIRO-Workshops/btp-workshop-ngs)

- [The Bioconductor 2018 Workshop Compilation](https://bioconductor.github.io/BiocWorkshops/), editors - Levi Waldron, Sean Davis, Marcel Ramos, Lori Shepherd, Martin Morgan. Task-oriented workshops covering a range of genomics/Bioconductor analyses. [GitHub](https://github.com/Bioconductor/BiocWorkshops)

- [CSAMA](https://github.com/Bioconductor/CSAMA) - Course material for CSAMA: Statistical Data Analysis for Genome Scale Biology, by Bioconductor team. Lectures, Labs

- [JHU EN.601.749: Computational Genomics: Applied Comparative Genomics](https://github.com/schatzlab/appliedgenomics) course, by Michael Schatz. Links to similar courses there. [Other courses by Michael](http://schatzlab.cshl.edu/teaching/)

- [JHU Data Science lab](https://jhudatascience.org/courses.html) - several data science courses, links to topic-specific GitHub repositories, other resources

- [CSE 549 - Introduction to Computational Biology](http://www3.cs.stonybrook.edu/~skiena/549/) by Steven Skiena. Includes video lectures. Another course: [CSE 519 - Data Science](https://www3.cs.stonybrook.edu/~skiena/519/)

- [DIYtranscriptomics.github.io](http://diytranscriptomics.com/) - Course material for the "do-it-yourself" RNA-seq course, by Daniel Beiting. [GitHub](https://github.com/DIYtranscriptomics/DIYtranscriptomics.github.io) 

- [Informatics for RNA-seq: A web resource for analysis on the cloud](https://github.com/griffithlab/rnaseq_tutorial) by Griffith lab. All aspects of RNA-seq analysis

- [RNA-seqlopedia](http://rnaseq.uoregon.edu/) - one long page of all steps of RNA-seq data analysis, from molecular biology to computational analysis, very detailed

- [BaRC Hot Topics](http://barc.wi.mit.edu/education/hot_topics/) - lecture slides and handouts on all genomics topics, from Unix to microarray, sequencing, genomics and statistics

- [Courses by Dr. Raghu Machiraju](http://web.cse.ohio-state.edu/~machiraju.1/teaching/). Topics: data visualization, biomedical informatics, computer graphics, linked from the homepage. [CSE5599-BMI7830](http://web.cse.ohio-state.edu/~machiraju.1/teaching/CSE5599-BMI7830/) - biomedical informatics. [CSE5544](http://web.cse.ohio-state.edu/~machiraju.1/teaching/CSE5544/) – Introduction to Data Visualization

- [alignment-and-variant-calling-tutorial](https://github.com/ekg/alignment-and-variant-calling-tutorial) - basic walk-throughs for alignment and variant calling from NGS sequencing data, PDF lecture, by Erik Garrison.

- [Oxford Nanopore sequencing tutorial](https://github.com/AbeelLab/integrated_bioinformatics) using procaryotic genomes. [Supplementary material](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007314#sec006) - walkthrough, tools, data, VM image available. <details>
    <summary>Paper</summary>
    Salazar, Alex N., Franklin L. Nobrega, Christine Anyansi, Cristian Aparicio-Maldonado, Ana Rita Costa, Anna C. Haagsma, Anwar Hiralal, et al. “[An Educational Guide for Nanopore Sequencing in the Classroom](https://doi.org/10.1371/journal.pcbi.1007314).” PLOS Computational Biology, (January 23, 2020)
</details>

## Videos

- [Machine Learning in Genomics - Fall 2019](https://www.youtube.com/playlist?list=PLypiXJdtIca6U5uQOCHjP9Op3gpa177fK) by Manolis Kellis. [Course website](http://stellar.mit.edu/S/course/6/fa19/6.047/)

- [Free online training in bioinformatics and biostatistics!](https://pickingupthetabb.wordpress.com/building-a-bioinformaticist/free-online-training-in-bioinformatics-and-biostatistics/) by David Tabb. Various topics beyond genomics

- [Bioconductor Workshop 2: RNA Seq and ChIP Seq Analysis](https://www.youtube.com/watch?v=J9LNYhhNhrk) - 6 hours workshop on RNA-seq and ChIP-seq technology and analysis by Levi Waldron and others

- [Differential Splicing Analysis with RNA-Seq: Current Applications, Approaches, & Limitations](https://www.youtube.com/watch?v=LknDQw08P5w) - 1 hour overview of differential splicing analysis

- [NHGRI_Genomics2016 - "Current Topics in Genome Analysis 2016" course](https://www.youtube.com/playlist?list=PL1ay9ko4A8skYqjhrA4INDZ7IHtebS0lY) - A lecture series covering contemporary areas in genomics and bioinformatics, [slides](https://www.genome.gov/12514288/current-topics-in-genome-analysis-2016-course-syllabus-handouts-and-videos/)

- [MIT_SysBiol2014 - "Foundations of Computational and Systems Biology"](https://www.youtube.com/watch?v=lJzybEXmIj0&list=PLUl4u3cNGP63uK-oWiLgO7LLJV6ZCWXac), [slides](https://ocw.mit.edu/courses/biology/7-91j-foundations-of-computational-and-systems-biology-spring-2014/lecture-slides/)

- [Regulatory Genomics and Epigenomics](https://www.youtube.com/playlist?list=PLgKuh-lKre10fYXnD-8ghi9b9xl83CKet) - series of genomcs-oriented talks by Simons Institute. 26 videos, \~30min each. 

- [Foundations of Data Science — Spring 2016](https://data-8.appspot.com/sp16/course) - course from UC Berkeley. Instructor: John DeNero. Video and slides. This course is accompanied by [other "connector" courses from UC Berkeley](https://data-8.appspot.com/sp16/modules/extra_tabs/render?index=3)

- [A Roadmap to the Living Genome](https://videocast.nih.gov/summary.asp?Live=12792&bhcp=1) by [John Stamatoyannopoulos](http://www.stamlab.org/). An overview of cell type-specific (epi)genomic landscapes, visualization and analysis techniques, association and enrichment of disease-associated variants in regulatory regions.

- [RNA-Seq Methods and Algorithms](https://www.youtube.com/watch?v=96yBPM8lEt8&list=PLfFNmoa-yUIb5cYG2R1zf5rtrQQKZvKwG) - short video course by Harold Pimentel, pseudoalignment, kallisto, sleuth, practical

- [Integrating ENCODE Data With Your Research: An Interactive Survey of ENCODE Tools and Resources](https://www.pathlms.com/ashg/courses/25101) - set of short videos about ENCODE data and functionality

## Bioinformatics core organization

- [Bioinfo-core.org](http://www.bioinfo-core.org) - the community of Bioinformatics Core FAcilities, ISMB workshops, slides, resources about Cores.

- [List of bioinformatics core facilities](http://www.bioinfo-core.org/index.php/BioWiki:Community_portal) at bioinfo-core.org

- Dragon, Julie A., Chris Gates, Shannan Ho Sui, John N. Hutchinson, R. Krishna Murthy Karuturi, Alper Kucukural, Shawn Polson, et al. “[Bioinformatics Core Survey Highlights the Challenges Facing Data Analysis Facilities](https://doi.org/10.7171/jbt.20-3102-005).” Journal of Biomolecular Techniques, June 2020 - Bioinformatics core facility considerations. Responsibilities, finansial models, software/data, reporting/accreditation, challenges and concerns, future.

- Kallioniemi, O., L. Wessels, and A. Valencia. “[On the Organization of Bioinformatics Core Services in Biology-Based Research Institutes](https://doi.org/10.1093/bioinformatics/btr125).” Bioinformatics, (May 15, 2011) - 1-page Bioinformatics core recommendations.

- Lewitter, Fran, Michael Rebhan, Brent Richter, and David Sexton. “[The Need for Centralization of Computational Biology Resources](https://doi.org/10.1371/journal.pcbi.1000372).” PLoS Computational Biology, (June 26, 2009) - Bioinformatics core as the center of computational resources, advantages and disadvantages. Questions to consider: IT and computational infrastructure, keeping up with science, training and education, funding model, hiring, evaluation, esternal affitiates with the core, outreach.

- Richter, Brent G., and David P. Sexton. “[Managing and Analyzing Next-Generation Sequence Data](https://doi.org/10.1371/journal.pcbi.1000369).” PLoS Computational Biology, (June 26, 2009) - Sequencing data computational analysis and storage, skills for data analysis (Unix, scripting, parallel computing, network, databases, biology/genomics, connect science with (novel) software solutions).

- Lewitter, Fran, and Michael Rebhan. “[Establishing a Successful Bioinformatics Core Facility Team](https://doi.org/10.1371/journal.pcbi.1000368).” PLoS Computational Biology, (June 26, 2009) - Considerations for successful bioinformatics core development (objectives, personnel, prioritization/time management, staying connected with research trends, outreach). [Slides from the ISMB 2008 BoF on best practices in running bioinformatics cores](http://www.bioinfo-core.org/index.php/ISMB_2008:_BoF_on_best_practices_in_running_bioinformatics_cores).


