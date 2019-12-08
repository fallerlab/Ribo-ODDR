# _Ribo-ODDR_

### _Ribo_-seq focused _O_ligo _D_esign pipeline to _D_eplete _R_ibosomal RNA fragments

###### v0.9 , peer review release
---

The software consists of two parts, ''**_Ribo-ODDR_**'' and ''**_Ribo-ODDR_:**_oligo-selector_''

**_Ribo-ODDR_ :** Based on given pilot sequencing data and rRNA sequences, _Ribo-ODDR_ pipeline designs oligos complementary to rRNA sequences and computes their depleting potentials in pilot Ribo-seq data. The tool outputs all oligo designs that pass user-defined filtering criterias, which are explained in later sections. Oligo designs are reported in FASTA, BED, GFF and CSV formats, and if needed, off-target information is outputted in TXT format.

**_Ribo-ODDR_: oligo-selector :** To make browsing through oligo designs easier, we also provide _Ribo-ODDR_:_oligo-selector_, an R-Shiny app.

---
#### Prerequisites for "**_Ribo-ODDR_**" pipeline
**_Ribo-ODDR_** does NOT pre-process Ribo-seq reads coming from pilot experiments. If adapter trimming and size selection marker cleaning is required for generated data, we suggest using [cutadapt](https://cutadapt.readthedocs.io) tool for such purposes.
* [**Python3**](https://www.python.org/) (>3.0)
  * With [**Biopython**](https://biopython.org/), [**pysam**](https://pysam.readthedocs.io) libraries installed
* [**TopHat2**](https://ccb.jhu.edu/software/tophat/index.shtml) & [**Bowtie2**](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/) & [**samtools**](http://www.htslib.org/doc/samtools.html) (optional) : The pipeline will use them to identify rRNA fragments, by aligning pre-processed Ribo-seq data to given rRNA sequences. However, one can use any other aligner and feed the alignment files (_BAM_ files) to the pipeline instead. Note that preprocessed reads from pilot experiments must be aligned to the same rRNA sequences that is passed to the _Ribo-ODDR_ pipeline and BAM files must be already indexed as (_XXX.bam.bai_).
* [**Vienna RNA Package**](https://www.tbi.univie.ac.at/RNA/) : (optional) RNA folding prediction tool. The pipeline will run the **RNAfold** program as a sub-process if possible.
* [**RIsearch2**](https://rth.dk/resources/risearch/) (mandatory for off-target predictions and "_cross-species optimization mode_") : RNA-RNA interaction prediction tool used for oligo-target interaction predictions.

#### How to run "**_Ribo-ODDR_**" pipeline?
asdasfasg
    asdasd
---
#### Prerequisites for "_**Ribo-ODDR: oligo-selector**_"
* [**R & RStudio**](https://rstudio.com/)
  * With [**shiny**](https://shiny.rstudio.com/), [**shinythemes**](https://rstudio.github.io/shinythemes/), [**DT**](https://rstudio.github.io/DT/), [**htmltools**](https://cran.r-project.org/web/packages/htmltools/index.html), [**ggplot2**](https://ggplot2.tidyverse.org/) libraries

#### How to run  "_**Ribo-ODDR: oligo-selector**_"
asfsafsaf
    asdasd
---
Ribo-ODDR was created by [Faller Lab]()
