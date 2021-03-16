# _Ribo-ODDR_

### Ribo-seq focused (O)ligo (D)esign pipeline to (D)eplete (R)ibosomal RNA fragments

###### v1.0 , Publication Release
---

The software consists of two parts, ''**_Ribo-ODDR_**'' and ''**_Ribo-ODDR_:**_oligo-selector_''

**_Ribo-ODDR_ :** Based on given pilot sequencing data and rRNA sequences, _Ribo-ODDR_ pipeline designs oligos complementary to rRNA sequences and computes their depleting potentials in pilot Ribo-seq data. The tool outputs all oligo designs that pass user-defined filtering criterias, which are explained in later sections. Oligo designs are reported in FASTA, BED, GFF and CSV formats, and if needed, off-target information is outputted in TXT format.

**_Ribo-ODDR_: oligo-selector :** To make browsing through oligo designs easier, we also provide _Ribo-ODDR_:_oligo-selector_, an R-Shiny app.

---
#### Prerequisites for "**_Ribo-ODDR_**" pipeline
**_Ribo-ODDR_** does NOT pre-process Ribo-seq reads coming from pilot experiments. If adapter trimming and size selection marker cleaning is required for generated data, we suggest using the [cutadapt](https://cutadapt.readthedocs.io) tool for such purposes.
* [**Python3**](https://www.python.org/) (>3.0)
  * With [**Biopython(v1.76)**](https://biopython.org/), [**pysam**](https://pysam.readthedocs.io) libraries installed
* [**RIsearch2**](https://rth.dk/resources/risearch/) : RNA-RNA interaction prediction tool used for oligo-target interaction predictions.
* [**TopHat2**](https://ccb.jhu.edu/software/tophat/index.shtml) & [**Bowtie2**](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/) & [**samtools**](http://www.htslib.org/doc/samtools.html) (optional) : The pipeline will use them to identify rRNA fragments, by aligning pre-processed Ribo-seq data to given rRNA sequences. However, one can use any other aligner and feed the alignment files (_BAM_ files) to the pipeline instead. Note that preprocessed reads from pilot experiments must be aligned to the same rRNA sequences that is passed to the _Ribo-ODDR_ pipeline and BAM files must be already indexed as (_XXX.bam.bai_).
* [**Vienna RNA Package - RNAfold**](https://www.tbi.univie.ac.at/RNA/) : (optional) RNA folding prediction tool. The pipeline will run the **RNAfold** program as a sub-process if possible.

#### How to run "**_Ribo-ODDR_**" pipeline?
There are two modes for Ribo-ODDR, "_cross species optimization_" and" _novel oligo design_". Both modes can be run using the same **Ribo-ODDR.py** script.

Please run the following code to get help on the arguments you can pass to this script.

    src/Ribo-ODDR.py -h

To optimize oligos from one organism for Ribo-seq in another organism, _cross species optimization_ mode should be used.  Before running this mode, please make sure the [RIsearch2](https://rth.dk/resources/risearch/) program is pre-installed. Below is an example code on how you can run _Ribo-ODDR_ in this mode (passing "_--RIsearch2_exe d_" argument assumes that RIsearch2 is available by default, 'risearch2.x' command being callable.).

######  _cross species optimization_
    src/Ribo-ODDR.py -r example_data/mouse_rRNAs.fa -o example_data/ -op example_data/human_oligos.fa --RIsearch2_exe d


######  _novel oligo design_
In a simple case, if you already have the bam files ready for sample data-rRNA alignments. You can run this mode as follows. Note that "_mouse_sample_1_" and "_mouse_sample_2_" will be used as sample labels. Please make sure the [RIsearch2](https://rth.dk/resources/risearch/) program is pre-installed.

    src/Ribo-ODDR.py -r example_data/mouse_rRNAs.fa -o example_data/ -ap ./ -as .bam -s mouse_sample_1 mouse_sample_2 --RIsearch2_exe d

If you want to use the trimmed-read files instead, you will have to use a code similar to the one below. The code below assumes that [TopHat2](https://ccb.jhu.edu/software/tophat/index.shtml), [Bowtie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/) and [samtools](http://www.htslib.org/doc/samtools.html) are already installed. You can also provide the path to their executables instead of passing "_d_" for default.

    src/Ribo-ODDR.py -r example_data/mouse_rRNAs.fa -o example_data/ -fp /PATH/PREFIX/TO/YOUR/FASTQ/FILES -as .fastq --bowtie2-build_exe d --samtools_exe d --tophat_exe -s mouse_sample_1 mouse_sample_2 --RIsearch2_exe d

######  _computing more features_
To compute self-folding characteristics, binding energy and off-target predictions for designed oligos you need the [RNAfold](https://www.tbi.univie.ac.at/RNA/) and [RIsearch2](https://rth.dk/resources/risearch/) pre-installed. When they are available by default, you can pass the "_--RIsearch2_exe d --RNAfold_exe d_" arguments or provide the path to their executables instead of "_d_". Note that, for off-target predictions, you also need to provide protein coding transcipt sequences by passing the "--OFFtargets /PATH/TO/FASTA" argument.

---
#### Prerequisites for "_**Ribo-ODDR: oligo-selector**_"
* [**R & RStudio**](https://rstudio.com/)
  * With [**shiny**](https://shiny.rstudio.com/), [**shinythemes**](https://rstudio.github.io/shinythemes/), [**DT**](https://rstudio.github.io/DT/), [**htmltools**](https://cran.r-project.org/web/packages/htmltools/index.html), [**ggplot2**](https://ggplot2.tidyverse.org/) libraries

#### How to run  "_**Ribo-ODDR: oligo-selector**_"
It is provided as an R script, which you can open in RStudio and run easily if you have the packages installed. The Shiny app self-contains the Instructions on how to navigate yourself through the app.

#### Publication link
If you use _Ribo-ODDR_ to optimize your Ribo-seq experiments, do not forget to cite our publication below.
[Ferhat Alkan, Joana Silva, Eric Pintó Barberà, William J Faller, Ribo-ODDR: Oligo design pipeline for experiment-specific rRNA depletion in ribo-seq, Bioinformatics, 2021, btab171.](https://doi.org/10.1093/bioinformatics/btab171)

---
Ribo-ODDR was created by [Faller Lab](https://www.fallerlab.com/)
