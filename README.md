<font size=20>__NPDtools 2.3.0 Manual__</font>

* [About NPDtools](#sec_intro)    
    * [Package content](#sec_intro_content)    
    * [Supported data types](#sec_intro_data)  
* [Installation](#sec_install)  
    * [Downloading NPDtools binaries for Linux](#sec_install_linux)  
    * [Downloading NPDtools binaries for macOS](#sec_install_mac)     
* [Running NPDtools](#sec_running)  
    * [Common command line options](#sec_run_opts)  
    * [Common output files](#sec_run_output)  
    * [Database search pipelines](#sec_run_db_pipelines)  
        * [Dereplicator](#sec_run_db_dereplicator)  
        * [VarQuest](#sec_run_db_varquest)  
        * [Dereplicator+](#sec_run_db_dereplicator+)  
    * [Metabologenomic pipelines](#sec_run_genomic_pipelines)  
        * [MetaMiner](#sec_run_genomic_metaminer)    
* [Citation](#sec_refs)  
* [Feedback and bug reports](#sec_feedback)  


<a name="sec_intro"></a>
# About NPDtools

NPDtools – Natural Product Discovery tools – is a toolkit containing various pipelines for 
_in silico_ analysis of natural product mass spectrometry data. 
This manual will help you to install and run NPDtools. 
The latest version of the manual is available online at <https://github.com/ablab/npdtools>. 
All projects news are at <http://cab.spbu.ru/software/npdtools/>.
 
NPDtools version 2.3.0 was released under the Apache 2.0 License on January 9, 2019 
and can be downloaded from <https://github.com/ablab/npdtools/releases>.
The software is developed in collaboration of [Saint Petersburg State University](http://cab.spbu.ru) (Russia), 
[University of California San Diego](http://cseweb.ucsd.edu/~ppevzner/) (CA, USA) 
and [Carnegie Mellon University](http://mohimanilab.cbd.cmu.edu) (PA, USA).

<a name="sec_intro_content"></a>
## Package content

The current version of NPDtools includes  
* **Dereplicator** — a tool for identification of peptidic natural products (PNPs) through database search of mass spectra
* **VarQuest** — a tool for modification-tolerant identification of novel variants of PNPs
* **Dereplicator+** — a tool for identification of metabolites (both peptidic and non-peptidic) through database search of mass spectra
* **MetaMiner** (former *RiPPquest*, *MetaRiPPquest*) — 
a tool implementing metabologenomics approach for the discovery of 
ribosomally synthesized and post-translationally modified peptides (RiPPs)

<a name="sec_intro_data"></a>
## Supported data

All pipelines in NPDtools work with liquid chromatography–tandem mass spectrometry data (LS-MS/MS). 
Spectra files must be centroided and be in an open spectrum format (**MGF**, **mzXML**, **mzML** or **mzData**). 
NPDtools natively supports [MGF](http://www.matrixscience.com/help/data_file_help.html) (Mascot Generic Format) 
and uses [msconvert](http://proteowizard.sourceforge.net/tools/msconvert.html ) 
utility from the [ProteoWizard package](http://proteowizard.sourceforge.net/index.html) 
to convert spectra in other formats to MGF. 

Database search pipelines (Dereplicator, VarQuest and Dereplicator+) also require 
a chemical structure database of known natural products. 
See [the corresponding section](#sec_run_db_pipelines) for the details on accepted data formats.

The metabologenomic pipelines (currently MetaMiner only) require either raw genome nucleotide sequences
or output of specific genome mining tools. 
See [the corresponding section](#sec_run_genomic_pipelines) for the details on accepted data formats.
 
<a name="sec_install"></a>
# Installation

NPDtools requires a 64-bit Linux system or macOS and Python 2.7 to be pre-installed on it. 
The MetaMiner pipeline also requires [GNU sed](https://www.gnu.org/software/sed/) 
to be present in the `PATH` environment variable as `sed` 
(this is always true for Linux systems but may require additional configurations on macOS since 
GNU sed is usually installed there as `gsed`).  

You can also use NPDtools pipelines online on the [GNPS platform](https://gnps.ucsd.edu/ProteoSAFe/static/gnps-theoretical.jsp). 
In this case, a registration is needed but it is quick and simple. 
The detailed documentation on registration, data upload and pipeline usage is available on the website.

Below are instructions on downloading and running the command line versions of the tools. 
We provide a separate package with NPDtools binaries for each OS type.
In case of successful installation the following files should be present in the `bin` directory:

-   `dereplicator.py` (main executable script for the Dereplicator pipeline)
-   `varquest.py` (main executable script for VarQuest)
-   `dereplicator+.py` (main executable script for Dereplicator+)
-   `metaminer.py` (main executable script for MetaMiner)
-   `dereplicate` (core binary for all database search pipelines)
-   `rippquest_ms` (core binary for the metabologenomic pipeline)
-   several auxiliary scripts and binaries (`npdtools_init.py`, `db_preprocessing.py`, 
`print_score`, `print_structure`, etc)

<a name="sec_install_linux"></a>
## 	Downloading NPDtools binaries for Linux

To download [NPDtools Linux binaries](https://github.com/ablab/npdtools/releases/download/npdtools-2.3.0/NPDtools-2.3.0-Linux.tar.gz) 
and extract them, go to the directory in which you wish NPDtools to be installed and run:

``` bash
    wget https://github.com/ablab/npdtools/releases/download/npdtools-2.3.0/NPDtools-2.3.0-Linux.tar.gz
    tar -xzf NPDtools-2.3.0-Linux.tar.gz
    cd NPDtools-2.3.0-Linux
```

We further refer to this directory as `<npdtools_installation_dir>`. 
All executables are located under `<npdtools_installation_dir>/bin/`, 
so consider adding this subdirectory to the `PATH` variable.
	
	
<a name="sec_install_mac"></a>
## 	Downloading NPDtools binaries for macOS
	
To download [NPDtools macOS binaries](https://github.com/ablab/npdtools/releases/download/npdtools-2.3.0/NPDtools-2.3.0-Darwin.tar.gz) 
and extract them, go to the directory in which you wish NPDtools to be installed and run:

``` bash
    curl https://github.com/ablab/npdtools/releases/download/npdtools-2.3.0/NPDtools-2.3.0-Darwin.tar.gz -o NPDtools-2.3.0-Darwin.tar.gz 
    tar -xzf NPDtools-2.3.0-Darwin.tar.gz
    cd NPDtools-2.3.0-Darwin
```

We further refer to this directory as `<npdtools_installation_dir>`. 
All executables are located under `<npdtools_installation_dir>/bin/`, 
so consider adding this subdirectory to the `PATH` variable.

<a name="sec_running"></a>
# Running NPDtools

To run any NPDtools pipeline, you need at least one LC-MS/MS file and 
either a chemical structure database (for Dereplicator, VarQuest, and Dereplicator+) or 
a genome information/sequence (for MetaMiner).   

Since all pipelines have many common command line options and similar output files, 
we first describe them altogether. Further, we explain specifics of each individual pipeline and 
give examples of running commands for each of them.

To run a NPDtools pipeline from the command line, type
``` bash
    <pipeline>.py [options] <spectra_file> [<spectra_file>] -o <output_dir>
```
where `<pipeline>` is one of `dereplicator`, `varquest`, `dereplicator+`, or `metaminer`; 
`<spectra_file>` could be either a path to a single spectra file or 
to a directory with multiple spectra files inside. In the latter case, 
NPDtools recursively walks through the directory and picks up all files with appropriate extensions 
(`.mgf`, `.mzML`, `.mzXML`, `.mzdata`, etc; case insensitive).
You can specify an unlimited number of input spectra files/directories, 
they will be processed independently (each one by a separate thread). 

Note that here and below we assume that NPDtools executable's directory (`<npdtools_installation_dir>/bin/`) 
is added to the `PATH` variable. Otherwise, you need to start your commands as follows
``` bash
    <npdtools_installation_dir>/bin/<pipeline>.py ...
```
Finally, if you have multiple Python versions installed on your system, 
you may need to explicitly specify the proper version before the running script (we support v.2.7 only):
``` bash
    python2.7 <npdtools_installation_dir>/bin/<pipeline>.py ...
```

To demonstrate working examples of each pipeline, 
we use small test datasets provided within the NPDtools package and also available online from 
<https://github.com/ablab/npdtools/tree/master/test_data>. 
We further assume that `test_data` is in the current working directory and give the corresponding relative paths. 
If you use test dataset from the installation package, 
you may need to specify the full paths as 
```bash
   <npdtools_installation_dir>/share/npdtools/test_data/
```

<a name="sec_run_opts"></a>
## Common command line options

### Basic options
`-o <output_dir> `  
    Specify the output directory. If the directory does not exist, it will be created automatically.
    **Required option**. 
 
`-m <mode>` (or `--mode <mode>`)  
    Running mode. Specify 'LL' for low-low mode (low resolution of both precursor (MS) and product (MS/MS) ions), 
    'HL' for high-low mode, or 'HH' for high-high mode. 
    We consider mass tolerances of ±0.5 Da for low resolution mass spectrometers and 
    ±0.02 Da for high resolution instruments (e.g. q-TOFs, q-Orbitrap).  
    User-defined accuracy thresholds can be specified using advanced options (see below), 
    in this case the running mode should be set to 'custom'. *The default value is 'custom'*. 

`-t <int>` (or `--threads <int>`)  
    Number of threads. *The default value is 50% of all available CPUs* but not less than 1. 
    If NPDtools fails to determine the number of CPUs, the maximum number of threads is set to *4*.

`-h` (or `--help`)                         
    Show all available options and exit.
 
### Advanced options

`--pm_thresh <float>`  
    Mass tolerance for precursor ion (also known as *parent ion*, *MS1*) in Daltons.
    This option is considered only if the running mode is 'custom' (see `-m` option above).
    *The default value is 0.02*.  

`--product_ion_thresh <float>`  
    Mass tolerance for product ion (also known as *fragment ion*, *MS2*, *MS/MS*) in Daltons.
    This option is considered only if the running mode is 'custom' (see `-m` option above).
    *The default value is 0.02*.  
    
`--ppm`  
    Mass tolerances are given in parts per million (ppm) rather than in absolute values in Daltons.
    Thus, this option modifies behaviour of `--pm_thresh` and `--product_ion_thresh` and also
    changes defaults for 'LL', 'HL' and 'HH' running modes (see `-m` option above). 
    *The default values are 500 ppm and 20 ppm for low and high resolution data,
    respectively*. Note that 20 ppm translates into 0.01 Da tolerance for the precursor/fragment 
    ion at m/z 500 Da/e, 0.02 Da at m/z 1000 Da/e, or 0.03 Da at m/z 1500 Da/e.
    *Currently not available in MetaMiner*.
    
`-e <int>` (or `--max-charge <int>`)     
    Max possible charge to consider for spectra without explicitly specified charge 
    (e.g. if an MGF file does not have `CHARGE=` field). *The default value is 2*. 

`--fdr`  
    Estimate False Discovery Rate (approximately doubles computation time). 

<a name="sec_"></a>
## Common output files
 
NPDtools pipelines store all their output files in `<output_dir>`, which is set by the user (`-o` option).
Auxiliary and intermediate files (mostly logs and configs) are located under `<output_dir>/work/` subdirectory.
The main result files are located directly under `<output_dir>`:

-   `summary.tsv` (summary numbers about the input data and identified natural products, elapsed time, and estimated FDR if available)
-   `all_matches.tsv` (the list of all identified compound–spectrum matches)
-   `significant_matches.tsv` (the list of the most reliable compound–spectrum matches; 
    each hit in the list passes a strict statistical significance threshold [P-value < 1e-10 by default])
-   `significant_unique_matches.tsv` (the list of the most reliable identified compounds (unique); 
    only the most statistically significant hit per each compound is listed)

If `--fdr` option is specified, NPDtools pipelines estimate FDR using Target–Decoy Approach. 
In this case, spectra are examined against both the target natural product database and 
an artificially generated decoy database of the same size. All hits in the decoy database are considered as false positives 
(see more details in [Elias & Gygi, 2007](https://www.ncbi.nlm.nih.gov/pubmed/17327847)). 
In this case, the output directory also contains:  
 
-   `all_decoy_matches.tsv` (the list of all identified compound–spectrum matches in the decoy database)
-   `significant_decoy_matches.tsv` (the list of the most reliable compound–spectrum matches in the decoy database)

### Format of output reports

All identifications are reported in plain text tab-separated value files (`.tsv`). 
Each file starts with a header line containing column descriptions. 
The rest lines represent compound–spectrum matches, so they include information about both 
the corresponding spectrum and the compound.

The column names present in all pipelines reports:
-  `SpecFile` (filepath of the spectra file)
-  `Scan` (scan number of the identified spectrum inside the spectra file)
-  `SpectrumMass` (mass of the spectrum in Daltons)
-  `Retention` (retention time of the spectrum in seconds)
-  `Charge` (charge of the spectrum)
-  `Score` (score of the compound–spectrum match)
-  `P-Value` (statistical significance of the compound–spectrum match; *currently NOT reported by Dereplicator+*)
-  `FDR` (estimated FDR at the corresponding P-Value level; 
*only if `--fdr` specified; only in `significant_` reports*)
-  `PeptideMass` (mass of the compound in Daltons)

Reports in database search pipelines (Dereplicator, VarQuest and Dereplicator+) also contain:
-  `Name` (name of the identified compound according to the database description file)
-  `LocalPeptideIdx` (index number of the identified compound in the database description file)
-  `SMILES` (the identified compound in SMILES format; 
*only if SMILES are available for the database; only in `significant_` reports*)
-  `LocalSpecIdx` (index number of the identified spectrum inside the spectra file; *usually differs from Scan*)
-  `Adduct` (adduct ion, e.g. 'M+H' or 'M+2H')
-  `VisualizationID` (auxiliary column for consistency with the online GNPS release; 
*only in Dereplicator and VarQuest pipelines; only in `significant_` reports*)

MetaMiner reports also contain:
-  `SeqFile` (filepath of the genome sequence file)
-  `Class` (class of the identified RiPP compound)
-  `FragmentSeq` (raw initial sequence of the identified compound)
-  `ModifiedSeq` (the sequence of the identified compound with all applied modifications 
[this sequence is actually matched with the spectrum])

<a name="sec_run_db_pipelines"></a>
## Database search pipelines

All database search pipelines (Dereplicator, VarQuest and Dereplicator+) **require** 
a chemical structure database of PNPs (Dereplicator and VarQuest) or general natural products (Dereplicator+).
Currently we support chemical structures only in [MDL MOL V3000](http://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf). 
Other formats (e.g. SMILES, sample peptide string, etc) could be converted into MOL V3000 format using 
[molconvert](https://docs.chemaxon.com/display/docs/Molecule+file+conversion+with+Molconverter) utility from the 
[ChemAxon Marvin package](https://chemaxon.com/products/marvin) (not included into NPDtools!). Our pipelines also require
custom database description file (`library.info`) listing all compounds filepaths, names, masses, numbers of amino acids, 
and metadata in a simple space-separated text format. You may find database example (both MOLs and description file) in 
the test data provided inside NPDtools package. You could also download a
[real database](http://cab.spbu.ru/files/VarQuest/pnpdatabase.tar.gz) with 5,021 PNPs used in the 
[VarQuest study](https://www.nature.com/articles/s41564-017-0094-2).


The following options are used to specify the database:  
`--db-path <dirpath>`     
  Relative or absolute path to the natural products database root directory. 
  Paths of individual compounds (in the description file) should be relative to this directory.
  For example, if a compound is located in `/path/to/db/mol_dir/compound.mol` and 
  the description file (`library.info`) contains `mol_dir/compound.mol` entry, you need to specify `--db-path /path/to/db/`.
  **Required option**.
  
`-l <filepath>` (or `--library-info <filepath>`)   
  Path to the database description file. *The default location is `library.info` 
  inside the database dir (provided with `--db-path` option)*. 
  **The description file is required: if the database dir does not contain `library.info`, this option is mandatory**. 
  
`-s <filepath>` (or `--smiles <filepath>`)   
  Path to the list of SMILES corresponding to the database compounds. 
  The number of entries and their order should be exactly the same to the order of entries 
  in the database description file. *The default location is `library.smiles` 
  inside the database dir (provided with `--db-path` option)*. The SMILES are used only for filling `significant_` reports,
  if not present, `SMILES` column in the reports will be missing (everything else will be working properly).
    
<a name="sec_run_db_dereplicator"></a>    
### Dereplicator

Dereplicator searches for exact identifications, that is it can match a known PNP of mass *m* 
with a mass spectrum of mass *M* only if *m* and *M* are within a small mass error threshold 
(set by `--pm_thresh` option, 0.02 Da by default).
The only exception is matching of chemical structures with isotopic shifts 
(see `--isotope` option below). In this case, a PNP and a spectrum can be matched 
if *m*, *m + c*, or *m + 2c* are within a mass threshold from *M* (*c* is the mass difference between 
Carbon-13 isotope and regular Carbon-12 which is *13.003355 - 12 = 1.003355 Da*).

The specific Dereplicator pipeline options are  
`-i <int>` (or `--isotope <int>`)     
    Maximum accepted isotopic shift (0, 1, or 2). *The default value is 0*.
    
`--fdr-limit <float>`  
    Maximum allowed FDR in percents for significant matches 
    (in 0.0-100.0 range). The hits above this threshold will go to `all_matches.tsv` report only.
     *The default value is 1.0 (that is up to 1%)*.
    
`--p-value-limit <float>`  
    Minimum allowed P-value for significant matches (in 0.0-1.0 range).
    The hits above this threshold will go to `all_matches.tsv` report only.
    *The default value is 1e-10*.
    
#### Usage example   
A sample run of Dereplicator may look like this:
```bash
    dereplicator.py test_data/dereplicator/ --db-path test_data/sample_database/ -o dereplicator_outdir
```
In this case, all spectra files in `test_data/dereplicator/` will be searched against 
a sample natural products database located in `test_data/sample_database/` 
(the description file is `test_data/sample_database/libary.info`) 
and the identification results will be saved in `dereplicator_outdir`. 
The search is performed with all default parameters, 
see the [corresponding subsection](#sec_run_opts) for the default values and available options.
See important notes on specifying paths of the running script and `test_data` in the 
[beginning of this section](#sec_running). 

If the run is finished correctly, you will see many *Putisolvins* identifications listed in 
`dereplicator_outdir/significant_matches.tsv`.

<a name="sec_run_db_varquest"></a>
### VarQuest

VarQuest enables search for mutated/modified PNP variants, that is it can match a known PNP of mass *m* 
with a mass spectrum of mass *M* if *m* and *M* are within a large *MaxMod* threshold 
(maximum allowed mutation/modification mass set by `--max-mod` option described below). By design,
VarQuest looks for a single mutation/modification. We assume that its actual mass is *mod = M - m* and
try to apply *mod* to each residue of the known PNP to find the most likely position of the modification.
The most likely novel PNP variant structure is scored against the experimental spectrum and the statistical 
significance of such score is reported.

The only specific VarQuest pipeline option is  
`--max-mod <float>`     
    Maximum allowed mutation/modification mass in Daltons. *The default value is 300.0*.   
    
VarQuest also accepts `--fdr-limit` and `--p-value-limit` options from Dereplicator pipeline 
(see few paragraphs above).

#### Usage example 
A sample run of VarQuest may look like this:
```bash
    varquest.py test_data/varquest/ --db-path test_data/sample_database/ -o varquest_outdir
```
In this case, all spectra files in `test_data/varquest/` will be searched in 
a modification-tolerant manner against 
a sample natural products database located in `test_data/sample_database/` 
(the description file is `test_data/sample_database/libary.info`) 
and the identification results will be saved in `varquest_outdir`. 
The search is performed with all default parameters, 
see the [corresponding subsection](#sec_run_opts) for the default values and available options.
See important notes on specifying paths of the running script and `test_data` in the 
[beginning of this section](#sec_running). 

If the run is finished correctly, you will see identifications of novel variants of *Surugamide B*,
*Venepeptide* and *Massetolide A* listed in `varquest_outdir/significant_matches.tsv`.
Note the mass differences between *PeptideMass* and *SpectrumMass* columns in the report, 
they correspond to the weights of modifications in the novel variants 
(comparing to the known PNPs from the database).
You may find more info about these three particular identifications in 
[Gurevich et al, 2018](https://www.nature.com/articles/s41564-017-0094-2).

<a name="sec_run_db_dereplicator+"></a>
### Dereplicator+

Dereplicator+ enables search for polyketides, lipids, terpenes, benzenoids, alkaloids, and other classes
of natural products (including PNPs, of course). Similarly to Dereplicator, this pipeline can identify only exact known compounds 
from the database. In contrast to PNP-focused pipelines (Dereplicator and VarQuest) which assume that 
mass spectrometers break only amide bonds (Nitrogen-to-Carbon), Dereplicator+ uses much more complicated
fragmentation model (general natural products may not include amide bonds at all). Dereplicator+ considers breakage 
of Carbon-to-Carbon (CC), Oxygen-to-Carbon (OC) and Nitrogen-to-Carbon (NC) bonds. This pipeline also considers
*bridge* and *2-cut* fragmentations of a metabolite. 
In the former case, breakage of a single bond (a bridge or a 1-cut) is enough to disconnect the whole structure into two fragments.
In the latter case, breakage of a pair of bonds (a 2-cut) is needed to disconnect the whole structure into two fragments.

To specify a particular fragmentation pattern, `--fragmentation_mode <mode>` option is used. 
The currently available modes are 'general_3_1_3', 'general_6_1_6', 'general_6_3_6', 'general_9_1_9', and 'amide_3_1_3'.
In each mode, the prefix word corresponds to use of CC, OC, and NC (*general*) or only NC (*amide*) bond breakages in 
the theoretical spectrum simulation. The numbers correspond to the maximum allowed number of *bridge*, *2-cut*, and *total* breakages, respectively.
For example, 'general_6_1_6' mode corresponds to up to six *bridge* and up to one *2-cut* breakages of any type (CC, OC, or NC) but 
no more than six breakages in total. Since the number of CC bonds in any chemical structure is huge, 
we additionally disallow use of CC bond breakage in any *2-cut* breakage and allow at max one CC breakage of type *bridge*.
*The default mode is 'general_6_1_6'*. 
 
Another important feature of Dereplicator+ is the ability to preprocess the metabolite database. 
This pipeline is much more time-consuming than Dereplicator/VarQuest, 
so a preprocessing of the database may save a lot of time when running many spectra against the same metabolite database.

The specific Dereplicator+ options related to the database preprocessing are  
`--preprocess`  
    Perform database preprocessing before dereplication.
    The preprocessed database files may be reused
    in consecutive runs against the same database and
    using the same fragmentation mode. The preprocessed files are saved in `<output_dir>/db_preproc/`. 
    There will be two files if you use `--fdr` (target and decoy databases) and one file otherwise (target database only).

`--preprocessed_ft <filepath>`  
    Path to the preprocessed database. You can preprocess your database using `--preprocess` (see above). 
    *The default location is `<library_info>.<fragmentation_mode>.ft.bin`, 
    e.g. `library.info.general_6_1_6.ft.bin` inside the database directory with `library.info` description file*. 
    
    
`--preprocessed_ftd <filepath>`  
    Path to the preprocessed decoy database (used for FDR computation only, so `--fdr` should be specified). 
    You can preprocess your database using `--preprocess --fdr` (see above). 
    *The default location is `<library_info>.<fragmentation_mode>.ftd.bin`, 
    e.g. `library.info.general_6_1_6.ftd.bin` inside the database directory with `library.info` description file*. 

Finally, Dereplicator+ introduces one more specific option:  
`--min-score <int>`  
    Minimum score for significant matches (a positive integer number).
    The hits below this threshold will go to `all_matches.tsv` report only.
    Note that the current version of Dereplicator+ does not compute P-values.
    *The default value is 12*.

Moreover, Dereplicator+ uses slightly different default values for some of the previously described options. 
The alternative (more strict) defaults are `--pm_thresh 0.005` and `--product_ion_thresh 0.01`.

#### Usage example 
A sample run of Dereplicator+ may look like this:
```bash
    dereplicator+.py test_data/dereplicator+/ --db-path test_data/sample_database/ -o dereplicator+_outdir
```
In this case, all spectra files in `test_data/dereplicator+/` will be searched against 
a sample natural products database located in `test_data/sample_database/` 
(the description file is `test_data/sample_database/libary.info`) 
and the identification results will be saved in `dereplicator+_outdir`. 
The search is performed with all default parameters, 
see the [corresponding subsection](#sec_run_opts) for the default values and available options.
See important notes on specifying paths of the running script and `test_data` in the 
[beginning of this section](#sec_running). 

If the run is finished correctly, you will see identifications of a PNP (Surugamide) and a polyketide (Chalcomycin) 
listed in `dereplicator+_outdir/significant_matches.tsv`.

<a name="sec_run_genomic_pipelines"></a>
## Metabologenomic pipelines

Metabologenomic pipelines combine metabolomic (mass spectra) and genomic data to identify novel metabolites and gene clusters encoding them.
Thus, the genomic data is **required** to run this type of workflows. The data could be either raw nucleotide sequences 
(a high-quality reference or a draft assembly in [FASTA format](https://www.ncbi.nlm.nih.gov/BLAST/fasta.shtml)) 
or specific genome mining tools' output (protein sequences of the translated gene clusters, see more info below).
The current version of NPDtools includes only one metabologenomic pipeline – *MetaMiner* – intended for identification of RiPPs.
We are also working on specialized tools for identification of NRPs and other classes of natural products (will be available in future releases).

  
The following options are used to specify the genomic sequences:  
`-s <path>` (or `--sequence <path>`)  
    Path to a sequnce file  or to a directory with multiple sequence files inside.
    In the latter case, NPDtools recursively walks through the directory and picks up all files 
    with appropriate extensions (by default: `.fna`, `.fasta`, or `.fa`; case insensitive). 
    You can specify an unlimited number of input sequence files/directories, 
    they will be processed independently (see also `--correspondence` option below).
    *By default, we assume that sequences are raw nucleotide input but this can be modified by specific options (see below)*.
    **At least one sequence file is required**.

`-C <filepath>` (or `--correspondence <filepath>`)  
    Path to a file describing correspondence between sequence and spectra files. 
    The file should be tab-separated and has two columns listing basenames of spectra and sequence filepaths.
    If not provided, the all-vs-all analysis will be performed.

<a name="sec_run_genomic_metaminer"></a>
### MetaMiner

MetaMiner can identify various classes of RiPPs by combining genome/metagenome mining with analysis of tandem mass spectra.
The tool either process raw genome nucleotide sequences with [HMMER](http://hmmer.org/) v.3.1 
(nucleotides are translated into proteins using all 6 possible frames) or works with output
of third-party genome mining tools ([BOA](https://github.com/nafizh/Boa)'s protein `.fasta` 
or [antiSMASH](https://antismash.secondarymetabolites.org)'s `.final.gbk`).

MetaMiner specific options are:  
`-a` (or `--antismash`)  
    Sequence files are antiSMASH output (`.final.gbk`). If not specified, the input files are expected to 
    be raw genome nucleotide sequences in FASTA format (see also `--boa` option). Tested with antiSMASH v.2 output.
       
`--boa`                   
    Sequence files are BOA output (protein `.fasta`). If not specified, the input files are expected to 
    be raw genome nucleotide sequences in FASTA format (see also `--antismash` option).
    
`-c <class>` (or `--class <class>`)  
    Class of RiPPs to look for. Valid choices are: 'formylated',
    'glycocin', 'lantibiotic', 'lap', 'lassopeptide', 'linaridin',
    'proteusin', 'cyanobactin', and 'methanobactin'. You can also specify 'all' to try all classes one by one.
    *The default value is 'lantibiotic'*.  

#### Usage example 
A sample run of MetaMiner may look like this:
```bash
    metaminer.py test_data/metaminer/ -s test_data/metaminer/ -o metaminer_outdir
```
In this case, all spectra files in `test_data/metaminer/` will be searched against 
all sequence files in the same directory. In this particular case, it is a search of `test_data/metaminer/AmfS.mgf` spectrum
against `test_data/metaminer/S.griseus_fragment.fasta` genome fragment. The search mode (considered RiPP class) is 'lantibiotic' (by default).
The identification results will be saved in `metaminer_outdir`. 
The search is performed with all default parameters, 
see the [corresponding subsection](#sec_run_opts) for the default values and available options.
See important notes on specifying paths of the running script and `test_data` in the 
[beginning of this section](#sec_running). 

If the run is finished correctly, you will see identification of a lantibiotic with "TGSQVSLLVCEYSSLSVVLCTP" original sequence 
and "T-18GS-18QVS-18LLVCEYS-18SLSVVLCTP" sequence after modifications in `dereplicator_outdir/significant_matches.tsv`.
The modifications "T-18" and "S-18" correspond to dehydrobutyrine and dehydroalanine, respectively.
These sequences correspond to AmfS peptide, you may read more about it in [Ueda et al, 2002](https://www.ncbi.nlm.nih.gov/pubmed/11844785).

<a name="sec_refs"></a>
# Citation
If you use NPDtools in your research, please cite the papers describing corresponding pipelines.

For Dereplicator please cite [Mohimani et al, Nature Chemical Biology, 2017](https://www.nature.com/articles/nchembio.2219).

For VarQuest please cite [Gurevich et al, Nature Microbiology, 2018](https://www.nature.com/articles/s41564-017-0094-2).

For Dereplicator+ please cite [Mohimani et al, Nature Communications, 2018](https://www.nature.com/articles/s41467-018-06082-8).

For MetaMiner please cite Cao et al, 
*MetaMiner: A Peptidogenomics Approach for the Discovery of Ribosomally Synthesized and Post-translationally Modified Peptides,* 
Submitted, *2019* (initial version of the paper is [available on bioRxiv](https://www.biorxiv.org/content/early/2017/12/03/227504)).

<a name="sec_feedback"></a>
# Feedback and bug reports

Your comments, bug reports, and suggestions are very welcomed. 
They will help us to further improve NPDtools.
You can leave them at [our GitHub repository tracker](https://github.com/ablab/npdtools/issues) 
or sent them via support e-mail: <npdtools.support@cab.spbu.ru>.  
