# MakeHub User Guide

Author and Contact Information {#authors-and-contact-information .unnumbered}
-------------------------------

Katharina J. Hoff, University of Greifswald, Institute for
Mathematics and Computer Science, Bioinformatics Group
(katharina.hoff@uni-greifswald.de)

Contents
========

-   [What is MakeHub?](#what-is-makehub)
-   [Installation](#installation)
    -   [Dependencies](#dependencies)
    -   [MakeHub](#makehub)
-   [Data preprocessing](#preprocessing)
    -   [BRAKER output](#braker_output)
    -   [MAKER output](#maker_output)
    -   [GeMoMa output](#gemoma_output)
-   [Running MakeHub](#running-makehub)
    -   [Creating a new hub](#creating-new-hub)
    -   [Adding tracks to existing hub](#adding_tracks)
    -   [Options explained](#options_explained)
-   [Output of MakeHub](#output-of-makehub)
-   [How to use MakeHub output with UCSC Genome Browser](#use_makehub)
-   [Bug reporting](#bug-reporting)
-   [Citing MakeHub](#citing-braker2-and-software-called-by-braker2)
-   [License](#license)


What is MakeHub? {#what-is-makehub}
================

MakeHub is a command line tool for the fully automatic generation
of track data hubs[^fn1] for visualizing genomes with the UCSC genome
browser[^fn2]. Track data hubs are data structures that contain all required information about a genome for visualizing with the UCSC genome browser.
Track data hubs need to be hosted on a publicly available webspace
(that might be user/password protected) for usage with the UCSC
genome browser.

MakeHub is implemented in Python3 and automatically executes tools
provided by UCSC for generation of assembly hubs
(<http://hgdownload.soe.ucsc.edu/admin/exe>) on Linux and MacOS X
x86_64 computers. For visualization of RNA-Seq alignment
data from BAM files, MakeHub uses Samtools [^fn3]. If installed,
the AUGUSTUS [^fn4] tool bam2wig is used to speed up BAM to
wig format conversion (<https://github.com/Gaius-Augustus/Augustus>),
which is otherwise performed without bam2wig.

MakeHub can either be used to create entirely new track data hubs,
or it can be used to add tracks to hubs that were previously created by MakeHub.

For display by the UCSC Genome Browser, track data hubs need to be
hosted on a publicly accessible web server.

Installation {#installation}
============

In the following, we provide instructions for installing MakeHub
on Linux or MacOS X with x86_64 architecture.

Dependencies {#dependencies}
------------

MakeHub is a Python3 script. It requires Python3, Biopython and
- in the case that BAM files are provided - samtools, gzip, sort,
and optionally the AUGUSTUS tool bam2hints. In the following, we
give instructions on where dependencies can be obtained, and
how they may be installed on Ubuntu Linux.

Python3 is available from <https://www.python.org/downloads/>, or
as package for many Unix systems.

For example, on Ubuntu, install Python3 with:

```
sudo apt install python3
```

We recommend to use pip for installing further python modules.
pip is available at <https://pypi.org/project/pip/>. It is also
avilable as package for many Unix systems.

For example, on ubuntu, install pip with:

```
sudo apt install python3-pip
```

Further, MakeHub uses Biopython (e.g. for parsing a genome file
in order to determine which parts of the genome have been masked
for repeats). Install biopython with pip as follows:

```
pip3 install biopython
```

MakeHub uses the following tools provided by UCSC at
<http://hgdownload.soe.ucsc.edu/admin/exe>:


* bedToBigBed
* genePredCheck
* faToTwoBit
* gtfToGenePred
* hgGcPercent
* ixIxx
* twoBitInfo
* wigToBigWig
* genePredToBed

You may download these binaries and make them available in your
$PATH. However, if you skip installing these tools, they will
be downloaded during MakeHub execution, automatically.

MakeHub uses Samtools for BAM file sorting and conversion. Samtools is avilable at <https://github.com/samtools/>. It is also avilable as
package with many linux distributions.

For example, on ubuntu, install samtools with:

```
sudo apt install samtools
```

MakeHub has been tested with Samtools 1.8-20-g4ff8062. It might
not be fully downward compatible with older versions.

MakeHub uses gzip for compressing wig files that were created
from BAM files. gzip is available at <https://ftp.gnu.org/gnu/gzip/>.
It often installed by default on Unix systems. If not, it is
usually available as a package.

If missing, on Ubuntu, install with:

```
sudo apt install gzip
```

MakeHub uses Unix sort. sort should be installed by default
on all Unix systems.

MakeHub can use the AUGUSTUS tool bam2wig, if that tool is
available in the $PATH. bam2wig is available as part of
AUGUSTUS at <https://github.com/Gaius-Augustus/Augustus>.
Please follow the compilation instructions in
Augustus/auxprogs/bam2wig/README.txt in case the default
make command fails.

MakeHub {#makehub}
-------

MakeHub is a python3 script named make_hub.py. It does not require
a particular installation procedure after download.

It can be executed either with

```
python3 make_hub.py
```

If you add make_hub.py to your $PATH (i.e. by adding the location
of make_hub.py at the bottom of your ~/.bashrc file similar to
```PATH=/path/to/MakeHub:$PATH```, followed by loading the 
~/.bashrc file in case you did not re-open a new bash session with
```source ~/.bashrc```)and make it executable (i.e.
with ```chmod u+x make_hub.py```), it can
be executed with

```
make_hub.py
```

from  any location on your computer.

Data preprocessing {#preprocessing}
==================

MakeHub accepts files in the following formats:

* genome file in FASTA format (simple FASTA headers without
  whitespaces or special characters); if the file is
  softmasked, a track with repeat information will 
  automatically be generated
* BAM file(s) with RNA-Seq to genome alignments
* gene prediction file(s) in GTF-format
* AUGUSTUS hints files in BRAKER-specific GFF hints format


BRAKER output {#braker_output}
-------------

No particular data preprocessing is required for using BRAKER output files. Please use the default GTF-files, not the optional GFF3 files for visualizing BRAKER gene predictions with MakeHub.


MAKER output {#maker_output}
------------

MAKER GFF3 files should be converted to GTF-format prior usage with MakeHub.
We recommend using GenomeTools for format conversion. GenomeTools are
available at <http://genometools.org/pub/>. GenomeTools are also
available as a package for some Unix distributions.

For example, install GenomeTools as follows on Ubuntu:

```
sudo apt install genometools
```

After the gff3_merges step of MAKER, which produces a file such as
maker.gff3, check the predictions for inconsistencies as follows:

```
gt  gff3  -force  -tidy  -o maker_by_gt.gff3  -retainids  -sort  makergff3  2>  errors_by_GenomeTools
cat errors_by_GenomeTools | grep -v gbunit | cut -f6 -d' ' | sort | uniq -c | wc -l
```

This will for example tell you that there are a number of issues 
- or no issues at all.

Convert the maker_by_gt.gff file to GTF-format as follows:

```
gt  gff3_to_gtf  -force  -o maker.gtf  maker_by_gt.gff3
```

Running this command may possibly tell you about features that have
been skipped during conversion. The resulting maker.gtf file can
serve as input for MakeHub.


Gemoma output {#gemoma_output}
-------------

No particular data preprocessing is required for using the Gemoma output file ```filtered_predictions.gff``` in GFF3-format. Please use the MakeHub
option ```--gemoma_filtered_predictions``` for passing this file.

Running MakeHub {#running-makehub}
===============

MakeHub can be used either to create new track data hubs, or to add
tracks to track data hubs that had previously been created.

Creating a new hub {#creating-new-hub}
------------------

The essential arguments for creating a new track data hub are:

* ```--email```/```-e``` EMAIL_ADDRESS - this should be the valid contact information of a hub creator.

* ```--genome```/```-g``` GENOME_IN_FASTA_FORMAT - this should be a 
  single species genome file in FASTA format. If the file contains
  softmasked repeats, a repeat masking track with softmasking 
  information will automatically be generated.
  
* ```--short_label```/```-l``` SHORT_LABEL - this should be a short
  name of the track data hub (without whitespaces). This label will
  for example be used as name for directories produced by
  MakeHub.

We strongly recommend the additional usage of

* ```--long_label```/```-L``` LONG_LABEL - this should be an
  appropriate description of the track data hub, e.g. containing
  the latin and english species names, a project name, or similar.
  Put the LONG_LABEL in quotation marks in case it contains
  whitespaces.

In the following, we show how to obtain example data and
how to generate a small example hub. We will use example data
provided by BRAKER:


```
mkdir example_hub
cd example_hub
wget https://github.com/Gaius-Augustus/BRAKER/raw/master/example/genome.fa
make_hub.py -l example_hub -L \
  "This is an example track data hub with data from BRAKER" \
  -g genome.fa -e katharina.hoff@uni-greifswald.de
```

The resulting track data hub is trivial an contains only a GC content
and repeat masking track.

Note that you can add numerous tracks already at the point of
initial hub creation by using the MakeHub options, e.g. for
gene prediction files. For the sake of illustration purposes,
we will here add a couple of more tracks to the now existing
example_hub in the next section.

Adding tracks to existing hub {#adding_tracks}
-----------------------------

If a hub already exists, you may add tracks to this existing hub
using the option ```--add_track```.


Options explained {#options_explained}
-----------------





Output of MakeHub {#output-of-makehub}
=================

make_hub.py produces a directory that is called identical to the
argument for option ```--short_label```/```-l```. Let's
assume the short label had been ```species```.

```species``` contains the following files:

* ```hub.txt``` - this file contains basic information about the 
  track data hub, for example, the short and long labels, a 
  reference to ```genomes.txt```, and contact information.
  
* ```genomes.txt``` - this file contains references to the 
  configuration files ```trackDb.txt``` and ```groups.txt```, as 
  well as for example a default browsing location.

Furthermore, ```species``` contains another directory ```species```
in which the files ```trackDb.txt``` and ```groups.txt```, as well
as all files that are required for browsing tracks, reside. The
number of files may differ depending on how many tracks have
actually been created.

How to use MakeHub output with UCSC Genome Browser {#use_makehub}
==================================================

Copy the complete hub folder (e.g. ```species```) to a publicly
accessible web server.

Go to <https://genome.ucsc.edu/index.html>, click on ```My Data```
-> ```Track Hubs``` -> ```My Hubs``` and add the link to your
publicly available ```hub.txt``` file  into the URL window.
Subsequently, click on ```Add Hub```.

Bug reporting {#bug-reporting}
=============

Before reporting bugs, please check that you are using the most recent
versions of MakeHub. Also, check the open and closed issues on
github at <https://github.com/Gaius-Augustus/MakeHub/issues>
for possible solutions to your problem.

Reporting bugs on github
------------------------

If you found a bug, please open an issue at
<https://github.com/Gaius-Augustus/MakeHub/issues> 
(or contact katharina.hoff@uni-greifswald.de).

Information worth mentioning in your bug report:

make_hub.py prints information about separate steps on
STDOUT. Please let us know at which step and with what
error message make_hub.py caused problems.


License {#license}
-------

All source code is under GNU public license 3.0 (see
<https://www.gnu.org/licenses/gpl-3.0.de.html>).



References
----------


[^fn1]: Raney BJ, Dreszer TR, Barber GP, Clawson H, Fujita PA, Wang T, Nguyen N, Paten B, Zweig AS, Karolchik D, Kent WJ. 2014.
“Track Data Hubs.” *Bioinformatics* 1;30(7):1003-5.

[^fn2]: Kent WJ, Sugnet CW, Furey TS, Roskin KM, Pringle TH, Zahler AM, Haussler D. 2002.
“UCSC Genome Browser.” *Genome Res.* 12(6):996-1006.

[^fn3]: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R. 2009.
“The sequence alignment/map format and SAMtools.” *Bioinformatics* 26(16):2078-2079.

[^fn4]: Stanke M, Diekhans M, Baertsch R, Haussler D. 2008.
“Using native and syntenically mapped cDNA alignments to improve de novo gene finding.” *Bioinformatics* 24(5):637-644






