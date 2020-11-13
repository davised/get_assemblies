get_assemblies
==============

.. contents:: **Table of Contents**
    :backlinks: none

Who might want this software?
-----------------------------

This software is intended to facilitate the use of the `NCBI assembly database.
<https://www.ncbi.nlm.nih.gov/assembly>`_

The intended audience is scientific researchers, computational biologists, and
bioinformaticians who need sequence data for comparative genomics projects.

Anyone using NCBI sequence data (for reference genomes/transcriptomes) can
benefit as well.

PhD/Masters students and undergraduates are especially encouraged to submit
issues if they are having trouble using this software.

This software is written in python but has a perl dependency. I'm happy to help
those who are having trouble with software installs.

Windows users, please install WSL to make use of this software. Using a Linux
distribution will make your life as a computational researcher significantly
easier.

Installation
------------

get_assemblies is distributed on `PyPI
<https://pypi.org/project/get-assemblies>`_ as a universal wheel and is
available on Linux/macOS and supports Python 3.7+. This software will work on
Windows using `WSL
<https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_.

get_assemblies depends on `NCBI Entrez Direct
<https://www.ncbi.nlm.nih.gov/books/NBK179288/>`_ which requires Perl. Perl is
installed by default on most \*nix systems. If edirect is not currently
installed, please run ``get_assemblies --dledirect`` to install.

.. code-block:: bash

    $ python3 -m pip install -U get-assemblies

Dependencies
------------

Python modules:

1. `rich <https://github.com/willmcgugan/rich>`_

See ``requirements.txt`` for more info.

External programs:

1. `NCBI Entrez Direct <https://www.ncbi.nlm.nih.gov/books/NBK179288/>`_

You can install external programs using the ``get_assemblies --dledirect``
command. These will be installed to ``${HOME}/edirect`` unless otherwise
specified.

Just tell me how to run it
--------------------------

.. code-block:: bash

    $ get_assemblies organism 'Pseudomonas fluorescens'

This will find all genomes tagged as 'Pseudomonas fluorescens' in NCBI's
database. By default, this command will only check to see how many genomes
fall into this category.

To download metadata for the genomes:

.. code-block:: bash

    $ get_assemblies organism 'Pseudomonas fluorescens' --function metadata

Check out the ``metadata.tab`` file that is created after running this command.
Generally you will want to select a subset from your search. One way to do this
is to select the lines that include the genomes of interest, and then saving
the assembly accessions to a file. You can either delete the lines that you
don't want, or use ``grep`` to pull out the lines that you want to keep. Then
you can use ``cut -f 14 > accs.txt`` to get the assembly accesions in a file.

.. code-block:: bash

    $ cat accs.txt | get_assemblies assembly_ids - --function genomes -o fna

This will download the nucleotide fasta files for your genomes of interest.

Overview
--------

This tool was written to make accessing genomic data from NCBI easier. The
output files are renamed such that each assembly has a Genus species strain in
the filename to make it easy to find the genomes that you're interested in. You
won't have to spend time renaming the files by hand.

This software is effectively a wrapper for the NCBI edirect tools that makes
getting genome files easier. If you are interested in starting a comparative
genomics project, this is the tool for you.

The software supports four types of input:

1. organism input, either taxonomy rank names (e.g. Genus species, Family) or
   taxids
2. assembly ids, either accessions or uids
3. nuccore ids (e.g. individual contig/chromosome names)
4. json input (e.g. the intermediate files - docsums - produced by this script)

Five file type outputs are supported:

1. Nucleotide genome sequence (fna)
2. Nucleotide coding sequence (ffn)
3. Amino acid coding sequence (faa)
4. General feature format (i.e. tab-delimited features) (gff)
5. GenBank format (gbk)

The program will attempt to find a unique prefix per genome assembly. This
prefix will be in the resulting filename. A metadata file that contains much
of the relevant information per genome will also be included. This file can
be included as a supplementary table for a manuscript in a comparative genomics
project.

If you need to make phylogenetic trees with these data, check out my other
python package, `automlsa2 <https://pypi.org/project/automlsa2/>`_.

More Examples
-------------

.. code-block:: bash

    $ get_assemblies organism 'Mycobacterium'
    2020-10-15 22:49:53,257 - INFO - Found 7522 genomes to download.
    2020-10-15 22:49:53,257 - INFO - Expect 37610MB to 52654MB of data.

.. code-block:: bash

    $ get_assemblies organism --type ID 167539 --function genomes -o gbk
    2020-10-15 23:10:13,822 - INFO - Found 1 genomes to download.
    2020-10-15 23:10:13,822 - INFO - Expect 5MB to 7MB of data pending the chosen file types for download.
    chunk: 1it [00:01,  1.21s/it]
    docsums: 100%|██████████████████████████████| 1/1 [00:00<00:00, 5146.39it/s]
    2020-10-15 23:10:16,262 - INFO - Downloading 1 files.
    100% [##################################################]           1M / 1M]
    2020-10-15 23:10:18,044 - INFO - P_marinus_CCMP1375_SS120.gbk successfully downloaded.
    download: 100%|███████████████████████████████| 1/1 [00:01<00:00,  1.78s/it]
    $ ls
    docsums0.json       metadata.tab
    get_assemblies.log  P_marinus_CCMP1375_SS120.gbk

.. code-block:: bash

    $ echo GCA_000269645.2 | get_assemblies assembly_ids -
    2020-10-15 23:18:04,107 - INFO - Found 1 genomes to download.
    2020-10-15 23:18:04,107 - INFO - Expect 5MB to 7MB of data pending the chosen file types for download.

Usage
-----

.. code-block:: bash

    $ get_assemblies -h
    usage: get_assemblies [-h] [--debug] [--version] [--dledirect [DLEDIRECT]] {organism,assembly_ids,nuccore_ids,json_input} ...

    Downloads assemblies & annotations from NCBI.

    positional arguments:
      {organism,assembly_ids,nuccore_ids,json_input}
                            Choose from this list of input types.
        organism            Valid NCBI organism or taxids.
        assembly_ids        Valid NCBI assembly IDs.
        nuccore_ids         Valid NCBI nucleotide accessions.
        json_input          Valid NCBI JSON docsums.

optional arguments:

-h, --help            show this help message and exit
--debug               Turn on debugging messages.
--version             show program's version number and exit
--dledirect <[DLEDIRECT]>
                      Download edirect to given location. [~/edirect]

Organism
^^^^^^^^

.. code-block:: bash

    $ get_assemblies organism -h
    usage: get_assemblies organism [-h] [--type {text,ID}] [--function {check,metadata,genomes} [{check,metadata,genomes} ...]]
                                   [--annotation] [--metadata_append] [--typestrain] [--keepmulti] [--force]
                                   [-f {abbr,full,strain}] [-o {fna,ffn,gff,gbk,faa,all} [{fna,ffn,gff,gbk,faa,all} ...]]
                                   [--edirect EDIRECT] [--debug]
                                   query

    positional arguments:
      query                 Valid NCBI organism text term or ID

optional arguments:

-h, --help            show this help message and exit
--type <{text,ID}>    Input is text term (default) or ID
--function <{check,metadata,genomes} [{check,metadata,genomes} ...]>
                      check counts, download metadata, or genomes. [check]
--annotation          Require annotation? False by default, True if gbk/faa/ffn requested
--metadata_append     Append to metadata, not overwrite.
--typestrain          Only download type strains.
--keepmulti           By default, genomes from large multi-isolatestudies are removed.
--force               Force download attempt of low-quality genomes.
-f <{abbr,full,strain}, --outformat {abbr,full,strain}>
                      Output file prefix. [full]
-o <{fna,ffn,gff,gbk,faa,all} [{fna,ffn,gff,gbk,faa,all} ...]>
                      Output file types.
--edirect EDIRECT     Path to edirect directory.
--debug               Turn on debugging messages.

Assembly IDs
^^^^^^^^^^^^

.. code-block:: bash

    $ get_assemblies assembly_ids -h
    usage: get_assemblies assembly_ids [-h] [--type {acc,uid}]
                                       [--function {check,metadata,genomes} [{check,metadata,genomes} ...]] [--annotation]
                                       [--metadata_append] [--typestrain] [--keepmulti] [--force] [-f {abbr,full,strain}]
                                       [-o {fna,ffn,gff,gbk,faa,all} [{fna,ffn,gff,gbk,faa,all} ...]] [--edirect EDIRECT]
                                       [--debug]
                                       infile

    positional arguments:
      infile                Input file with NCBI assembly IDs; "-" for stdin

optional arguments:

-h, --help            show this help message and exit
--type <{acc,uid}>    Input is Accession (default) or ID
--function <{check,metadata,genomes} [{check,metadata,genomes} ...]>
                      check counts, download metadata, or genomes. [check]
--annotation          Require annotation? False by default, True if gbk/faa/ffn requested
--metadata_append     Append to metadata, not overwrite.
--typestrain          Only download type strains.
--keepmulti           By default, genomes from large multi-isolatestudies are removed.
--force               Force download attempt of low-quality genomes.
-f <{abbr,full,strain}, --outformat {abbr,full,strain}>
                      Output file prefix. [full]
-o <{fna,ffn,gff,gbk,faa,all} [{fna,ffn,gff,gbk,faa,all} ...]>
                      Output file types.
--edirect EDIRECT     Path to edirect directory.
--debug               Turn on debugging messages.

Nucleotide IDs
^^^^^^^^^^^^^^

.. code-block:: bash

    $ get_assemblies nuccore_ids -h
    usage: get_assemblies nuccore_ids [-h] [--function {check,metadata,genomes} [{check,metadata,genomes} ...]] [--annotation]
                                      [--metadata_append] [--typestrain] [--keepmulti] [--force] [-f {abbr,full,strain}]
                                      [-o {fna,ffn,gff,gbk,faa,all} [{fna,ffn,gff,gbk,faa,all} ...]] [--edirect EDIRECT] [--debug]
                                      infile

    positional arguments:
      infile                Input file with NCBI nuccore IDs; "-" for stdin

optional arguments:

-h, --help            show this help message and exit
--function <{check,metadata,genomes} [{check,metadata,genomes} ...]>
                      check counts, download metadata, or genomes. [check]
--annotation          Require annotation? False by default, True if gbk/faa/ffn requested
--metadata_append     Append to metadata, not overwrite.
--typestrain          Only download type strains.
--keepmulti           By default, genomes from large multi-isolatestudies are removed.
--force               Force download attempt of low-quality genomes.
-f <{abbr,full,strain}, --outformat {abbr,full,strain}>
                      Output file prefix. [full]
-o <{fna,ffn,gff,gbk,faa,all} [{fna,ffn,gff,gbk,faa,all} ...]>
                      Output file types.
--edirect EDIRECT     Path to edirect directory.
--debug               Turn on debugging messages.

NCBI JSON Docsum input
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    $ get_assemblies json_input -h
    usage: get_assemblies json_input [-h] [--function {metadata,genomes} [{metadata,genomes} ...]] [--annotation]
                                     [--metadata_append] [--typestrain] [--keepmulti] [--force] [-f {abbr,full,strain}]
                                     [-o {fna,ffn,gff,gbk,faa,all} [{fna,ffn,gff,gbk,faa,all} ...]] [--edirect EDIRECT] [--debug]
                                     jsonfile [jsonfile ...]

    positional arguments:
      jsonfile              Input JSON file with docsums; "-" for stdin

optional arguments:

-h, --help            show this help message and exit
--function <{metadata,genomes} [{metadata,genomes} ...]>
                      Download metadata and/or genomes. [metadata]
--annotation          Require annotation? False by default, True if gbk/faa/ffn requested
--metadata_append     Append to metadata, not overwrite.
--typestrain          Only download type strains.
--keepmulti           By default, genomes from large multi-isolatestudies are removed.
--force               Force download attempt of low-quality genomes.
-f <{abbr,full,strain}, --outformat {abbr,full,strain}>
                      Output file prefix. [full]
-o <{fna,ffn,gff,gbk,faa,all} [{fna,ffn,gff,gbk,faa,all} ...]>
                      Output file types.
--edirect EDIRECT     Path to edirect directory.
--debug               Turn on debugging messages.


Bugs
----

Viruses are currently not handled well, if at all. Look elsewhere to download
those.

Contributing
------------

Feel free to submit bug reports or pull requests so we can improve this
software. Undoubtedly there will be some erroneous prefixes generated out
there, and I'd like to fix them.

Author Contact
--------------

`Ed Davis <mailto:ed@cgrb.oregonstate.edu>`_

Acknowledgments
----------------

Special thanks for helping me test the software and get the python code packaged:

* `Alex Weisberg <https://github.com/alexweisberg>`_
* `Shawn O'Neil <https://github.com/oneilsh>`_

Also, thanks to these groups for supporting me through my scientific career:

* `OSU Chang Lab <https://github.com/osuchanglab>`_
* `Center for Genome Research and Biocomputing @ OSU <https://cgrb.oregonstate.edu>`_

License
-------

get_assemblies is distributed under the terms listed in the ``LICENSE`` file.
The software is free for non-commercial use.

Copyrights
----------

Copyright (c) 2020 Oregon State University

All Rights Reserved.
