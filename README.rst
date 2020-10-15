get_assemblies
==============

.. contents:: **Table of Contents**
    :backlinks: none

Installation
------------

get_assemblies is distributed on `PyPI <https://pypi.org>`_ as a universal
wheel and is available on Linux/macOS and supports Python 3.7+. This software
will work on Windows using `WSL <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_.

get_assemblies depends on `NCBI Entrez Direct <https://www.ncbi.nlm.nih.gov/books/NBK179288/>`_
which requires Perl. Perl is installed by default on most \*nix systems. If
edirect is not currently installed, please run ``get_assemblies --dledirect``
to install it.

.. code-block:: bash

    $ python3 -m pip install -U get-assemblies

Dependencies
------------

Python modules:

1. wget
2. tqdm

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

^ This will find all genomes tagged as 'Pseudomonas fluorescens' in NCBI's
database. By default, this command will only check to see how many genomes
fall into this category. To download metadata for the genomes:

.. code-block:: bash

    $ get_assemblies organism 'Pseudomonas fluorescens' --function metadata

Check out the ``metadata.tab`` file that is created after running this command.

Overview
--------

This tool was written to make accessing genomic data from NCBI easier. Although
downloading a single genome is relatively easy from the NCBI website is
relatively easy, how can we make it easy to download many genomes and have
those files named in a way that is clear what the contents are? Other tools
have been written to download genome files from NCBI, but getting them renamed
in a way that makes sense to a human has been a missing component for a while.

This software is essentially a wrapper for the NCBI edirect tools that makes
getting genome files easier. If you are interested in starting a comparative
genomics project, this is the tool for you.

If you need to make phylogenetic trees with this data, check out my other
python package, `automlsa2 <https://pypi.org/project/automlsa2/>`_.

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
