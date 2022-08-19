NPDtools
Natural Product Discovery tools
Version: see VERSION.txt
License: see LICENSE.txt

Developed in Saint Petersburg State University, St. Petersburg, Russia
Developed in University of California San Diego, La Jolla, CA, USA
Developed in Carnegie Mellon University, Pittsburgh, PA, USA

The package contains:
* DEREPLICATOR: in silico identification of peptidic natural products through database search of mass spectra
* VarQuest: modification-tolerant identification of novel variants of peptidic antibiotics and other natural products
* MetaMiner: a peptidogenomics approach for the discovery of ribosomally synthesized and post-translationally modified peptides
* Dereplicator+: indentification of metabolites through database search of mass spectra
* NPS: scoring and evaluating the statistical significance of peptidic natural product–spectrum matches (a part of Dereplicator and VarQuest)


Usage examples:
    ./bin/dereplicator.py \
            --db-path share/npdtools/test_data/sample_database/ \
            -o dereplicator_test_output_dir \
            share/npdtools/test_data/dereplicator/

    ./bin/varquest.py \
            --db-path share/npdtools/test_data/sample_database/ \
            -o varquest_test_output_dir \
            share/npdtools/test_data/varquest/

    ./bin/metaminer.py \
            -o metaminer_test_output_dir \
            share/npdtools/test_data/metaminer/ \
            -s share/npdtools/test_data/metaminer/

    ./bin/dereplicator+.py \
            --db-path share/npdtools/test_data/sample_database/ \
            -o dereplicator+_test_output_dir \
            share/npdtools/test_data/dereplicator+/

For the full list of available options please run
    ./bin/dereplicator.py --help
    ./bin/varquest.py --help
    ./bin/metaminer.py --help
    ./bin/dereplicator+.py --help


Output:
    significant_matches.tsv          list of the most reliable metabolite-spectrum matches
    significant_unique_matches.tsv   list of the most reliable unique metabolite identifications
    summary.tsv                      short summary report


System requirements:
- Linux or macOS
- Python 2.6-2.7 or Python 3.3+


Contacts and additional information:
    npdtools.support@cab.spbu.ru
    http://cab.spbu.ru/software/dereplicator
    http://cab.spbu.ru/software/varquest
    http://cab.spbu.ru/software/metaminer
    http://cab.spbu.ru/software/dereplicator-plus
    http://cab.spbu.ru/software/nps


References: see LICENSE.txt
