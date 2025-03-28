# RESCRIPt

[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.3891931.svg)](https://doi.org/10.5281/zenodo.3891931)
 ![CI](https://github.com/bokulich-lab/RESCRIPt/actions/workflows/ci-dev.yaml/badge.svg)
 [![DOI](https://img.shields.io/badge/DOI-10.1371/journal.pcbi.1009581-B31B1B)](http://dx.doi.org/10.1371/journal.pcbi.1009581)
 <p align="left"><img src="logo.png" height="150" /></p>

RESCRIPt is a python 3 package to support a variety of operations for managing and curating reference sequence databases, DNA/RNA sequence data, and taxonomic data. See citation information below for a full benchmark and description.

## Install via QIIME 2.

RESCRIPt is now installed as part of [QIIME 2](https://qiime2.org/) and [MOSHPIT](https://moshpit.readthedocs.io), and you can install as outlined in the [QIIME 2 Documentation](https://docs.qiime2.org/), for:
 -  "Amplicon Distrubtion" version `2024.2` and later.
 -  "Metagenome Distribution" (MOSHPIT) version `2023.9` and later. 


### Prior versions of RESCRIPt / QIIME 2.

For details on how to install prior and minimal versions of RESCRIPt for the various versions of QIIME 2, please see the `install-prior-versions.md` document.

### Read help documentation

To view a help menu for using RESCRIPt via the QIIME 2 CLI:
```
qiime dev refresh-cache
qiime rescript --help
```

## Using RESCRIPt

These tutorials demonstrate some of the basic functionality of RESCRIPt, via the q2CLI (QIIME 2 command-line interface):
- [General Overview and working with SILVA data](https://forum.qiime2.org/t/rescript-sequence-reference-database-management-tutorial/15494)
- [Getting sequences and taxonomy with get-ncbi-data](https://forum.qiime2.org/t/using-rescript-to-compile-an-sequence-databases-and-taxonomy-classifiers-from-ncbi-genbank/15947)
- [Building a COI database with BOLD sequences](https://forum.qiime2.org/t/building-a-coi-database-from-bold-references/16129)
- [Building a COI database with NCBI sequences](https://forum.qiime2.org/t/building-a-coi-database-from-ncbi-references/16500)
- [Using RESCRIPt's 'extract-seq-segments' to extract reference sequences without PCR primer pairs](https://forum.qiime2.org/t/using-rescripts-extract-seq-segments-to-extract-reference-sequences-without-pcr-primer-pairs/23618)
- [How to train a GTDB SSU classifier using RESCRIPt](https://forum.qiime2.org/t/how-to-train-a-gtdb-ssu-classifier-using-rescript/25725)
- [Constructing an RDP classifier](https://forum.qiime2.org/t/importing-sequence-data-with-lower-case-nucleotide-characters-constructing-an-rdp-classifier-as-an-example/25158)
- [How to train a UNITE classifier using RESCRIPt](https://forum.qiime2.org/t/how-to-train-a-unite-classifier-using-rescript/28285)

Examples of visualizations produced by RESCRIPt actions can be found in this [Visualization Gallery](https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494#heading--seventeenth-header). Other code examples can be found [here](https://github.com/bokulich-lab/db-benchmarks-2020).

## Getting Help

Problem? Suggestion? Technical errors and user support requests can be filed on the [QIIME 2 Forum](https://forum.qiime2.org/).


## Citation

If you use RESCRIPt in your research, please cite the following:

Michael S Robeson II, Devon R O'Rourke, Benjamin D Kaehler, Michal Ziemski, Matthew R Dillon, Jeffrey T Foster, Nicholas A Bokulich. (2021) *RESCRIPt: Reproducible sequence taxonomy reference database management*. PLoS Computational Biology 17 (11): e1009581. doi: [10.1371/journal.pcbi.1009581](http://dx.doi.org/10.1371/journal.pcbi.1009581).


## License

RESCRIPt is released under a BSD-3-Clause license. See LICENSE for more details.

However, other resources accessible via RESCRIPt are released under different licenses, as detailed below.

**If using the SILVA database** (*e.g.*, with `get-silva-data`): Versions are released under different licenses. Refer to the [current SILVA release license information](https://www.arb-silva.de/silva-license-information/) for more details. [How to cite SILVA](https://www.arb-silva.de/contact/).

**If using NCBI Genbank data** (*e.g.*, with `get-ncbi-data`): See the [NCBI disclaimer and copyright notice](https://www.ncbi.nlm.nih.gov/home/about/policies/) for more details. [How to cite NCBI](https://support.nlm.nih.gov/knowledgebase/article/KA-03391/en-us).

**If using GTDB data** (*e.g.*, with `get-gtdb-data`): See the [GTDB "about" page](https://gtdb.ecogenomic.org/about) for more details. [How to cite GTDB](https://gtdb.ecogenomic.org/about).

**If using UNITE data** (*e.g.*, with `get-unite-data`): See the [UNITE citation page](https://unite.ut.ee/cite.php) for more details. [How to cite UNITE](https://unite.ut.ee/cite.php). 
