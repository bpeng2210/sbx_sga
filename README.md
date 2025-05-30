<img src="https://github.com/sunbeam-labs/sunbeam/blob/main/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_sga

<!-- Badges start -->
[![Tests](https://github.com/sunbeam-labs/sbx_sga/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_sga/actions/workflows/tests.yml)
![Condabot](https://img.shields.io/badge/condabot-active-purple)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_sga)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_sga/)
<!-- Badges end -->

## Introduction

sbx_sga (Single Genome Assembly) is a [sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for isolate QC, assembly, and classification. This pipeline uses [Mash](https://github.com/marbl/mash) for quality control, [Shovill](https://github.com/tseemann/shovill) for bacterial isolate assembly, [CheckM2](https://github.com/chklovski/CheckM2) and [QUAST](https://github.com/ablab/quast) for assembly QC, [MLST](https://github.com/tseemann/mlst) for typing, [Bakta](https://github.com/oschwengers/bakta) for annotation, [abriTAMR](https://github.com/MDU-PHL/abritamr) for AMR profiling, and [Sylph](https://github.com/bluenote-1577/sylph) for taxonomic classification.

## Config

  - mash_ref: the reference file for running Mash (should be a file ending in `.msh`)
  - checkm_ref: the diamond database for running CheckM2 (should be a file ending in `.dmnd`)
  - bakta_ref: the bakta reference database (should be a directory similar to `.../bakta_db/db/`)
    
## Docs

More [docs](https://sunbeam.readthedocs.io/en/stable/extensions.html).