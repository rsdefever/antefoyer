antefoyer
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/antefoyer.svg?branch=master)](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/antefoyer)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/REPLACE_WITH_APPVEYOR_LINK/branch/master?svg=true)](https://ci.appveyor.com/project/REPLACE_WITH_OWNER_ACCOUNT/antefoyer/branch/master)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/antefoyer/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/antefoyer/branch/master)

This repository provides support for the General Amber Force Field (GAFF) in foyer. Foyer is a part of the Molecular Simulation Design Framework (MoSDeF). A few relevant links are provided below:

- [MoSDeF](https://mosdef.org)
- [Foyer GitHub](https://github.com/mosdef-hub/foyer)
- [Foyer paper](https://doi.org/10.1016/j.commatsci.2019.05.026) 
- [GAFF paper](https://doi.org/10.1002/jcc.20035)
- [Antechamber paper](https://doi.org/10.1016/j.jmgm.2005.12.005)

There are two primary components to this repository: 

- Foyer-compatible XML file for GAFF
- Antechamber wrapper to integrate antechamber atomtyping and AM1-BCC charge assignment directly into MoSDeF workflows.

## Installation
The following instructions assume you already have foyer installed. If you do not please install foyer __before__ proceeding with these instructions.

    git clone https://github.com/rsdefever/antefoyer
    cd antefoyer/.
    pip install . 

The `pip install` command should be issued from within the same `conda` environment as `foyer` is installed. 

## Dependencies

- foyer >= 0.7.2
- ambertools (for the `antefoyer.ante_atomtyping` and `antefoyer.ante_charges` functions)

## Usage
There are two primary workflows for atomtyping with `antefoyer`. The first uses the traditional `foyer` approach with SMARTS strings defined in `antefoyer/xml/gaff.xml`. 

    # Assuming we also have mbuild installed for this example
    import foyer
    import mbuild as mb
    
    # Build/import a molecule
    molecule = mb.load('CCC',smarts=True)
    
    # Load GAFF
    GAFF = foyer.forcefields.load_GAFF()
    
    # Use foyer to parameterize
    parmaeterized_molecule = GAFF.apply(molecule.to_parmed())
    
    # Save to whatever file format you desire
    
Note that `antefoyer` was _never_ imported. If `antefoyer` is properly installed, the GAFF forcefield will become available under `foyer.forcefields.load_GAFF()` via `entrypoints`. 

The second workflow uses the antechamber wrapper. This requires `antechamber` to be installed an accesible in your `PATH`. `Antechamber` can be installed from conda: `conda install -c ambermd ambertools`. __THIS WORKFLOW IS UNDER CONSTRUCTION AND NOT FULLY SUPPORTED BY FOYER__.


    # Assuming we also have mbuild installed for this example
    import foyer
    import mbuild as mb
    
    # For this we need to import antefoyer
    import antefoyer
    
    # Build/import a molecule
    molecule = mb.load('CCC',smarts=True)
    
    # Load GAFF
    GAFF = foyer.forcefields.load_GAFF()
    
    # Identify GAFF atomtypes with antechamber
    typed_molecule = antefoyer.ante_atomtyping(molecule,'gaff')
    
    # Use foyer to parameterize
    parmaeterized_molecule = GAFF.apply(typed_molecule, run_atomtyping=False)
    
    # We can also apply AM1-BCC charges from antechamber
    molecule_with_charges = antefoyer.ante_charges(parameterized_molecule,'bcc')
    
    # Save to whatever file format you desire

## Details of SMARTS strings development and testing
Further details of the SMARTS strings development and testing are provided in a separate repository, [GAFF-Foyer](https://github.com/rsdefever/GAFF-foyer).

## Acknowledgements
 
Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.


This material is based upon work supported by the National Science Foundation under Grant #1835874

