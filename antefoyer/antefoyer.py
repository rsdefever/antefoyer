from __future__ import division

import os
import sys
import warnings

import parmed as pmd
import networkx as nx

from distutils.spawn import find_executable
from subprocess import PIPE, Popen
from antefoyer.utils.tempdir import temporary_directory
from antefoyer.utils.tempdir import temporary_cd

from foyer.exceptions import FoyerError
from foyer.utils.io import import_, has_mbuild

ANTECHAMBER = find_executable("antechamber")


def ante_atomtyping(molecule, atype_style):
    """Perform atomtyping by calling antechamber

    Parameters
    ----------
    molecule : parmed.Structure or mbuild.Compound
        Molecular structure to perform atomtyping on
    atype_style : str
        Style of atomtyping. Options include 'gaff', 'gaff2',
        'amber', 'bcc', 'sybyl'.

    Returns
    -------
    typed_molecule : parmed.Structure
        The molecule with antechamber atomtyping applied
    """
    _check_antechamber(ANTECHAMBER)

    # Check valid atomtype name
    supported_atomtypes = ["gaff", "gaff2", "amber", "bcc", "sybyl"]
    if atype_style not in supported_atomtypes:
        raise FoyerError(
            "Unsupported atomtyping style requested. "
            "Please select from {}".format(supported_atomtypes)
        )

    # Check for parmed.Structure. Convert from mbuild.Compound if possible
    molecule = _check_structure(molecule)
    # Confirm single connected molecule
    _check_single_molecule(molecule)

    # Get current directory to write any error logs
    workdir = os.getcwd()
    # Work within a temporary directory
    # to clean up after antechamber
    with temporary_directory() as tmpdir:
        with temporary_cd(tmpdir):
            # Save the existing molecule to file
            _write_pdb(molecule, "ante_in.pdb")
            # Call antechamber
            command = (
                "antechamber -i ante_in.pdb -fi pdb "
                "-o ante_out.mol2 -fo mol2 " + "-at " + atype_style + " " + "-s 2"
            )

            proc = Popen(
                command, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True
            )

            out, err = proc.communicate()

            # Error handling here
            if "Fatal Error" in err or proc.returncode != 0:
                _antechamber_error(out, err, workdir)

            # Now read in the mol2 file with atomtyping
            typed_molecule = pmd.load_file("ante_out.mol2", structure=True)

    # Foyer requires that the atom type info is stored under atom.id
    for atom in typed_molecule:
        atom.id = atom.type

    # And return it
    return typed_molecule


def ante_charges(
    molecule, charge_style, net_charge=0.0, multiplicity=1, charge_tol=0.005
):
    """Calculates partial charges by calling antechamber

    Parameters
    ----------
    molecule : parmed.Structure or mbuild.Compound
        Molecular structure to perform atomtyping on
    charge_style : str
        Style of partial charges calculation. Options include
        'bcc', 'gas', and 'mul'. See antechamber documentation
        by running 'antechamber -L' for details.
    net_charge : float, optional, default=0.0
        Net charge of the molecule
    multiplicity : int, optional, default=1
        Spin multiplicity, 2S + 1

    Returns
    -------
    molecule : parmed.Structure
        The molecule with charges applied
    """
    _check_antechamber(ANTECHAMBER)

    # Check valid atomtype name
    supported_chargetypes = ["bcc", "gas", "mul"]
    if charge_style not in supported_chargetypes:
        raise FoyerError(
            "Unsupported charge style requested. "
            "Please select from {}".format(supported_chargetypes)
        )

    # Check for parmed.Structure. Convert from mbuild.Compound if possible
    molecule = _check_structure(molecule)
    # Confirm single connected molecule
    _check_single_molecule(molecule)

    # Get current directory to write any error logs
    workdir = os.getcwd()
    # Work within a temporary directory
    # to clean up after antechamber
    with temporary_directory() as tmpdir:
        with temporary_cd(tmpdir):
            # Save the existing molecule to file
            _write_pdb(molecule, "ante_in.pdb")
            # Call antechamber
            command = (
                "antechamber -i ante_in.pdb -fi pdb "
                "-o ante_out.mol2 -fo mol2 "
                + "-c "
                + charge_style
                + " "
                + "-nc "
                + str(net_charge)
                + " "
                + "-m "
                + str(multiplicity)
                + " "
                + "-s 2"
            )

            proc = Popen(
                command, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True
            )

            out, err = proc.communicate()

            # Error handling here
            if "Fatal Error" in err or proc.returncode != 0:
                _antechamber_error(out, err, workdir)

            # Now read in the mol2 file with atomtyping
            charges = pmd.load_file("ante_out.mol2", structure=True)

            charge_list = [atom.charge for atom in charges]
            if abs(net_charge - sum(charge_list)) > charge_tol:
                raise ValueError(
                    "The sum of charges defined by antechamber"
                    " is {}, which differs from the desired net charge"
                    " of {} by a value greater than {}".format(
                        sum(charge_list), net_charge, charge_tol
                    )
                )
            elif abs(net_charge - sum(charge_list)) < charge_tol:
                charge_delta = (net_charge - sum(charge_list)) / len(charge_list)
                for atom in charges:
                    atom.charge += charge_delta

    # Combine charge information with existing molecule structure
    assert len(molecule.atoms) == len(charges.atoms)
    for atom_idx in range(len(molecule.atoms)):
        assert molecule.atoms[atom_idx].element == charges.atoms[atom_idx].element
        molecule.atoms[atom_idx].charge = charges.atoms[atom_idx].charge
    return molecule


def _write_pdb(molecule, filename):
    """Write a pdb file with CONECT records."""

    # Check that we have a box. If not, define one.
    if molecule.box is None:
        molecule.box = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]

    assert len(molecule.box) == 6, "Invalid box object."

    # Generate CONECT records
    conect = {atom.idx: [] for atom in molecule.atoms}
    for atom in molecule.atoms:
        atidx = atom.idx
        for bond in molecule.bonds:
            if bond.atom1.idx == atidx:
                conect[atidx].append(bond.atom2.idx)
            elif bond.atom2.idx == atidx:
                conect[atidx].append(bond.atom1.idx)

    with open(filename, "w") as pdb:
        pdb.write("REMARK 1   Created by antefoyer\n")
        pdb.write(
            "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}"
            " {:9s}{:3d}\n".format(
                molecule.box[0],
                molecule.box[1],
                molecule.box[2],
                molecule.box[3],
                molecule.box[4],
                molecule.box[5],
                molecule.space_group,
                1,
            )
        )
        for atom in molecule.atoms:
            pdb.write(
                "ATOM  {:5d} {:4s} RES A{:4d}    "
                "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}"
                "          {:>2s}  \n".format(
                    atom.idx + 1,
                    atom.name,
                    0,
                    atom.xx,
                    atom.xy,
                    atom.xz,
                    1.0,
                    0.0,
                    atom.element_name,
                )
            )
        for atidx, atomlist in conect.items():
            pdb.write("CONECT{:5d}".format(atidx + 1))
            for at2idx in atomlist:
                pdb.write("{:5d}".format(at2idx + 1))
            pdb.write("\n")


def _check_structure(molecule):
    """ Confirm that input is parmed.Structure. Convert
    from mbuild.Compound to parmed.Structure if possible.
    """
    if not isinstance(molecule, pmd.Structure) and has_mbuild:
        mb = import_("mbuild")
        if isinstance(molecule, mb.Compound):
            molecule = molecule.to_parmed()

    if not isinstance(molecule, pmd.Structure):
        raise FoyerError(
            "Unknown molecule format: {}\n"
            "Supported formats are: "
            '"parmed.Structure" and '
            '"mbuild.Compound"'.format(molecule)
        )

    return molecule


def _check_single_molecule(molecule):
    """ Confirms that the parmed structure represents a single
    connect molecule with connectivity info present.
    """
    graph_nodes = []
    graph_edges = []
    for atom in molecule.atoms:
        graph_nodes.append(atom.idx)
    for bond in molecule.bonds:
        graph_edges.append([bond.atom1.idx, bond.atom2.idx])
    bond_graph = nx.Graph()
    bond_graph.add_edges_from(graph_edges)
    bond_graph.add_nodes_from(graph_nodes)
    if not nx.is_connected(bond_graph):
        raise FoyerError(
            "Antechamber requires connectivity information and "
            "only supports single molecules (i.e., all atoms "
            "in the molecule are connected by bonds."
        )


def _antechamber_error(out, err, workdir):
    """Log antechamber output to file. """
    with open(workdir + "/ante_errorlog.txt", "w") as log_file:
        log_file.write("STDOUT:\n\n")
        log_file.write(out)
        log_file.write("STDERR:\n\n")
        log_file.write(err)
    raise RuntimeError("Antechamber failed. See 'ante_errorlog.txt'")


def _check_antechamber(ANTECHAMBER):
    if not ANTECHAMBER:
        msg = (
            "Antechamber not found. Please ensure that antechamber "
            "is available. It can be installed via conda, "
            "'conda install -c conda-forge ambertools'"
        )
        raise IOError(msg)
