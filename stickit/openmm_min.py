from typing import Any, Mapping, Optional

import openmm
from openff.units.openmm import from_openmm
from openmm import unit
from openmm.app.simulation import Simulation
from rdkit import Chem
from rdkit.Chem import Mol as ROMol

from .utils import pathlib_which


def _is_valid_forcefield(ff_name: str) -> bool:
    from openff.toolkit import get_available_force_fields
    return ff_name in get_available_force_fields()


def have_antechamber(deep: bool = False, timeout: float = 5.0) -> bool:
    """
    Return True if 'antechamber' is on PATH.
    If deep=True, also try launching 'antechamber -h' to confirm it runs.
    """
    path = pathlib_which("antechamber")
    if path is None:
        return False
    if not deep:
        return True
    try:
        cp = subprocess.run([path, "-h"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
        # Treat 'it started' as success; help may exit 0 or 1 depending on build
        return cp.returncode in (0, 1)
    except (OSError, subprocess.TimeoutExpired):
        return False


def _to_openmm_context(mol: ROMol, ff_name) -> Simulation:
    from openff.toolkit.topology import Molecule as OFFMol
    from openff.interchange import Interchange
    # from openff.toolkit.topology import Topology as OFFTop
    from openff.toolkit.typing.engines.smirnoff import ForceField as OFFForceField
    from openff.toolkit.utils import RDKitToolkitWrapper, AmberToolsToolkitWrapper

    offmol = OFFMol.from_rdkit(mol, allow_undefined_stereo=True)
    offtop = offmol.to_topology()
    if not _is_valid_forcefield(ff_name):
        raise ValueError(f"Invalid OpenFF forcefield name: '{ff_name}'")
    off_ff = OFFForceField(ff_name, load_plugins=True)

    # Check if AmberTools is available.  If not, fallback to other charge models
    charge_kwargs = {}
    if have_antechamber():
        AmberToolsToolkitWrapper().assign_partial_charges(offmol, partial_charge_method="am1bcc")
        charge_kwargs["charge_from_molecules"] = [offmol]
    elif offmol.partial_charges is not None:
        charge_kwargs["charge_from_molecules"] = [offmol]
    else:
        RDKitToolkitWrapper().assign_partial_charges(offmol, partial_charge_method="gasteiger")
        charge_kwargs["charge_from_molecules"] = [offmol]

    interchange = Interchange.from_smirnoff(off_ff, offtop, **charge_kwargs)
    # Use a simple integrator since we're just minimizing
    integrator = openmm.VerletIntegrator(1)
    simulation = interchange.to_openmm_simulation(integrator)
    return simulation


def minimize_openmm(
        mol: ROMol,
        cfg: Mapping[str, Any],
        conf_ids: Optional[list[int]] = None
) -> tuple[ROMol, dict[int, tuple[unit.Quantity, unit.Quantity]]]:
    from openff.toolkit import Molecule
    ff = cfg['minimization']['forcefield_spec']
    gradient = cfg['minimization']['gradient'] * unit.kilocalorie_per_mole / unit.angstrom
    max_steps = cfg['minimization']['max_steps']

    # Process a subset of conformers if requested
    sim_mol = Molecule.from_rdkit(mol)
    if conf_ids is not None:
        sim_mol.clear_conformers()
        conf_ids = [conf.GetId() for conf in mol.GetConformers()]
        for conf_id in conf_ids:
            temp_mol = Molecule.from_rdkit(Chem.Mol(mol, confId=conf_id))
            sim_mol.add_conformer(temp_mol.conformers[0])

    # Unit conversion for the minimization tolerance
    kJ_per_mol_nm = unit.kilojoules_per_mole / (unit.nano * unit.meter)
    tol = gradient.in_units_of(kJ_per_mol_nm)

    # Create a new molecule to hold the optimized conformers
    minimized_mol = Molecule.from_rdkit(mol)
    minimized_mol.conformers.clear()

    # Keep track of the energies
    initial_energies = []
    energies = []
    simulation = _to_openmm_context(mol, ff)
    for conformer in sim_mol.conformers:
        simulation.context.setPositions(conformer.to_openmm())
        initial_energies.append(simulation.context.getState(getEnergy=True).getPotentialEnergy())
        simulation.minimizeEnergy(tolerance=tol, maxIterations=max_steps)
        min_state = simulation.context.getState(getEnergy=True, getPositions=True)
        energies.append(min_state.getPotentialEnergy())
        minimized_mol.add_conformer(from_openmm(min_state.getPositions()))
        result_mol = minimized_mol.to_rdkit()

    initial_energies = [e.in_units_of(unit.kilocalorie_per_mole)._value for e in initial_energies]
    energies = [e.in_units_of(unit.kilocalorie_per_mole)._value for e in energies]
    conformer_energies = dict(enumerate(zip(initial_energies, energies), start=0))

    return result_mol, conformer_energies
