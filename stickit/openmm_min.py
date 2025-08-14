from typing import Dict, List
from rdkit import Chem
import numpy as np

def _to_openmm_context(mol, conf_id, ff_name):
    from openff.toolkit.topology import Molecule as OFFMol, Topology as OFFTop
    from openff.toolkit.typing.engines.smirnoff import ForceField as OFFForceField
    import openmm as mm
    from openmm import unit
    from openmm.app import Simulation, PDBFile

    offmol = OFFMol.from_rdkit(mol, allow_undefined_stereo=True)
    offtop = OFFTop.from_molecules([offmol])
    off_ff = OFFForceField(ff_name)
    system = off_ff.create_openmm_system(offtop)

    pdb_str = Chem.MolToPDBBlock(mol, confId=conf_id)
    pdb = PDBFile.fromPDBString(pdb_str)

    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    platform = mm.Platform.getPlatformByName("CPU")
    sim = Simulation(pdb.topology, system, integrator, platform)

    conf = mol.GetConformer(conf_id)
    xyz = np.array(conf.GetPositions(), dtype=float) / 10.0  # Å → nm
    sim.context.setPositions(xyz * unit.nanometer)
    return sim

def minimize_openmm(mol, conf_ids: List[int], cfg) -> Dict[int, float]:
    import openmm as mm
    from openmm import unit

    ff = cfg['minimization']['openff_forcefield']
    tol = cfg['minimization']['tol_kj_per_mol_nm'] * unit.kilojoule / unit.mole / unit.nanometer
    max_steps = cfg['minimization']['max_steps']

    energies = {}
    for cid in conf_ids:
        sim = _to_openmm_context(mol, cid, ff)
        mm.LocalEnergyMinimizer.minimize(sim.context, tolerance=tol, maxIterations=max_steps)
        state = sim.context.getState(getEnergy=True)
        e = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        energies[cid] = float(e)
    return energies

