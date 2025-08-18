from typing import Iterable, List, Optional, Tuple

from rdkit import Chem

from .conformers import make_conformers
from .data import ConformerRecord, STIC, STICKey, STICSet
from .ionization import enum_ionization_states_dimorphite, ionization_key, ti_enables_new_tautomer_rules
from .mopac import mopac_refine_and_prune
from .openmm_min import minimize_openmm
from .stereo import enumerate_stereo_filtered, stereo_key_for
from .tautomer import enum_tautomers, tautomer_key
from .utils import canonical_parent_key, dump_outputs, parallel_map, smiles_iter_from_file


def _single_parent(smiles: str, name: str, cfg) -> Optional[STICSet]:
    try:
        parent = Chem.MolFromSmiles(smiles)
        if Chem.rdmolops.NeedsHs(parent):
            parent = Chem.rdmolops.AddHs(parent)
    except:
        return None

    parent_key = canonical_parent_key(parent)
    sticset = STICSet(parent_key=parent_key)

    # T
    tautomers = enum_tautomers(parent, cfg)

    # TI fixed-point loop
    seen_ti = set()
    frontier: List[Tuple[Chem.Mol, Optional[str]]] = [(t, None) for t in tautomers]
    for _ in range(cfg['chem']['max_ti_iters']):
        next_frontier: List[Tuple[Chem.Mol, Optional[str]]] = []
        for mol_T, _ in frontier:
            t_key = tautomer_key(mol_T)
            # I (Dimorphite-DL on the tautomer SMILES)
            smi_T = Chem.MolToSmiles(mol_T, isomericSmiles=False)
            ion_smiles = enum_ionization_states_dimorphite(smi_T, cfg['chem']['ph'], cfg['chem']['ph_delta'],
                                                           max_variants=256)
            for smi_TI in ion_smiles:
                mol_TI = Chem.MolFromSmiles(smi_TI)
                if mol_TI is None:
                    continue
                if Chem.rdmolops.NeedsHs(mol_TI):
                    mol_TI = Chem.rdmolops.AddHs(mol_TI)
                Chem.SanitizeMol(mol_TI)
                i_key = ionization_key(mol_TI)

                ti_key = (t_key, i_key)
                if ti_key in seen_ti:
                    continue
                seen_ti.add(ti_key)

                # charge-gated tautomerization
                if ti_enables_new_tautomer_rules(mol_TI, cfg):
                    for t2 in enum_tautomers(mol_TI, cfg, charge_aware=True):
                        next_frontier.append((t2, i_key))

                # S
                for mol_STI in enumerate_stereo_filtered(mol_TI, cfg):
                    s_key = stereo_key_for(mol_STI)
                    key = STICKey(parent_key, t_key, i_key, s_key)

                    # C (Gypsum or RDKit) + OpenMM minimization
                    mol3d, conf_ids = make_conformers(mol_STI, cfg)
                    mol3d, e_kcal = minimize_openmm(mol3d, cfg, conf_ids)
                    conformers = [ConformerRecord(conf_id=cid, enthalpy={"OpenMM": e_kcal[cid][1]}) for cid in conf_ids]

                    stic = STIC(key=key, mol=mol3d, conformers=conformers,
                                annotations={"name": name, "pH": str(cfg['chem']['ph'])})
                    sticset.stics.append(stic)
        frontier = next_frontier
        if not frontier:
            break

    # Optional MOPAC (auto-prunes negative frequencies)
    if cfg['mopac'].get('enable', False):
        mopac_refine_and_prune(sticset, cfg)

    return sticset


def stic_generation(
        molecules: Optional[Iterable[tuple]] = None,  # iterable of (smiles, name)
        config: Optional[dict] = None
) -> List[STICSet]:
    cfg = config or {}
    if molecules is None:
        molecules = smiles_iter_from_file(cfg['io']['input'])
    work = list(molecules)
    results = parallel_map(lambda pair: _single_parent(pair[0], pair[1], cfg), work, cfg['parallel'])
    sticsets = [r for r in results if r is not None]
    dump_outputs(sticsets, cfg)
    return sticsets
