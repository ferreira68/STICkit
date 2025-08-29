from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Union

from rdkit import Chem


def apply_ph_from_config(
    config: Dict[str, Any], gypsum_opts: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Update `gypsum_opts` with pH-related options from `config['chem']`, if present.
    Mirrors the original snippet's behavior: only updates when 'chem' exists/truthy.

    Returns the (possibly) updated `gypsum_opts`.
    """
    if cfg_ops := config.get("chem"):
        base_ph = cfg_ops.get("ph", 7.4)
        ph_delta = cfg_ops.get("ph_delta", 1.0)
        ph_tol = cfg_ops.get("ph_tol", 0.75)
        gypsum_opts.update(
            min_ph=(base_ph - ph_delta),
            max_ph=(base_ph + ph_delta),
            pka_precision=ph_tol,
        )
    return gypsum_opts


def enum_tautomers(
    mols: Union[Chem.Mol, Sequence[Chem.Mol]],
    config: Optional[Dict[str, Any]] = None,
    charge_aware: bool = True,
) -> List[Chem.Mol]:
    """
    Enumerate tautomeric and stereoisomeric states using Gypsum-DL's SMILES stage.
    No file I/O and no 3D generation. Returns RDKit molecules per input name.

    Parameters
    ----------
    mols : Chem.Mol or sequence of Chem.Mol
        Input molecules. Names are taken from the RDKit property "_Name" when present.
    config : dict, optional
        Gypsum-DL options to override defaults (see defaults below). Notable keys:
          - "job_manager": "multiprocessing" | "serial" | "mpi"   (default: "multiprocessing")
          - "num_processors": int (default: -1, like CLI = use all)
          - "max_variants_per_compound": int (default: 8)
          - "thoroughness": int (default: 3)
          - "min_ph", "max_ph", "pka_precision"
          - "let_tautomers_change_chirality", "use_durrant_lab_filters"
          - "skip_making_tautomers", "skip_enumerate_chiral_mol", "skip_enumerate_double_bonds"
    charge_aware : bool, optional
        Whether to perform ionization state enumeration (default: True).
    """
    # Lazy imports from Gypsum-DL
    from gypsum_dl.MolContainer import MolContainer
    from gypsum_dl.Steps.SMILES.PrepareSmiles import prepare_smiles
    from gypsum_dl.Parallelizer import Parallelizer

    # Normalize inputs to a list
    mol_list = [mols] if isinstance(mols, Chem.Mol) else list(mols)

    # Build a default configuration (matches CLI semantics where practical)
    gypsum_opts: Dict[str, Any] = {
        "job_manager": "multiprocessing",
        "num_processors": -1,  # like CLI: -1 => use all cores
        "max_variants_per_compound": 10,
        "thoroughness": 16,
        "min_ph": 6.4,
        "max_ph": 8.4,
        "pka_precision": 0.75,
        "let_tautomers_change_chirality": True,
        "use_durrant_lab_filters": True,
        "skip_adding_hydrogen": False,
        "skip_making_tautomers": False,
        "skip_enumerate_chiral_mol": False,
        "skip_enumerate_double_bonds": False,
        "skip_optimize_geometry": True,
        # Keep anything downstream from trying 3D/output if you expand usage later
        "2d_output_only": True,
        "separate_output_files": False,
        "add_pdb_output": False,
        "add_html_output": False,
    }

    # Set non-default options from config
    if not charge_aware:
        gypsum_opts["skip_adding_hydrogen"] = True
    if config:
        gypsum_opts = apply_ph_from_config(config, gypsum_opts)

    # Parallelizer (Gypsum-DL uses this object in the SMILES stage too)
    gypsum_opts["Parallelizer"] = Parallelizer(
        gypsum_opts["job_manager"], gypsum_opts["num_processors"], True
    )

    # Convert RDKit mols to MolContainers (Gypsum-DL SMILES stage input is SMILES)
    def _name(i: int, m: Chem.Mol) -> str:
        return m.GetProp("_Name").strip() if m.HasProp("_Name") else f"mol{i + 1:02d}"

    containers: List[MolContainer] = []
    for i, m in enumerate(mol_list):
        smi = Chem.MolToSmiles(m, isomericSmiles=True)
        containers.append(MolContainer(smi, _name(i, m), i, {}))

    # Run ONLY the SMILES-stage enumeration (tautomers/stereo/ionization)
    prepare_smiles(containers, gypsum_opts)

    # Collect RDKit molecules from the Gypsum containers
    # out: Dict[str, List[Chem.Mol]] = {}
    for cont in containers:
        variants: List[Chem.Mol] = []
        for obj in getattr(cont, "mols", []):
            mol = getattr(obj, "rdkit_mol", None) or getattr(obj, "mol", None)
            if mol is None and hasattr(obj, "get_rdkit_mol"):
                mol = obj.get_rdkit_mol()
            if mol is None and hasattr(obj, "smiles"):
                mol = Chem.MolFromSmiles(obj.smiles)
            if mol is not None:
                variants.append(Chem.Mol(mol))  # copy to detach from Gypsum internals

    for idx, mol in enumerate(variants, start=1):
        mol.SetProp("_Name", f"{cont.name}-T{idx:02d}")

    return variants


def tautomer_key(mol):
    return Chem.MolToSmiles(Chem.AddHs(mol), isomericSmiles=False)
