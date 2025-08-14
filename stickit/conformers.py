from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess, tempfile
from pathlib import Path

def _rdkit_confs(mol, num_confs, rmsd_thresh):
    m = Chem.AddHs(mol, addCoords=True)
    ps = AllChem.ETKDGv3()
    ps.pruneRmsThresh = rmsd_thresh
    ids = AllChem.EmbedMultipleConfs(m, numConfs=num_confs, params=ps)
    AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=200)
    return m, list(ids)

def _gypsum_confs(mol, num_confs, rmsd_thresh):
    smi = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    with tempfile.TemporaryDirectory() as td:
        inp = Path(td)/"inp.smi"
        outd = Path(td)/"out"
        inp.write_text(f"{smi}\tSTIC\n")
        cmd = [
            "gypsum_dl","--source",str(inp),"--output_folder",str(outd),
            "--add_hydrogens","True",
            "--enumerate_protomers","False","--enumerate_tautomers","False","--enumerate_enantiomers","False",
            "--max_variants_per_compound","1",
            "--max_conf",str(num_confs),
            "--rmsd_thresh",str(rmsd_thresh),
            "--precise_bond_geo","True","--thoroughness","3"
        ]
        subprocess.run(cmd, check=True)
        sdf = list(outd.rglob("*.sdf"))
        suppl = [Chem.SDMolSupplier(str(p), removeHs=False) for p in sdf]
        mols = [m for sup in suppl for m in sup if m]
        if not mols:
            raise RuntimeError("Gypsum produced no molecules.")
        base = Chem.Mol(mols[0])
        ids = [c.GetId() for c in base.GetConformers()]
        for m2 in mols[1:]:
            for c in m2.GetConformers():
                base.AddConformer(c, assignId=True)
                ids.append(c.GetId())
        return base, ids

def make_conformers(mol, cfg):
    engine = cfg['conformers']['engine']
    if engine == "gypsum":
        return _gypsum_confs(mol, cfg['conformers']['num_confs'], cfg['conformers']['rmsd_thresh'])
    else:
        return _rdkit_confs(mol, cfg['conformers']['num_confs'], cfg['conformers']['rmsd_thresh'])

