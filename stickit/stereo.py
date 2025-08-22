from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions


def _exclude_labile_centers(mol, allow_pyramidal_N=False):
    mask = set()
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 7:
            if a.GetTotalDegree() == 4:  # ammonium (stable)
                continue
            if not allow_pyramidal_N:
                mask.add(a.GetIdx())
    return mask


def enumerate_stereo_filtered(mol, cfg):
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    mask = _exclude_labile_centers(mol, cfg['chem']['allow_pyramidal_N'])
    opts = StereoEnumerationOptions(tryEmbedding=False,
                                    maxIsomers=cfg['chem']['max_stereoisomers'],
                                    onlyUnassigned=True)
    for m in EnumerateStereoisomers(mol, options=opts):
        bad = False
        chiral_atoms = [a.GetIdx() for a in m.GetAtoms() if a.HasProp('_CIPCode')]
        if any(i in mask for i in chiral_atoms):
            bad = True
        if not bad:
            yield m


def stereo_key_for(mol):
    return Chem.MolToSmiles(mol, isomericSmiles=True)
