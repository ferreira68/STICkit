from pathlib import Path
from rdkit import Chem
import subprocess, tempfile
from .data import STICSet

def _net_charge_from_formal(mol):
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())

def _write_mopac(mol, conf_id, charge, method, tempK, path):
    xyz = Chem.MolToXYZBlock(mol, confId=conf_id).splitlines()[2:]
    keywords = f"{method} OPT FREQ THERMO({int(round(tempK))} K) CHARGE={charge}"
    with open(path, "w") as f:
        f.write(keywords+"\nSTIC conformer\n\n")
        for line in xyz:
            if not line.strip(): continue
            sym, x, y, z = line.split()
            f.write(f"{sym} {x} 1 {y} 1 {z} 1\n")

def _parse_out(out_path):
    d = {"n_imag": None, "G_kcal": None, "HOF_kcal": None}
    with open(out_path) as f:
        for s in f:
            if "NEGATIVE EIGENVALUES" in s:
                d["n_imag"] = int(s.split()[-1])
            elif "GIBBS FREE ENERGY" in s and "KCAL/MOL" in s:
                d["G_kcal"] = float(s.split()[4])
            elif "HEAT OF FORMATION" in s and "KCAL/MOL" in s:
                d["HOF_kcal"] = float(s.split()[4])
    return d

def mopac_refine_and_prune(sticset: STICSet, cfg):
    method = cfg['mopac']['method']
    tempK = cfg['mopac']['temperature_K']
    dG_keep = cfg['mopac']['keep_within_dG_kcal']
    drop_imag = cfg['mopac'].get('delete_imaginary', True)

    for stic in sticset.stics:
        q = _net_charge_from_formal(stic.mol)
        results = []
        with tempfile.TemporaryDirectory() as td:
            for c in stic.conformers:
                mop = Path(td)/f"stic_{c.conf_id}.mop"
                _write_mopac(stic.mol, c.conf_id, q, method, tempK, mop)
                subprocess.run(["mopac", str(mop)], check=True)
                out = mop.with_suffix(".out")
                d = _parse_out(out)
                c.n_imag_freq = d["n_imag"]
                c.method_energy["MOPAC_HOF"] = d["HOF_kcal"] if d["HOF_kcal"] is not None else None
                c.free_energy = d["G_kcal"]
                results.append(c)

        # prune: any negative frequencies → drop; then ΔG window
        keep = []
        finite = [c for c in results if (c.free_energy is not None)]
        if drop_imag:
            finite = [c for c in finite if (c.n_imag_freq or 0) == 0]
        if finite:
            gmin = min(c.free_energy for c in finite)
            keep = [c for c in finite if c.free_energy - gmin <= dG_keep]
        stic.conformers = keep if keep else results

