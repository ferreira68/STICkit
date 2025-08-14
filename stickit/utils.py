import yaml, json, os
from rdkit import Chem
from concurrent.futures import ProcessPoolExecutor, as_completed

def load_config(path): 
    with open(path) as f: 
        return yaml.safe_load(f)

def smiles_iter_from_file(path):
    with open(path) as f:
        for line in f:
            if not line.strip(): 
                continue
            parts = line.strip().split()
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else smi
            yield smi, name

def canonical_parent_key(mol):
    m = Chem.Mol(mol)
    [a.ClearProp('_CIPCode') for a in m.GetAtoms() if a.HasProp('_CIPCode')]
    return Chem.MolToSmiles(m, isomericSmiles=False)

def parallel_map(func, items, parallel_cfg):
    backend = (parallel_cfg or {}).get("backend", "none")
    if backend == "ray":
        import ray
        if not ray.is_initialized():
            ray.init(ignore_reinit_error=True, num_cpus=parallel_cfg.get("num_workers", None))
        @ray.remote
        def _wrap(x): 
            return func(x)
        return ray.get([_wrap.remote(x) for x in items])

    elif backend == "multiprocessing":
        n = parallel_cfg.get("num_workers", os.cpu_count() or 1)
        out = []
        with ProcessPoolExecutor(max_workers=n) as ex:
            futs = [ex.submit(func, it) for it in items]
            for f in as_completed(futs):
                out.append(f.result())
        return out
    else:
        return [func(x) for x in items]

def dump_outputs(sticsets, cfg):
    outdir = cfg['io']['output_dir']
    os.makedirs(outdir, exist_ok=True)
    if cfg['io'].get('write_json', True):
        payload = []
        for s in sticsets:
            payload.append({
                "parent_key": s.parent_key,
                "stics": [{
                    "key": {
                        "parent_key": st.key.parent_key,
                        "tautomer_key": st.key.tautomer_key,
                        "ion_key": st.key.ion_key,
                        "stereo_key": st.key.stereo_key,
                    },
                    "n_conformers": len(st.conformers),
                    "annotations": st.annotations
                } for st in s.stics]
            })
        with open(os.path.join(outdir, "summary.json"), "w") as f:
            json.dump(payload, f, indent=2)

    if cfg['io'].get('write_sdf', True):
        from rdkit.Chem import SDWriter
        sdf = os.path.join(outdir, "stickit_conformers.sdf")
        w = SDWriter(sdf)
        for s in sticsets:
            for st in s.stics:
                for c in st.conformers:
                    m = Chem.Mol(st.mol)
                    m.SetProp("_Name", st.annotations.get("name", s.parent_key))
                    w.write(m, confId=c.conf_id)
        w.close()

