from stickit.pipeline import stic_generation

def test_import_and_call():
    # just verify it imports and handles empty list
    out = stic_generation(molecules=[], config={
        "io":{"output_dir":"./_out","write_json":False,"write_sdf":False},
        "chem":{"ph":7.4,"dph":2.0,"max_ti_iters":1,"max_stereoisomers":8,"include_atropisomers":True,"allow_pyramidal_N":False},
        "conformers":{"engine":"rdkit","num_confs":2,"rmsd_thresh":0.5},
        "minimization":{"engine":"openmm","openff_forcefield":"openff-2.0.0.offxml","implicit_solvent":None,"max_steps":10,"tol_kj_per_mol_nm":10.0},
        "mopac":{"enable":False},
        "parallel":{"backend":"none"}
    })
    assert out == []

