import re
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Union, cast

from rdkit import Chem

from .data import MopacCalcResult, STICSet, ThermoRow


def _resolve_mopac_exe(cfg: Dict[str, Any], fpath_mopac=None) -> str:
    """
    Return the first executable MOPAC path found, or None if not found.
    Precedence:
        1) explicit fpath_mopac argument
        2) cfg['mopac']['mopac_exe']
        3) environment variables STICKIT_MOPAC_EXE or MOPAC_EXE
        4) PATH lookup: 'mopac', 'MOPAC', 'mopac7', 'mopac2016', 'mopac2012'
        5) historical default '/opt/mopac/bin/mopac'
    """

    candidates: list[str] = []

    if fpath_mopac:
        candidates.append(fpath_mopac)

    mopac_opts = (cfg or {}).get("mopac", {})
    if isinstance(mopac_opts, dict) and mopac_opts.get("mopac_exe"):
        candidates.append(mopac_opts["mopac_exe"])

    for env_key in ("STICKIT_MOPAC_EXE", "MOPAC_EXE"):
        val = os.environ.get(env_key)
        if val:
            candidates.append(val)

    for name in ("mopac", "MOPAC", "mopac7", "mopac2016", "mopac2012"):
        found = shutil.which(name)
        if found:
            candidates.append(found)

    candidates.append("/opt/mopac/bin/mopac")

    for c in candidates:
        p = Path(c).expanduser()
        if p.exists() and os.access(str(p), os.X_OK):
            return str(p)
    return ""


def _have_mopac(cfg: Dict[str, Any], fpath_mopac=None) -> bool:
    """Return True if MOPAC is installed and available for execution."""
    return bool(_resolve_mopac_exe(cfg, fpath_mopac))


def _net_charge_from_formal(mol):
    """Calculates the net formal charge of a molecule by summing the formal charges of its atoms."""
    return sum(atom.GetFormalCharge() for atom in mol.GetAtoms())


def _write_mopac(mol, conf_id, charge, method, tempK, path):
    xyz = Chem.MolToXYZBlock(mol, confId=conf_id).splitlines()[2:]
    mol_name = f"{mol.GetProp('_Name')} : " if mol.HasProp("_Name") else ""
    preopt_keywords = (
        f"{method} CHARGE={charge} EPS=78.4 RSOLV=1.4 NOMM DIIS +\n"
        f"GNORM=2.0 LET GEO-OK\n"
        f"{mol_name}Geometry pre-optimization\n"
    )

    opt_keywords = (
        f"{method} CHARGE={charge} EPS=78.4 RSOLV=1.4 NOMM DIIS +\n"
        f"PRECISE DDMIN=0.0 GNORM=0.001 LET OLDGEO\n"
        f"{mol_name}Geometry optimization\n"
    )

    freq_keywords = (
        f"{method} CHARGE={charge} EPS=78.4 RSOLV=1.4 OLDGEO +\n"
        f"THERMO PRECISE LET\n{mol_name}Free energy calculation\n"
    )

    with open(path, "w") as f:
        f.write(preopt_keywords)
        f.write("\n")
        f.writelines("\n".join(xyz))
        f.write("\n\n")
        f.write(opt_keywords)
        f.write("\n")
        f.write(freq_keywords)


def _split_calculations(text: str) -> List[str]:
    # New calc starts at the MOPAC banner; version varies.
    return [
        blk.strip()
        for blk in re.split(
            r"(?=^\s*\*{5,}\s*$\n^\s*\*\*.*MOPAC v[0-9.]+.*$\n)", text, flags=re.M
        )
        if blk.strip()
    ]


def _parse_method(block: str) -> Optional[str]:
    m = re.search(r"\n\s*([A-Z0-9+\-]+(?:-[A-Z0-9]+)?)\s+CALCULATION RESULTS", block)
    if m:
        return m.group(1)
    m2 = re.search(r"\n\s*([A-Z0-9+\-]+)\s+CALCULATION\s", block)
    return m2.group(1) if m2 else None


def _parse_version(block: str) -> Optional[str]:
    m = re.search(r"MOPAC v([0-9.]+)", block)
    return m.group(1) if m else None


def _parse_results_block(block: str) -> MopacCalcResult | None:
    res = MopacCalcResult()
    res.method = _parse_method(block)
    res.version = _parse_version(block)
    res.scf_achieved = bool(re.search(r"SCF FIELD WAS ACHIEVED", block))
    res.job_ended_normally = bool(re.search(r"JOB ENDED NORMALLY", block))

    m = re.search(r"CHARGE ON SYSTEM\s*=\s*([+-]?\d+)", block)
    if m:
        res.charge = int(m.group(1))

    m = re.search(r"\bEPS\s*=\s*([0-9.]+)", block)
    if m:
        res.eps = float(m.group(1))

    # ---- FINAL HEAT OF FORMATION: take the **last** occurrence in this block ----
    # Accept "HEAT OF FORMATION = ... KCAL/MOL", "KCALS/MOLE", with/without "FINAL".
    hof_matches = re.findall(
        r"(?:FINAL\s+)?HEAT OF FORMATION\s*=\s*([-\d.]+)\s*KCAL(?:S)?/MOL(?:E)?",
        block,
        flags=re.I,
    )
    if hof_matches:
        res.heat_of_formation_kcal_mol = float(hof_matches[-1])
    # -----------------------------------------------------------------------------

    m = re.search(r"GRADIENT NORM\s*=\s*([0-9.]+)", block)
    if m:
        res.gradient_norm = float(m.group(1))

    m = re.search(r"IONIZATION POTENTIAL\s*=\s*([-\d.]+)\s*EV", block)
    if m:
        res.ionization_potential_ev = float(m.group(1))

    m = re.search(r"HOMO LUMO ENERGIES.*?=\s*([-\d.]+)\s+([-\d.]+)", block)
    if m:
        res.homo_lumo_ev = [float(m.group(1)), float(m.group(2))]

    m = re.search(
        r"MOLECULAR WEIGHT\s*=\s*([0-9.]+)\s+POINT GROUP:\s*([A-Za-z0-9]+)", block
    )
    if m:
        res.molecular_weight = float(m.group(1))
        res.point_group = m.group(2)

    # Rotational constants / principal moments (grab if present)
    m = re.search(
        r"ROTATIONAL CONSTANTS.*?\n\s*A\s*=\s*([-\d.]+)\s+B\s*=\s*([-\d.]+)\s+C\s*=\s*([-\d.]+)",
        block,
        flags=re.S | re.I,
    )
    if m:
        res.rotational_constants_cm_inv = {
            "A": float(m.group(1)),
            "B": float(m.group(2)),
            "C": float(m.group(3)),
        }

    m = re.search(
        r"PRINCIPAL MOMENTS OF INERTIA.*?\n\s*A\s*=\s*([-\d.]+)\s+B\s*=\s*([-\d.]+)\s+C\s*=\s*([-\d.]+)",
        block,
        flags=re.S | re.I,
    )
    if m:
        res.principal_moi_1e_minus40_g_cm2 = {
            "A": float(m.group(1)),
            "B": float(m.group(2)),
            "C": float(m.group(3)),
        }

    # Vibrational frequencies (many decks list individual lines with "FREQUENCY ...")
    if re.findall(r"\bFREQUENC(?:Y|IES)\b.*?([0-9.]+)", block):
        # Extract numeric tails; tolerant of extra text on same line
        if nums := re.findall(r"\bFREQUENCY\s+([0-9.]+)", block):
            res.vibrational_frequencies_cm_inv = [float(x) for x in nums]

    # ZERO POINT ENERGY: only if vibrational data present; accept with or **without** '='.
    if res.vibrational_frequencies_cm_inv:
        if m := re.search(
            r"ZERO\s+POINT\s+ENERGY(?:\s*=\s*|\s+)\s*([-\d.]+)\s*"
            r"(KCAL(?:S)?/MOL(?:E)?|KJ/MOL|EV)",
            block,
            flags=re.I,
        ):
            zpe_value = float(m[1])
            zpe_unit = str(m[2])
            res.zero_point_energy = {"value": zpe_value, "unit": zpe_unit}
        else:
            return None
    # Thermodynamics table → ΔG
    thermos: List[ThermoRow] = []
    pat = re.compile(
        r"^\s*([0-9]+(?:\.[0-9]+)?)\s+VIB\..*?(?:\n.*?)*?^\s*TOT\.\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)",
        flags=re.M | re.S,
    )
    for m in pat.finditer(block):
        T = float(m.group(1))
        hof = float(m.group(2))  # KCAL/MOL (from TOT row)
        H_cal = float(m.group(3))  # CAL/MOL
        Cp = float(m.group(4))  # CAL/K/MOL
        S = float(m.group(5))  # CAL/K/MOL
        G_kcal = (H_cal - T * S) / 1000.0
        thermos.append(
            ThermoRow(
                temperature_K=T,
                hof_kcal_mol=hof,
                enthalpy_cal_mol=H_cal,
                heat_capacity_cal_K_mol=Cp,
                entropy_cal_K_mol=S,
                gibbs_free_energy_kcal_mol=G_kcal,
            )
        )
    if thermos:
        res.thermo = thermos

    return res


def parse_mopac_output(path_or_text: Union[str, Path]) -> Dict[str, Any]:
    """
    Parse a MOPAC output deck (possibly with multiple calculations).
    Returns: {"calculations": [ { per-calc fields... }, ... ]}
    A calculation is “valid” if both are present:
      - "SCF FIELD WAS ACHIEVED"
      - "JOB ENDED NORMALLY"
    """
    if isinstance(path_or_text, (str, Path)) and Path(path_or_text).exists():
        text = Path(path_or_text).read_text(errors="replace")
    else:
        text = str(path_or_text)

    calcs: List[Dict[str, Any]] = []
    for blk in _split_calculations(text):
        if re.search(r"CALCULATION RESULTS|CALCULATION\s", blk):
            rb = _parse_results_block(blk)
            if rb is not None:
                calcs.append(rb.to_dict())
    return {"calculations": calcs}


def mopac_refine_and_prune(sticset: STICSet, cfg):
    mopac_opts = cfg.get("mopac", {})
    method = mopac_opts.get("method", "PM6-D3H4")
    tempK = mopac_opts.get("temperature", 298.15)
    dG_keep = mopac_opts.get("energy_window", 0.5)
    drop_imag = mopac_opts.get("delete_imaginary", True)
    mopac_exe = _resolve_mopac_exe(cfg)

    if not _have_mopac(cfg):
        raise RuntimeError("MOPAC not installed")

    for stic in sticset.stics:
        q = _net_charge_from_formal(stic.mol)
        conf_results = []
        with tempfile.TemporaryDirectory() as td:
            for conf in stic.conformers:
                input_deck = Path(td) / f"stic_{conf.conf_id}.mop"
                _write_mopac(stic.mol, conf.conf_id, q, method, tempK, input_deck)
                subprocess.run([mopac_exe, str(input_deck)], check=True)
                output_deck = input_deck.with_suffix(".out")
                parsed = parse_mopac_output(output_deck)

                calcs_any = parsed.get("calculations")
                if not isinstance(calcs_any, list) or not calcs_any:
                    continue
                results = cast(Dict[str, Any], calcs_any[-1])

                freqs_any = results.get("vibrational_frequencies_cm_inv")
                freqs: List[float] = (
                    [float(x) for x in freqs_any] if isinstance(freqs_any, list) else []
                )
                conf.n_imag_freq = sum(freq < 0 for freq in freqs)

                method_name = str(results.get("method") or "MOPAC")
                hof_any = results.get("heat_of_formation_kcal_mol")
                hof_val = (
                    float(hof_any)
                    if isinstance(hof_any, (int, float))
                    else float("nan")
                )
                energy_record: Dict[str, float] = {method_name: hof_val}
                conf.enthalpy = energy_record

                thermo_any = results.get("thermo")
                fe: Optional[float] = None
                if isinstance(thermo_any, list):
                    for r in thermo_any:
                        if (
                            isinstance(r, dict)
                            and r.get("temperature_K") == 298.0
                            and isinstance(
                                r.get("gibbs_free_energy_kcal_mol"), (int, float)
                            )
                        ):
                            fe = float(r["gibbs_free_energy_kcal_mol"])
                            break
                conf.free_energy = fe
                conf_results.append(conf)

        # prune: any negative frequencies → drop; then ΔG window
        keep = []
        finite = [c for c in conf_results if (c.free_energy is not None)]
        if drop_imag:
            finite = [c for c in finite if (c.n_imag_freq or 0) == 0]
        if finite:
            fes: List[float] = [cast(float, c.free_energy) for c in finite]
            gmin = min(fes)
            keep = [c for c in finite if cast(float, c.free_energy) - gmin <= dG_keep]
        stic.conformers = keep or conf_results
