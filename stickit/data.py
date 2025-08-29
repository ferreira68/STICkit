from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional

from rdkit import Chem


@dataclass(frozen=True)
class STICKey:
    parent_key: str
    tautomer_key: str
    ion_key: str
    stereo_key: str


@dataclass
class ConformerRecord:
    conf_id: int
    enthalpy: Dict[str, float] = field(default_factory=dict)
    free_energy: Optional[float] = None
    n_imag_freq: Optional[int] = None
    tags: Dict[str, str] = field(default_factory=dict)


@dataclass
class STIC:
    key: STICKey
    mol: Chem.Mol
    conformers: List[ConformerRecord]
    annotations: Dict[str, str] = field(default_factory=dict)


@dataclass
class STICSet:
    parent_key: str
    stics: List[STIC] = field(default_factory=list)

    def all_conformers(self):
        for s in self.stics:
            for c in s.conformers:
                yield s, c


@dataclass
class ThermoRow:
    temperature_K: float
    hof_kcal_mol: Optional[float] = None  # from TOT row (KCAL/MOL)
    enthalpy_cal_mol: Optional[float] = None  # CAL/MOL
    heat_capacity_cal_K_mol: Optional[float] = None
    entropy_cal_K_mol: Optional[float] = None  # CAL/K/MOL
    gibbs_free_energy_kcal_mol: Optional[float] = None  # ΔG = H − T·S (kcal/mol)


@dataclass
class MopacCalcResult:
    method: Optional[str] = None  # e.g., PM6-D3H4
    version: Optional[str] = None  # e.g., 22.0.6
    scf_achieved: bool = False  # “SCF FIELD WAS ACHIEVED”
    job_ended_normally: bool = False  # “JOB ENDED NORMALLY”
    charge: Optional[int] = None
    eps: Optional[float] = None
    heat_of_formation_kcal_mol: Optional[float] = None  # final HOF in this block
    gradient_norm: Optional[float] = None
    ionization_potential_ev: Optional[float] = None
    homo_lumo_ev: Optional[List[float]] = None  # [HOMO, LUMO]
    molecular_weight: Optional[float] = None
    point_group: Optional[str] = None
    rotational_constants_cm_inv: Optional[Dict[str, float]] = None  # A,B,C
    principal_moi_1e_minus40_g_cm2: Optional[Dict[str, float]] = None  # A,B,C
    zero_point_energy: Optional[Dict[str, float | str]] = None  # only if freqs present
    vibrational_frequencies_cm_inv: Optional[List[float]] = None
    thermo: Optional[List[ThermoRow]] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
