from dataclasses import dataclass, field
from typing import List, Optional, Dict
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
    method_energy: Dict[str, float] = field(default_factory=dict)
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

