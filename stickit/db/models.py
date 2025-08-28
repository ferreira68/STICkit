from __future__ import annotations
from sqlalchemy import String, Integer, Float
from sqlalchemy.orm import Mapped, mapped_column
from .session import Base

class Molecule(Base):
    __tablename__ = "molecules"
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    compound_id: Mapped[str] = mapped_column(String(64), index=True, unique=True)
    smiles: Mapped[str] = mapped_column(String(4096))
    inchi: Mapped[str] = mapped_column(String(512), nullable=True)
    mw: Mapped[float] = mapped_column(Float, nullable=True)
