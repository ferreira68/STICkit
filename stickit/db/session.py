from __future__ import annotations
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, DeclarativeBase
from ..config.settings import Settings

class Base(DeclarativeBase): pass

_settings = Settings.load()
engine = create_engine(_settings.db.url, future=True)
SessionLocal = sessionmaker(bind=engine, expire_on_commit=False)

def init_db() -> None:
    from . import models as _  # ensure models imported
    Base.metadata.create_all(bind=engine)
