"""SQLite database schema, connection management, and query helpers.

Uses SQLAlchemy 2.0 ORM with DeclarativeBase. Schema matches Paper 2's
galaxy_dynamics.db (galaxies, radial_profiles, model_fits).

Paper 3 adds model_name='constrained_rt' rows to model_fits.
"""

from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd
from sqlalchemy import (
    Float, Integer, String, DateTime, Boolean, ForeignKey,
    create_engine, delete,
)
from sqlalchemy.orm import (
    DeclarativeBase, Mapped, Session, mapped_column,
    relationship, sessionmaker,
)

from src.utils import get_db_path, setup_logger

logger = setup_logger(__name__)


# ---------------------------------------------------------------------------
# ORM Models
# ---------------------------------------------------------------------------

class Base(DeclarativeBase):
    pass


class Galaxy(Base):
    __tablename__ = "galaxies"

    galaxy_id: Mapped[str] = mapped_column(String, primary_key=True)
    distance_mpc: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    inclination: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    luminosity_band_36: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    quality_flag: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    r_disk_kpc: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    sb_disk: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    v_flat: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    data_source: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    profiles: Mapped[list["RadialProfile"]] = relationship(
        back_populates="galaxy", cascade="all, delete-orphan"
    )
    fits: Mapped[list["ModelFit"]] = relationship(
        back_populates="galaxy", cascade="all, delete-orphan"
    )


class RadialProfile(Base):
    __tablename__ = "radial_profiles"

    profile_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    galaxy_id: Mapped[str] = mapped_column(ForeignKey("galaxies.galaxy_id"))
    radius_kpc: Mapped[float] = mapped_column(Float)
    v_obs: Mapped[float] = mapped_column(Float)
    v_err: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    v_gas: Mapped[float] = mapped_column(Float)
    v_disk: Mapped[float] = mapped_column(Float)
    v_bulge: Mapped[float] = mapped_column(Float)
    v_baryon_total: Mapped[Optional[float]] = mapped_column(Float, nullable=True)

    galaxy: Mapped["Galaxy"] = relationship(back_populates="profiles")


class ModelFit(Base):
    """One row per (galaxy_id, model_name) fit.

    model_name values:
        'rational_taper'  -- free RT (k=2), param1=omega, param2=Rt
        'constrained_rt'  -- constrained RT (k=1), param1=omega, param2=Rt (derived)
        'nfw'             -- NFW halo (k=2), param1=c, param2=V_200
        'mond_fixed'      -- Fixed MOND (k=0), param1=a0
        'mond_free'       -- Free MOND (k=1), param1=a0_best
    """
    __tablename__ = "model_fits"

    fit_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    galaxy_id: Mapped[str] = mapped_column(ForeignKey("galaxies.galaxy_id"))
    model_name: Mapped[str] = mapped_column(String)
    n_params: Mapped[int] = mapped_column(Integer)
    chi_squared: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    reduced_chi_squared: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    bic: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    residuals_rmse: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    n_points: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    converged: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True)
    flag_v_obs_lt_v_bary: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True)
    method_version: Mapped[str] = mapped_column(String, default="v1")
    upsilon_disk: Mapped[float] = mapped_column(Float, default=0.5)
    upsilon_bulge: Mapped[float] = mapped_column(Float, default=0.7)
    param1: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    param1_err: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    param2: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    param2_err: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    timestamp: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)

    galaxy: Mapped["Galaxy"] = relationship(back_populates="fits")


# ---------------------------------------------------------------------------
# Engine / Session helpers
# ---------------------------------------------------------------------------

def get_engine(db_path: str | None = None):
    """Create a SQLAlchemy engine for the galaxy_dynamics database."""
    path = db_path or str(get_db_path())
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    return create_engine(f"sqlite:///{path}", echo=False)


def init_db(db_path: str | None = None):
    """Create all tables if they don't exist. Returns the engine."""
    engine = get_engine(db_path)
    Base.metadata.create_all(engine)
    return engine


def get_session(engine=None) -> Session:
    """Create and return a new SQLAlchemy Session."""
    if engine is None:
        engine = get_engine()
    return sessionmaker(bind=engine)()


# ---------------------------------------------------------------------------
# Query helpers
# ---------------------------------------------------------------------------

def query_profiles_as_dataframe(session: Session, galaxy_id: str) -> pd.DataFrame:
    """Retrieve all radial profiles for a galaxy, sorted by radius."""
    profiles = (
        session.query(RadialProfile)
        .filter(RadialProfile.galaxy_id == galaxy_id)
        .order_by(RadialProfile.radius_kpc)
        .all()
    )
    if not profiles:
        return pd.DataFrame()

    return pd.DataFrame([
        {
            "radius_kpc": p.radius_kpc,
            "v_obs": p.v_obs,
            "v_err": p.v_err,
            "v_gas": p.v_gas,
            "v_disk": p.v_disk,
            "v_bulge": p.v_bulge,
            "v_baryon_total": p.v_baryon_total,
        }
        for p in profiles
    ])


def query_fits_as_dataframe(
    session: Session,
    galaxy_id: Optional[str] = None,
    model_name: Optional[str] = None,
) -> pd.DataFrame:
    """Retrieve model fit results as a DataFrame."""
    query = session.query(ModelFit)
    if galaxy_id is not None:
        query = query.filter(ModelFit.galaxy_id == galaxy_id)
    if model_name is not None:
        query = query.filter(ModelFit.model_name == model_name)

    fits = query.all()
    if not fits:
        return pd.DataFrame()

    return pd.DataFrame([
        {
            "galaxy_id": f.galaxy_id,
            "model_name": f.model_name,
            "n_params": f.n_params,
            "chi_squared": f.chi_squared,
            "reduced_chi_squared": f.reduced_chi_squared,
            "bic": f.bic,
            "residuals_rmse": f.residuals_rmse,
            "n_points": f.n_points,
            "converged": f.converged,
            "param1": f.param1,
            "param1_err": f.param1_err,
            "param2": f.param2,
            "param2_err": f.param2_err,
            "method_version": f.method_version,
        }
        for f in fits
    ])


def insert_model_fit(session: Session, fit_result: dict) -> ModelFit:
    """Insert a new model fit result."""
    fit = ModelFit(**fit_result)
    session.add(fit)
    session.commit()
    return fit


def delete_fits(session: Session, galaxy_id: str = None, model_name: str = None):
    """Delete model fits matching the given filters."""
    stmt = delete(ModelFit)
    if galaxy_id is not None:
        stmt = stmt.where(ModelFit.galaxy_id == galaxy_id)
    if model_name is not None:
        stmt = stmt.where(ModelFit.model_name == model_name)
    result = session.execute(stmt)
    session.commit()
    return result.rowcount
