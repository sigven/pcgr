import os
from sqlmodel import Session, create_engine
from core.kb_models import SQLModel

DATABASE_URL = os.getenv("DATABASE_URL", "postgresql://localhost/dgg_rules")
engine = create_engine(DATABASE_URL)

def create_tables():
    """Create all database tables."""
    SQLModel.metadata.create_all(engine)

def get_session():
    """Database session dependency for FastAPI."""
    with Session(engine) as session:
        yield session
