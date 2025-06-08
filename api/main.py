from fastapi import FastAPI
from contextlib import asynccontextmanager
from api.routes import router
from api.database import create_tables

@asynccontextmanager
async def lifespan(app: FastAPI):
    create_tables()
    yield

app = FastAPI(
    title="DGG Variant Interpretation Engine",
    description="Clinical variant classification and interpretation API",
    version="1.0.0",
    lifespan=lifespan
)

app.include_router(router, prefix="/v1")
