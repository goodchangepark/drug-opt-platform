from datetime import datetime
from typing import Optional

from pydantic import BaseModel, ConfigDict, Field


class ProjectCreate(BaseModel):
    name: str = Field(min_length=1, max_length=200)
    target: str = ""
    indication: str = ""
    mechanism_modality: str = ""
    description: str = ""


class ProjectUpdate(BaseModel):
    name: Optional[str] = None
    target: Optional[str] = None
    indication: Optional[str] = None
    mechanism_modality: Optional[str] = None
    description: Optional[str] = None


class ProjectOut(ProjectCreate):
    model_config = ConfigDict(from_attributes=True)
    id: int
    created_at: datetime
    updated_at: datetime
    compound_count: int = 0


class CompoundCreate(BaseModel):
    compound_id: str = Field(min_length=1, max_length=50)
    name: str = ""
    smiles: str = Field(alias="smiles", min_length=1)
    notes: str = ""


class CompoundUpdate(BaseModel):
    name: Optional[str] = None
    notes: Optional[str] = None
    smiles: Optional[str] = None
    change_note: str = "Structure update"


class CompoundOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)
    row_id: int
    project_id: int
    compound_id: str
    name: str
    notes: str
    current_version: int
    created_at: datetime
    updated_at: datetime
    version: dict
