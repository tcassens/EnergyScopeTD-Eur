"""Common helpers shared across conversion_energyscope Snakemake scripts."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Mapping, Optional

import pandas as pd

# Ensure that the parent directory of the given path exists
def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


# Load a CSV file into a DataFrame with the first column as index
def load_indexed_csv(
    path: Path,
    *,
    non_numeric_columns: Optional[Iterable[str]] = None,
    coerce_numeric: bool = True,
    fillna_value: float = 0.0,
) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    if frame.empty:
        frame.index.name = frame.index.name or "country"
        return frame

    frame.index = frame.index.astype(str).str.strip()
    frame.index.name = frame.index.name or "country"
    if coerce_numeric:
        numeric_cols = [c for c in frame.columns if non_numeric_columns is None or c not in non_numeric_columns]
        cols = [c for c in numeric_cols if c in frame.columns]
        if cols:
            frame.loc[:, cols] = (
                frame.loc[:, cols]
                .apply(pd.to_numeric, errors="coerce")
                .fillna(fillna_value)
            )
    return frame


# Load a template configuration into a DataFrame, ensuring required metadata exists
def load_template(
    records: Optional[Iterable[Mapping[str, Any]]] = None,
    *,
    required_column: str = "resource",
) -> pd.DataFrame:
    template_records = list(records) if records is not None else []
    template = pd.DataFrame(template_records)
    if template.empty:
        raise ValueError("Template configuration must not be empty")
    if required_column not in template.columns:
        raise ValueError(f"Template must contain a '{required_column}' column")
    return template.reset_index(drop=True)


# Clean strings in a Series by stripping whitespace and replacing empty strings with NA
def clean_strings(series: pd.Series) -> pd.Series:
    cleaned = series.astype("string").str.strip()
    return cleaned.replace({"": pd.NA, pd.NA: pd.NA})

# Check if the code is a valid two-letter country code
def is_valid_country(code: str) -> bool: 
    code = str(code).strip()
    return len(code) == 2 and code.isalpha()

# Extract country code from a given string
def country_key(value: str) -> str: # Extract first two characters as country code
    value = str(value).strip()
    if not value:
        return ""
    return value[:2]

# Define the public API of this module
__all__ = ["ensure_parent", "load_indexed_csv", "load_template", "clean_strings"]
