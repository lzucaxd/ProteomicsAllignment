"""Raw shared matrix method — no harmonization baseline."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from harmonize.methods.base import MethodInterface, MethodResult
from harmonize.utils.io import intersect_genes
from harmonize.preprocessing.metadata import build_feature_meta

logger = logging.getLogger(__name__)


class RawMethod(MethodInterface):
    name = "raw"
    display_name = "Raw"

    def run(
        self,
        cptac_matrix: pd.DataFrame,
        ccle_matrix: pd.DataFrame,
        sample_meta: pd.DataFrame,
        config: dict[str, Any],
        outdir: str | Path,
    ) -> MethodResult:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        shared = intersect_genes(cptac_matrix, ccle_matrix)
        combined = pd.concat(
            [cptac_matrix.loc[shared], ccle_matrix.loc[shared]], axis=1
        )

        valid_sids = set(sample_meta["sample_id"])
        combined = combined[[c for c in combined.columns if c in valid_sids]]

        feature_meta = build_feature_meta(shared)

        notes = (
            f"Raw shared matrix: {len(shared)} genes, {len(combined.columns)} samples.\n"
            "No cross-domain transformation applied.\n"
            "This is the unaligned baseline."
        )
        logger.info("Raw: %d genes x %d samples", len(shared), len(combined.columns))

        result = MethodResult(
            matrix=combined,
            sample_meta=sample_meta,
            feature_meta=feature_meta,
            method_name=self.name,
            display_name=self.display_name,
            notes=notes,
            transforms_values=False,
            value_scale="log2_abundance",
        )
        result.save(outdir)
        return result
