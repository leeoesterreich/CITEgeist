"""Single source of truth for the breast cancer cohort hand marker dictionary.

``BREAST_PROFILE`` is the v41 dictionary. It is the v40 ``UNIFIED_PROFILE`` plus
one change: ``PAX5-1`` added to the ``B_Cells`` Major markers. The 7 types and
their order are the v40 ones.

A Neutrophils type (CEACAM8-1) was tried and dropped: on raw antibody UMIs
CEACAM8-1 does not exceed its matched isotype control (mouse IgG1 kappa) in any of
the 12 specimens, and neutrophil RNA genes show no signal.

``BREAST_PROFILE_V40`` is the unchanged v40 dictionary. It is kept for the
internal diagnostic and for protocols that do not declare
``inputs.dictionary_version == "v41"``.

Markers carry the ``-1`` suffix used by the real patient data.
"""

BREAST_PROFILE_V40 = {
    "Cancer_Cells": {
        "Major": ["EPCAM-1"],
        "Minor": ["SDC1-1", "KRT5-1"],
    },
    "Macrophages": {
        "Major": ["CD68-1"],
        "Minor": ["CD14-1"],
    },
    "CD4_T_Cells": {
        "Major": ["CD3E-1", "CD4-1"],
    },
    "CD8_T_Cells": {
        "Major": ["CD3E-1", "CD8A-1"],
    },
    "B_Cells": {
        "Major": ["MS4A1-1", "CD19-1"],
    },
    "Endothelial": {
        "Major": ["PECAM1-1"],
    },
    "Fibroblasts": {
        "Major": ["ACTA2-1"],
    },
}

BREAST_PROFILE = {
    "Cancer_Cells": {
        "Major": ["EPCAM-1"],
        "Minor": ["SDC1-1", "KRT5-1"],
    },
    "Macrophages": {
        "Major": ["CD68-1"],
        "Minor": ["CD14-1"],
    },
    "CD4_T_Cells": {
        "Major": ["CD3E-1", "CD4-1"],
    },
    "CD8_T_Cells": {
        "Major": ["CD3E-1", "CD8A-1"],
    },
    "B_Cells": {
        "Major": ["MS4A1-1", "CD19-1", "PAX5-1"],
    },
    "Endothelial": {
        "Major": ["PECAM1-1"],
    },
    "Fibroblasts": {
        "Major": ["ACTA2-1"],
    },
}
