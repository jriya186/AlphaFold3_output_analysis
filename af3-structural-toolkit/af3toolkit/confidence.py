"""
Pick the best-scoring AF3 model from a folder of outputs, by highest
iPTM then highest pTM among ties.
"""

import os
import glob
import json


def get_best_models_by_iptm_ptm(folder_path: str) -> list:
    """
    Scan folder_path recursively for AF3 *summary_confidences*.json files,
    and return the model(s) with the highest iPTM (ties broken by pTM).

    Returns
    -------
    list of dicts: [{"file": ..., "iptm": ..., "ptm": ...}, ...]
    Empty list if no valid models are found.
    """
    pattern = os.path.join(folder_path, "**", "*summary_confidences*.json")
    all_files = glob.glob(pattern, recursive=True)

    candidates = []
    for file in all_files:
        try:
            with open(file, "r") as f:
                data = json.load(f)
            iptm = data.get("iptm")
            ptm = data.get("ptm")
            if iptm is not None and ptm is not None:
                candidates.append({"file": file, "iptm": iptm, "ptm": ptm})
        except Exception as e:
            print(f"Error reading {file}: {e}")

    if not candidates:
        print("No valid models with both iPTM and pTM found.")
        return []

    max_iptm = max(c["iptm"] for c in candidates)
    top_iptm_models = [c for c in candidates if c["iptm"] == max_iptm]

    max_ptm = max(c["ptm"] for c in top_iptm_models)
    best_models = [c for c in top_iptm_models if c["ptm"] == max_ptm]

    return best_models
