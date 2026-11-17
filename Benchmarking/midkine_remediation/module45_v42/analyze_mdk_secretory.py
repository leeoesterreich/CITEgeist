"""Read-only summary of the Module 4/4b/5 run: MDK and secretory-gene loadings per program, 4b pairs of MDK programs, Module 5 groups.
    python analyze_mdk_secretory.py <run_dir>
Secretory set = SECRETORY of src/a12_10_real.py (the v42 MDK-secretory analysis). Rank = 1 is the largest W loading among all genes of that program.
"""
import pickle, sys
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
SEC = ["HSP90B1", "HSPA5", "CALR", "CANX", "PDIA4", "PDIA6", "SEC23A", "SEC61B", "ATF6", "MAN1A1", "XBP1"]
run = Path(sys.argv[1])
SPECS = sorted(p.stem for p in (run / "per_specimen").glob("*.pkl"))
rows, pairs = [], []
for s in SPECS:
    rec = pickle.loads((run / "per_specimen" / f"{s}.pkl").read_bytes())
    for t, r in rec["m4"].results_by_anchor.items():
        g = list(r.gene_names)
        gi = {x: i for i, x in enumerate(g)}
        W = np.asarray(r.W)
        for k, p in enumerate(r.programs):
            col = W[:, k]
            order = np.argsort(-col)
            rank = np.empty(len(col), int); rank[order] = np.arange(1, len(col) + 1)
            d = {"specimen": s, "cell_type": t, "program": k, "moran_i": p.spatial_moran_i, "moran_p": p.spatial_moran_pvalue,
                 "coherent": bool(p.spatial_moran_i > 0.15 and p.spatial_moran_pvalue < 0.01)}
            for name in ["MDK"] + SEC:
                if name in gi:
                    d[f"{name}_rank"] = int(rank[gi[name]]); d[f"{name}_w"] = float(col[gi[name]])
                else:
                    d[f"{name}_rank"] = np.nan; d[f"{name}_w"] = np.nan
            d["n_sec_top50"] = int(sum(1 for x in SEC if x in gi and rank[gi[x]] <= 50))
            d["n_sec_top200"] = int(sum(1 for x in SEC if x in gi and rank[gi[x]] <= 200))
            d["top10"] = ",".join(p.top_genes[:10])
            rows.append(d)
    for q in rec["m4b"].all_pairs:
        pairs.append({"specimen": s, "a1": q.anchor1, "p1": q.program1_id, "a2": q.anchor2, "p2": q.program2_id, "I": q.bivariate_morans_i,
                      "p": q.bivariate_pvalue, "pearson": q.pearson_correlation, "n": q.n_spots_used, "type": q.relationship_type})
df = pd.DataFrame(rows); df.to_csv(run / "mdk_secretory_loadings.csv", index=False)
pr = pd.DataFrame(pairs); pr.to_csv(run / "module4b_all_pairs.csv", index=False)
print("programs per type per specimen:"); print(df.groupby(["specimen", "cell_type"]).size().unstack().to_string())
print("coherent per specimen (all/cancer):"); print(df.groupby("specimen").coherent.sum().to_string()); print(df[df.cell_type == "Cancer_Cells"].groupby("specimen").coherent.sum().to_string())
print("\nCANCER programs with MDK in top50 or any secretory in top50:")
c = df[df.cell_type == "Cancer_Cells"]
print(c[(c.MDK_rank <= 50) | (c.n_sec_top50 > 0)][["specimen", "program", "moran_i", "moran_p", "MDK_rank", "n_sec_top50", "n_sec_top200", "HSPA5_rank", "PDIA6_rank", "CALR_rank", "SEC61B_rank"]].to_string())
print("\nMDK rank <=200 in any type:"); print(df[df.MDK_rank <= 200][["specimen", "cell_type", "program", "MDK_rank", "MDK_w", "moran_i", "n_sec_top50"]].to_string())
print("\nBest MDK rank per specimen (cancer):"); print(c.loc[c.groupby("specimen").MDK_rank.idxmin()][["specimen", "program", "MDK_rank", "MDK_w", "n_sec_top50", "n_sec_top200"]].to_string())
print("\nBest n_sec_top200 per specimen (cancer):"); print(c.loc[c.groupby("specimen").n_sec_top200.idxmax()][["specimen", "program", "n_sec_top50", "n_sec_top200", "MDK_rank"]].to_string())
print("\n4b: min p, #p<0.05, #FDR-significant per specimen"); print(pr.groupby("specimen").agg(minp=("p", "min"), n_p05=("p", lambda x: (x < .05).sum()), n_tested=("p", "size")).to_string())
