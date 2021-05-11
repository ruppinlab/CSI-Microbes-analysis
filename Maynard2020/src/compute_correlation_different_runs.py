import pandas as pd
from scipy.stats import spearmanr
import statistics

out = []
for i in snakemake.input:
    out.append(pd.read_csv(i, sep="\t"))

auc_out = []
p_out = []
for i in range(0, len(out)):
    for j in range(i, len(out)):
        r = spearmanr(out[i].sort_index()["summary.AUC"], out[j].sort_index()["summary.AUC"])
        auc_out.append(r[0])
        r = spearmanr(out[i].sort_index()["p.value"], out[j].sort_index()["p.value"])
        p_out.append(r[0])

print(min(auc_out))
print(max(auc_out))
print(statistics.mean(auc_out))

print(min(p_out))
print(max(p_out))
print(statistics.mean(p_out))
