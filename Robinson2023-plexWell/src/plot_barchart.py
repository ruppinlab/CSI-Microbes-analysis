import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

top_df = pd.read_csv("output/P1/summed_microbe_reads.tsv", sep="\t")
sns.barplot(data=top_df, x="name", y="reads")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
