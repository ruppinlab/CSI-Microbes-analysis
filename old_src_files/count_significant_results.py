import pandas as pd

sig_events = 0
events = 0
for input in snakemake.input:
    print(input)
    df = pd.read_csv(input, sep="\t", index_col=0)
    print(df)
    # if the microbe doesn't show up, we count it as not finding it to be DA
    try:
        if df.loc[snakemake.wildcards["microbe"]].FDR < .05:
            sig_events += 1
    except:
        pass
    events +=1
try:
    df = pd.DataFrame(data=[sig_events/events], columns=["sig_percent"], index=[snakemake.wildcards["ncells"]])
except:
    df = pd.DataFrame(data=[sig_events/events], columns=["sig_percent"], index=[snakemake.wildcards["nreads"]])
df.to_csv(snakemake.output[0], sep="\t")
