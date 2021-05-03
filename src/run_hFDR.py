import pandas as pd
import numpy as np
from anytree import Node, RenderTree
from anytree.exporter import DotExporter
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.conversion import localconverter

structSSI = importr("structSSI")
hFDR = ro.r["hFDR.adjust"]


# when no children, it is a leaf node and you just return the node and the p-value
# when you have children, you call build_tree on each child
# build_tree returns a node (the root node of a tree)
def build_tree(node_name, el, pval_dict, auc_dict, n_reads_dict):
    pvalue = pval_dict.get(node_name)
    auc = auc_dict.get(node_name)
    n_reads = n_reads_dict.get(node_name)
    node = Node(node_name, pvalue=pvalue, metric=auc, n_reads=n_reads, removed_children=[])
    children = el.loc[el.parent == node_name]
    child_nodes = []
    for index, child in children.iterrows():
        child_node = build_tree(child["child"], el, pval_dict, auc_dict, n_reads_dict)
        child_node.parent = node
    return node

# the function prunes nodes that were not evaluated due to low expression
# in addition, it removes single child nodes
def update_tree(node):
    n_children = len(node.children)
    if n_children == 0:
        return node
    if n_children == 1:
        #print("removing {}".format(node.children[0].name))
        child = update_tree(node.children[0])
        node.children = child.children
        # to do, add removed children as node attribute
        node.removed_children.extend(child.removed_children)
        node.removed_children.append(child.name)
        return node
    if n_children > 1:
        updated_children = []
        for child in node.children:
            updated_children.append(update_tree(child))
        node.children = updated_children
        children_with_pvalues = []
        children_without_pvalues = []
        for child in node.children:
            if child.pvalue is None:
                children_without_pvalues.append(child)
            else:
                children_with_pvalues.append(child)
        # if >= 2 child has non-na p-value - drop all children with na p-value from the tree
        if len(children_with_pvalues) > 1:
            node.children = children_with_pvalues
            return node
        # if == 1 child has non-na p-value - pick one child with na p-value and set p-value=1
        if len(children_with_pvalues) == 1:
            added_child = children_without_pvalues[0]
            added_child.pvalue = 1
            children_with_pvalues.append(added_child)
            node.children = children_with_pvalues
            return node
        # if 0 children have non-na p-value - drop all children from the tree
        if len(children_with_pvalues) == 0:
            node.children = []
            return node

# this returns a dictionary of names to p-values from the tree
def get_pvalues_from_tree(node):
    node_dict = {node.name: node.pvalue}
    for child in node.children:
        node_dict.update(get_pvalues_from_tree(child))
        #print(node_dict)
    return node_dict

# this returns a dictionary of names to p-values from the tree
def get_auc_dict_from_tree(node):
    node_dict = {node.name: node.metric}
    for child in node.children:
        node_dict.update(get_auc_dict_from_tree(child))
    return node_dict

def get_nreads_dict_from_tree(node):
    node_dict = {node.name: node.n_reads}
    for child in node.children:
        node_dict.update(get_nreads_dict_from_tree(child))
    return node_dict

# this returns a dataframe representing an edgelist of the tree
def get_edgelist_from_tree(node):
    el = pd.DataFrame(columns=["parent", "child"])
    for child in node.children:
        el = el.append({"parent": node.name, "child": child.name}, ignore_index=True)
        el = el.append(get_edgelist_from_tree(child), ignore_index=True)
    return el

def get_removed_children_from_tree(node):
    removed_child_df = pd.DataFrame(columns=["parent", "child"])
    for removed_child in node.removed_children:
        removed_child_df = removed_child_df.append({"parent": node.name, "child": removed_child}, ignore_index=True)
    for child in node.children:
        removed_child_df = removed_child_df.append(get_removed_children_from_tree(child), ignore_index=True)
    return removed_child_df

def run_hFDR_adjust(name, df, el):
    # build dictionary of p-values
    df.index = df.index.astype("str")
    pval_dict = df["p.value"].to_dict()
    auc_dict = df["summary.{}".format(metric)].to_dict()
    n_reads_dict = df["n.reads"].to_dict()
    # build tree
    tree = build_tree(str(name), el, pval_dict, auc_dict, n_reads_dict)
    tree = update_tree(tree)
    pvalues = get_pvalues_from_tree(tree)
    aucs = get_auc_dict_from_tree(tree)
    n_reads = get_nreads_dict_from_tree(tree)
    removed_children = get_removed_children_from_tree(tree)
    print(pvalues)
    # construct R object to pass to hFDR.adjust
    vec = ro.FloatVector(pvalues.values())
    vec.names = ro.StrVector(pvalues.keys())

    # get updated edgelist
    el = get_edgelist_from_tree(tree)
    l = el["parent"].tolist()
    l.extend(el["child"].tolist())
    #print(tree)
    #print(el)
    #print(l)
    if el.empty: # edge case where m is empty, which means this is a single node (not a tree)
        print("edge case: l is empty")
        if tree.pvalue < .001:
            adjsig = "***"
        elif tree.pvalue < .01:
            adjsig = "**"
        elif tree.pvalue < .05:
            adjsig = "*"
        else:
            adjsig = "-"
        pd_from_r_df = pd.DataFrame(index=[tree.name], data={"unadjp": tree.pvalue, "adjp": tree.pvalue, "adj.significance": adjsig})
        #print(pd_from_r_df)
    else:
        v = ro.StrVector(l)
        m = ro.r['matrix'](v, ncol = 2)
        # this is where we call the R function hFDR
        out = hFDR(vec, m, alpha = 0.05)
        adj_pvals = out.slots['p.vals']
        with localconverter(ro.default_converter + pandas2ri.converter):
            pd_from_r_df = ro.conversion.rpy2py(adj_pvals)
    #print(pd.Series(data=aucs))
    #print("printing expected output")
    #print(pd_from_r_df)
    out = pd_from_r_df.merge(pd.DataFrame(data={"metric": aucs, "n.reads": n_reads}), left_index=True, right_index=True)
    # now add back in the removed children
    removed_children_list = []
    removed_children_list.append(out)
    for index, row in removed_children.iterrows():
        #print(row)
        parent_row = out.loc[row["parent"]]
        #print(parent_row)
        new_df = pd.DataFrame(index=[str(row["child"])], columns=["unadjp", "adjp", "adj.significance", "metric", "n.reads"], data=[[parent_row["unadjp"], parent_row["adjp"], parent_row["adj.significance"], parent_row["metric"], parent_row["n.reads"]]])
        removed_children_list.append(new_df)
    return pd.concat(removed_children_list)


el = pd.read_csv(snakemake.input["edgelist"], sep="\t", dtype={"parent": "str", "child": "str"})

taxID_map_df = pd.read_csv(snakemake.input["tax_id_map"], sep="\t", dtype={"tax_id": "str"})
# print(taxID_map_df)
class_df = pd.read_csv(snakemake.input["class_markers"], sep="\t", index_col=0)
order_df = pd.read_csv(snakemake.input["order_markers"], sep="\t", index_col=0)
family_df = pd.read_csv(snakemake.input["family_markers"], sep="\t", index_col=0)
genus_df = pd.read_csv(snakemake.input["genus_markers"], sep="\t", index_col=0)
species_df = pd.read_csv(snakemake.input["species_markers"], sep="\t", index_col=0)

df = pd.concat([class_df, order_df, family_df, genus_df, species_df])
#nreads_df = df["n.reads"]
metric = snakemake.params["metric"]

# when there are no microbial reads found above the background
if df.empty:
    df = pd.DataFrame(columns=["name", "taxa_level", "unadjp", "adjp", "summary.{}".format(metric), "n.reads", "tax_id"])
    df.to_csv(snakemake.output[0], sep="\t", index=False)
    quit()

output = []
# loop through the class p-values and run hFDR if they have p-value < .05
for index, class_row in class_df.iterrows():
    print(class_row.name)
    if class_row["p.value"] < .05:
        #print(class_row.name)
        #print(df)
        #print(el)
        out = run_hFDR_adjust(class_row.name, df, el)
        out = out.rename(columns={"metric": "summary.{}".format(metric)})
        #print(out)
        output.append(out)
    else:
        new_df = pd.DataFrame(index=[str(class_row.name)], columns=["unadjp", "adjp", "adj.significance", "summary.{}".format(metric), "n.reads"], data=[[class_row["p.value"], class_row["p.value"], "_", class_row["summary.{}".format(metric)], class_row["n.reads"]]])
        print(new_df)
        output.append(new_df)

n_hypotheses = class_df.shape[0]
df = pd.concat(output)
print(df)
df.adjp = df.adjp.apply(lambda x: np.nan if np.isnan(x) else min(1, x * n_hypotheses))
df = df.merge(taxID_map_df, left_index=True, right_on="tax_id")
# if metric == "AUC":
#     df["summary.AUC"] = df["summary.AUC"].fillna(.5)
# elif metric == "logFC":
#     df["summary.LogFC"] = df["summary.logFC"].fillna(0)
# elif metric == "odds.ratio":
#     df["summary.odds.ratio"] = df["summary.odds.ratio"].fillna(1)
# remove "artifacts"
df = df.loc[~df["summary.{}".format(metric)].isna()]
df["n.reads"] = df["n.reads"].astype("int64")
df[["name", "taxa_level", "unadjp", "adjp", "summary.{}".format(metric), "n.reads", "tax_id"]].sort_values(by=["adjp", "unadjp", "n.reads"]).to_csv(snakemake.output[0], sep="\t", index=False, na_rep="NA")

# DotExporter(tree).to_picture("output/plots/Actinobacteria_tree.png")
