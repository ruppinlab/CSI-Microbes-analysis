import pandas as pd
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
def build_tree(node_name, el, pval_dict, auc_dict):
    pvalue = pval_dict.get(node_name)
    auc = auc_dict.get(node_name)
    node = Node(node_name, pvalue=pvalue, auc=auc, removed_children=[])
    children = el.loc[el.parent == node_name]
    child_nodes = []
    for index, child in children.iterrows():
        child_node = build_tree(child["child"], el, pval_dict, auc_dict)
        child_node.parent = node
    return node

# the function prunes nodes that were not evaluated due to low expression
# in addition, it removes single child nodes
def update_tree(node):
    n_children = len(node.children)
    if n_children == 0:
        return node
    if n_children == 1:
        print("removing {}".format(node.children[0].name))
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
    node_dict = {node.name: node.auc}
    for child in node.children:
        node_dict.update(get_auc_dict_from_tree(child))
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
    auc_dict = df["summary.AUC"].to_dict()
    # build tree
    tree = build_tree(str(name), el, pval_dict, auc_dict)
    tree = update_tree(tree)
    pvalues = get_pvalues_from_tree(tree)
    aucs = get_auc_dict_from_tree(tree)
    removed_children = get_removed_children_from_tree(tree)

    # construct R object to pass to hFDR.adjust
    vec = ro.FloatVector(pvalues.values())
    vec.names = ro.StrVector(pvalues.keys())

    # get updated edgelist
    el = get_edgelist_from_tree(tree)
    l = el["parent"].tolist()
    l.extend(el["child"].tolist())
    v = ro.StrVector(l)
    m = ro.r['matrix'](v, ncol = 2)

    out = hFDR(vec, m, alpha = 0.05)
    adj_pvals = out.slots['p.vals']
    with localconverter(ro.default_converter + pandas2ri.converter):
        pd_from_r_df = ro.conversion.rpy2py(adj_pvals)
    #print(pd.Series(data=aucs))
    out = pd_from_r_df.merge(pd.Series(data=aucs).to_frame(name="AUC"), left_index=True, right_index=True)
    # now add back in the removed children
    removed_children_list = []
    removed_children_list.append(out)
    for index, row in removed_children.iterrows():
        print(row)
        parent_row = out.loc[row["parent"]]
        new_df = pd.DataFrame(index=[str(row["child"])], columns=["unadjp", "adjp", "adj.significance", "AUC"], data=[[parent_row["unadjp"], parent_row["adjp"], parent_row["adj.significance"], parent_row["AUC"]]])
        removed_children_list.append(new_df)
    return pd.concat(removed_children_list)


el = pd.read_csv(snakemake.input["edgelist"], sep="\t", dtype={"parent": "str", "child": "str"})

taxID_map_df = pd.read_csv(snakemake.input["tax_id_map"], sep="\t", dtype={"tax_id": "str"})
#print(taxID_map_df)
class_df = pd.read_csv(snakemake.input["class_markers"], sep="\t")
order_df = pd.read_csv(snakemake.input["order_markers"], sep="\t")
family_df = pd.read_csv(snakemake.input["family_markers"], sep="\t")
genus_df = pd.read_csv(snakemake.input["genus_markers"], sep="\t")
species_df = pd.read_csv(snakemake.input["species_markers"], sep="\t")

df = pd.concat([class_df, order_df, family_df, genus_df, species_df])

output = []
# loop through the class p-values and run hFDR if they have p-value < .05
for index, class_row in class_df.iterrows():
    if class_row["p.value"] < .05:
        #print(class_row.name)
        out = run_hFDR_adjust(class_row.name, df, el)
        print(out)
        output.append(out)
    else:
        new_df = pd.DataFrame(index=[str(class_row.name)], columns=["unadjp", "adjp", "adj.significance", "AUC"], data=[[class_row["p.value"], class_row["p.value"], "_", class_row["summary.AUC"]]])
        output.append(new_df)

n_hypotheses = class_df.shape[0]
df = pd.concat(output)
df.adjp = df.adjp.apply(lambda x: min(1, x * n_hypotheses))
print(df)
df = df.merge(taxID_map_df, left_index=True, right_on="tax_id")
df[["name", "taxa_level", "unadjp", "adjp", "AUC", "tax_id"]].sort_values(by=["adjp", "unadjp"]).to_csv(snakemake.output[0], sep="\t", index=False)

# DotExporter(tree).to_picture("output/plots/Actinobacteria_tree.png")
