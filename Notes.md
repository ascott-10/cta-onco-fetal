# Everyday github

git add .
git commit -m "update"
git push

# 4/1
Get marker genes
Filter for epithlieal and mesenchymal cells

# Generation of pseudo-bulk profiles
https://decoupler.readthedocs.io/en/v1.9.2/notebooks/pseudobulk.html

The current best practice to correct for this is using a pseudo-bulk approach (Squair J.W., et al 2021), which involves the following steps:

Subsetting the cell type(s) of interest to perform DEA.

Extracting their raw integer counts.

Summing their counts per gene into a single profile if they pass quality control.

Performing DEA if at least two biological replicates per condition are available (more replicates are recommended).