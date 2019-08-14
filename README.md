# Seurat wrapper for kidney
This is a very simple wrapper for a integrative analysis with two single-sample conditions (e.g. disease and normal).

# Usage
Example of single run:
```
NPC_ANCHOR=20
NPC_CLUST=20
NPC_RES=0.5

Rscript ./src/integrated_seurat.R --CASE_INPUT "./data/sc/disease/filtered_feature_bc_matrix/" \
	--CASE_SNAME "disease" \
	--CONTROL_INPUT "./data/sc/control/filtered_feature_bc_matrix/" \
	--CONTROL_SNAME "control" \
	--OUTDIR "./results/Integrated_disease_control_ANCHOR${NPC_ANCHOR}_CLUST${NPC_CLUST}_RES${RES}/" \
	--NPC_ANCHOR ${NPC_ANCHOR} --NPC_CLUSTERING ${NPC_CLUST} --RES ${RES};
```

Example of multiple runnings tuning parameters:
```
for NPC_ANCHOR in 5 10 15 20 25 30 35 40 45 50;do
	for NPC_CLUST in 5 10 15 20 25 30 35 40 45 50;do
		for RES in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0;do
			Rscript ./src/integrated_seurat.R --CASE_INPUT "./data/sc/disease/filtered_feature_bc_matrix/" \
			--CASE_SNAME "disease" \
			--CONTROL_INPUT "./data/sc/control/filtered_feature_bc_matrix/" \
			--CONTROL_SNAME "control" \
			--OUTDIR "./results/Integrated_disease_control_ANCHOR${NPC_ANCHOR}_CLUST${NPC_CLUST}_RES${RES}/" \
			--NPC_ANCHOR ${NPC_ANCHOR} --NPC_CLUSTERING ${NPC_CLUST} --RES ${RES};
		done
	done

done
```
