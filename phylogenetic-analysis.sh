#!/bin/bash
META=2021-04-19-08-36_B117-globDK.csv
SEQS=2021-04-19-08-36_B117-globDK-aligned.fasta

# Build initial tree.
mkdir -p init-tree

singularity --silent exec \
 -B /srv/rbd:/srv/rbd \
 -B ${PWD}:${PWD} \
 /srv/rbd/thecontainer/covid19_latest.sif \
 bash -c "source activate nextstrain; iqtree -ninit 2 -me 0.05 -nt 40 -s $SEQS -m GTR+F+R3 --prefix init-tree/out"
 
# filtering and pruning --------------------------------------------------------
# make the dates.
awk -F',' 'BEGIN {print "strain,date"} NR > 1 && $8 == "B.1.1.7" {print $1","$3}' $META > B117_dates.csv

mkdir -p filtering

# Root-to-tip regression pruning.
Rscript --vanilla root-to-tip.R \
  --tree init-tree/out.treefile \
  --dates B117_dates.csv \
  --iqr 3 --name B117_rtt \
  --outdir ./filtering \
  --threads 10
  
# long branch pruning.
Rscript --vanilla prune-long-branches.R \
  --tree filtering/B117_rtt.filtered.nwk \
  --outdir ./filtering --name B117_brlen \
  --percentile .99
  
# Combine the removed into one file for downstream analysis.
cat filtering/*filtered > to-remove.txt

# subset and find optimal model.--------------------------------------
# find sequences to retain.
comm -23 <(grep ">" $SEQS | sed 's/>//' | sort) <(sort to-remove.txt) > to-keep.txt

seqtk subseq $SEQS to-keep.txt > masked-filtered.fasta
  
# quick tree.
mkdir -p quick-tree

singularity --silent exec \
 -B /srv/rbd:/srv/rbd \
 -B ${PWD}:${PWD} \
 /srv/rbd/thecontainer/covid19_latest.sif \
 bash -c "source activate nextstrain; iqtree -ninit 2 -me 0.05 -nt 30 -s masked-filtered.fasta -m GTR+F+R3 --prefix quick-tree/out"

#Find the optimal model (already done previously!).
mkdir model-finder

singularity --silent exec \
 -B /srv/rbd:/srv/rbd \
 -B ${PWD}:${PWD} \
 /srv/rbd/thecontainer/covid19_latest.sif \
 bash -c "source activate nextstrain; iqtree -ninit 2 -me 0.05 -nt 20 -s masked-filtered.fasta -m MFP --prefix model-finder/out"
 
# Best model was GTR+F+R3

# BOOTSTRAP -------------------------------------------------------------
mkdir -p bootstrap_1000iter

singularity --silent exec \
 -B /srv/rbd:/srv/rbd \
 -B ${PWD}:${PWD} \
 /srv/rbd/thecontainer/covid19_latest.sif \
 bash -c "source activate nextstrain; iqtree -ninit 2 -me 0.05 -nt 40 -s masked-filtered.fasta -m GTR+F+R3 --prefix bootstrap_1000iter/out -B 1000 -alrt 1000"

# We need a unique identifier for each node to identify it later.
Rscript -e "args <- commandArgs(trailingOnly=T); library(ape); tree <- read.tree(args[1]); tree[['node.label']] <- paste0(tree[['node.label']],'-',seq_along(tree[['node.label']])); write.tree(tree,'bootstrap_1000iter/out.renamed.contree')" bootstrap_1000iter/out.contree 
 
# time-adjusted tree.----------------------------------------------------------

# Make a date-file.
awk -F',' 'NR > 1 && $8 == "B.1.1.7" {a[$1]=$3} END {print length(a); for (i in a) print i"\t"a[i]}' $META > date-file_B117.tsv

# Make the outgroup file.
awk -F',' 'NR > 1 && $1 == "MN908947.3" {a[$1]=$3} END {print length(a); for (i in a) print i"\t"a[i]}' $META > date-file_outgroup.tsv

# Run lsd2
mkdir -p timetree

lsd2 \
 -i bootstrap_1000iter/out.renamed.contree \
 -d date-file_B117.tsv \
 -g date-file_outgroup.tsv \
 -G \
 -s 29903 \
 -e 4 \
 -w <(echo 0.00056) \
 -l 0.00001 \
 -o timetree/out

# Convert nexus tree to nwk.
Rscript -e "args <- commandArgs(trailingOnly=T); library(ape); tree <- read.nexus(args[1]); write.tree(tree,'timetree/out.date.nwk')" timetree/out.date.nexus

# Pastml analysis.-------------------------------------------------------------
TIMETREE=timetree/out.date.nwk

conda activate pastmlenv_v2

### Group by country. ###
byCountry=pastml/by_country
mkdir -p $byCountry

pastml \
  --tree $TIMETREE \
  --data $META \
  --columns country_region \
  --data_sep , \
  --html_compressed $byCountry/B117_map.html \
  --pajek $byCountry/B117_map.net \
  --threads 10  \
  --resolve_polytomies \
  --tip_size_threshold 10000 \
  --out_data $byCountry/assignment.tsv \
  --work_dir $byCountry \
  --verbose

### Group by DK-vs-nonDK. ###
byNonDK=pastml/by_nonDK
mkdir -p $byNonDK

pastml \
  --tree $TIMETREE \
  --data $META \
  --columns nonDK_region \
  --data_sep , \
  --html_compressed $byNonDK/B117_map.html \
  --pajek $byNonDK/B117_map.net \
  --threads 10  \
  --resolve_polytomies \
  --tip_size_threshold 10000 \
  --out_data $byNonDK/assignment.tsv \
  --work_dir $byNonDK \
  --verbose

# Quantify number of state-changes.
python3 /srv/rbd/tym/software/misc_code/calculate_changes.py \
  --tree $byNonDK/named.tree_out.date.nwk \
  --acr $byNonDK/assignment.tsv \
  --columns nonDK_region \
  --out_log $byNonDK/state_changes.tsv
