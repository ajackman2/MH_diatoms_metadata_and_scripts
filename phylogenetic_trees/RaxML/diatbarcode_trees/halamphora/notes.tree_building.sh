#notes for sequence prep and tree building for thalassiophysales backbone tree
#using eukref procedure, modified for vsearch
#author: Evan Morien
#last modified: Dec 10 2024
#working directory: fucus:/home/evan.morien/projects/parfreylab/diatbarcode_trees/halamphora

#OPTIONAL: pull sequences from PR2 & SILVA/GenBank if there aren't enough backbone sequences to make a decent tree #we are okay for this clade
#seqkit grep -r -n -p '.*Pseudoalteromonas.*' ~/projects/taxonomyDBs/SILVAv138.2/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta> thalassiophysales.tree_input.fasta 
#sed -i '/^[^>]/s/U/T/g' thalassiophysales.tree_input.fasta

#### 18S ####
#remove sequences: HM805021 (somehow made through the filters in even though it's ITS2 not 18S), JX896689 (long branching)	

#removing all illegal characters in fasta headrs per raxml's documentation (";", ":", ",", ")", "(", "]", "[", "'"), replace with - character
sed -i '/^>/ s/;/-/g' *.fasta
sed -i '/^>/ s/[\(\)]/-/g' *.fasta
sed -i '/^[actg|ACTG]/ s/-//g' *.fasta #fussing around because vsearch doesn't like - characters in the sequences

#combine backbone and target clade sequences, there are not too many and the tree should be OK
cat thalassiophysales_backbone.plus_additional_seqs.18S.fasta thalassiophysales_backbone.18S.fasta halamphora.18S.fasta amphora.18S.fasta > thalassiophysales.tree_input.fasta

#sorting, clustering
vsearch -sortbylength thalassiophysales.tree_input.fasta --output thalassiophysales.tree_input.sorted.fasta -minseqlength 200 -notrunclabels
vsearch -cluster_smallmem thalassiophysales.tree_input.sorted.fasta -id 0.99 -centroids thalassiophysales.tree_input.clustered.fasta -uc thalassiophysales.tree_input.clusters -notrunclabels

#align, trim, clean up
mafft --thread 14 --reorder --auto thalassiophysales.tree_input.clustered.fasta > thalassiophysales.tree_input.aligned.fasta #align with mafft
sed -i '/^[-nywkmsractg|NACTG]/ s/n/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "n" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/y/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "y" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/w/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "y" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/k/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "y" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/m/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "y" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/s/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "y" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/r/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "r" characters, but apparently will accept "N"??
trimal -in thalassiophysales.tree_input.aligned.fasta -out thalassiophysales.tree_input.trimal.fasta -gt 0.3 -st 0.001 #trim with trimal
sed -i '/^>/ s/_.*$//' thalassiophysales.tree_input.trimal.fasta #clean up alignment fasta headers so they just have the accessions, not sequence lengths (we need them to match the IDs in the tree output for the next steps)

raxmlHPC-PTHREADS -T 14 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n thalassiophysales.tree_input -s thalassiophysales.tree_input.trimal.fasta -w /parfreylab/morien/diatbarcode_trees/halamphora/18S

#blast ASVs to reference sequences, filter ASVs to 96% similarity or better
mkdir blast_filtering
sed '/^>/ s/_.*$//' thalassiophysales.tree_input.clustered.fasta > blast_filtering/ref.fasta
#make blastdb from reference seqs
makeblastdb -in blast_filtering/ref.fasta -out blast_filtering/18S_thalassiophysales.blastDB -parse_seqids -dbtype nucl -title 18S_thalassiophysales.blastDB
#blast ASVs to reference sequences
blastn -task megablast -num_threads 22 -evalue 1e-5 -max_target_seqs 10 -perc_identity 98 -qcov_hsp_perc 50 -db blast_filtering/18S_thalassiophysales.blastDB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_18S_asv.fasta  -out blast_filtering/ASV_sequences.blast.out

#filter original ASV set by the set of ASVs that get high identity hits to the clade of interest
cut -f1 blast_filtering/ASV_sequences.blast.out | sort | uniq > blast_filtering/to_keep
seqkit grep -r -n -f blast_filtering/to_keep /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_18S_asv.fasta -o blast_filtering/subset.fasta

#remove following sequqences due to long branching
#asv606_Bacillariophyceae
#asv326_Bacillariophyceae
#asv489_Bacillariophyceae


#align short sequences to template alignment
mothur > align.seqs(candidate=blast_filtering/subset.fasta, template=thalassiophysales.tree_input.trimal.fasta)

#place short sequences
~/programs/epa-ng/bin/epa-ng --ref-msa thalassiophysales.tree_input.trimal.fasta --tree /home/evan.morien/projects/parfreylab/diatbarcode_trees/halamphora/18S/RAxML_bestTree.thalassiophysales.tree_input --query blast_filtering/subset.align --model GTR+G --redo
#convert jplace output to newick format
gappa examine graft --allow-file-overwriting  --jplace-path epa_result.jplace


#### rbcL ####
#removing all illegal characters in fasta headrs per raxml's documentation (";", ":", ",", ")", "(", "]", "[", "'"), replace with - character
sed -i '/^>/ s/;/-/g' *.fasta
sed -i '/^>/ s/[\(\)]/-/g' *.fasta
sed -i '/^[actg|ACTG]/ s/-//g' *.fasta #fussing around because vsearch doesn't like - characters in the sequences

#combine backbone and target clade sequences, there are not too many and the tree should be OK
cat thalassiophysales_backbone.rbcL.fasta halamphora.rbcL.fasta amphora.rbcL.fasta > thalassiophysales.tree_input.fasta

#sorting, clustering
vsearch -sortbylength thalassiophysales.tree_input.fasta --output thalassiophysales.tree_input.sorted.fasta -minseqlength 200 -notrunclabels
vsearch -cluster_smallmem thalassiophysales.tree_input.sorted.fasta -id 0.99 -centroids thalassiophysales.tree_input.clustered.fasta -uc thalassiophysales.tree_input.clusters -notrunclabels

#align, trim, clean up
mafft --thread 14 --reorder --auto thalassiophysales.tree_input.clustered.fasta > thalassiophysales.tree_input.aligned.fasta #align with mafft
sed -i '/^[-nywkmsractg|NACTG]/ s/n/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "n" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/y/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "y" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/w/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "y" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/k/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "y" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/m/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "y" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/s/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "y" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/r/N/g' thalassiophysales.tree_input.aligned.fasta #fussing around because trimal doesn't like "r" characters, but apparently will accept "N"??
trimal -in thalassiophysales.tree_input.aligned.fasta -out thalassiophysales.tree_input.trimal.fasta -gt 0.3 -st 0.001 #trim with trimal
sed -i '/^>/ s/_.*$//' thalassiophysales.tree_input.trimal.fasta #clean up alignment fasta headers so they just have the accessions, not sequence lengths (we need them to match the IDs in the tree output for the next steps)

raxmlHPC-PTHREADS -T 14 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n thalassiophysales.tree_input -s thalassiophysales.tree_input.trimal.fasta -w /parfreylab/morien/diatbarcode_trees/halamphora/rbcL

#blast ASVs to reference sequences, filter ASVs to 96% similarity or better
mkdir blast_filtering
sed '/^>/ s/_.*$//' thalassiophysales.tree_input.clustered.fasta > blast_filtering/ref.fasta
#make blastdb from reference seqs
makeblastdb -in blast_filtering/ref.fasta -out blast_filtering/rbcL_thalassiophysales.blastDB -parse_seqids -dbtype nucl -title rbcL_thalassiophysales.blastDB
#blast ASVs to reference sequences
blastn -task megablast -num_threads 22 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db blast_filtering/rbcL_thalassiophysales.blastDB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_rbcl_asv.fasta  -out blast_filtering/ASV_sequences.blast.out
#filter original ASV set by the set of ASVs that get high identity hits to the clade of interest
cut -f1 blast_filtering/ASV_sequences.blast.out | sort | uniq > blast_filtering/to_keep
seqkit grep -r -n -f blast_filtering/to_keep /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_rbcl_asv.fasta -o blast_filtering/subset.fasta

#align short sequences to template alignment
mothur > align.seqs(candidate=blast_filtering/subset.fasta, template=thalassiophysales.tree_input.trimal.fasta)

#place short sequences
~/programs/epa-ng/bin/epa-ng --ref-msa thalassiophysales.tree_input.trimal.fasta --tree /home/evan.morien/projects/parfreylab/diatbarcode_trees/halamphora/rbcL/RAxML_bestTree.thalassiophysales.tree_input --query blast_filtering/subset.align --model GTR+G --redo
#convert jplace output to newick format
gappa examine graft --allow-file-overwriting --jplace-path epa_result.jplace