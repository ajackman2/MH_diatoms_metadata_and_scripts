#notes for sequence prep and tree building for rhabdonema backbone tree
#using eukref procedure, modified for vsearch
#author: Evan Morien
#last modified: Dec 10 2024
#working directory: fucus:/home/evan.morien/projects/parfreylab/diatbarcode_trees/rhabdonema

#OPTIONAL: pull sequences from PR2 & SILVA/GenBank if there aren't enough backbone sequences to make a decent tree #we are okay for this clade
#seqkit grep -r -n -p '.*Pseudoalteromonas.*' ~/projects/taxonomyDBs/SILVAv138.2/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta> rhabdonema.tree_input.fasta 
#sed -i '/^[^>]/s/U/T/g' rhabdonema.tree_input.fasta

#### 18S ####
#removing all illegal characters in fasta headrs per raxml's documentation (";", ":", ",", ")", "(", "]", "[", "'"), replace with - character
sed -i '/^>/ s/;/-/g' *.fasta
sed -i '/^>/ s/[\(\)]/-/g' *.fasta
sed -i '/^[actg|ACTG]/ s/-//g' *.fasta #fussing around because vsearch doesn't like - characters in the sequences

#combine backbone and target clade sequences, there are not too many and the tree should be OK
cat rhabdonematales_backbone.18S.fasta rhabdonema.18S.fasta > rhabdonema.tree_input.fasta

#sorting, clustering
vsearch -sortbylength rhabdonema.tree_input.fasta --output rhabdonema.tree_input.sorted.fasta -minseqlength 200 -notrunclabels
vsearch -cluster_smallmem rhabdonema.tree_input.sorted.fasta -id 0.99 -centroids rhabdonema.tree_input.clustered.fasta -uc rhabdonema.tree_input.clusters -notrunclabels

#align, trim, clean up
mafft --thread 14 --reorder --auto rhabdonema.tree_input.clustered.fasta > rhabdonema.tree_input.aligned.fasta #align with mafft
sed -i '/^[-nywkmsractg|NACTG]/ s/n/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like "n" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/y/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/w/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/k/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/m/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/s/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/r/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
trimal -in rhabdonema.tree_input.aligned.fasta -out rhabdonema.tree_input.trimal.fasta -gt 0.3 -st 0.001 #trim with trimal
sed -i '/^>/ s/_.*$//' rhabdonema.tree_input.trimal.fasta #clean up alignment fasta headers so they just have the accessions, not sequence lengths (we need them to match the IDs in the tree output for the next steps)

raxmlHPC-PTHREADS -T 14 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n rhabdonema.tree_input -s rhabdonema.tree_input.trimal.fasta -w /parfreylab/morien/diatbarcode_trees/rhabdonema/18S

#blast ASVs to reference sequences, filter ASVs to 96% similarity or better
mkdir blast_filtering
sed '/^>/ s/_.*$//' rhabdonema.tree_input.clustered.fasta > blast_filtering/ref.fasta
#make blastdb from reference seqs
makeblastdb -in blast_filtering/ref.fasta -out blast_filtering/18S_rhabdonema.blastDB -parse_seqids -dbtype nucl -title 18S_rhabdonema.blastDB
#blast ASVs to reference sequences
blastn -task megablast -num_threads 22 -evalue 1e-5 -max_target_seqs 10 -perc_identity 91 -qcov_hsp_perc 50 -db blast_filtering/18S_rhabdonema.blastDB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_18S_asv.fasta  -out blast_filtering/ASV_sequences.blast.out

#filter original ASV set by the set of ASVs that get high identity hits to the clade of interest
cut -f1 blast_filtering/ASV_sequences.blast.out | sort | uniq > blast_filtering/to_keep
seqkit grep -r -n -f blast_filtering/to_keep /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_18S_asv.fasta -o blast_filtering/subset.fasta

#align short sequences to template alignment
mothur > align.seqs(candidate=blast_filtering/subset.fasta, template=rhabdonema.tree_input.trimal.fasta)

#place short sequences
~/programs/epa-ng/bin/epa-ng --ref-msa rhabdonema.tree_input.trimal.fasta --tree /home/evan.morien/projects/parfreylab/diatbarcode_trees/rhabdonema/18S/RAxML_bestTree.rhabdonema.tree_input --query blast_filtering/subset.align --model GTR+G --redo
#convert jplace output to newick format
gappa examine graft --allow-file-overwriting  --jplace-path epa_result.jplace


#### rbcL ####
#removing all illegal characters in fasta headrs per raxml's documentation (";", ":", ",", ")", "(", "]", "[", "'"), replace with - character
sed -i '/^>/ s/;/-/g' *.fasta
sed -i '/^>/ s/[\(\)]/-/g' *.fasta
sed -i '/^[actg|ACTG]/ s/-//g' *.fasta #fussing around because vsearch doesn't like - characters in the sequences

#combine backbone and target clade sequences, there are not too many and the tree should be OK
cat rhabdonematales_backbone.rbcL.fasta rhabdonema.rbcL.fasta > rhabdonema.tree_input.fasta

#sorting, clustering
vsearch -sortbylength rhabdonema.tree_input.fasta --output rhabdonema.tree_input.sorted.fasta -minseqlength 200 -notrunclabels
vsearch -cluster_smallmem rhabdonema.tree_input.sorted.fasta -id 0.99 -centroids rhabdonema.tree_input.clustered.fasta -uc rhabdonema.tree_input.clusters -notrunclabels

#align, trim, clean up
mafft --thread 14 --reorder --auto rhabdonema.tree_input.clustered.fasta > rhabdonema.tree_input.aligned.fasta #align with mafft
sed -i '/^[-nywkmsractg|NACTG]/ s/n/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like "n" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/y/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/w/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/k/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/m/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/s/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/r/N/g' rhabdonema.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
trimal -in rhabdonema.tree_input.aligned.fasta -out rhabdonema.tree_input.trimal.fasta -gt 0.3 -st 0.001 #trim with trimal
sed -i '/^>/ s/_.*$//' rhabdonema.tree_input.trimal.fasta #clean up alignment fasta headers so they just have the accessions, not sequence lengths (we need them to match the IDs in the tree output for the next steps)

raxmlHPC-PTHREADS -T 14 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n rhabdonema.tree_input -s rhabdonema.tree_input.trimal.fasta -w /parfreylab/morien/diatbarcode_trees/rhabdonema/rbcL


#blast ASVs to reference sequences, filter ASVs to 96% similarity or better
mkdir blast_filtering
sed '/^>/ s/_.*$//' rhabdonema.tree_input.clustered.fasta > blast_filtering/ref.fasta
#make blastdb from reference seqs
makeblastdb -in blast_filtering/ref.fasta -out blast_filtering/rbcL_rhabdonema.blastDB -parse_seqids -dbtype nucl -title rbcL_rhabdonema.blastDB
#blast ASVs to reference sequences
blastn -task megablast -num_threads 22 -evalue 1e-5 -max_target_seqs 10 -perc_identity 91 -qcov_hsp_perc 50 -db blast_filtering/rbcL_rhabdonema.blastDB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_rbcl_asv.fasta  -out blast_filtering/ASV_sequences.blast.out
#filter original ASV set by the set of ASVs that get high identity hits to the clade of interest
cut -f1 blast_filtering/ASV_sequences.blast.out | sort | uniq > blast_filtering/to_keep
seqkit grep -r -n -f blast_filtering/to_keep /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_rbcl_asv.fasta -o blast_filtering/subset.fasta

#align short sequences to template alignment
mothur > align.seqs(candidate=blast_filtering/subset.fasta, template=rhabdonema.tree_input.trimal.fasta)

#place short sequences
~/programs/epa-ng/bin/epa-ng --ref-msa rhabdonema.tree_input.trimal.fasta --tree /home/evan.morien/projects/parfreylab/diatbarcode_trees/rhabdonema/rbcL/RAxML_bestTree.rhabdonema.tree_input --query blast_filtering/subset.align --model GTR+G --redo
#convert jplace output to newick format
gappa examine graft --allow-file-overwriting --jplace-path epa_result.jplace