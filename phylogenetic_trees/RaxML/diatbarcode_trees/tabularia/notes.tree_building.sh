#notes for sequence prep and tree building for tabularia backbone tree
#using eukref procedure, modified for vsearch
#author: Evan Morien
#last modified: Dec 10 2024
#working directory: fucus:/home/evan.morien/projects/parfreylab/diatbarcode_trees/tabularia

#OPTIONAL: pull sequences from PR2 & SILVA/GenBank if there aren't enough backbone sequences to make a decent tree #we are okay for this clade
#seqkit grep -r -n -p '.*Pseudoalteromonas.*' ~/projects/taxonomyDBs/SILVAv138.2/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta> tabularia.tree_input.fasta 
#sed -i '/^[^>]/s/U/T/g' tabularia.tree_input.fasta

#### 18S ####
#removing all illegal characters in fasta headrs per raxml's documentation (";", ":", ",", ")", "(", "]", "[", "'"), replace with - character
sed -i '/^>/ s/;/-/g' *.fasta
sed -i '/^>/ s/[\(\)]/-/g' *.fasta
sed -i '/^[actg|ACTG]/ s/-//g' *.fasta #fussing around because vsearch doesn't like - characters in the sequences

#combine backbone and target clade sequences, there are not too many and the tree should be OK
cat licmophorales_backbone.18S.fasta tabularia.18S.fasta > tabularia.tree_input.fasta

#sorting, clustering
vsearch -sortbylength tabularia.tree_input.fasta --output tabularia.tree_input.sorted.fasta -minseqlength 200 -notrunclabels
vsearch -cluster_smallmem tabularia.tree_input.sorted.fasta -id 0.99 -centroids tabularia.tree_input.clustered.fasta -uc tabularia.tree_input.clusters -notrunclabels

#align, trim, clean up
mafft --thread 14 --reorder --auto tabularia.tree_input.clustered.fasta > tabularia.tree_input.aligned.fasta #align with mafft
sed -i '/^[-nywkmsractg|NACTG]/ s/n/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like "n" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/y/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/w/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/k/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/m/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/s/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/r/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
trimal -in tabularia.tree_input.aligned.fasta -out tabularia.tree_input.trimal.fasta -gt 0.3 -st 0.001 #trim with trimal
sed -i '/^>/ s/_.*$//' tabularia.tree_input.trimal.fasta #clean up alignment fasta headers so they just have the accessions, not sequence lengths (we need them to match the IDs in the tree output for the next steps)

raxmlHPC-PTHREADS -T 14 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n tabularia.tree_input -s tabularia.tree_input.trimal.fasta -w /parfreylab/morien/diatbarcode_trees/tabularia/18S

#blast ASVs to reference sequences, filter ASVs to 96% similarity or better
mkdir blast_filtering
sed '/^>/ s/_.*$//' tabularia.tree_input.clustered.fasta > blast_filtering/ref.fasta
#make blastdb from reference seqs
makeblastdb -in blast_filtering/ref.fasta -out blast_filtering/18S_tabularia.blastDB -parse_seqids -dbtype nucl -title 18S_tabularia.blastDB
#blast ASVs to reference sequences
blastn -task megablast -num_threads 22 -evalue 1e-5 -max_target_seqs 10 -perc_identity 97 -qcov_hsp_perc 50 -db blast_filtering/18S_tabularia.blastDB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_18S_asv.fasta  -out blast_filtering/ASV_sequences.blast.out

#filter original ASV set by the set of ASVs that get high identity hits to the clade of interest
cut -f1 blast_filtering/ASV_sequences.blast.out | sort | uniq > blast_filtering/to_keep
seqkit grep -r -n -f blast_filtering/to_keep /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_18S_asv.fasta -o blast_filtering/subset.fasta

#align short sequences to template alignment
mothur > align.seqs(candidate=blast_filtering/subset.fasta, template=tabularia.tree_input.trimal.fasta)

#place short sequences
~/programs/epa-ng/bin/epa-ng --ref-msa tabularia.tree_input.trimal.fasta --tree /home/evan.morien/projects/parfreylab/diatbarcode_trees/tabularia/18S/RAxML_bestTree.tabularia.tree_input --query blast_filtering/subset.align --model GTR+G --redo
#convert jplace output to newick format
gappa examine graft --allow-file-overwriting  --jplace-path epa_result.jplace



#### rbcL ####
#removing all illegal characters in fasta headrs per raxml's documentation (";", ":", ",", ")", "(", "]", "[", "'"), replace with - character
sed -i '/^>/ s/;/-/g' *.fasta
sed -i '/^>/ s/[\(\)]/-/g' *.fasta
sed -i '/^[actg|ACTG]/ s/-//g' *.fasta #fussing around because vsearch doesn't like - characters in the sequences

#combine backbone and target clade sequences, there are not too many and the tree should be OK
cat licmophorales_backbone.rbcL.fasta tabularia.rbcL.fasta > tabularia.tree_input.fasta

#sorting, clustering
vsearch -sortbylength tabularia.tree_input.fasta --output tabularia.tree_input.sorted.fasta -minseqlength 200 -notrunclabels
vsearch -cluster_smallmem tabularia.tree_input.sorted.fasta -id 0.99 -centroids tabularia.tree_input.clustered.fasta -uc tabularia.tree_input.clusters -notrunclabels

#align, trim, clean up
mafft --thread 14 --reorder --auto tabularia.tree_input.clustered.fasta > tabularia.tree_input.aligned.fasta #align with mafft
sed -i '/^[-nywkmsractg|NACTG]/ s/n/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like "n" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/y/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/w/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/k/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/m/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/s/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/r/N/g' tabularia.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
trimal -in tabularia.tree_input.aligned.fasta -out tabularia.tree_input.trimal.fasta -gt 0.3 -st 0.001 #trim with trimal
sed -i '/^>/ s/_.*$//' tabularia.tree_input.trimal.fasta #clean up alignment fasta headers so they just have the accessions, not sequence lengths (we need them to match the IDs in the tree output for the next steps)

raxmlHPC-PTHREADS -T 14 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n tabularia.tree_input -s tabularia.tree_input.trimal.fasta -w /parfreylab/morien/diatbarcode_trees/tabularia/rbcL

#blast ASVs to reference sequences, filter ASVs to 96% similarity or better
mkdir blast_filtering
sed '/^>/ s/_.*$//' tabularia.tree_input.clustered.fasta > blast_filtering/ref.fasta
#make blastdb from reference seqs
makeblastdb -in blast_filtering/ref.fasta -out blast_filtering/rbcL_tabularia.blastDB -parse_seqids -dbtype nucl -title rbcL_tabularia.blastDB
#blast ASVs to reference sequences
blastn -task megablast -num_threads 22 -evalue 1e-5 -max_target_seqs 10 -perc_identity 98 -qcov_hsp_perc 50 -db blast_filtering/rbcL_tabularia.blastDB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_rbcl_asv.fasta  -out blast_filtering/ASV_sequences.blast.out
#filter original ASV set by the set of ASVs that get high identity hits to the clade of interest
cut -f1 blast_filtering/ASV_sequences.blast.out | sort | uniq > blast_filtering/to_keep
seqkit grep -r -n -f blast_filtering/to_keep /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_rbcl_asv.fasta -o blast_filtering/subset.fasta

#align short sequences to template alignment
mothur > align.seqs(candidate=blast_filtering/subset.fasta, template=tabularia.tree_input.trimal.fasta)

#place short sequences
~/programs/epa-ng/bin/epa-ng --ref-msa tabularia.tree_input.trimal.fasta --tree /home/evan.morien/projects/parfreylab/diatbarcode_trees/tabularia/rbcL/RAxML_bestTree.tabularia.tree_input --query blast_filtering/subset.align --model GTR+G --redo
#convert jplace output to newick format
gappa examine graft --allow-file-overwriting --jplace-path epa_result.jplace