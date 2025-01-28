#notes for sequence prep and tree building for navicuales backbone tree
#using eukref procedure, modified for vsearch
#author: Evan Morien
#last modified: Dec 10 2024
#working directory: fucus:/home/evan.morien/projects/parfreylab/diatbarcode_trees/navicuales

#OPTIONAL: pull sequences from PR2 & SILVA/GenBank if there aren't enough backbone sequences to make a decent tree #we are okay for this clade
#seqkit grep -r -n -p '.*Pseudoalteromonas.*' ~/projects/taxonomyDBs/SILVAv138.2/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta> navicuales.tree_input.fasta 
#sed -i '/^[^>]/s/U/T/g' navicuales.tree_input.fasta

#### 18S ####
#this set of sequences clusters together. all likely ITS2, rather than 18S as some are labeled
grep HM805024 navicuales.tree_input.fasta
grep HM805025 navicuales.tree_input.fasta
grep HM805026 navicuales.tree_input.fasta
grep HM805027 navicuales.tree_input.fasta
grep HM805028 navicuales.tree_input.fasta
grep KC017468 navicuales.tree_input.fasta
grep KC017471 navicuales.tree_input.fasta
grep KC329503 navicuales.tree_input.fasta
grep KF015200 navicuales.tree_input.fasta
grep KT247427 navicuales.tree_input.fasta
grep KT247431 navicuales.tree_input.fasta
grep KT247438 navicuales.tree_input.fasta
grep KT247443 navicuales.tree_input.fasta
grep KT390088 navicuales.tree_input.fasta
# >HM805024|ITS2_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Naviculales-Sellaphoraceae-Fistulifera-Fistulifera saprophila
# >HM805025|ITS2_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Bacillariales-Bacillariaceae-Nitzschia-Nitzschia fonticola
# >HM805026|ITS2_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Naviculales-Naviculaceae-Navicula-Navicula gregaria
# >HM805027|ITS2_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Bacillariales-Bacillariaceae-Nitzschia-Nitzschia microcephala
# >HM805028|ITS2_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Naviculales-Naviculaceae-Navicula-Navicula perminuta
# >KC017468|18S_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Bacillariales-Bacillariaceae-Pseudo-nitzschia-Pseudo-nitzschia hasleana
# >KC017471|18S_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Bacillariales-Bacillariaceae-Pseudo-nitzschia-Pseudo-nitzschia pungens
# >KC329503|18S_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Bacillariales-Bacillariaceae-Pseudo-nitzschia-Pseudo-nitzschia fraudulenta
# >KF015200|18S_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Bacillariales-Bacillariaceae-Pseudo-nitzschia-Pseudo-nitzschia delicatissima
# >KT247427|18S_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Bacillariales-Bacillariaceae-Pseudo-nitzschia-Pseudo-nitzschia delicatissima
# >KT247431|18S_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Bacillariales-Bacillariaceae-Pseudo-nitzschia-Pseudo-nitzschia pungens
# >KT247438|18S_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Bacillariales-Bacillariaceae-Pseudo-nitzschia-Pseudo-nitzschia calliantha
# >KT247443|18S_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Bacillariales-Bacillariaceae-Pseudo-nitzschia-Pseudo-nitzschia multistriata
# >KT390088|18S_Eukaryota-Chromista-Chromobiota-Bacillariophyta-Bacillariophyceae-Bacillariales-Bacillariaceae-Nitzschia-Nitzschia cf. pusilla

#removing all illegal characters in fasta headrs per raxml's documentation (";", ":", ",", ")", "(", "]", "[", "'"), replace with - character
sed -i '/^>/ s/;/-/g' *.fasta
sed -i '/^>/ s/[\(\)]/-/g' *.fasta
sed -i '/^[actg|ACTG]/ s/-//g' *.fasta #fussing around because vsearch doesn't like - characters in the sequences

#clustering backbone tree at 97% identity
vsearch -sortbylength navicuales_backbone.18S.fasta --output navicuales_backbone.18S.sorted.fasta -minseqlength 200 -notrunclabels
vsearch -cluster_smallmem navicuales_backbone.18S.sorted.fasta -id 0.97 -centroids navicuales_backbone.18S.clustered.fasta -uc navicuales_backbone.tree_input.clusters -notrunclabels

#combine backbone and target clade sequences, there are not too many and the tree should be OK
cat navicuales_backbone.18S.clustered.fasta berkeleya.18S.fasta gyrosigma.18S.fasta plagiotropis.18S.fasta pleurosigma.18S.fasta fistulifera.18S.fasta haslea.18S.fasta navicula.18S.fasta plagiotropis.additional_seqs.18S.fasta proschkinia.18S.fasta > navicuales.tree_input.fasta

#sorting, clustering
vsearch -sortbylength navicuales.tree_input.fasta --output navicuales.tree_input.sorted.fasta -minseqlength 200 -notrunclabels
vsearch -cluster_smallmem navicuales.tree_input.sorted.fasta -id 0.99 -centroids navicuales.tree_input.clustered.fasta -uc navicuales.tree_input.clusters -notrunclabels

#align, trim, clean up
mafft --thread 14 --reorder --auto navicuales.tree_input.clustered.fasta > navicuales.tree_input.aligned.fasta #align with mafft
sed -i '/^[-nywkmsractg|NACTG]/ s/n/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like "n" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/y/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/w/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/k/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/m/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/s/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/r/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
trimal -in navicuales.tree_input.aligned.fasta -out navicuales.tree_input.trimal.fasta -gt 0.3 -st 0.001 #trim with trimal
sed -i '/^>/ s/_.*$//' navicuales.tree_input.trimal.fasta #clean up alignment fasta headers so they just have the accessions, not sequence lengths (we need them to match the IDs in the tree output for the next steps)

raxmlHPC-PTHREADS -T 14 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n navicuales.tree_input -s navicuales.tree_input.trimal.fasta -w /parfreylab/morien/diatbarcode_trees/navicuales/18S

#blast ASVs to reference sequences, filter ASVs to 96% similarity or better
mkdir blast_filtering
sed '/^>/ s/_.*$//' navicuales.tree_input.clustered.fasta > blast_filtering/ref.fasta
#make blastdb from reference seqs
makeblastdb -in blast_filtering/ref.fasta -out blast_filtering/18S_navicuales.blastDB -parse_seqids -dbtype nucl -title 18S_navicuales.blastDB
#blast ASVs to reference sequences
blastn -task megablast -num_threads 22 -evalue 1e-5 -max_target_seqs 10 -perc_identity 98 -qcov_hsp_perc 50 -db blast_filtering/18S_navicuales.blastDB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_18S_asv.fasta  -out blast_filtering/ASV_sequences.blast.out

#filter original ASV set by the set of ASVs that get high identity hits to the clade of interest
cut -f1 blast_filtering/ASV_sequences.blast.out | sort | uniq > blast_filtering/to_keep
seqkit grep -r -n -f blast_filtering/to_keep /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_18S_asv.fasta -o blast_filtering/subset.fasta

#align short sequences to template alignment
mothur > align.seqs(candidate=blast_filtering/subset.fasta, template=navicuales.tree_input.trimal.fasta)


#place short sequences
~/programs/epa-ng/bin/epa-ng --ref-msa navicuales.tree_input.trimal.fasta --tree /home/evan.morien/projects/parfreylab/diatbarcode_trees/navicuales/18S/RAxML_bestTree.navicuales.tree_input --query blast_filtering/subset.align --model GTR+G --redo
#convert jplace output to newick format
gappa examine graft --jplace-path epa_result.jplace




#### rbcL ####
#removing all illegal characters in fasta headrs per raxml's documentation (";", ":", ",", ")", "(", "]", "[", "'"), replace with - character
sed -i '/^>/ s/;/-/g' *.fasta
sed -i '/^>/ s/[\(\)]/-/g' *.fasta
sed -i '/^[actg|ACTG]/ s/-//g' *.fasta #fussing around because vsearch doesn't like - characters in the sequences

#clustering backbone tree at 97% identity
vsearch -sortbylength navicuales_backbone.rbcL.fasta --output navicuales_backbone.rbcL.sorted.fasta -minseqlength 200 -notrunclabels
vsearch -cluster_smallmem navicuales_backbone.rbcL.sorted.fasta -id 0.97 -centroids navicuales_backbone.rbcL.clustered.fasta -uc navicuales_backbone.tree_input.clusters -notrunclabels

#combine backbone and target clade sequences, there are not too many and the tree should be OK
cat navicuales_backbone.rbcL.clustered.fasta berkeleya.rbcL.fasta gyrosigma.rbcL.fasta plagiotropis.rbcL.fasta pleurosigma.rbcL.fasta fistulifera.rbcL.fasta haslea.rbcL.fasta navicula.rbcL.fasta pleurosigma.additional_seqs.rbcL.fasta proschkinia.rbcL.fasta > navicuales.tree_input.fasta

#sorting, clustering
vsearch -sortbylength navicuales.tree_input.fasta --output navicuales.tree_input.sorted.fasta -minseqlength 200 -notrunclabels
vsearch -cluster_smallmem navicuales.tree_input.sorted.fasta -id 0.99 -centroids navicuales.tree_input.clustered.fasta -uc navicuales.tree_input.clusters -notrunclabels

#align, trim, clean up
mafft --thread 14 --reorder --auto navicuales.tree_input.clustered.fasta > navicuales.tree_input.aligned.fasta #align with mafft
sed -i '/^[-nywkmsractg|NACTG]/ s/n/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like "n" characters, but apparently will accept "N"??
sed -i '/^[-nywkmsractg|NACTG]/ s/y/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/w/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/k/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/m/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/s/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
sed -i '/^[-nywkmsractg|NACTG]/ s/r/N/g' navicuales.tree_input.aligned.fasta #fussing around because trimal doesn't like degenerate base characters
trimal -in navicuales.tree_input.aligned.fasta -out navicuales.tree_input.trimal.fasta -gt 0.3 -st 0.001 #trim with trimal
sed -i '/^>/ s/_.*$//' navicuales.tree_input.trimal.fasta #clean up alignment fasta headers so they just have the accessions, not sequence lengths (we need them to match the IDs in the tree output for the next steps)

raxmlHPC-PTHREADS -T 14 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n navicuales.tree_input -s navicuales.tree_input.trimal.fasta -w /parfreylab/morien/diatbarcode_trees/navicuales/rbcL

#blast ASVs to reference sequences, filter ASVs to 96% similarity or better
mkdir blast_filtering
sed '/^>/ s/_.*$//' navicuales.tree_input.clustered.fasta > blast_filtering/ref.fasta
#make blastdb from reference seqs
makeblastdb -in blast_filtering/ref.fasta -out blast_filtering/rbcL_navicuales.blastDB -parse_seqids -dbtype nucl -title rbcL_navicuales.blastDB
#blast ASVs to reference sequences
blastn -task megablast -num_threads 22 -evalue 1e-5 -max_target_seqs 10 -perc_identity 98 -qcov_hsp_perc 50 -db blast_filtering/rbcL_navicuales.blastDB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query /home/evan.morien/projects/parfreylab/diatbarcode_trees/  -out blast_filtering/ASV_sequences.blast.out

#filter original ASV set by the set of ASVs that get high identity hits to the clade of interest
cut -f1 blast_filtering/ASV_sequences.blast.out | sort | uniq > blast_filtering/to_keep
seqkit grep -r -n -f blast_filtering/to_keep /home/evan.morien/projects/parfreylab/diatbarcode_trees/2024-11-06_MH_diatoms_rbcl_asv.fasta -o blast_filtering/subset.fasta

#align short sequences to template alignment
mothur > align.seqs(candidate=blast_filtering/subset.fasta, template=navicuales.tree_input.trimal.fasta)

#place short sequences
~/programs/epa-ng/bin/epa-ng --ref-msa navicuales.tree_input.trimal.fasta --tree /home/evan.morien/projects/parfreylab/diatbarcode_trees/navicuales/rbcL/RAxML_bestTree.navicuales.tree_input --query blast_filtering/subset.align --model GTR+G --redo
#convert jplace output to newick format
gappa examine graft --jplace-path epa_result.jplace