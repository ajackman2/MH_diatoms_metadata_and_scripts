#dietbarcode database trees for parfreylab
#last modified: nov 28, 2024

#initial prep for database done (downloading it from website https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/TOMBYZ)
#manual cleaning (lots of newlines where there shouldn't be, mostly in the attribution fields of the tsv, but that's it. most fixes possible with regex in text editor, a handful (<5 needed manual fix))
#final starting dataset: /home/evan.morien/projects/parfreylab/diatbarcode_trees/diatbarcode.all.modified.txt

#first split into 18S and rbcL sequence sets
grep -i 18S diatbarcode.all.modified.txt > diatbarcode.all.filt.18S.txt
grep -i rbcl diatbarcode.all.modified.txt > diatbarcode.all.filt.rbcL.txt

#first, figure out what fields to retain, and just keep those. we want a unique name for each sequence, plus the sequence itself, plus the full taxonomy
#retain column 2 (accession), 4 (sequence), 5(amplified region), 33 (species), 34 (genus), 35(family), 36(order), 37(class), 38(phylum), 39(subkingdom), 40(kingdom), 41(empire)
awk -F "\t" '{ print $2"\t"$5"\t"$41"\t"$40"\t"$39"\t"$38"\t"$37"\t"$36"\t"$35"\t"$34"\t"$33"\t"$4}' diatbarcode.all.filt.18S.txt > diatbarcode.all.modified.18S.txt
awk -F "\t" '{ print $2"\t"$5"\t"$41"\t"$40"\t"$39"\t"$38"\t"$37"\t"$36"\t"$35"\t"$34"\t"$33"\t"$4}' diatbarcode.all.filt.rbcL.txt > diatbarcode.all.modified.rbcL.txt

# eight trees are needed
#1. Achnanthales backbone tree that should include representation of the genera
#	Cocconeis, Planothidium, Achnanthidium ,Achnanthes
grep -wi "Achnanthales" diatbarcode.all.modified.18S.txt | grep -viw "Cocconeis" | grep -viw "Planothidium" | grep -viw "Achnanthidium" | grep -viw "Achnanthes" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > achnanthales_backbone.18S.fasta
grep -wi "Cocconeis" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > cocconeis.18S.fasta
grep -wi "Planothidium" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > planothidium.18S.fasta
grep -wi "Achnanthidium" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > achnanthidium.18S.fasta
grep -wi "Achnanthes" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > achnanthes.18S.fasta

#2. Halamphora and Amphora with outgroups from other members of Thalassiophysales
grep -wi "Thalassiophysales" diatbarcode.all.modified.18S.txt | grep -viw "Halamphora" | grep -viw "Amphora" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > thalassiophysales_backbone.18S.fasta
grep -wi "Halamphora" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > halamphora.18S.fasta
grep -wi "Amphora" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > amphora.18S.fasta

#3. Attheya with outgroups from other members of Chaetocerotales
grep -wi "Chaetocerotales" diatbarcode.all.modified.18S.txt | grep -viw "Attheya" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > chaetocerotales_backbone.18S.fasta
grep -wi "Attheya" diatbarcode.all.modified.18S.txt| awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > attheya.18S.fasta

#4. Navicuales backbone tree that should include representation of the genera
#	Navicula, Gyrosigma Plagiotropis Pleurosigma Proschkinia Fallacia Berkeleya Haslea Diploneis Fistulifera
grep -wi "Bacillariaceae" diatbarcode.all.modified.18S.txt | grep -viw "Navicula" | grep -viw "Gyrosigma" | grep -viw "Plagiotropis" | grep -viw "Pleurosigma" | grep -viw "Proschkinia" | grep -viw "Proschkinia" | grep -viw "Fallacia" | grep -viw "Berkeleya" | grep -viw "Haslea" | grep -viw "Diploneis" | grep -viw "Fistulifera" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > navicuales_backbone.18S.fasta
grep -wi "Navicula" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > navicula.18S.fasta
grep -wi "Gyrosigma" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > gyrosigma.18S.fasta
grep -wi "Plagiotropis" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > plagiotropis.18S.fasta
grep -wi "Pleurosigma" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > pleurosigma.18S.fasta
grep -wi "Proschkinia" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > proschkinia.18S.fasta
grep -wi "Berkeleya" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > berkeleya.18S.fasta
grep -wi "Haslea" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > haslea.18S.fasta
grep -wi "Fistulifera" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > fistulifera.18S.fasta

#5. Rhabdonema with outgroups from other members of Rhabdonematales
grep -wi "Rhabdonema" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > rhabdonema.18S.fasta
grep -wi "Rhabdonematales" diatbarcode.all.modified.18S.txt | grep -viw "Rhabdonema" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > rhabdonematales_backbone.18S.fasta

#6. Tabularia with outgroups from other members of licmophorales
grep -wi "Licmophorales" diatbarcode.all.modified.18S.txt | grep -viw "Tabularia" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > licmophorales_backbone.18S.fasta
grep -wi "Tabularia" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > tabularia.18S.fasta

#7. Bacillariaceae  backbone tree that should include representation of the genera, Nitzschia and Tryblionella
grep -wi "Bacillariaceae" diatbarcode.all.modified.18S.txt | grep -viw "Nitzschia" | grep -viw "Tryblionella" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > bacillariaceae_backbone.18S.fasta
grep -wi "Nitzschia" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > bacillariaceae.nitzschia.18S.fasta
grep -wi "Tryblionella" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > bacillariaceae.tryblionella.18S.fasta

#8. Stauroneidaceae backbone tree
grep -wi "Stauroneidaceae" diatbarcode.all.modified.18S.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > stauroneidaceae.18S.fasta

#must be completed twice, above for 18S and here for rbcL
# eight trees are needed
#1. Achnanthales backbone tree that should include representation of the genera
#	Cocconeis, Planothidium, Achnanthidium ,Achnanthes
grep -wi "Achnanthales" diatbarcode.all.modified.rbcL.txt | grep -viw "Cocconeis" | grep -viw "Planothidium" | grep -viw "Achnanthidium" | grep -viw "Achnanthes" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > achnanthales_backbone.rbcL.fasta
grep -wi "Cocconeis" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > cocconeis.rbcL.fasta
grep -wi "Planothidium" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > planothidium.rbcL.fasta
grep -wi "Achnanthidium" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > achnanthidium.rbcL.fasta
grep -wi "Achnanthes" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > achnanthes.rbcL.fasta

#2. Halamphora and Amphora with outgroups from other members of Thalassiophysales
grep -wi "Thalassiophysales" diatbarcode.all.modified.rbcL.txt | grep -viw "Halamphora" | grep -viw "Amphora" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > thalassiophysales_backbone.rbcL.fasta
grep -wi "Halamphora" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > halamphora.rbcL.fasta
grep -wi "Amphora" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > amphora.rbcL.fasta

#3. Attheya with outgroups from other members of Chaetocerotales
grep -wi "Chaetocerotales" diatbarcode.all.modified.rbcL.txt | grep -viw "Attheya" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > chaetocerotales_backbone.rbcL.fasta
grep -wi "Attheya" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > attheya.rbcL.fasta

#4. Navicuales backbone tree that should include representation of the genera
#	Navicula, Gyrosigma Plagiotropis Pleurosigma Proschkinia Fallacia Berkeleya Haslea Diploneis Fistulifera
grep -wi "Bacillariaceae" diatbarcode.all.modified.rbcL.txt | grep -viw "Navicula" | grep -viw "Gyrosigma" | grep -viw "Plagiotropis" | grep -viw "Pleurosigma" | grep -viw "Proschkinia" | grep -viw "Proschkinia" | grep -viw "Fallacia" | grep -viw "Berkeleya" | grep -viw "Haslea" | grep -viw "Diploneis" | grep -viw "Fistulifera" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > navicuales_backbone.rbcL.fasta
grep -wi "Navicula" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > navicula.rbcL.fasta
grep -wi "Gyrosigma" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > gyrosigma.rbcL.fasta
grep -wi "Plagiotropis" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > plagiotropis.rbcL.fasta
grep -wi "Pleurosigma" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > pleurosigma.rbcL.fasta
grep -wi "Proschkinia" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > proschkinia.rbcL.fasta
grep -wi "Berkeleya" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > berkeleya.rbcL.fasta
grep -wi "Haslea" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > haslea.rbcL.fasta
grep -wi "Fistulifera" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > fistulifera.rbcL.fasta

#5. Rhabdonema with outgroups from other members of Rhabdonematales
grep -wi "Rhabdonema" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > rhabdonema.rbcL.fasta
grep -wi "Rhabdonematales" diatbarcode.all.modified.rbcL.txt | grep -viw "Rhabdonema" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > rhabdonematales_backbone.rbcL.fasta

#6. Tabularia with outgroups from other members of Licmophorales
grep -wi "Licmophorales" diatbarcode.all.modified.rbcL.txt | grep -viw "Tabularia" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > licmophorales_backbone.rbcL.fasta
grep -wi "Tabularia" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > tabularia.rbcL.fasta

#7. Bacillariaceae  backbone tree that should include representation of the genera, Nitzschia and Tryblionella
grep -wi "Bacillariaceae" diatbarcode.all.modified.rbcL.txt | grep -viw "Nitzschia" | grep -viw "Tryblionella" | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > bacillariaceae_backbone.rbcL.fasta
grep -wi "Nitzschia" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > bacillariaceae.nitzschia.rbcL.fasta
grep -wi "Tryblionella" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > bacillariaceae.tryblionella.rbcL.fasta

#8. Stauroneidaceae backbone tree
grep -wi "Stauroneidaceae" diatbarcode.all.modified.rbcL.txt | awk -F "\t" '{ print ">"$1"|"$2"_"$3";"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11"\n"$12}' > stauroneidaceae.rbcL.fasta

#cleanup (make separate folders for each tree and sort files into proper places)
mkdir achnanthales/18S
mkdir achnanthales/rbcL
mv achnanthales_backbone.rbcL.fasta cocconeis.rbcL.fasta planothidium.rbcL.fasta achnanthidium.rbcL.fasta achnanthes.rbcL.fasta achnanthales/rbcL
mv achnanthales_backbone.18S.fasta cocconeis.18S.fasta planothidium.18S.fasta achnanthidium.18S.fasta achnanthes.18S.fasta achnanthales/18S

mkdir halamphora/18S
mkdir halamphora/rbcL
mv thalassiophysales_backbone.rbcL.fasta halamphora.rbcL.fasta amphora.rbcL.fasta halamphora/rbcL
mv thalassiophysales_backbone.18S.fasta halamphora.18S.fasta amphora.18S.fasta halamphora/18S

mkdir attheya/18S
mkdir attheya/rbcL
mv chaetocerotales_backbone.rbcL.fasta attheya.rbcL.fasta attheya/rbcL
mv chaetocerotales_backbone.18S.fasta attheya.18S.fasta attheya/18S

mkdir navicuales/18S
mkdir navicuales/rbcL
mv navicuales_backbone.rbcL.fasta navicula.rbcL.fasta gyrosigma.rbcL.fasta plagiotropis.rbcL.fasta pleurosigma.rbcL.fasta proschkinia.rbcL.fasta berkeleya.rbcL.fasta haslea.rbcL.fasta fistulifera.rbcL.fasta navicuales/rbcL
mv navicuales_backbone.18S.fasta navicula.18S.fasta gyrosigma.18S.fasta plagiotropis.18S.fasta pleurosigma.18S.fasta proschkinia.18S.fasta berkeleya.18S.fasta haslea.18S.fasta fistulifera.18S.fasta navicuales/18S

mkdir rhabdonema/18S
mkdir rhabdonema/rbcL
mv rhabdonema.rbcL.fasta rhabdonematales_backbone.rbcL.fasta rhabdonema/rbcL
mv rhabdonema.18S.fasta rhabdonematales_backbone.18S.fasta rhabdonema/18S

mkdir tabularia/18S
mkdir tabularia/rbcL
mv licmophorales_backbone.rbcL.fasta tabularia.rbcL.fasta tabularia/rbcL
mv licmophorales_backbone.18S.fasta tabularia.18S.fasta tabularia/18S

mkdir bacillariaceae/18S
mkdir bacillariaceae/rbcL
mv bacillariaceae_backbone.rbcL.fasta bacillariaceae.nitzschia.rbcL.fasta bacillariaceae.tryblionella.rbcL.fasta bacillariaceae/rbcL
mv bacillariaceae_backbone.18S.fasta bacillariaceae.nitzschia.18S.fasta bacillariaceae.tryblionella.18S.fasta bacillariaceae/18S

mkdir stauroneidaceae/18S
mkdir stauroneidaceae/rbcL
mv stauroneidaceae.rbcL.fasta stauroneidaceae/rbcL
mv stauroneidaceae.18S.fasta stauroneidaceae/18S

#### clades for which there are no target sequences, or not enough backbone sequences need sequences from public databases ####

# Clades with no/very few backbone sequences in non-target taxa:
# Thalassiophysales (backbone for Amphora, Halamphora), 18S & rbcL
# Rhabdonematales (18S & rbcL)

# Clades with no/very few sequences in general:
# plagiotropis (18S & rbcL)
# proschkinia (18S)

# used NCBI Genbank to supplement above, search parameters ("18S" or "rbcL" limited to clade of interest, and sequences longer than 200bp, shorter than 2000bp, then hand-pruned in case of errant, non-target sequences appearing in results)