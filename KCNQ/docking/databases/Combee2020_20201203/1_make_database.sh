#!/bin/sh

# This makes a database for a small set of compounds given as smiles


rm -f fort.*
rm -f name.txt
rm -f outputhex.log
rm -f outputwat.log
rm -f ring_count
rm -f temp.*



SUBSTANCES_FNAME=substances.smi

rm -f substances.mol2
cat raw_substances/*mol2 >> substances.mol2
obabel -imol2 -osmi substances.mol2 -O substances.smi

time $DOCKBASE/ligand/generate/build_smiles_ligand.sh $SUBSTANCES_FNAME
echo "$(pwd)/$(basename $SUBSTANCES_FNAME .smi).db2.gz" > "database.sdi"

# used for enrichment
awk '{print substr($2, $2<length($2)-16+1?$2:length($2)-16+1, 16)}' substances.smi > substances.name
