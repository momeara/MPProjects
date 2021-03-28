#!/bin/sh

# This makes a database for a small set of compounds given as smiles


rm -f fort.*
rm -f name.txt
rm -f outputhex.log
rm -f outputwat.log
rm -f ring_count
rm -f temp.*

SUBSTANCES_FNAME=substances.smi
SUBSTANCES=$(cut -d' ' -f2 substances.smi)


rm -f substances.mol2
for substance in $SUBSTANCES; do rm -rf $substance; done
rm -rf *-swap.db2



cat raw_substances/*mol2 >> substances.mol2
obabel -imol2 -osmi substances.mol2 -O substances.smi

time $DOCKBASE/ligand/generate/build_smiles_ligand.sh $SUBSTANCES_FNAME
rm substances.db2.gz


# re run omega conformation with additional sampling
unset OMEGA_MAX_CONFS
unset OMEGA_ENERGY_WINDOW

for substance in $SUBSTANCES
do
  echo "use omega to build extensive conformations for substance: $substance"
  pushd $substance/0
  rm output.*.*
  rm ring_count
  time python ../../../../scripts/omega_db2.e.py \
    --fixrms 0.01 \
    --rms 0.001 \
    --ewindow 100 \
    --maxconfs 10000 \
    --zero_ring_maxconfs 10000 \
    output.mol2
  popd
done
export OMEGA_MAX_CONFS
export OMEGA_ENERGY_WINDOW


for substance in $SUBSTANCES
do
echo "build conformations as rigid elements into .db2 files for substance $substance"
rm -rf ${substance}_rigid
rm -rf ${substance}_rigid.db2.gz
bash ../../scripts/build_smiles_ligand_rigid.sh \
  ${substance}/0 \
  ${substance}_rigid
done


for substance in $SUBSTANCES
do
  echo "For substance $substance, $(zgrep '^E$' ${substance}_rigid.db2.gz | wc -l) conformations generated"
done


rm -rf database.sdi
for substance in $SUBSTANCES
do
  find $PWD/${substance}_rigid -name "*db2.gz" >> database.sdi
done

# used for enrichment
awk '{print substr($2, $2<length($2)-16+1?$2:length($2)-16+1, 16)}' substances.smi > substances.name
