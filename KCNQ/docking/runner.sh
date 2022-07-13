

for docking_run in docking_runs/KCNQ2_model_Z21_*; do
    pushd ${docking_run}
    bash 1_run.sh
done

for docking_run in docking_runs/KCNQ2_model_Z21_*; do
  pushd ${docking_run}  
  echo ${docking_run}
  rm -rf extract_all.*
  rm -rf pose*
  rm -rf dock_statistics*
  ls results/ | grep -v joblist | sed "s#^#results/#" > dirlist
  python ~/opt/DOCK/analysis/extract_all_blazing_fast.py dirlist extract_all.txt 10
  python ~/opt/DOCK/analysis/getposes_blazing_faster.py '' extract_all.sort.uniq.txt 500 poses.mol2 test.mol2.gz 
#  source ${DOCK_TEMPLATE}/scripts/dock_statistics.sh
#  Rscript ${DOCK_TEMPLATE}/scripts/analysis/gather_pose_features.R --verbose
  popd
done

# I messed up and assumed ZA was Z21 but it is actually Z16
mv structures/KCNQ2_model_Z21_axial_20220309 structures/KCNQ2_model_Z16_axial_20220309
mv structures/KCNQ2_model_Z21_equatorial_20220309 structures/KCNQ2_model_Z16_equatorial_20220309
mv prepared_structures/KCNQ2_model_Z21_axial_20220309 prepared_structures/KCNQ2_model_Z16_axial_20220309
mv prepared_structures/KCNQ2_model_Z21_axial_tarted_20220309 prepared_structures/KCNQ2_model_Z16_axial_tarted_20220309
mv prepared_structures/KCNQ2_model_Z21_equatorial_20220309 prepared_structures/KCNQ2_model_Z16_equatorial_20220309
sed -i "s/Z21/Z16/g" prepared_structures/KCNQ2_model_Z16_axial_20220309/1_prepare_structure.sh
sed -i "s/Z21/Z16/g" prepared_structures/KCNQ2_model_Z16_axial_tarted_20220309/1_prepare_structure.sh
sed -i "s/Z21/Z16/g" prepared_structures/KCNQ2_model_Z16_equatorial_20220309/1_prepare_structure.sh

for docking_run in docking_runs/KCNQ2_model_Z21_*; do
    new_docking_run=$(echo ${docking_run} | sed "s/Z21/Z16/g")
    echo ${new_docking_run}
    mv ${docking_run} ${new_docking_run}
    sed -i "s/Z21/Z16/g" ${new_docking_run}/1_run.sh
done


# now run Z21 like Z16
for Z16_structure in structures/*Z16*; do
    Z21_structure = $(echo ${Z16_structure} | sed "s/Z16/Z21/g")
    echo ${Z21_structure}
    cp -r ${Z16_structure} ${Z21_structure}
done
# edit to for Z21


pushd prepared_structures/KCNQ2_model_Z21_axial_20220309
bash 1_prepare_structure.sh
popd

pushd prepared_structures/KCNQ2_model_Z21_axial_20220309
bash 1_prepare_structure.sh
popd

for Z16_prepared_structure in prepared_structures/*Z16*; do
   Z21_prepared_structure=$(echo ${Z16_prepared_structure} | sed "s/Z16/Z21/g")
   mkdir ${Z21_prepared_structure}
   #cp ${Z16_prepared_structure}/1_prepare_structure.sh ${Z21_prepared_structure}/1_prepare_structure.sh
   #sed -i "s/Z16/Z21/g" ${Z21_prepared_structure}/1_prepare_structure.sh
   #sed -i "s/20220309/20220310/g" ${Z21_prepared_structure}/1_prepare_structure.sh
   bash ${Z21_prepared_structure}/1_prepare_structure.sh
   # make max_atoms 300 in the INDOCK
done

for Z16_docking_run in docking_runs/KCNQ2_model_Z16*; do
   Z21_docking_run=$(echo ${Z16_docking_run} | sed "s/Z16/Z21/g" | sed "s/20220309/20220310/g" )
#   mkdir ${Z21_docking_run}
#   cp ${Z16_docking_run}/1_run.sh ${Z21_docking_run}/1_run.sh
#   sed -i "s/Z16/Z21/g" ${Z21_docking_run}/1_run.sh
#   sed -i "s/20220309/20220310/g" ${Z21_docking_run}/1_run.sh
   pushd ${Z21_docking_run}    
   bash 1_run.sh
   popd
done


for Z21_docking_run in docking_runs/KCNQ2_model_Z21*; do
  pushd ${Z21_docking_run}
    echo "Collecint dock results ..."
    ls results/ | grep -v joblist | sed "s#^#results/#" > dirlist
    python ~/opt/DOCK/analysis/extract_all_blazing_fast.py dirlist extract_all.txt 10
    python ~/opt/DOCK/analysis/getposes_blazing_faster.py '' extract_all.sort.uniq.txt 500 poses.mol2 test.mol2.gz 
    source ${DOCK_TEMPLATE}/scripts/dock_statistics.sh  
    Rscript ${DOCK_TEMPLATE}/scripts/analysis/gather_pose_features.R --verbose
  popd
done

# looks like some of the compounds are big enough that they are
# reaching all the way over to the voltage sensor domain, so try
# cropping the binding site see if that helps descrimination


mkdir prepared_structures/KCNQ2_model_crop_Z21_equatorial_20220405
cp \
  prepared_structures/KCNQ2_model_Z21_equatorial_20220310/1_prepare_structure.sh \
  prepared_structures/KCNQ2_model_crop_Z21_equatorial_20220405

pushd   prepared_structures/KCNQ2_model_crop_Z21_equatorial_20220405
  sed -i "s/KCNQ2_model_Z21_equatorial_20220310/KCNQ2_model_crop_Z21_equatorial_20220405/" 1_prepare_structure.sh
  bash 1_prepare_structure.sh
popd

# wait till it's finished
pushd   prepared_structures/KCNQ2_model_crop_Z21_equatorial_20220405
  bash 2_modify_indock.sh
popd

for Z16_docking_run in docking_runs/KCNQ2_model_Z16_equatorial*; do
   Z21_docking_run=$(echo ${Z16_docking_run} | sed "s/model_Z16/model_crop_Z21/g" | sed "s/20220309/20220405/g" )
   rm -rf ${Z21_docking_run}
   mkdir ${Z21_docking_run}
   cp ${Z16_docking_run}/1_run.sh ${Z21_docking_run}/1_run.sh
   sed -i "s/model_Z16/model_crop_Z21/g" ${Z21_docking_run}/1_run.sh
   sed -i "s/20220309/20220405/g" ${Z21_docking_run}/1_run.sh
   pushd ${Z21_docking_run}    
   bash 1_run.sh
   popd
done


for Z21_docking_run in docking_runs/KCNQ2_model_crop_Z21*; do
    pushd ${Z21_docking_run}
    rm -rf dock_statistics.* extract_all* pose*
    echo "Collecint dock results ..."
    ls results/ | grep -v joblist | sed "s#^#results/#" > dirlist
    python ~/opt/DOCK/analysis/extract_all_blazing_fast.py dirlist extract_all.txt 10
    python ~/opt/DOCK/analysis/getposes_blazing_faster.py '' extract_all.sort.uniq.txt 500 poses.mol2 test.mol2.gz 
    source ${DOCK_TEMPLATE}/scripts/dock_statistics.sh  
    Rscript ${DOCK_TEMPLATE}/scripts/analysis/gather_pose_features.R --verbose
  popd
done


###
# Try building in both ZK21 and Retigabine into the ZK21 model

old_structure_id=KCNQ2_model_Z21_equatorial_20220310
structure_id=
mkdir prepared_structures/${structure_id}
cp \
  prepared_structures/${old_structure_id}/1_prepare_structure.sh \
  prepared_structures/${structure_id}/

pushd   prepared_structures/${structure_id}
  sed -i "s/${old_structure_id}/${structure_id}/" 1_prepare_structure.sh
  bash 1_prepare_structure.sh
popd

pushd   prepared_structures/${structure_id}
  bash 2_modify_indock.sh
popd

for old_docking_run in docking_runs/KCNQ2_model_Z16_equatorial*; do
   docking_run=$(echo ${old_docking_run} | sed "s/model_Z16_equatorial/model_crop_Z21_equatorial_and_ret_sph/g" | sed "s/20220309/20220405/g" )
   rm -rf ${docking_run}
   mkdir ${docking_run}
   cp ${old_docking_run}/1_run.sh ${docking_run}/1_run.sh
   sed -i "s/model_Z16_equatorial/model_crop_Z21_equatorial_and_ret_sph/g" ${docking_run}/1_run.sh
   sed -i "s/20220309/20220405/g" ${docking_run}/1_run.sh
   pushd ${docking_run}    
   bash 1_run.sh
   popd
done


for docking_run in docking_runs/KCNQ2_model_crop_Z21_equatorial_and_ret_sph*; do
    pushd ${docking_run}
    rm -rf dock_statistics.* extract_all* pose*
    echo "Collecint dock results ..."
    ls results/ | grep -v joblist | sed "s#^#results/#" > dirlist
    python ~/opt/DOCK/analysis/extract_all_blazing_fast.py dirlist extract_all.txt 10
    python ~/opt/DOCK/analysis/getposes_blazing_faster.py '' extract_all.sort.uniq.txt 500 poses.mol2 test.mol2.gz 
    source ${DOCK_TEMPLATE}/scripts/dock_statistics.sh  
    Rscript ${DOCK_TEMPLATE}/scripts/analysis/gather_pose_features.R --verbose
  popd
done


###
# Try docking all the project compounds into ZK21_equatorial

old_docking_run=KCNQ2_model_crop_Z21_equatorial_20220405,Hernandez_rigid_2022_20220112,,20220405
docking_run=KCNQ2_model_crop_Z21_equatorial_20220405,project_20210913,,20220405
mkdir docking_runs/${docking_run}
cp docking_runs/${old_docking_run}/1_run.sh \
   docking_runs/${docking_run}/

pushd docking_runs/${docking_run}
  #sed -i "s/Hernandez_rigid_2022_20220112/project_20210913/g" 1_run.sh
  bash 1_run.sh
popd
# maybe the terminal benzene in the acrylamide series fits into the slot?


# Try docking into the ZK21 OMEGA model
for old_docking_run in docking_runs/KCNQ2_model_Z16_equatorial*; do
   docking_run=$(echo ${old_docking_run} | sed "s/KCNQ2_model_Z16_equatorial_20220309/KCNQ2_7CR2_docked_Z21_20220217/g")
   rm -rf ${docking_run}
   mkdir ${docking_run}
   cp ${old_docking_run}/1_run.sh ${docking_run}/1_run.sh
   sed -i "s/KCNQ2_model_Z16_equatorial_20220309/KCNQ2_7CR2_docked_Z21_20220217/g" ${docking_run}/1_run.sh
   pushd ${docking_run}    
   bash 1_run.sh
   popd
done


for docking_run in docking_runs/KCNQ2_7CR2_docked_Z21_20220217*; do
    pushd ${docking_run}
    rm -rf dock_statistics.* extract_all* pose*
    echo "Collecint dock results ..."
    ls results/ | grep -v joblist | sed "s#^#results/#" > dirlist
    python ~/opt/DOCK/analysis/extract_all_blazing_fast.py dirlist extract_all.txt 10
    python ~/opt/DOCK/analysis/getposes_blazing_faster.py '' extract_all.sort.uniq.txt 500 poses.mol2 test.mol2.gz 
    source ${DOCK_TEMPLATE}/scripts/dock_statistics.sh  
    Rscript ${DOCK_TEMPLATE}/scripts/analysis/gather_pose_features.R --verbose
  popd
done

# Try increase the dock sampling for the ZK21 omega
for old_docking_run in docking_runs/KCNQ2_model_Z16_equatorial*; do
   docking_run=$(echo ${old_docking_run} | sed "s/KCNQ2_model_Z16_equatorial_20220309/KCNQ2_7CR2_docked_Z21_extrasample_20220217/g")
   rm -rf ${docking_run}
   mkdir ${docking_run}
   cp ${old_docking_run}/1_run.sh ${docking_run}/1_run.sh
   sed -i "s/KCNQ2_model_Z16_equatorial_20220309/KCNQ2_7CR2_docked_Z21_extrasample_20220217/g" ${docking_run}/1_run.sh
   pushd ${docking_run}    
   bash 1_run.sh
   popd
done


for docking_run in docking_runs/KCNQ2_7CR2_docked_Z21_extrasample_20220217*; do
    pushd ${docking_run}
    rm -rf dock_statistics.* extract_all* pose*
    echo "Collecint dock results ..."
    ls results/ | grep -v joblist | sed "s#^#results/#" > dirlist
    python ~/opt/DOCK/analysis/extract_all_blazing_fast.py dirlist extract_all.txt 10
    python ~/opt/DOCK/analysis/getposes_blazing_faster.py '' extract_all.sort.uniq.txt 500 poses.mol2 test.mol2.gz 
    source ${DOCK_TEMPLATE}/scripts/dock_statistics.sh  
    Rscript ${DOCK_TEMPLATE}/scripts/analysis/gather_pose_features.R --verbose
  popd
done


# Try increase the dock sampling and use monte carlo for the ZK21 omega
for old_docking_run in docking_runs/KCNQ2_model_Z16_equatorial*; do
   docking_run=$(echo ${old_docking_run} | sed "s/KCNQ2_model_Z16_equatorial_20220309/KCNQ2_7CR2_docked_Z21_extrasample_MC_20220217/g")
   rm -rf ${docking_run}
   mkdir ${docking_run}
   cp ${old_docking_run}/1_run.sh ${docking_run}/1_run.sh
   sed -i "s/KCNQ2_model_Z16_equatorial_20220309/KCNQ2_7CR2_docked_Z21_extrasample_MC_20220217/g" ${docking_run}/1_run.sh
   pushd ${docking_run}    
   bash 1_run.sh
   popd
done


for docking_run in docking_runs/KCNQ2_7CR2_docked_Z21_extrasample_MC_20220217*; do
    pushd ${docking_run}
    rm -rf dock_statistics.* extract_all* pose*
    echo "Collecint dock results ..."
    ls results/ | grep -v joblist | sed "s#^#results/#" > dirlist
    python ~/opt/DOCK/analysis/extract_all_blazing_fast.py dirlist extract_all.txt 10
    python ~/opt/DOCK/analysis/getposes_blazing_faster.py '' extract_all.sort.uniq.txt 500 poses.mol2 test.mol2.gz 
    source ${DOCK_TEMPLATE}/scripts/dock_statistics.sh  
    Rscript ${DOCK_TEMPLATE}/scripts/analysis/gather_pose_features.R --verbose
  popd
done
