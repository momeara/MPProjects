

import pandas as pd
from rdkit import Chem

substances = pd.read_csv("raw_data/substances_20210323.tsv", sep = "\t")

rdkit_smiles_all = []
for i in range(substances.shape[0]):
    substance_name = substances['substance_name'][i]
    input_smiles = substances['substance_smiles'][i]
    try:
        m = Chem.MolFromSmiles(substances['substance_smiles'][i])
        Chem.SanitizeMol(m)
        rdkit_smiles_isomeric = Chem.MolToSmiles(m)
        rdkit_smiles = Chem.MolToSmiles(m, isomericSmiles=False)
    except:
        print(f"failed to read/write smiles for substance {substance_name}")
        print(f"  input: {input_smiles}")
        rdkit_smiles_all.append(None)
        continue

    if input_smiles == rdkit_smiles_isomeric:
        rdkit_smiles_all.append(rdkit_smiles_isomeric)
        continue
    elif input_smiles == rdkit_smiles:
        rdkit_smiles_all.append(rdkit_smiles)
        continue
    else:
        print(f"substance_name: {substance_name}")
        print(f"   input:  {input_smiles}")
        print(f"   rdkit i:{rdkit_smiles_isomeric}")
        print(f"   rdkit:  {rdkit_smiles}")
        rdkit_smiles_all.append(rdkit_smiles_isomeric)

substances['substance_smiles_rdkit'] = rdkit_smiles_all
substances.to_csv(
    path_or_buf = "intermediate_data/substances_sanitized_20210323.tsv",
    sep = "\t",
    index = False)
