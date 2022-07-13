



# Install hugging face transformers package
#
#    pip install transformers[torch]
#

# Install selfies for robust chemical representation
# https://github.com/aspuru-guzik-group/selfies
#
#   pip install selfies
#



# https://huggingface.co/ncfrey/ChemGPT-4.7M
# https://huggingface.co/ncfrey/ChemGPT-19M
# https://huggingface.co/ncfrey/ChemGPT-1.2B

import torch
from transformers import pipeline
from transformers import AutoTokenizer, AutoModelForCausalLM
import selfies
import pyarrow.parquet
import pyarrow as pa

import pandas as pd

tokenizer = AutoTokenizer.from_pretrained("ncfrey/ChemGPT-4.7M")
model = AutoModelForCausalLM.from_pretrained("ncfrey/ChemGPT-4.7M")


def embed_substance(substance_smiles, model, tokenizer):
    try:
        substance_selfies = selfies.encoder(substance_smiles)
        substance_tokens = tokenizer(substance_selfies)
        with torch.no_grad():
            substance_embedding = list(
                model.forward(
                    input_ids = torch.tensor(substance_tokens.input_ids)).values())[0].mean(0)
        return substance_embedding.numpy()
    except:
        print(f"Failed to embed smiles '{substance_smiles}'")
        return None


substances = pd.read_csv("intermediate_data/substances_sanitized_20210323.tsv", sep = "\t")
substance_ids = substances["substance_name"]

substance_embeddings = []
substance_ids_generated = []
for index, substance_smiles in enumerate(substances["substance_smiles_rdkit"]):
    if index % 20 == 0:
        print(f"Processing molecule {i}")

    try:
        substance_embedding = embed_substance(
            substance_smiles = substance_smiles,
            model = model,
            tokenizer = tokenizer)
    except:
        continue

    substance_embeddings.append[substance_embedding]
    substance_ids_generated.append(substance_ids[index])
    
substance_embeddings = np.array(substance_embeddings)
    
feature_columns = pd.DataFrame({
    'feature' : [ f"feature_{i}" for i in range(substance_embeddings.shape[1])],
    'transform' : ["identity" for i in range(substance_embeddings.shape[1])]})
feature_columns.to_csv(
    "../intermediate_data/project_substances_ChemGPT-4.7M/fingerprint_feature_columns.tsv", sep = "\t")

fingerprints_df = pd.DataFrame(substance_embeddings, columns = feature_columns['feature'])
if substance_ids.dtype == np.dtype('O'):
    substance_ids = substance_ids.astype(str)

fingerprints_df.insert(loc=0, column="substance_id", value=substance_ids)
print(f"Writing fingerprints DataFrame to disk '../intermediate_data/project_substances_ChemGPT-4.7M/fingerprints.parquet")
pa.parquet.write_table(
    table=pa.Table.from_pandas(fingerprints_df),
    where = f"../intermediate_data/project_substances_ChemGPT-4.7M/fingerprints.parquet")

