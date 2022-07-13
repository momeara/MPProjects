

import pandas as pd
import numpy as np
import pyarrow.parquet
import pyarrow as pa


data_path = "/nfs/turbo/umms-jzsexton/MPalign/covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet"

data = pa.parquet.read_table(
    source = data_path).to_pandas()


for column_name in data.column_names:
    print(column_name)

graph = []
for image in data:
    for i in rows of :
        for j in rows of data:
            if distance(cell_i, cell_j) < threshold:
               add edge(cell_i, cell_j) to graph
