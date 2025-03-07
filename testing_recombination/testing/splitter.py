# import pandas as pd
# import pyarrow as pa  # Importing pyarrow directly
# import pyarrow.parquet as pq
# import os
# import multiprocessing
#
# def get_slurm_processors():
#     """
#     Determines the number of available processors in a SLURM environment.
#     If SLURM_CPUS_PER_TASK is set, use that; otherwise, fall back to multiprocessing count.
#     """
#     slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
#     if slurm_cpus:
#         try:
#             return int(slurm_cpus)
#         except ValueError:
#             pass  # If it's not an integer, ignore and fall back to default
#
#     return multiprocessing.cpu_count()
#
# # Set the number of processors
# nproc = get_slurm_processors()
# print(f"Using {nproc} processors for splitting.")
#
#
#
# def split_parquet(input_parquet, output_prefix, nproc):
#     """Split a Parquet file into nproc smaller Parquet files."""
#     table = pq.read_table(input_parquet)
#     df = table.to_pandas()
#     split_dfs = [df.iloc[i::nproc] for i in range(nproc)]
#
#     for i, split_df in enumerate(split_dfs, 1):
#         pq.write_table(pa.Table.from_pandas(split_df), f"{output_prefix}_{i}.parquet")  # ✅ Use pa.Table
#         print(f"Written: {output_prefix}_{i}.parquet")
#
#
# if __name__ == "__main__":
#     nproc = get_slurm_processors()  # Set max CPUs to 4
#     print(f"Using {nproc} processors for splitting.")
#
#     split_parquet("ring.parquet", "ring", nproc)
#
#     # Load all ring files and sidechain into memory
#     ring_data = {i: pq.read_table(f"ring_{i}.parquet").to_pandas() for i in range(1, nproc + 1)}
#     sidechain_data = pq.read_table("sidechain.parquet").to_pandas()
#
#     print("All data loaded into memory.")

import sys
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

def split_parquet(input_parquet, output_prefix, nproc):
    """Split a Parquet file into nproc smaller Parquet files."""
    table = pq.read_table(input_parquet)
    df = table.to_pandas()

    split_dfs = [df.iloc[i::nproc] for i in range(nproc)]
    for i, split_df in enumerate(split_dfs, 1):
        pq.write_table(pa.Table.from_pandas(split_df), f"{output_prefix}_{i}.parquet")
        print(f"✅ Written: {output_prefix}_{i}.parquet")

if __name__ == "__main__":
    # Ensure we get the correct number of processors from build_db.py
    if len(sys.argv) > 1:
        try:
            nproc = int(sys.argv[1])  # Read from command-line argument
        except ValueError:
            print("⚠️ Invalid processor count. Defaulting to 1.")
            nproc = 1
    else:
        print("⚠️ No processor count provided. Defaulting to 1.")
        nproc = 1  # Default to 1 processor if not specified

    print(f"🚀 Using {nproc} processors for splitting.")

    split_parquet("ring.parquet", "ring", nproc)

    # Load all ring files and sidechain into memory
    ring_data = {i: pq.read_table(f"ring_{i}.parquet").to_pandas() for i in range(1, nproc + 1)}
    sidechain_data = pq.read_table("sidechain.parquet").to_pandas()

    print("✅ All data loaded into memory.")
