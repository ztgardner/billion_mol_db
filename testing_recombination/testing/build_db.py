
# TODO
# Create option to have detailed logs and debug in build_db.py
# Add functionality for builder.py to run on all rows without specifying the number
# Add functionality to evenly spread the ring distribution with splitter. It would be nice to be able to calculate how many outputs it will have and have a target size, like 1 million / proc

# import subprocess
# import multiprocessing
# import os
# import sys
# from pathlib import Path
# import glob
#
# # Define constants
# RING_FILE = "ring.parquet"
# SIDECHAIN_FILE = "sidechain.parquet"
# OUTPUT_DIR = Path("output")
# OUTPUT_DIR.mkdir(exist_ok=True)
#
# def get_available_cores():
#     """Returns the number of available CPU cores."""
#     return multiprocessing.cpu_count()
#
# def parse_user_requested_cores():
#     """Parses the number of cores requested by the user from command-line arguments."""
#     if len(sys.argv) > 1:
#         try:
#             requested_cores = int(sys.argv[1])
#             return max(1, requested_cores)  # Ensure at least 1 core is requested
#         except ValueError:
#             print("⚠️ Invalid core count provided. Using default settings.")
#     return None  # No user input, return None
#
# def determine_processor_count():
#     """Determines the number of processors to use, based on user input and system availability."""
#     available_cores = get_available_cores()
#     requested_cores = parse_user_requested_cores()
#
#     if requested_cores is not None:
#         if requested_cores > available_cores:
#             print(f"⚠️ {requested_cores} cores requested, but only {available_cores} available.")
#             print(f"▶️ Running with {available_cores} processors instead.")
#             return available_cores
#         return requested_cores
#     return available_cores  # Default to using all available cores
#
# def run_splitter(num_processors):
#     """Runs splitter.py with the correct number of processors."""
#     print("\n🚀 Running splitter.py...")
#     subprocess.run([sys.executable, "splitter.py", str(num_processors)], check=True)
#
#     # Count actual number of created files
#     ring_files = sorted(glob.glob("ring_*.parquet"))  # Get list of files
#     detected_processors = len(ring_files)
#
#     print(f"\n🔢 Detected {detected_processors} ring files. Using {num_processors} processors.")
#     return min(num_processors, detected_processors)  # Ensure we don't exceed actual files
#
# def run_combinations(procnumber):
#     """Runs combinations.py for a specific processor."""
#     ring_file = f"ring_{procnumber}.parquet"
#     if not Path(ring_file).exists():
#         print(f"⚠️ Skipping processor {procnumber}: {ring_file} does not exist.")
#         return
#     print(f"▶️ Running combinations.py for processor {procnumber}...")
#     subprocess.run([sys.executable, "combinations.py", str(procnumber)], check=True)
#
# def run_builder(procnumber):
#     """Runs builder.py for a specific processor."""
#     combinations_file = f"output/combinations_{procnumber}.parquet"
#     if not Path(combinations_file).exists():
#         print(f"⚠️ Skipping processor {procnumber}: {combinations_file} does not exist.")
#         return
#     print(f"▶️ Running builder.py for processor {procnumber}...")
#     subprocess.run([sys.executable, "builder.py", str(procnumber)], check=True)
#
# def run_parallel_processing(script_func, num_processors):
#     """Runs a given function in parallel for all available processor datasets."""
#     with multiprocessing.Pool(num_processors) as pool:
#         pool.map(script_func, range(1, num_processors + 1))
#
# def run_appender():
#     """Runs appender.py to merge all combination files into a final parquet file."""
#     print("\n🚀 Running appender.py to combine all parquet files...")
#     subprocess.run([sys.executable, "appender.py"], check=True)
#
# if __name__ == "__main__":
#     print("==========================================")
#     print("🔥 Starting full workflow 🔥")
#     print("==========================================\n")
#
#     # Determine the number of processors to use
#     NPROC = determine_processor_count()
#
#     # Step 1: Run the splitter with user-specified processor count
#     NPROC = run_splitter(NPROC)
#
#     # Step 2: Run combinations.py in parallel for all processors
#     print("\n🚀 Running combinations.py in parallel...")
#     run_parallel_processing(run_combinations, NPROC)
#
#     # Step 3: Run builder.py in parallel for all processors
#     print("\n🚀 Running builder.py in parallel...")
#     run_parallel_processing(run_builder, NPROC)
#
#     # Step 4: Run appender.py to merge results
#     run_appender()
#
#     print("\n✅ Workflow completed successfully!")

import subprocess
import multiprocessing
import os
import sys
import time
from pathlib import Path
import glob

# Define constants
RING_FILE = "ring.parquet"
SIDECHAIN_FILE = "sidechain.parquet"
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

# Store runtimes
runtimes = {}


def get_available_cores():
    """Returns the number of available CPU cores."""
    return multiprocessing.cpu_count()

def parse_user_requested_cores():
    """Parses the number of cores requested by the user from command-line arguments."""
    if len(sys.argv) > 1:
        try:
            requested_cores = int(sys.argv[1])
            return max(1, requested_cores)  # Ensure at least 1 core is requested
        except ValueError:
            print("⚠️ Invalid core count provided. Using default settings.")
    return None  # No user input, return None

def determine_processor_count():
    """Determines the number of processors to use, based on user input and system availability."""
    available_cores = get_available_cores()
    requested_cores = parse_user_requested_cores()

    if requested_cores is not None:
        if requested_cores > available_cores:
            print(f"⚠️ {requested_cores} cores requested, but only {available_cores} available.")
            print(f"▶️ Running with {available_cores} processors instead.")
            return available_cores
        return requested_cores
    return available_cores  # Default to using all available cores

def time_function(func, *args, **kwargs):
    """Measures the execution time of a function while passing all arguments correctly."""
    start_time = time.time()
    result = func(*args, **kwargs)  # Pass *args and **kwargs to the function
    end_time = time.time()
    return result, round(end_time - start_time, 3)

def run_splitter(num_processors):
    """Runs splitter.py with the correct number of processors and measures execution time."""
    print("\n🚀 Running splitter.py...")
    _, runtime = time_function(subprocess.run, [sys.executable, "splitter.py", str(num_processors)], check=True)
    runtimes["splitter.py"] = runtime

    # Count actual number of created files
    ring_files = sorted(glob.glob("ring_*.parquet"))  # Get list of files
    detected_processors = len(ring_files)

    print(f"\n🔢 Detected {detected_processors} ring files. Using {num_processors} processors.")
    return min(num_processors, detected_processors)  # Ensure we don't exceed actual files

def run_combinations(procnumber):
    """Runs combinations.py for a specific processor and returns execution time."""
    ring_file = f"ring_{procnumber}.parquet"
    if not Path(ring_file).exists():
        print(f"⚠️ Skipping processor {procnumber}: {ring_file} does not exist.")
        return procnumber, 0  # Return processor number and 0 time

    print(f"▶️ Running combinations.py for processor {procnumber}...")
    _, runtime = time_function(subprocess.run, [sys.executable, "combinations.py", str(procnumber)], **{"check": True})

    return procnumber, runtime  # Return processor number and runtime

def run_builder(procnumber):
    """Runs builder.py for a specific processor and returns execution time."""
    combinations_file = f"output/combinations_{procnumber}.parquet"
    if not Path(combinations_file).exists():
        print(f"⚠️ Skipping processor {procnumber}: {combinations_file} does not exist.")
        return procnumber, 0  # Return processor number and 0 time

    print(f"▶️ Running builder.py for processor {procnumber}...")
    _, runtime = time_function(subprocess.run, [sys.executable, "builder.py", str(procnumber)], **{"check": True})

    return procnumber, runtime  # Return processor number and runtime

def run_appender():
    """Runs appender.py to merge all combination files into a final parquet file and measures execution time."""
    print("\n🚀 Running appender.py to combine all parquet files...")
    _, runtime = time_function(subprocess.run, [sys.executable, "appender.py"], check=True)
    runtimes["appender.py"] = runtime

def run_parallel_processing(script_func, num_processors, script_name):
    """Runs a given function in parallel for all available processor datasets and captures runtimes."""
    with multiprocessing.Pool(num_processors) as pool:
        results = pool.starmap(script_func, [(i,) for i in range(1, num_processors + 1)])  # ✅ Uses starmap

    # Store runtimes
    script_runtime = {}
    for proc_id, runtime in results:
        script_runtime[proc_id] = runtime

    runtimes[script_name] = script_runtime  # ✅ Store per-processor times


def print_runtimes():
    """Prints the summary of script runtimes."""
    print("\n==========================================")
    print("🕒 **Runtimes Summary**")
    print("==========================================\n")

    for script, time_taken in runtimes.items():
        if isinstance(time_taken, dict):
            # For multi-processor scripts
            total_time = sum(time_taken.values())
            print(f"{script} --- {total_time:.3f} s")
            for proc, proc_time in time_taken.items():
                print(f"\tproc {proc} --- {proc_time:.3f} s")
        else:
            print(f"{script} --- {time_taken:.3f} s")


if __name__ == "__main__":
    print("==========================================")
    print("🔥 Starting full workflow 🔥")
    print("==========================================\n")

    overall_start_time = time.time()  # Start overall timer

    # Determine the number of processors to use
    NPROC = determine_processor_count()

    # Step 1: Run the splitter with user-specified processor count
    NPROC = run_splitter(NPROC)

    # Step 2: Run combinations.py in parallel for all processors
    print("\n🚀 Running combinations.py in parallel...")
    run_parallel_processing(run_combinations, NPROC, "combinations.py")  # ✅ Fixed

    # Step 3: Run builder.py in parallel for all processors
    print("\n🚀 Running builder.py in parallel...")
    run_parallel_processing(run_builder, NPROC, "builder.py")  # ✅ Fixed

    # Step 4: Run appender.py to merge results
    run_appender()

    overall_end_time = time.time()  # End overall timer
    runtimes["Total Runtime"] = round(overall_end_time - overall_start_time, 3)

    # Print runtime summary
    print_runtimes()

    print("\n✅ Workflow completed successfully!")
