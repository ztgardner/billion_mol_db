import subprocess
import multiprocessing
import os
import sys
from pathlib import Path
import glob

# Define constants
RING_FILE = "ring.parquet"
SIDECHAIN_FILE = "sidechain.parquet"
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)


def run_splitter():
    """Runs splitter.py to split ring.parquet into multiple processor-specific files."""
    print("\n🚀 Running splitter.py...")
    subprocess.run([sys.executable, "splitter.py", RING_FILE, SIDECHAIN_FILE], check=True)

    # Count actual number of created files
    ring_files = sorted(glob.glob("ring_*.parquet"))  # Get list of files
    num_processors = len(ring_files)

    print(f"\n🔢 Detected {num_processors} processors based on available ring_*.parquet files.")
    return num_processors


def run_combinations(procnumber):
    """Runs combinations.py for a specific processor."""
    ring_file = f"ring_{procnumber}.parquet"
    if not Path(ring_file).exists():
        print(f"⚠️ Skipping processor {procnumber}: {ring_file} does not exist.")
        return
    print(f"▶️ Running combinations.py for processor {procnumber}...")
    subprocess.run([sys.executable, "combinations.py", str(procnumber)], check=True)


def run_builder(procnumber):
    """Runs builder.py for a specific processor."""
    combinations_file = f"output/combinations_{procnumber}.parquet"
    if not Path(combinations_file).exists():
        print(f"⚠️ Skipping processor {procnumber}: {combinations_file} does not exist.")
        return
    print(f"▶️ Running builder.py for processor {procnumber}...")
    subprocess.run([sys.executable, "builder.py", str(procnumber)], check=True)


def run_parallel_processing(script_func, num_processors):
    """Runs a given function in parallel for all available processor datasets."""
    with multiprocessing.Pool(num_processors) as pool:
        pool.map(script_func, range(1, num_processors + 1))


def run_appender():
    """Runs appender.py to merge all combination files into a final parquet file."""
    print("\n🚀 Running appender.py to combine all parquet files...")
    subprocess.run([sys.executable, "appender.py"], check=True)


if __name__ == "__main__":
    print("==========================================")
    print("🔥 Starting full workflow 🔥")
    print("==========================================\n")

    # Step 1: Run the splitter
    NPROC = run_splitter()

    # Step 2: Run combinations.py in parallel for all processors
    print("\n🚀 Running combinations.py in parallel...")
    run_parallel_processing(run_combinations, NPROC)

    # Step 3: Run builder.py in parallel for all processors
    print("\n🚀 Running builder.py in parallel...")
    run_parallel_processing(run_builder, NPROC)

    # Step 4: Run appender.py to merge results
    run_appender()

    print("\n✅ Workflow completed successfully!")
