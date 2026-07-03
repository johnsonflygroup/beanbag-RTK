import numpy as np
import tifffile
import h5py
import os
import pandas as pd
from skimage.measure import regionprops_table, label

# Set file paths
dapi_path = r"D:\UNI\POLYPLOIDY_ANALYSIS\DF.WT\tube1.slide10\gut4\DF.WT.tube1.LHS4.6_DAPI.tif"
prob_path = r"D:\UNI\POLYPLOIDY_ANALYSIS\DF.WT\tube1.slide10\gut4\DF.WT.tube1.LHS4.6_DAPI_Probabilities.h5"

# Load data
dapi_stack = tifffile.imread(dapi_path)

with h5py.File(prob_path, 'r') as f:
    probs = f['exported_data'][:]

# Threshold and label nuclei
threshold = 0.6 #adjust based on QC
nuclei_mask = probs[..., 0] > threshold  # Channel 0 assumed to be nuclei
labeled_nuclei = label(nuclei_mask)

# Measure properties 
props = regionprops_table(
    labeled_nuclei,
    intensity_image=dapi_stack,
    properties=('label', 'area', 'mean_intensity')
)

df = pd.DataFrame(props)

# Calculate integrated intensity 
df["integrated_intensity"] = df["mean_intensity"] * df["area"]

# Filter out small regions 
min_area = 30  # in pixels
before = len(df)
df_filtered = df[df["area"] >= min_area].copy()
after = len(df_filtered)

print(f"Original nuclei count: {before}")
print(f"Nuclei kept (area ≥ {min_area} pixels): {after}")
print(f"Nuclei excluded: {before - after}")

# Output CSV 
output_csv = dapi_path.replace(".tif", "_nuclei_measurements.csv")
df_filtered.to_csv(output_csv, index=False)
print(f"Saved measurements to: {output_csv}")
