import os
import numpy as np
import tifffile
import h5py
import pandas as pd
from skimage.measure import label, regionprops_table

# Paths
input_dir = r"D:/UNI/POLYPLOIDY_ALL_DAPI"
output_dir = os.path.join(input_dir, "measurements")
os.makedirs(output_dir, exist_ok=True)

# Parameters
threshold = 0.6
min_size = 30  # in pixels

# Process each file
for file in os.listdir(input_dir):
    if file.endswith(".tif"):
        base_name = os.path.splitext(file)[0]
        dapi_path = os.path.join(input_dir, file)
        prob_path = os.path.join(input_dir, f"{base_name}_Probabilities.h5")
        
        if not os.path.exists(prob_path):
            print(f"Skipping {file} – no matching probability map found.")
            continue

        try:
            # Load DAPI image and probability map
            dapi = tifffile.imread(dapi_path)

            with h5py.File(prob_path, "r") as f:
                prob = f["exported_data"][:]
                nuclei_channel = prob[..., 0]  # Assuming channel 0 is nuclei

            # Threshold and label
            mask = nuclei_channel > threshold
            labeled = label(mask)

            # Filter small objects
            region_props = regionprops_table(
                labeled,
                intensity_image=dapi,
                properties=("label", "area", "mean_intensity", "intensity_image")
            )

            # Convert to DataFrame
            df = pd.DataFrame(region_props)
            df["integrated_intensity"] = df["area"] * df["mean_intensity"]

            # Remove small objects
            df = df[df["area"] >= min_size]

            # Save
            csv_path = os.path.join(output_dir, f"{base_name}.csv")
            df.to_csv(csv_path, index=False)
            print(f"Saved: {csv_path}")

        except Exception as e:
            print(f"Error processing {file}: {e}")
