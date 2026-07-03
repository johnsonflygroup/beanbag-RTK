import os
import pandas as pd
import re

# Define paths
measurements_folder = "D:/UNI/POLYPLOIDY_DAPI_ANALYSIS/filtered_measurements"
output_csv = "D:/UNI/POLYPLOIDY_DAPI_ANALYSIS/summary_with_metadata.csv"

# Pattern to extract metadata from filename
filename_pattern = re.compile(
    r'(?P<genotype>[A-Za-z0-9.]+)\.tube(?P<tube>\d+)\.(?P<side>LHS|RHS)(?P<gut>\d+)\.(?P<image>\d+)_DAPI'
)

# Initialize summary list
summary_data = []

# Process each CSV file
for filename in os.listdir(measurements_folder):
    if filename.endswith(".csv"):
        filepath = os.path.join(measurements_folder, filename)
        try:
            df = pd.read_csv(filepath)
            n_nuclei = len(df)
            total_intensity = df["integrated_intensity"].sum()
            mean_intensity = total_intensity / n_nuclei if n_nuclei > 0 else 0

            name_without_ext = filename.replace(".csv", "")
            match = filename_pattern.match(name_without_ext)
            if match:
                metadata = match.groupdict()
                summary_data.append({
                    "filename": name_without_ext,
                    "genotype": metadata["genotype"],
                    "tube": int(metadata["tube"]),
                    "side": metadata["side"],
                    "gut_number": int(metadata["gut"]),
                    "image_number": int(metadata["image"]),
                    "n_nuclei": n_nuclei,
                    "total_integrated_intensity": total_intensity,
                    "mean_intensity_per_image": mean_intensity
                })
            else:
                print(f"Filename didn't match pattern: {filename}")

        except Exception as e:
            print(f"Error processing {filename}: {e}")

# Save the summary CSV
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(output_csv, index=False)

print("Summary file saved to:", output_csv)
