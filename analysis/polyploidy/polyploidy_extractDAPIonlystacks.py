import os
import tifffile
import numpy as np

# Path to your raw data
input_root = r"D:/UNI/POLYPLOIDY"
# Path to save DAPI-only outputs
output_root = r"D:/UNI/POLYPLOIDY_ANALYSIS"

# Helper function to ensure output folders exist
def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

# Walk through every folder
for root, dirs, files in os.walk(input_root):
    for file in files:
        if file.lower().endswith(".tif") or file.lower().endswith(".tiff"):
            filepath = os.path.join(root, file)
            print(f"Processing: {filepath}")
            try:
                # Load image
                stack = tifffile.imread(filepath)

                # If image is 5D: (T, Z, C, Y, X) or (Z, C, Y, X)
                if stack.ndim == 5:
                    # Assume format (T, Z, C, Y, X)
                    stack = stack[0]  # Take T=0

                if stack.ndim == 4:
                    # Format: (Z, C, Y, X)
                    z, c, y, x = stack.shape
                    dapi_stack = stack[:, 0, :, :]  # Channel 0 assumed to be DAPI
                elif stack.ndim == 3:
                    print(f"Skipping: {filepath} – already single-channel?")
                    continue
                else:
                    print(f"Skipping: {filepath} – unknown format")
                    continue

                # Determine relative path
                rel_path = os.path.relpath(root, input_root)
                output_dir = os.path.join(output_root, rel_path)
                ensure_dir(output_dir)

                # Save DAPI-only stack
                output_file = os.path.join(output_dir, f"{os.path.splitext(file)[0]}_DAPI.tif")
                tifffile.imwrite(output_file, dapi_stack.astype(np.uint16))
            except Exception as e:
                print(f"Error processing {filepath}: {e}")

