import h5py
import numpy as np
import tifffile
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops

# File paths
dapi_path = r"D:\UNI\POLYPLOIDY_ANALYSIS\DF.WT\tube1.slide10\gut4\DF.WT.tube1.LHS4.6_DAPI.tif"
probability_path = r"D:\UNI\POLYPLOIDY_ANALYSIS\DF.WT\tube1.slide10\gut4\DF.WT.tube1.LHS4.6_DAPI_Probabilities.h5"

# Load DAPI image
dapi_stack = tifffile.imread(dapi_path)

# Load probability map (Channel 0 = nuclei)
with h5py.File(probability_path, "r") as f:
    probability_data = f["exported_data"][:]  # shape: (Z, Y, X, C)
nuclei_prob = probability_data[..., 0]

# Apply threshold to get binary mask
binary_mask = nuclei_prob > 0.6 # adjust

# Label connected components in 3D
labeled_mask = label(binary_mask)

# Measure DAPI intensity in each labeled nucleus
props = regionprops(labeled_mask, intensity_image=dapi_stack)
measurements = [{"label": p.label, "mean_intensity": p.mean_intensity, "area": p.area} for p in props]

# Preview: overlay a mid-Z slice for visual inspection
mid_z = dapi_stack.shape[0] // 2

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
axes[0].imshow(dapi_stack[mid_z], cmap="gray")
axes[0].set_title("DAPI Slice")
axes[1].imshow(nuclei_prob[mid_z], cmap="gray")
axes[1].set_title("Probability Map (Channel 0)")
axes[2].imshow(binary_mask[mid_z], cmap="gray")
axes[2].set_title("Thresholded Mask @ 0.1") #adjust label to match threshold
for ax in axes:
    ax.axis("off")
plt.tight_layout()
plt.show()

# Print a few measurement results
for m in measurements[:5]:
    print(m)
