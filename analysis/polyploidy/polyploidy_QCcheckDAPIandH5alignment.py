import tifffile
import h5py
import os

# Path to files
dapi_path = r"D:\UNI\POLYPLOIDY_ANALYSIS\DF.WT\tube1.slide10\gut4\DF.WT.tube1.LHS4.6_DAPI.tif"
prob_path = r"D:\UNI\POLYPLOIDY_ANALYSIS\DF.WT\tube1.slide10\gut4\DF.WT.tube1.LHS4.6_DAPI_Probabilities.h5"

# Load the DAPI image
dapi_image = tifffile.imread(dapi_path)
print(f"DAPI image shape: {dapi_image.shape}")

# Load the probability map
with h5py.File(prob_path, 'r') as f:
    # You may need to inspect keys; often it's 'exported_data'
    print(f"HDF5 keys: {list(f.keys())}")
    prob_map = f['exported_data'][()]
    print(f"Probability map shape: {prob_map.shape}")
