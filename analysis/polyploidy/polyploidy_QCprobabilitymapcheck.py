import h5py
import tifffile
import matplotlib.pyplot as plt

#.h5 probability map file
probability_path = r"D:\UNI\POLYPLOIDY_ANALYSIS\DF.WT\tube1.slide10\gut4\DF.WT.tube1.LHS4.6_DAPI_Probabilities.h5"

with h5py.File(probability_path, 'r') as f:
    data = f['exported_data'][()]
    print(f"Probability map shape: {data.shape}")  # Should be (Z, Y, X, C)

    z_mid = data.shape[0] // 2  # Choose middle Z-slice for inspection

    # Plot each channel of the middle slice
    fig, axes = plt.subplots(1, data.shape[-1], figsize=(10, 5))
    for i in range(data.shape[-1]):
        axes[i].imshow(data[z_mid, :, :, i], cmap='gray')
        axes[i].set_title(f'Channel {i}')
        axes[i].axis('off')

    plt.tight_layout()
    plt.show()
