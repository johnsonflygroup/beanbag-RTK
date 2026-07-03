import os, glob, re
import numpy as np
import h5py
import tifffile
import matplotlib.pyplot as plt
import pandas as pd
from skimage import filters, measure, morphology, segmentation, color
from skimage.measure import marching_cubes, mesh_surface_area
from scipy import ndimage as ndi
from skimage.segmentation import watershed
from skimage.feature import peak_local_max

# --- CONFIG ---
fas3_dir = r"D:\UNI\HUB_CELLS\EXPERIMENTAL\Hub cells\SplitChannels\Fas3H5"
dapi_dir = r"D:\UNI\HUB_CELLS\EXPERIMENTAL\Hub cells\SplitChannels\DAPIH5"
output_dir = r"D:\UNI\HUB_CELLS\EXPERIMENTAL\Hub cells\Analysisfinal"
os.makedirs(output_dir, exist_ok=True)

fas3_threshold = 0.25
dapi_threshold = 0.5
GENOTYPES = ["t2a.df", "t2a.wt", "df.wt", "wt"]

# --- Helper: parse genotype ---
def parse_genotype(name):
    lname = name.lower()
    for gt in GENOTYPES:
        if gt in lname:
            return gt
    return "unknown"

# --- Helper: QC panels ---
def save_panel(images, titles, outpath, suptitle):
    cols = len(images)
    fig, axes = plt.subplots(1, cols, figsize=(4*cols, 4))
    for ax, img, title in zip(axes, images, titles):
        ax.imshow(img, cmap="gray")
        ax.set_title(title)
        ax.axis("off")
    plt.suptitle(suptitle)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

# --- Hub metrics including surface area & sphericity ---
def extract_hub_metrics(mask, label_id):
    props = measure.regionprops(mask)[label_id - 1]
    metrics = {
        "hub_label": label_id,
        "hub_volume_voxels": props.area,
        "hub_centroid_z": props.centroid[0],
        "hub_centroid_y": props.centroid[1],
        "hub_centroid_x": props.centroid[2],
        "hub_bbox": props.bbox,
        "hub_equivalent_diameter": props.equivalent_diameter,
        "hub_extent": props.extent,
        "hub_solidity": props.solidity if props.solidity != np.inf else np.nan,
        "hub_major_axis_length": props.major_axis_length,
        "hub_minor_axis_length": props.minor_axis_length,
    }

    # Check if object spans multiple Z-slices
    zmin, zmax = props.bbox[0], props.bbox[3]
    if zmax - zmin < 2:
        print(f"Skipping hub {label_id} in {core}: flat in Z (z-range = {zmax - zmin})")
        metrics["hub_surface_area"] = np.nan
        metrics["hub_sphericity"] = np.nan
        return metrics

    try:
        verts, faces, _, _ = marching_cubes(mask == label_id, level=0.5)
        surface_area = mesh_surface_area(verts, faces)
        volume = props.area
        sphericity = (np.pi**(1/3) * (6*volume)**(2/3)) / surface_area
        metrics["hub_surface_area"] = surface_area
        metrics["hub_sphericity"] = sphericity
    except Exception:
        metrics["hub_surface_area"] = np.nan
        metrics["hub_sphericity"] = np.nan

    return metrics


# --- Utility: strip for pairing ---
def strip_core(name):
    name = os.path.basename(name)
    name = re.sub(r"^C[12]-", "", name)
    name = name.replace(".tif_Fas3_Probabilities.h5", "")
    name = name.replace(".h5", "")
    return name

# --- Main loop ---
fas3_files = sorted(glob.glob(os.path.join(fas3_dir, "*_Fas3_Probabilities.h5")))
dapi_files = sorted(glob.glob(os.path.join(dapi_dir, "*.h5")))
rows = []

for fas3_path in fas3_files:
    base = os.path.basename(fas3_path)
    core = strip_core(base)
    genotype = parse_genotype(base)

    # Match DAPI
    dapi_matches = [d for d in dapi_files if strip_core(d) == core]
    if not dapi_matches:
        print("No DAPI match for", base)
        continue
    dapi_path = dapi_matches[0]

    # Load Fas3 probability map
    try:
        with h5py.File(fas3_path, "r") as f:
            fas3_probs = f["exported_data"][:]
        fas3_prob = fas3_probs[:, :, :, 0]  
    except Exception as e:
        print("Error loading Fas3:", base, e)
        continue

    # Load DAPI probability map
    try:
        with h5py.File(dapi_path, "r") as f:
            dapi_probs = f["exported_data"][:]
        dapi_prob = dapi_probs[:, :, :, 0]
    except Exception as e:
        print("Error loading DAPI:", dapi_path, e)
        continue

    fas3_prob = fas3_prob.astype(np.float32)
    dapi_prob = dapi_prob.astype(np.float32)
    fas3_max = fas3_prob.max()
    print(f"{core} Fas3 prob min/max:", fas3_prob.min(), fas3_max)
    if fas3_max > 0:
        fas3_prob /= fas3_max
    else:
        print(f"Warning: Fas3 prob max is 0 in {core}")

    dapi_max = dapi_prob.max()
    print(f"{core} DAPI prob min/max:", dapi_prob.min(), dapi_max)
    if dapi_max > 0:
        dapi_prob /= dapi_max
    else:
        print(f"Warning: DAPI prob max is 0 in {core}")


    fas3_mask = (fas3_prob > fas3_threshold).astype(np.uint8)
    fas3_mask = morphology.remove_small_objects(fas3_mask.astype(bool), min_size=100)
    fas3_mask = morphology.remove_small_holes(fas3_mask, area_threshold=100)
    fas3_mask = fas3_mask.astype(np.uint8)
    print(f"{core} Hub mask unique values:", np.unique(fas3_mask))

    dapi_mask = (dapi_prob > dapi_threshold).astype(np.uint8)
    dapi_mask = morphology.remove_small_objects(dapi_mask.astype(bool), min_size=50)
    dapi_mask = morphology.remove_small_holes(dapi_mask, area_threshold=50)
    dapi_mask = dapi_mask.astype(np.uint8)
    print(f"{core} DAPI mask unique values:", np.unique(dapi_mask))

    # QC panels
    print("Fas3 prob shape:", fas3_prob.shape)
    print("Fas3 mask shape:", fas3_mask.shape)
    save_panel([fas3_prob.max(axis=0), fas3_mask.max(axis=0)],
               ["Fas3 Prob (MIP)", f"Mask>{fas3_threshold}"],
               os.path.join(output_dir, f"{core}_Fas3_QC.png"),
               f"{genotype} Fas3 QC")
    
    print("DAPI prob shape:", dapi_prob.shape)
    print("DAPI mask shape:", dapi_mask.shape)
    save_panel([dapi_prob.max(axis=0), dapi_mask.max(axis=0)],
               ["DAPI Prob (MIP)", f"Mask>{dapi_threshold}"],
               os.path.join(output_dir, f"{core}_DAPI_QC.png"),
               f"{genotype} DAPI QC")

    # Label hub objects
    hub_labels = measure.label(fas3_mask)
    hub_count = hub_labels.max()
    if hub_count == 0:
        print("No hubs detected in", base)
        continue

    if hub_count > 1:
        mip = color.label2rgb(hub_labels.max(axis=0), bg_label=0)
        fig, ax = plt.subplots(figsize=(6,6))
        ax.imshow(mip)
        for region in measure.regionprops(hub_labels.max(axis=0)):
            y, x = region.centroid
            ax.text(x, y, str(region.label), color='white', fontsize=8, ha='center')
        ax.set_title(f"{genotype} Hub Object Map")
        ax.axis("off")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{core}_HubObjects_Labeled.png"), dpi=200)
        plt.close()


    # Segment nuclei
    distance = ndi.distance_transform_edt(dapi_mask)
    coords = peak_local_max(distance, labels=dapi_mask, footprint=np.ones((7,7,7)))
    if len(coords) > 1000:
        print(f"Warning: {len(coords)} markers found — too many?")
    markers = np.zeros_like(dapi_mask, dtype=int)
    for i, (z,y,x) in enumerate(coords, start=1):
        markers[z,y,x] = i
    nuclei_labels = watershed(-distance, markers, mask=dapi_mask)
    nuclei_props = measure.regionprops(nuclei_labels)
    def is_compact(nuc):
        return nuc.solidity > 0.75 and nuc.extent > 0.4
    nuclei_props = [n for n in nuclei_props if n.area >= 50 and is_compact(n)]
    volumes = [n.area for n in nuclei_props]
    print(f"{core} nucleus volume stats: min={min(volumes)}, max={max(volumes)}, median={np.median(volumes)}")
    sizes = [n.area for n in nuclei_props]
    print(f"{core} nucleus size stats: min={min(sizes)}, max={max(sizes)}, median={np.median(sizes)}")


    # Loop through hub objects
    for label_id in range(1, hub_count + 1):
        hub_mask = (hub_labels == label_id).astype(np.uint8)
        hub_mets = extract_hub_metrics(hub_labels, label_id)
        if hub_mets is None:
            continue

        nuclei_in_hub = 0
        for nuc in nuclei_props:
            coords = nuc.coords  # all voxel coordinates of this nucleus
            if np.any(hub_mask[tuple(coords.T)] > 0):
                nuclei_in_hub += 1

                print(f"Nucleus at ({z},{y},{x}) is inside hub {label_id}")
        print(f"{core} hub {label_id}: nuclei total = {len(nuclei_props)}, nuclei in hub = {nuclei_in_hub}")

        if nuclei_in_hub > 0:
            row = {
                "file": base,
                "core": core,
                "genotype": genotype,
                "hub_label": label_id,
                "nuclei_total": len(nuclei_props),
                "nuclei_in_hub": nuclei_in_hub,
            }
            row.update(hub_mets)
            rows.append(row)

# Save CSV
df = pd.DataFrame(rows)
df.to_csv(os.path.join(output_dir, "hub_analysis.csv"), index=False)
print(f"Done. Saved {len(rows)} rows to CSV.")

