import nibabel as nib
import numpy as np
import os

from nimare.decode import discrete
from nimare.extract import fetch_neurosynth, fetch_neuroquery, download_abstracts
from nimare.io import convert_neurosynth_to_dataset

import nibabel as nib
import numpy as np
from nilearn.plotting import plot_roi

from nimare.dataset import Dataset
from nimare.decode import discrete
from nimare.utils import get_resource_path
from nimare.decode.continuous import CorrelationDecoder
from nimare.meta.cbma import mkda


# Load the parcellation
img = nib.load("stuff_for_revisions/schaefer_neurosynth/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.nii.gz")
atlas_ras = nib.as_closest_canonical(img) 
data = atlas_ras.get_fdata().astype(int)
affine = atlas_ras.affine

# Loop parcels
for i in range(1, 1001):
    mask = (data == i).astype(np.uint8)
    out_img = nib.Nifti1Image(mask, affine)
    nib.save(out_img, f"stuff_for_revisions/schaefer_neurosynth/parcel_{i}.nii.gz")
    

out_dir = os.path.abspath("stuff_for_revisions/neurosynth")
os.makedirs(out_dir, exist_ok=True)

files = nimare.extract.fetch_neurosynth(
    data_dir=out_dir,  # version 0.0.10 switched to data directory
    version="7",
    overwrite=False,
    source="abstract",
    vocab="LDA50",  # Note the difference here
)


neurosynth_db = files[0]

neurosynth_dset = nimare.io.convert_neurosynth_to_dataset(
    coordinates_file=neurosynth_db["coordinates"],
    metadata_file=neurosynth_db["metadata"],
    annotations_files=neurosynth_db["features"],
)

neurosynth_dset.save(os.path.join(out_dir, "neurosynth_dataset.pkl.gz"))
neurosynth_dset = download_abstracts(neurosynth_dset, "j.rittmo@gmail.com")
neurosynth_dset.save(os.path.join(out_dir, "neurosynth_dataset_with_abstracts.pkl.gz"))

####
# Neuroquery
####

files = fetch_neuroquery(
    data_dir=out_dir,
    version="1",
    overwrite=False,
    source="combined",
    vocab="neuroquery6308",
    type="tfidf",
)


neuroquery_db = files[0]

# Note that the conversion function says "neurosynth".
# This is just for backwards compatibility.
neuroquery_dset = convert_neurosynth_to_dataset(
    coordinates_file=neuroquery_db["coordinates"],
    metadata_file=neuroquery_db["metadata"],
    annotations_files=neuroquery_db["features"],
)


neuroquery_dset.save(os.path.join(out_dir, "neuroquery_dataset.pkl.gz"))


neuroquery_dset = download_abstracts(neuroquery_dset, "j.rittmo@gmail.com")
neuroquery_dset.save(os.path.join(out_dir, "neuroquery_dataset_with_abstracts.pkl.gz"))

neuroquery_dset = Dataset.load("neuroquery_dataset_with_abstracts.pkl.gz")

parcel_path = os.path.join(parcel_dir, f"parcel_1.nii.gz")
parc = nib.load(parcel_path)

decoder = ROIAssociationDecoder(masker=parc)
decoder.fit(neuroquery_dset)

decoded_df = decoder.transform()

x = decoded_df.sort_values(by="r", ascending=False)
x = x.reset_index().rename(columns={"index": "feature"})
x = x.filter(r>0.025)
#ids = neurosynth_dset.get_studies_by_mask(parc)

    

#dset = Dataset(os.path.join(get_resource_path(), "neurosynth_laird_studies.json"))
# dset.annotations.head(5)
# g1 = nib.load("stuff_for_revisions/gradient1.nii.gz")
# 
# decoder = CorrelationDecoder(
#     frequency_threshold=0.001,
#     meta_estimator=mkda.MKDAChi2,
#     target_image='z_desc-association',
# )
# 
# decoder.fit(dset)
# decoding_results = decoder.transform(g1)

# atlas_img = nib.load("stuff_for_revisions/schaefer_neurosynth/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.nii.gz")
# atlas_ras = nib.as_closest_canonical(atlas_img) 
# decoder = ROIAssociationDecoder(masker=atlas_ras)
# decoder.fit(neurosynth_dset)



parcel_dir = "stuff_for_revisions/schaefer_neurosynth"
all_results = []

for i in range(1, 1001):  # Schaefer 1000 parcels
    parcel_path = os.path.join(parcel_dir, f"parcel_{i}.nii.gz")
    parc = nib.load(parcel_path)

    # run decoder
    decoder = ROIAssociationDecoder(masker=parc)
    decoder.fit(neurosynth_dset)

    # extract decoded values
    decoded_df = decoder.transform()

    # add parcel ID column
    decoded_df["parcel"] = i
    
    # store results
    all_results.append(decoded_df)
    print(i)


results_df = pd.concat(all_results, ignore_index=True)
results_df.to_csv("schaefer1000_neurosynth_terms.csv", index=False)




parc = nib.load("stuff_for_revisions/schaefer_neurosynth/parcel_1.nii.gz")
decoder = discrete.ROIAssociationDecoder(parc)
decoder.fit(neurosynth_dset)

decoded_df = decoder.transform()
x = decoded_df.sort_values(by="r", ascending=False)

##### Parallell


import os
import nibabel as nib
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
from nimare.decode.discrete import ROIAssociationDecoder

# --- recommended to avoid CPU over-subscription (MKL/OpenMP) ---
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

parcel_dir = "stuff_for_revisions/schaefer_neurosynth"
parcel_indices = range(1, 1001)  # Schaefer-1000

def decode_one_parcel(i):
    """Decode a single parcel and return a DataFrame with results, or None on failure."""
    try:
        parcel_path = os.path.join(parcel_dir, f"parcel_{i}.nii.gz")
        parc = nib.load(parcel_path)

        # one decoder per process to avoid shared-state issues
        dec = ROIAssociationDecoder(masker=parc)
        dec.fit(neurosynth_dset)

        # NiMARE API differences: try transform() first, then .results fallback
        try:
            df = dec.transform()
        except Exception:
            # older versions store results on the decoder; convert to DataFrame if needed
            res = getattr(dec, "results", None)
            if res is None:
                return None
            if isinstance(res, pd.DataFrame):
                df = res.copy()
            else:
                # res could be dict-like
                df = pd.DataFrame(res).reset_index().rename(columns={"index": "feature"})

        df["parcel"] = i
        return df
    except Exception as e:
        # optional: log the error
        # print(f"Parcel {i} failed: {e}")
        return None

# choose #cores (leave 1 for the OS)
n_jobs = max(1, (os.cpu_count() or 2) - 1)

# run in parallel with a progress bar
dfs = Parallel(n_jobs=n_jobs, backend="loky")(
    delayed(decode_one_parcel)(i) for i in tqdm(parcel_indices, desc="Decoding parcels")
)


# drop failures and concatenate
all_results = [d for d in dfs if d is not None]


fixed_results = []
for df in all_results:
    if df.index.name or not df.reset_index().equals(df):
        df = df.reset_index().rename(columns={"index": "feature"})
    fixed_results.append(df)

results_df = pd.concat(fixed_results, ignore_index=True)

# save
results_df.to_csv("schaefer1000_neurosynth_terms.csv", index=False)

