import numpy as np
import matplotlib.pyplot as plt
from validphys.loader import Loader
from validphys.convolution import central_predictions
from validphys.covmats import covmat_from_systematics, sqrt_covmat

# Initialize Loader
l = Loader()

# Datasets
DATASETS = ["ATLAS_2JET_13TEV_DIF_MJJ-Y"]

# PDF
PDF = l.check_pdf(name="NNPDF40_nnlo_as_01180")

# Initialize counters and containers
ndata = 0
sqrtcm = {}

print("Loading datasets...")

for name in DATASETS:
    # Load dataset with new-style PineAPPL theory
    ds = l.check_dataset(name=name, theoryid=41000000)

    # Load commondata
    cd = ds.commondata.load()
    ndata += cd.ndata

    # Load cut mask
    mask = ds.cuts.load()
    print(f"Dataset {name}: mask shape = {mask.shape}")
    
    # Compute central predictions and apply mask
    pred = np.array(central_predictions(dataset=ds, pdf=PDF)).flatten()[mask]
    print(f"Predictions for {name}:\n", pred)

    # Compute covariance matrix and apply mask
    covmat = covmat_from_systematics(
        loaded_commondata_with_cuts=cd,
        dataset_input=None,
        use_weights_in_covmat=False
    )
    covmat = covmat[np.ix_(mask, mask)]
    sqrtcm[name] = sqrt_covmat(covariance_matrix=covmat)

print(f"Total data points: {ndata}")
print("Square-root covariance matrices computed for all datasets.")

