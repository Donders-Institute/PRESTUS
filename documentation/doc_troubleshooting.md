# Troubleshooting
This document contains fixes to known usage issues. If you find and solve an issue not listed here, please add it with an elaborate explanation. If you stumble upon an error in the working of the pipeline, or have suggestions for improvement, please [open an issue on GitHub](https://github.com/Donders-Institute/PRESTUS/issues).

## Simulation errors
- An error occurs where the script can't find the `makeBowl.m` function.
    - This means that the path to k-Wave has not been defined properly.

## Segmentation errors
- When SimNIBS has not yet run the segmentation (found in the m2m folder), the pipeline will only run the segmentation.
    - Only after segmentation has been completed can the pipeline be restarted to run the acoustic simulation.
- When running bilateral simulations simultaneously when segmentation has not yet been completed, both will try to start a segmentation run.
    - One of these will produce an `segment_error` file in the `batch_job_logs` folder. Simply wait for the segmentation to complete and run both again.
- A common segmentation error is `ValueError: The qform and sform of do not match. Please run charm with the --forceqform option`. 
    - This can be solved by pasting `use_forceqform: 1` into your config file.
