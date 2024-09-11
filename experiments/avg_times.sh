#!/bin/bash

# TODO : could be sorted by modules/environment if need be:
pipeline_path=/software/projects/director2183/msok/blink_pipeline/gpu/pipeline/

ls avg_time??????.fits > fits_list_time
echo "sbatch ${pipeline_path}/experiments/avg_images.sh fits_list_time avg.fits rms.fits"
sbatch ${pipeline_path}/experiments/avg_images.sh fits_list_time avg.fits rms.fits
