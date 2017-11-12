** Prerequest

"LevelSetMethods" for image process (http://vision.ece.ucsb.edu/~sumengen/).
"InterX" to obtain intersection between curves and straight lines.

3D reconstruction

1. Open matlab and excute "image_processing.m". This will generate the boundary for each Z-stack frame. Make sure the loop upper bound to be the same as the image source, for example, set 35 for "wt_C001Z035T055.tif".
2. Excute "boundary_identification.m", this will differentiate inner and outer boundary and estimate the central axis.

**Please feel free to contact huyc@tsinghua.edu.cns
