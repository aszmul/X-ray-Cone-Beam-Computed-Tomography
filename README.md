# X-ray Cone Beam Computed Tomography

## THIS REPOSITORY IS ARCHIVED

X-ray Cone-beam Computed tomography (CBCT) is a technique for imaging cross-sections of an object using a series of X-ray measurements taken from
different angles around the object. CBCT uses a cone shaped X-ray beam.
The projection data acquired at the 2D detector array, which is then reconstructed into a 3D image volume using the FDK reconstruction
algorithm.

This repository contains computer program to simulate CBCT.


```
├── Astra-toolbox
│   ├── img
│   │   ├── output_10_1.png
│   │   ├── output_13_0.png
│   │   ├── output_6_0.png
│   │   └── output_8_0.png
│   ├── Phantom3DLibrary.dat
│   ├── Project.py
│   └── README.md
├── LICENSE
├── MATLAB
│   ├── main
│   │   ├── back_project.m
│   │   ├── phantom.m
│   │   ├── projections.m
│   │   ├── ramp_filter.m
│   │   ├── shepp_logan.m
│   │   └── ye_yu_wang.m
│   ├── parameter.mat
│   ├── xray_cone_beam_reconstruction_filtered.mlx
│   └── xray_cone_beam_reconstruction.mlx
├── Python
│   ├── back_projection.py
│   ├── phantom_features.py
│   ├── phatom_const.py
│   ├── projection.py
│   ├── __pycache__
│   └── ramp.py
└── README.md
```
