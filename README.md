# Surface Generator

This repository consists of CAE files, python scripts used to generate abaqus input files for the paper "Micro-scale deterministic asperity contact FEM simulation"

The folder "code" contains the python classes and methods, which can read point cloud data of measured surface from any standard optical microscope and convert it into a meshed abaqus _part_. Additionally, they can also be used to create artifical surfaces like flat plane or a hemi-spherical asperity of specified radius. An example python script _example.py_ in the same folder shows how these classes and methods could be used in a program.

The folder "paper" contains all the abaqus input files used as per section.