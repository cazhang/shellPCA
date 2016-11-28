# shellPCA: code of ShellPCA paper
The code implements a simple demo to show how to use ShellPCA given a bunch of meshes (.ply).

# Dependency
The code is implemented based on quocMesh library http://numod.ins.uni-bonn.de/software/quocmesh/
To run the demo, you need to follow instructions there first. 
Once installation is done, the script could be added and compiled.

# Averaging
The spheres demo uses linear average although an elastic average is expected to perform better
You can find the script to compute elastic average under quocMesh repositories easily

# Usage
1. determine the number of data to use;
2. provide bendWeight and memWeight, default values are 1.0;
3. run and visualize eigen modes

# Reference
@inproceedings{ZhHeRu15,
  author = {Zhang, Chao and Heeren, Behrend and Rumpf, Martin and Smith, William},
  title = {Shell {PCA}: statistical shape modelling in shell space},
  booktitle = {Proc. of IEEE International Conference on Computer Vision},
  year = {2015},
  pdf = {http://numod.ins.uni-bonn.de/research/papers/public/ZhHeRu15.pdf 1}
}

