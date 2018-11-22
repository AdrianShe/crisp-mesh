# crisp-mesh	
Crisp Meshes for Implicit Functions
Torus example copied from https://github.com/gaschler/bounding-mesh/blob/master/examples/torus/torus_decimated.off 
Letter C created from AutoDesk TinkerCAD

# How to use marching_cubes_root_finding

Copy the .h and .cpp files to the location where the original marching_cubes implementation exists under IGL.

Import using marching_cubes_root_finding.h

API spec:

Example -

igl::copyleft::marching_cubes_root_finding(double radius_of_influence, double cutoff, Point cloud input P, Grid points x, Grid sizes nx, ny, nz, Function Type based on different papers: igl::copyleft::LOCAL_IMPLICIT_FUNCTION_DEFAULT, Returns Vertices: V, Returns Faces: F);


Leave Cutoff above to some <0 value to automatically pick up based on the heuristics mentioned in our HW writeup. That is I calculate the field value at all the point cloud points using the radius of influence, and average them to get the cutoff.

If a cutoff > 0 is specified, then that will be used.

You may include the original marching cubes and this in the same cpp file.

