=================================================================================================

"An Algorithm for Triangulating Multiple 3D Polygons"

Publication: Eurographics Symposium on Geometry Processing (SGP) 2013
Author: Ming Zou (mingzou@wustl.edu), Tao Ju, Nathan Carr
Version: 1.0 
Last update: June 28, 2013


I. Introduction==================================================================================

TMP (Triangulating Multiple Polygon) is based on the algorithm presented in paper "An Algorithm for Triangulating Multiple 3D Polygons" (SGP 2013). The algorithm reads in a set of 3D non-intersecting polygons (.curve file, see bellow) and generates a triangulation mesh (.obj file) bounded by those polygons. The shape of this triangulation can be controlled by specifying different metrics; the output surface is the optimal triangulation with a minimum metric cost. This algorithm is suitable for various surfacing application, like hole filling and lofting 3D sketches. Additionally, the algorithm can take in user defined boundary normals (.normal file, see bellow) to better control the shape of the output surface.

+ Choices of metric

There are currently 4 types of metrics
1. minimizing the total area of the mesh triangles
2. minimizing the sum of perimeter of each of the mesh triangle
3. minimizing the average of dihedral angle between each pair of adjacent mesh triangles 
4. minimizing the worst of dihedral angle between each pair of adjacent mesh triangles



II. How to run===================================================================================

[ Usage ]

TMP.exe <curveName> <useDT> <useMinSet> <areaWeight> <edgeWeight> <dihedralWeight> <useMinMaxDihedral> <saveObj> <useNormal>
-----------------------------------------------------------------------------------------------
ARGUMENT             |VALUE           |DESCRIPTION                                    |EXAMPLE
-----------------------------------------------------------------------------------------------
1) curveName         |string          |the name of the .curve file                    |monkey
2) useDT             |{0,1}           |search in Delaunay triangle space (1)          |  1
                     |                |or in all triangle set (0)                     |
3) useMinSet         |{0,1}           |use minimal set to speed up the algorithm (1)  |  1
                     |                |or not (0)                                     |
4) areaWeight        |rational number |use area metric                                |  0
5) edgeWeight        |rational number |use perimeter metric                           |  0
6) dihedralWeight    |rational number |use average dihedral metric                    |  1
7) useMinMaxDihedral |{0,1}           |minimizing the worst of dihedral angle between |  0
                     |                |each pair of adjacent mesh triangles (1)       |
                     |                |or use other metrics (0)                       |
8) saveObj           |{0,1}           |save surface (1) or not (0)                    |  1
9) useNormal         |{0,1}           |include in normal file (1) or not (0)          |  0
-----------------------------------------------------------------------------------------------

Notes: 
a. normal file should have the same name of the curve file. e.g. "monkey.curve" + "monkey.normal"
b. currently, the algorithm takes at most 6 polygons as input
c. search in all triangle set (useDT=0) could take a considerable time and memory

[ Example]

1. Read in monkey saddle curve (monkey.curve) to compute a triangulation (monkey.obj), and minimize the average of dihedral. Use Delaunay triangles; use minimal set to speed up; and do not use additional normals.
e.g. TMP.exe monkey 1 1 0 0 1 0 1 0 

2. Read in monkey saddle curve (monkey.curve) to compute a triangulation (monkey.obj), and minimize the average of dihedral. Use Delaunay triangles; use minimal set to speed up; and use additional normals (monkey.normal).
e.g. TMP.exe monkey 1 1 0 0 1 0 1 1



III. File Format=================================================================================

1. .curve file
-----------------------------------------------------------------------------------------------
FILE FORMAT                   |DESCRIPTION
-----------------------------------------------------------------------------------------------
# of polygons                 |current the number is limited to 1~6
# of total vertices           |total number of vertices of all input polygons
# of vertices in polygon_1    |number of vertices in polygon_1
v[x] v[y] v[z]                |the xyz coordinates of each vertex in polygon_1
                              |separated by blank, one vertex per line
# of vertices in polygon_2    |number of vertices in polygon_2
v[x] v[y] v[z]                |the xyz coordinates of each vertex in polygon_2
                              |separated by blank, one vertex per line
…                             |more polygon info
# of vertices in polygon_k    |number of vertices in polygon_k
v[x] v[y] v[z]                |the xyz coordinates of each vertex in polygon_k
                              |separated by blank, one vertex per line
-----------------------------------------------------------------------------------------------

2. .normal file
The same format as .curve, except that the coordinates are for the normals. 
Each normal is defined on an input edge. For example, the ith normal in a polygon is actually the normal defined on edge (v_i, v_i+1) in this polygon.


=================================================================================================

Tested on Window7 64bit.
External Lib: tetgen1.4.3; boost_1_53_0

---

Updated by Yotam Gingold to compile on macOS. Light modernization removed boost dependency. Compile with CMake:

```
mkdir build
cd build
cmake ..
make
```
