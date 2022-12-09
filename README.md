# voronoi-on-amr

The labeling the adaptive grid with ids of the nearest points using radially expanding bubble method. This method has O(N) complexity and can be used to make Voronoi meta-grid over very large AMR-grids. At each iteration only the cells on the bubbles surface are proceeded. So the number operation in theory can be as small as number of cells in underlying grid. The video shows 2d-slice of 3d-cells one side of which are touching xy-plane, note the process is full 3d.
