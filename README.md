# DROPS: surface Navier-Stokes solver

Please see `BUILD.md` for building instructions.

Surface Navier-Stokes solver source: [src/surfnavierstokes/surfnavierstokes.cpp](https://github.com/56th/drops/blob/surfaceNSE_06/03/2021/src/surfnavierstokes/surfnavierstokes.cpp)

Example input files to run Kelvin-Hemholtz instability simulation on a sphere and torus are provided, see `sphere.json` and `torus.json`.

To run the torus simulation, move `torus.json` to `param/surfnavierstokes/No_Bnd_Condition.json`. Modify input parameters if necessary and run `[YOUR BUILD DIR]/surfnavierstokes/surfnavierstokes` executable.
Output files such as solver statistics, .vtu files for visualization in Paraview/Visit etc. will be written to `output/torus` directory (location can be changed in .json input file).
