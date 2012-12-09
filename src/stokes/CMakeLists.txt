set(HOME stokes)
libs(instatstokes2phase)
libs_ser(integrTime stokesCoeff)

target_link_libraries_par(stokes-instatstokes2phase num-discretize)

exec_ser(sdrops geom-boundary geom-builder geom-deformation stokes-sdrops geom-simplex geom-multigrid misc-scopetimer misc-progressaccu misc-params num-unknowns geom-topo num-fe misc-problem num-interfacePatch misc-utils num-discretize out-output geom-principallattice geom-reftetracut)

exec_ser(sdropsP2 geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid misc-scopetimer misc-progressaccu misc-params num-unknowns geom-topo num-fe misc-problem num-interfacePatch misc-utils out-output num-discretize num-MGsolver misc-params out-ensightOut out-vtkOut stokes-integrTime geom-geomselect misc-funcmap misc-vectorFunctions stokes-stokesCoeff geom-principallattice geom-reftetracut num-oseenprecond)

exec_ser(errorestimator stokes-errorestimator geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid misc-scopetimer misc-progressaccu misc-params num-unknowns geom-topo num-fe misc-problem  num-interfacePatch misc-utils num-discretize num-MGsolver misc-params stokes-integrTime  geom-geomselect out-output misc-funcmap misc-vectorFunctions  stokes-stokesCoeff geom-principallattice geom-reftetracut num-oseenprecond)
