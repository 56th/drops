#CMakeLists.txt for stokes package


#LIBRARIES
add_library(stokes-instatstokes2phase SHARED instatstokes2phase)
add_library(stokes-integrTime SHARED integrTime)

target_link_libraries(stokes-instatstokes2phase num-discretize)

if(NOT MPI)
  add_library(stokes-sdrops SHARED sdrops)
endif(NOT MPI)

#EXECUTABLES
if(NOT MPI)
        add_library(stokes-stokesCoeff SHARED stokesCoeff)
        add_executable(sdrops.exe sdrops)
        target_link_libraries(sdrops.exe geom-boundary geom-builder stokes-sdrops geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem num-interfacePatch misc-utils num-discretize out-output geom-principallattice geom-reftetracut)
        
        add_library(stokes-sdropsP2 SHARED sdropsP2)
        add_executable(sdropsP2.exe sdropsP2)
        target_link_libraries(sdropsP2.exe stokes-sdropsP2 geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem num-interfacePatch misc-utils out-output num-discretize num-MGsolver misc-params out-ensightOut out-vtkOut stokes-integrTime geom-geomselect misc-funcmap misc-vectorFunctions stokes-stokesCoeff geom-principallattice geom-reftetracut num-oseenprecond)

        add_library(stokes-errorestimator SHARED errorestimator)
        add_executable(errorestimator.exe errorestimator)
        target_link_libraries(errorestimator.exe stokes-errorestimator geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-fe misc-problem num-interfacePatch misc-utils num-discretize num-MGsolver misc-params stokes-integrTime geom-geomselect out-output misc-funcmap misc-vectorFunctions stokes-stokesCoeff geom-principallattice geom-reftetracut num-oseenprecond)
    
endif(NOT MPI)
