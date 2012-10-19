#CMakeLists.txt for poisson package


#LIBRARIES

add_library(poisson-ale SHARED ale)
add_library(poisson-cdrdrops SHARED cdrdrops)
add_library(poisson-poissonCoeff SHARED poissonCoeff)
add_library(poisson-poisson SHARED poisson)
add_library(poisson-transport2phase SHARED transport2phase)

#EXECUTABLES
if(NOT MPI)

        add_executable(cdrdrops.exe cdrdrops)
        target_link_libraries(cdrdrops.exe poisson-poissonCoeff poisson-ale poisson-cdrdrops geom-boundary geom-builder out-ensightOut geom-simplex geom-multigrid num-unknowns geom-topo poisson-poisson misc-problem misc-utils out-output num-fe num-discretize num-interfacePatch geom-geomselect misc-params out-vtkOut misc-funcmap misc-scalarFunctions misc-vectorFunctions)

        #add_executable(ipdrops.exe ipdrops)
        #target_link_libraries(ipdrops.exe poisson-poissonCoeff geom-boundary geom-builder geom-simplex geom-multigrid num-unknowns geom-topo num-interfacePatch poisson-poisson misc-problem out-output misc-utils)

        #add_executable(ipfilm.exe ipfilm)
        #target_link_libraries(ipfilm.exe poisson-poissonCoeff geom-boundary geom-builder num-interfacePatch geom-simplex geom-multigrid num-unknowns geom-topo poisson-poisson misc-problem out-output misc-utils)

endif(NOT MPI)
