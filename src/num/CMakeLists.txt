set(HOME num)

libs(bndData discretize fe hypre interfacePatch MGsolver oseenprecond quadrature renumber unknowns spacetime_geom spacetime_quad spacetime_map stokespardiso)

if (MKL_HOME STREQUAL "")
else (MKL_HOME STREQUAL "")
     target_link_libraries_ser(num-stokespardiso -fopenmp -Wl,--start-group ${MKL_HOME}/lib/intel64/libmkl_intel_ilp64.a ${MKL_HOME}/lib/intel64/libmkl_gnu_thread.a ${MKL_HOME}/lib/intel64/libmkl_core.a -Wl,--end-group -ldl -lpthread -lm) 
endif (MKL_HOME STREQUAL "")


target_link_libraries_par(num-unknowns misc-problem)
target_link_libraries_par(num-interfacePatch parallel-exchange num-fe)
target_link_libraries_par(num-quadrature num-discretize)

target_link_libraries(num-bndData misc-params)
target_link_libraries(num-discretize num-fe)
target_link_libraries(num-interfacePatch num-fe misc-utils)
target_link_libraries(num-interfacePatch num-discretize)
target_link_libraries(num-spacetime_quad num-spacetime_geom)

target_link_libraries_par(num-spacetime_quad num-spacetime_geom)

add_my_custom_targets(num)
