set(HOME num)
libs(bndData discretize fe hypre interfacePatch MGsolver oseenprecond quadrature renumber unknowns)

target_link_libraries_par(num-unknowns misc-problem)
target_link_libraries_par(num-interfacePatch parallel-exchange num-fe)
target_link_libraries_par(num-quadrature num-discretize)

target_link_libraries(num-bndData misc-params)
target_link_libraries(num-discretize num-fe)
target_link_libraries(num-interfacePatch num-fe misc-utils)
target_link_libraries(num-interfacePatch num-discretize)
