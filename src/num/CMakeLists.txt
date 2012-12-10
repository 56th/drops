set(HOME num)
libs(discretize fe hypre interfacePatch MGsolver oseenprecond quadrature renumber unknowns)

target_link_libraries_par(num-unknowns misc-problem)
target_link_libraries_par(num-interfacePatch parallel-exchange num-fe)
target_link_libraries_par(num-quadrature num-discretize)

target_link_libraries(num-discretize num-fe)
target_link_libraries(num-interfacePatch misc-problem num-fe misc-utils)
