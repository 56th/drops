set(HOME num)
libs(discretize fe hypre interfacePatch MGsolver oseenprecond quadrature renumber unknowns)

if(MPI)
	target_link_libraries(num-unknowns misc-problem)
	target_link_libraries(num-interfacePatch parallel-exchange num-fe)
	target_link_libraries(num-quadrature num-discretize)
endif(MPI)
