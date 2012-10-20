set(HOME stokes)
libs(instatstokes2phase)
libs_ser(errorestimator integrTime sdrops sdropsP2 stokesCoeff)

if(${MPI})
	target_link_libraries(stokes-instatstokes2phase num-discretize)
endif(${MPI})
