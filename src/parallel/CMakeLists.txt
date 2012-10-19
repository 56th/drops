#CMakeLists.txt for parallel package


#LIBRARIES
if(MPI)
	#~ add_library(parallel-decompose SHARED decompose)
	#~ add_library(parallel-exchange SHARED exchange)
	#~ add_library(parallel-loadbal SHARED loadbal)
	#~ add_library(parallel-logger SHARED logger)
	#~ add_library(parallel-migrateunknowns SHARED migrateunknowns)
	#~ add_library(parallel-parallel SHARED parallel)
	#~ add_library(parallel-parmultigrid SHARED parmultigrid)
	#~ add_library(parallel-partime SHARED partime)
	#~ #add_library(parallel-partitioner SHARED partitioner)
	
	set(home parallel)
	set(libs decompose exchange loadbal logger migrateunknowns parallel parmultigrid partime partitioner)
	foreach(lib ${libs})
	    add_library(${home}-${lib} SHARED ${lib})
	    target_link_libraries(${home}-${lib} ${MPI_CXX_LIBRARIES})
	endforeach(lib)
endif(MPI)
#EXECUTABLES
