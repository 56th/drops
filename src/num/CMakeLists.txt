#CMakeLists.txt for num package


#LIBRARIES

add_library(num-unknowns SHARED unknowns)
add_library(num-interfacePatch SHARED interfacePatch)
add_library(num-discretize SHARED discretize)
add_library(num-fe SHARED fe)
add_library(num-hypre SHARED hypre)
add_library(num-MGsolver SHARED MGsolver)
add_library(num-oseenprecond SHARED oseenprecond)
add_library(num-oseensolver SHARED oseensolver)
add_library(num-quadrature SHARED quadrature)
add_library(num-renumber SHARED renumber)

target_link_libraries(num-discretize num-fe)
target_link_libraries(num-interfacePatch num-fe)

#EXECUTABLES
