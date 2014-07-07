set(HOME DiST)
libs_par(DiST mpistream remotedata)

target_link_libraries_par(DiST-remotedata DiST-DiST DiST-mpistream)
target_link_libraries_par(DiST-DiST geom-simplex num-unknowns)

add_my_custom_targets(DiST)
