set(HOME parallel)

libs_par(decompose exchange loadbal migrateunknowns parallel parmultigrid partime)

target_link_libraries_par(parallel-parmultigrid geom-multigrid)

add_my_custom_targets(parallel)
