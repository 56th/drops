set(HOME parallel)

libs_par(decompose exchange loadbal migrateunknowns parallel parmultigrid partime)

target_link_libraries_par(parallel-parmultigrid geom-multigrid)
