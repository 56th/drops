set(HOME navstokes)
libs(instatnavstokes2phase)

exec_ser(nsdrops geom-boundary geom-builder geom-deformation misc-scopetimer misc-progressaccu geom-simplex geom-multigrid num-unknowns geom-topo num-fe num-interfacePatch misc-problem misc-utils out-output num-discretize geom-principallattice geom-reftetracut stokes-slipBndOnePhase)

exec_ser(nsdrops_begehung geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid num-unknowns geom-topo num-fe num-interfacePatch misc-problem misc-utils out-output num-discretize geom-principallattice geom-reftetracut stokes-slipBndOnePhase)

exec_ser(insdrops geom-boundary geom-builder geom-deformation geom-simplex geom-multigrid num-unknowns geom-topo num-fe num-interfacePatch misc-problem misc-utils out-output num-discretize out-ensightOut geom-principallattice geom-reftetracut stokes-slipBndOnePhase)

exec_ser(insadrops geom-boundary geom-builder geom-deformation misc-scopetimer misc-progressaccu geom-simplex geom-multigrid num-unknowns geom-topo num-fe num-interfacePatch misc-problem misc-utils out-output num-discretize out-ensightOut geom-principallattice geom-reftetracut stokes-slipBndOnePhase)

add_my_custom_targets(navstokes)
