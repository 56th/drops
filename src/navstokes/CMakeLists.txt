set(HOME navstokes)
libs(instatnavstokes2phase)

exec_ser(nsdrops geom-boundary geom-builder geom-deformation navstokes-nsdrops misc-scopetimer misc-progressaccu geom-simplex geom-multigrid num-unknowns geom-topo num-fe num-interfacePatch misc-problem misc-utils out-output num-discretize geom-principallattice geom-reftetracut)

exec_ser(nsdrops_begehung geom-boundary geom-builder geom-deformation navstokes-nsdrops_begehung geom-simplex geom-multigrid num-unknowns geom-topo num-fe num-interfacePatch misc-problem misc-utils out-output num-discretize geom-principallattice geom-reftetracut)

exec_ser(insdrops geom-boundary geom-builder geom-deformation navstokes-insdrops geom-simplex geom-multigrid num-unknowns geom-topo num-fe num-interfacePatch misc-problem misc-utils out-output num-discretize out-ensightOut geom-principallattice geom-reftetracut)

exec_ser(insadrops geom-boundary geom-builder geom-deformation navstokes-insadrops misc-scopetimer misc-progressaccu geom-simplex geom-multigrid num-unknowns geom-topo num-fe num-interfacePatch misc-problem misc-utils out-output num-discretize out-ensightOut geom-principallattice geom-reftetracut)
