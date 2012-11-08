set(HOME misc)

libs(funcmap params problem scalarFunctions utils vectorFunctions)
libs_par(scopetimer)

target_link_libraries(misc-problem num-interfacePatch)
exec_ser(compensight misc-utils)
