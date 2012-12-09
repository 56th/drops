set(HOME misc)

libs(funcmap params problem scalarFunctions utils vectorFunctions scopetimer progressaccu)

target_link_libraries(misc-problem num-interfacePatch)
exec_ser(compensight misc-utils)
