set(HOME misc)

libs(funcmap params problem utils scopetimer progressaccu dynamicload)

target_link_libraries(misc-problem num-interfacePatch)

if(NOT WIN32)
    target_link_libraries(misc-dynamicload dl)
endif(NOT WIN32)

exec_ser(compensight misc-utils)

add_library(misc-scalarFunctions SHARED scalarFunctions)
add_library(misc-vectorFunctions SHARED vectorFunctions)
add_library(misc-csgFunctions SHARED csgFunctions)

if(NOT WIN32)
    target_link_libraries(misc-scalarFunctions)
    target_link_libraries(misc-vectorFunctions)
    target_link_libraries(misc-csgFunctions misc-params geom-csg)
endif(NOT WIN32)

add_my_custom_targets(misc)
