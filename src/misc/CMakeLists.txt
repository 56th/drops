set(HOME misc)

libs(funcmap params problem utils scopetimer progressaccu dynamicload)

target_link_libraries(misc-problem num-interfacePatch)

if(NOT WIN32)
    target_link_libraries(misc-dynamicload dl)
endif(NOT WIN32)

exec_ser(compensight misc-utils)

add_library(misc-scalarFunctions SHARED scalarFunctions)
add_library(misc-vectorFunctions SHARED vectorFunctions)

#remove suffix from shared libraries for inter-operability of osses.
#set_property(TARGET misc-scalarFunctions misc-vectorFunctions PROPERTY SUFFIX "")

if(NOT WIN32)
    target_link_libraries(misc-scalarFunctions)
    target_link_libraries(misc-vectorFunctions)
endif(NOT WIN32)
