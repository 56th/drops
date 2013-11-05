set(HOME out)

libs(ensightOut output vtkOut)

target_link_libraries(out-vtkOut misc-utils)

add_my_custom_targets(out)
