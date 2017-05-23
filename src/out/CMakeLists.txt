set(HOME out)

libs(ensightOut output vtkOut)

target_link_libraries(out-vtkOut misc-utils misc-base64)

add_my_custom_targets(out)
