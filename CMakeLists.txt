# save file with path to simulation sources 
set(SIM_SRC_DIR "sim_src_dir.txt") 
execute_process(COMMAND echo ${CMAKE_CURRENT_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}  
    OUTPUT_FILE "${CMAKE_CURRENT_BINARY_DIR}/${SIM_SRC_DIR}")

  add_simulation(crack_face_forces crack_face_forces.cc)
