project(count)

set(SOURCE_FILES
        main.cpp helper.h)

# if all packages could be found
if(Boost_FOUND)

    # compile binary
    add_executable(${PROJECT_NAME} ${SOURCE_FILES})

    # link against libraries
    target_link_libraries(${PROJECT_NAME} ${MARATHON_LIBRARIES})

    # install target
    install(TARGETS ${PROJECT_NAME} DESTINATION bin)

else(Boost_FOUND)

    message("-- Warning: Can not build subproject '${PROJECT_NAME}' due to missing dependencies.")

endif(Boost_FOUND)