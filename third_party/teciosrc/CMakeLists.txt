if("${HAMERS_USE_TECIO}")
    set(BaseFlags "-DLINUX -DLINUX64 -DTP_PROJECT_USES_BOOST -DBOOST_ALL_NO_LIB -DMAKEARCHIVE -DNO_THIRD_PARTY_LIBS -DNO_ASSERTS -fmessage-length=0 -fPIC -fvisibility=hidden -w")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BaseFlags}")
    
    FILE(GLOB mainFiles "*.cpp" "*.h")
    LIST(REMOVE_ITEM mainFiles "${CMAKE_CURRENT_SOURCE_DIR}/szcombine.cpp")
    
    add_library(tecio STATIC ${mainFiles})
    
    add_executable(szcombine "szcombine.cpp")
    target_link_libraries(szcombine tecio pthread)
    
    # Had to add this so building tecio care package within loci-chem works...
    set_target_properties(tecio PROPERTIES LINKER_LANGUAGE CXX)
    set_target_properties(szcombine PROPERTIES LINKER_LANGUAGE CXX)
endif()