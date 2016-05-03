set(CONFIGS
    #"1024, 6200, uint64_t"
	"8192, 124, uint64_t"
    #"1024, 124, uint64_t"
    #"1UL<<14, 620, uint64_t"
    #"1024, 62, uint64_t"
    #"1024, 37200, uint64_t"
    "1024, 60, uint32_t"
    #"1024, 3000, uint32_t"
    "512, 30, uint32_t"
	#"512, 14, uint16_t"
	#"256, 14, uint16_t"
    "8, 60, uint32_t"
)

find_package(CEREAL)
if(CEREAL_FOUND)
MESSAGE(STATUS "Found CEREAL: ${CEREAL_INCLUDE_DIR}")
include_directories(${CEREAL_INCLUDE_DIR})
endif()

#always double check during tests
add_definitions(-DCHECK_STRICTMOD)

foreach(CONFIG ${CONFIGS})
    string(REGEX REPLACE ", " "_" CONFIG_ID ${CONFIG})

    foreach(VERSION op func)

        file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/nfllib_demo${CONFIG_ID}_${VERSION}.cpp"
            "#define CONFIG ${CONFIG}\n#include \"${CMAKE_CURRENT_SOURCE_DIR}/nfllib_demo_main_${VERSION}.cpp\"\n")
        add_executable(nfllib_demo${CONFIG_ID}_${VERSION} ${CMAKE_CURRENT_BINARY_DIR}/nfllib_demo${CONFIG_ID}_${VERSION}.cpp)
        target_link_libraries(nfllib_demo${CONFIG_ID}_${VERSION} nfllib)
        add_test(run_nfllib_demo${CONFIG_ID}_${VERSION} nfllib_demo${CONFIG_ID}_${VERSION})

    endforeach()

    foreach(OP nfl_add nfl_sub nfl_mul nfl_eq nfl_neq nfl_stream poly_p poly_set poly_mpz)
        file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${OP}${CONFIG_ID}_${VERSION}.cpp"
            "#define CONFIG ${CONFIG}\n#include \"${CMAKE_CURRENT_SOURCE_DIR}/${OP}.cpp\"\n")
        add_executable(${OP}${CONFIG_ID}_${VERSION} ${CMAKE_CURRENT_BINARY_DIR}/${OP}${CONFIG_ID}_${VERSION}.cpp)
        target_link_libraries(${OP}${CONFIG_ID}_${VERSION} nfllib)
        add_test(run_${OP}${CONFIG_ID}_${VERSION} ${OP}${CONFIG_ID}_${VERSION})
    endforeach()

    if(CEREAL_FOUND)
        foreach(OP poly_serialize)
            file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${OP}${CONFIG_ID}_${VERSION}.cpp"
                "#define CONFIG ${CONFIG}\n#include \"${CMAKE_CURRENT_SOURCE_DIR}/${OP}.cpp\"\n")
            add_executable(${OP}${CONFIG_ID}_${VERSION} ${CMAKE_CURRENT_BINARY_DIR}/${OP}${CONFIG_ID}_${VERSION}.cpp)
            target_link_libraries(${OP}${CONFIG_ID}_${VERSION} nfllib)
            add_test(run_${OP}${CONFIG_ID}_${VERSION} ${OP}${CONFIG_ID}_${VERSION})
        endforeach() 
    endif()

endforeach()

add_executable(prng_demo prng_demo_main.cpp)
target_link_libraries(prng_demo nfllib)
add_test(run_prng_demo prng_demo)

add_executable(ntt_perfs ntt_perfs.cpp)
target_link_libraries(ntt_perfs nfllib)
