cmake_minimum_required(VERSION 3.16)

project(
    GasAccess
    VERSION 0.1.0
    DESCRIPTION "Geometric gas-accessibility analysis on Cartesian voxel grids"
    LANGUAGES C CXX
)

add_library(gasaccess
    src/accessibility_query.cpp
    src/atom_voxelizer.cpp
    src/exterior_classifier.cpp
    src/gasaccess_c.cpp
    src/gas_grid.cpp
)
add_library(GasAccess::gasaccess ALIAS gasaccess)

target_compile_features(gasaccess PUBLIC cxx_std_17)
target_include_directories(gasaccess
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
set_target_properties(gasaccess PROPERTIES
    CXX_EXTENSIONS OFF
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(gasaccess PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
        -Wsign-conversion
    )
elseif(MSVC)
    target_compile_options(gasaccess PRIVATE /W4)
endif()

include(CTest)

if(BUILD_TESTING)
    add_executable(gasaccess_grid_tests
        tests/test_gas_grid.cpp
    )
    target_link_libraries(gasaccess_grid_tests PRIVATE GasAccess::gasaccess)

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(gasaccess_grid_tests PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(gasaccess_grid_tests PRIVATE /W4)
    endif()

    add_test(NAME gasaccess_grid_tests COMMAND gasaccess_grid_tests)

    add_executable(gasaccess_atom_voxelizer_tests
        tests/test_atom_voxelizer.cpp
    )
    target_link_libraries(gasaccess_atom_voxelizer_tests PRIVATE GasAccess::gasaccess)

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(gasaccess_atom_voxelizer_tests PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(gasaccess_atom_voxelizer_tests PRIVATE /W4)
    endif()

    add_test(NAME gasaccess_atom_voxelizer_tests COMMAND gasaccess_atom_voxelizer_tests)

    add_executable(gasaccess_exterior_classifier_tests
        tests/test_exterior_classifier.cpp
    )
    target_link_libraries(gasaccess_exterior_classifier_tests
        PRIVATE GasAccess::gasaccess
    )

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(gasaccess_exterior_classifier_tests PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(gasaccess_exterior_classifier_tests PRIVATE /W4)
    endif()

    add_test(
        NAME gasaccess_exterior_classifier_tests
        COMMAND gasaccess_exterior_classifier_tests
    )

    add_executable(gasaccess_accessibility_query_tests
        tests/test_accessibility_query.cpp
    )
    target_link_libraries(gasaccess_accessibility_query_tests
        PRIVATE GasAccess::gasaccess
    )

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(gasaccess_accessibility_query_tests PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )

        if(UNIX AND NOT APPLE AND CMAKE_SIZEOF_VOID_P EQUAL 8)
            target_compile_definitions(gasaccess_accessibility_query_tests
                PRIVATE GASACCESS_TEST_WRAP_ALLOCATIONS=1
            )
            target_link_options(gasaccess_accessibility_query_tests PRIVATE
                -Wl,--wrap=_Znwm
                -Wl,--wrap=_Znam
            )
        endif()
    elseif(MSVC)
        target_compile_options(gasaccess_accessibility_query_tests PRIVATE /W4)
    endif()

    add_test(
        NAME gasaccess_accessibility_query_tests
        COMMAND gasaccess_accessibility_query_tests
    )

    add_executable(gasaccess_c_api_tests
        tests/test_c_api.c
    )
    target_link_libraries(gasaccess_c_api_tests PRIVATE GasAccess::gasaccess)
    set_target_properties(gasaccess_c_api_tests PROPERTIES
        C_STANDARD 99
        C_STANDARD_REQUIRED ON
        C_EXTENSIONS OFF
    )

    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(gasaccess_c_api_tests PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(gasaccess_c_api_tests PRIVATE /W4)
    endif()

    add_test(NAME gasaccess_c_api_tests COMMAND gasaccess_c_api_tests)
endif()
