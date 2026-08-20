cmake_minimum_required(VERSION 3.16)

project(
    GasAccess
    VERSION 0.1.0
    DESCRIPTION "Geometric gas-accessibility analysis on Cartesian voxel grids"
    LANGUAGES C CXX
)

add_library(gasaccess
    src/accessibility_query.cpp
    src/affected_region_repair.cpp
    src/atom_voxelizer.cpp
    src/deposition_updater.cpp
    src/exterior_classifier.cpp
    src/gasaccess_c.cpp
    src/gas_grid.cpp
    src/local_topology_filter.cpp
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

option(GASACCESS_ENABLE_MPI "Build the optional MPI gas-grid layer" OFF)
set(GASACCESS_SPPARKS_SOURCE_DIR "" CACHE PATH
    "Optional SPPARKS src directory used for adapter compatibility compilation"
)

if(GASACCESS_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)

    add_library(gasaccess_mpi
        src/distributed_deposition_updater.cpp
        src/distributed_exterior_classifier.cpp
        src/mpi_gas_grid.cpp
    )
    add_library(GasAccess::gasaccess_mpi ALIAS gasaccess_mpi)
    target_compile_features(gasaccess_mpi PUBLIC cxx_std_17)
    target_include_directories(gasaccess_mpi
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>
    )
    target_link_libraries(gasaccess_mpi
        PUBLIC
            GasAccess::gasaccess
            MPI::MPI_CXX
    )
    set_target_properties(gasaccess_mpi PROPERTIES
        CXX_EXTENSIONS OFF
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
    )

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(gasaccess_mpi PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(gasaccess_mpi PRIVATE /W4)
    endif()

    add_library(gasaccess_mock_spparks STATIC
        tools/mock_spparks.cpp
        tools/mock_spparks_fixture.cpp
    )
    target_compile_features(gasaccess_mock_spparks PUBLIC cxx_std_17)
    target_include_directories(gasaccess_mock_spparks
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/tools
    )
    target_link_libraries(gasaccess_mock_spparks
        PUBLIC GasAccess::gasaccess_mpi
    )
    set_target_properties(gasaccess_mock_spparks PROPERTIES CXX_EXTENSIONS OFF)

    add_executable(gasaccess_spparks_mock_driver
        tools/spparks_mock_driver.cpp
    )
    target_link_libraries(gasaccess_spparks_mock_driver
        PRIVATE gasaccess_mock_spparks
    )

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(gasaccess_mock_spparks PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
        target_compile_options(gasaccess_spparks_mock_driver PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(gasaccess_mock_spparks PRIVATE /W4)
        target_compile_options(gasaccess_spparks_mock_driver PRIVATE /W4)
    endif()
endif()

add_executable(gasaccess_reference_driver
    tools/reference_driver.cpp
)
target_link_libraries(gasaccess_reference_driver PRIVATE GasAccess::gasaccess)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(gasaccess_reference_driver PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
        -Wsign-conversion
    )
elseif(MSVC)
    target_compile_options(gasaccess_reference_driver PRIVATE /W4)
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

    add_executable(gasaccess_kmc_query_integration_tests
        tests/test_kmc_query_integration.cpp
    )
    target_link_libraries(gasaccess_kmc_query_integration_tests
        PRIVATE GasAccess::gasaccess
    )

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(gasaccess_kmc_query_integration_tests PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(gasaccess_kmc_query_integration_tests PRIVATE /W4)
    endif()

    add_test(
        NAME gasaccess_kmc_query_integration_tests
        COMMAND gasaccess_kmc_query_integration_tests
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

    add_executable(gasaccess_deposition_updater_tests
        tests/test_deposition_updater.cpp
    )
    target_link_libraries(gasaccess_deposition_updater_tests
        PRIVATE GasAccess::gasaccess
    )

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(gasaccess_deposition_updater_tests PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(gasaccess_deposition_updater_tests PRIVATE /W4)
    endif()

    add_test(
        NAME gasaccess_deposition_updater_tests
        COMMAND gasaccess_deposition_updater_tests
    )

    add_executable(gasaccess_local_topology_filter_tests
        tests/test_local_topology_filter.cpp
    )
    target_link_libraries(gasaccess_local_topology_filter_tests
        PRIVATE GasAccess::gasaccess
    )

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(gasaccess_local_topology_filter_tests PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(gasaccess_local_topology_filter_tests PRIVATE /W4)
    endif()

    add_test(
        NAME gasaccess_local_topology_filter_tests
        COMMAND gasaccess_local_topology_filter_tests
    )

    add_executable(gasaccess_affected_region_repair_tests
        tests/test_affected_region_repair.cpp
    )
    target_link_libraries(gasaccess_affected_region_repair_tests
        PRIVATE GasAccess::gasaccess
    )

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(gasaccess_affected_region_repair_tests PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
    elseif(MSVC)
        target_compile_options(gasaccess_affected_region_repair_tests PRIVATE /W4)
    endif()

    add_test(
        NAME gasaccess_affected_region_repair_tests
        COMMAND gasaccess_affected_region_repair_tests
    )

    if(GASACCESS_ENABLE_MPI)
        add_executable(gasaccess_mpi_grid_tests
            tests/test_mpi_gas_grid.cpp
        )
        target_link_libraries(gasaccess_mpi_grid_tests
            PRIVATE GasAccess::gasaccess_mpi
        )

        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(gasaccess_mpi_grid_tests PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
            )
        elseif(MSVC)
            target_compile_options(gasaccess_mpi_grid_tests PRIVATE /W4)
        endif()

        foreach(mpi_rank_count IN ITEMS 1 2 4)
            add_test(
                NAME gasaccess_mpi_grid_tests_${mpi_rank_count}_ranks
                COMMAND
                    ${MPIEXEC_EXECUTABLE}
                    ${MPIEXEC_NUMPROC_FLAG} ${mpi_rank_count}
                    ${MPIEXEC_PREFLAGS}
                    $<TARGET_FILE:gasaccess_mpi_grid_tests>
                    ${MPIEXEC_POSTFLAGS}
            )
            set_tests_properties(
                gasaccess_mpi_grid_tests_${mpi_rank_count}_ranks
                PROPERTIES
                    LABELS "mpi"
                    PROCESSORS ${mpi_rank_count}
                    TIMEOUT 60
            )
        endforeach()

        add_executable(gasaccess_spparks_mock_integration_tests
            tests/test_spparks_mock_integration.cpp
        )
        target_link_libraries(gasaccess_spparks_mock_integration_tests
            PRIVATE gasaccess_mock_spparks
        )

        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(gasaccess_spparks_mock_integration_tests PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
            )
        elseif(MSVC)
            target_compile_options(gasaccess_spparks_mock_integration_tests
                PRIVATE /W4
            )
        endif()

        foreach(mpi_rank_count IN ITEMS 1 2 4)
            add_test(
                NAME gasaccess_spparks_mock_tests_${mpi_rank_count}_ranks
                COMMAND
                    ${MPIEXEC_EXECUTABLE}
                    ${MPIEXEC_NUMPROC_FLAG} ${mpi_rank_count}
                    ${MPIEXEC_PREFLAGS}
                    $<TARGET_FILE:gasaccess_spparks_mock_integration_tests>
                    ${MPIEXEC_POSTFLAGS}
            )
            set_tests_properties(
                gasaccess_spparks_mock_tests_${mpi_rank_count}_ranks
                PROPERTIES
                    LABELS "mpi;spparks"
                    PROCESSORS ${mpi_rank_count}
                    TIMEOUT 120
            )
        endforeach()

        foreach(mpi_rank_count IN ITEMS 1 2 4)
            add_test(
                NAME gasaccess_spparks_million_${mpi_rank_count}_ranks
                COMMAND
                    ${MPIEXEC_EXECUTABLE}
                    ${MPIEXEC_NUMPROC_FLAG} ${mpi_rank_count}
                    ${MPIEXEC_PREFLAGS}
                    $<TARGET_FILE:gasaccess_spparks_mock_driver>
                    --scenario million-slab
                    ${MPIEXEC_POSTFLAGS}
            )
            set_tests_properties(
                gasaccess_spparks_million_${mpi_rank_count}_ranks
                PROPERTIES
                    LABELS "benchmark;spparks"
                    PROCESSORS ${mpi_rank_count}
                    TIMEOUT 180
                    PASS_REGULAR_EXPRESSION "acceptance=pass"
            )
        endforeach()

        if(GASACCESS_SPPARKS_SOURCE_DIR
           AND EXISTS "${GASACCESS_SPPARKS_SOURCE_DIR}/domain.h"
           AND EXISTS "${GASACCESS_SPPARKS_SOURCE_DIR}/app.h")
            add_executable(gasaccess_spparks_header_compatibility
                tests/test_spparks_header_compatibility.cpp
            )
            target_include_directories(gasaccess_spparks_header_compatibility
                SYSTEM PRIVATE ${GASACCESS_SPPARKS_SOURCE_DIR}
            )
            target_link_libraries(gasaccess_spparks_header_compatibility
                PRIVATE GasAccess::gasaccess_mpi
            )
            if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
                target_compile_options(gasaccess_spparks_header_compatibility PRIVATE
                    -Wall
                    -Wextra
                    -Wpedantic
                    -Wconversion
                    -Wsign-conversion
                )
            elseif(MSVC)
                target_compile_options(gasaccess_spparks_header_compatibility
                    PRIVATE /W4
                )
            endif()
            add_test(
                NAME gasaccess_spparks_header_compatibility
                COMMAND gasaccess_spparks_header_compatibility
            )
            set_tests_properties(gasaccess_spparks_header_compatibility
                PROPERTIES LABELS "spparks"
            )
        endif()

        add_executable(gasaccess_distributed_updater_tests
            tests/test_distributed_deposition_updater.cpp
        )
        target_link_libraries(gasaccess_distributed_updater_tests
            PRIVATE GasAccess::gasaccess_mpi
        )

        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(gasaccess_distributed_updater_tests PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
            )
        elseif(MSVC)
            target_compile_options(gasaccess_distributed_updater_tests
                PRIVATE /W4
            )
        endif()

        foreach(mpi_rank_count IN ITEMS 1 2 4)
            add_test(
                NAME gasaccess_distributed_updater_tests_${mpi_rank_count}_ranks
                COMMAND
                    ${MPIEXEC_EXECUTABLE}
                    ${MPIEXEC_NUMPROC_FLAG} ${mpi_rank_count}
                    ${MPIEXEC_PREFLAGS}
                    $<TARGET_FILE:gasaccess_distributed_updater_tests>
                    ${MPIEXEC_POSTFLAGS}
            )
            set_tests_properties(
                gasaccess_distributed_updater_tests_${mpi_rank_count}_ranks
                PROPERTIES
                    LABELS "mpi"
                    PROCESSORS ${mpi_rank_count}
                    TIMEOUT 120
            )
        endforeach()

        add_executable(gasaccess_distributed_classifier_tests
            tests/test_distributed_exterior_classifier.cpp
        )
        target_link_libraries(gasaccess_distributed_classifier_tests
            PRIVATE GasAccess::gasaccess_mpi
        )

        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(gasaccess_distributed_classifier_tests PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
            )
        elseif(MSVC)
            target_compile_options(gasaccess_distributed_classifier_tests
                PRIVATE /W4
            )
        endif()

        foreach(mpi_rank_count IN ITEMS 1 2 4)
            add_test(
                NAME gasaccess_distributed_classifier_tests_${mpi_rank_count}_ranks
                COMMAND
                    ${MPIEXEC_EXECUTABLE}
                    ${MPIEXEC_NUMPROC_FLAG} ${mpi_rank_count}
                    ${MPIEXEC_PREFLAGS}
                    $<TARGET_FILE:gasaccess_distributed_classifier_tests>
                    ${MPIEXEC_POSTFLAGS}
            )
            set_tests_properties(
                gasaccess_distributed_classifier_tests_${mpi_rank_count}_ranks
                PROPERTIES
                    LABELS "mpi"
                    PROCESSORS ${mpi_rank_count}
                    TIMEOUT 60
            )
        endforeach()
    endif()

    add_test(
        NAME gasaccess_reference_open_trench
        COMMAND gasaccess_reference_driver
            --scenario open-trench
            --nx 12 --ny 8 --nz 10
            --query-count 2000
            --update-count 2
    )
    add_test(
        NAME gasaccess_reference_sealed_trench
        COMMAND gasaccess_reference_driver
            --scenario sealed-trench
            --nx 12 --ny 8 --nz 10
            --query-count 2000
            --update-count 2
    )
    add_test(
        NAME gasaccess_reference_sealed_trench_periodic_z
        COMMAND gasaccess_reference_driver
            --scenario sealed-trench
            --nx 12 --ny 8 --nz 10
            --query-count 2000
            --update-count 2
            --periodic-z
    )
    add_test(
        NAME gasaccess_reference_pinch_off
        COMMAND gasaccess_reference_driver
            --scenario pinch-off
            --nx 12 --ny 8 --nz 10
            --query-count 2000
            --update-count 16
    )
    add_test(
        NAME gasaccess_reference_bulk_small
        COMMAND gasaccess_reference_driver
            --scenario bulk
            --nx 12 --ny 8 --nz 10
            --atom-count 200
            --query-count 2000
            --update-count 2
            --seed 12345
    )
    add_test(
        NAME gasaccess_reference_bulk_periodic
        COMMAND gasaccess_reference_driver
            --scenario bulk
            --nx 32 --ny 16 --nz 24
            --atom-count 5000
            --query-count 10000
            --update-count 2
            --seed 12345
            --periodic-x --periodic-y --periodic-z
    )
    set_tests_properties(
        gasaccess_reference_open_trench
        gasaccess_reference_sealed_trench
        gasaccess_reference_sealed_trench_periodic_z
        gasaccess_reference_pinch_off
        gasaccess_reference_bulk_small
        gasaccess_reference_bulk_periodic
        PROPERTIES LABELS "reference"
    )
    set_tests_properties(gasaccess_reference_open_trench PROPERTIES
        PASS_REGULAR_EXPRESSION
            "initial_state_checksum=0x438e38acbd5c56b5"
    )
    set_tests_properties(gasaccess_reference_sealed_trench PROPERTIES
        PASS_REGULAR_EXPRESSION
            "initial_state_checksum=0xd3d951011d1548c5"
    )
    set_tests_properties(
        gasaccess_reference_sealed_trench_periodic_z PROPERTIES
        PASS_REGULAR_EXPRESSION
            "initial_state_checksum=0xf28e91515bcd9585"
    )
    set_tests_properties(gasaccess_reference_pinch_off PROPERTIES
        PASS_REGULAR_EXPRESSION
            "final_state_checksum=0xd3d951011d1548c5"
    )
    set_tests_properties(gasaccess_reference_bulk_small PROPERTIES
        PASS_REGULAR_EXPRESSION
            "initial_state_checksum=0x4800ebb6f44af30c"
    )
    set_tests_properties(gasaccess_reference_bulk_periodic PROPERTIES
        PASS_REGULAR_EXPRESSION
            "initial_state_checksum=0x8bb2ec06c6c68d92"
    )
endif()
