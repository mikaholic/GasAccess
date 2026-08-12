#include "gasaccess/gasaccess_c.h"

#include <stdio.h>
#include <string.h>

static int failure_count = 0;

#define REQUIRE(condition)                                                        \
    do {                                                                          \
        if (!(condition)) {                                                       \
            fprintf(                                                              \
                stderr,                                                           \
                "%s:%d: check failed: %s\n",                                    \
                __FILE__,                                                         \
                __LINE__,                                                         \
                #condition);                                                      \
            ++failure_count;                                                      \
            return;                                                               \
        }                                                                         \
    } while (0)

static ga_grid_spec make_grid_spec(uint64_t x, uint64_t y, uint64_t z)
{
    ga_grid_spec grid_spec;
    memset(&grid_spec, 0, sizeof(grid_spec));
    grid_spec.spacing = 1.0;
    grid_spec.dimensions.x = x;
    grid_spec.dimensions.y = y;
    grid_spec.dimensions.z = z;
    return grid_spec;
}

static ga_voxel_id get_voxel_id(ga_grid* grid, int64_t x, int64_t y, int64_t z)
{
    ga_voxel_coord voxel_coord;
    ga_voxel_id voxel_id = 0;
    voxel_coord.x = x;
    voxel_coord.y = y;
    voxel_coord.z = z;
    if (ga_grid_get_voxel_id(grid, voxel_coord, &voxel_id) != GA_STATUS_SUCCESS) {
        fprintf(stderr, "could not get voxel id: %s\n", ga_last_error_message());
        ++failure_count;
    }
    return voxel_id;
}

static void test_c_grid_classification_and_queries(void)
{
    ga_grid_spec grid_spec = make_grid_spec(3, 1, 3);
    ga_grid* grid = NULL;
    ga_voxel_id solid_ids[5];
    ga_voxel_id target_id;
    ga_voxel_id source_id;
    ga_classification_summary summary;
    ga_gas_state gas_state = GA_GAS_STATE_UNCLASSIFIED;
    uint8_t is_accessible = 0;
    size_t index;

    grid_spec.reservoir_faces.z_high = 1;
    REQUIRE(ga_grid_create(&grid_spec, &grid) == GA_STATUS_SUCCESS);
    REQUIRE(grid != NULL);

    solid_ids[0] = get_voxel_id(grid, 0, 0, 0);
    solid_ids[1] = get_voxel_id(grid, 2, 0, 0);
    solid_ids[2] = get_voxel_id(grid, 0, 0, 1);
    solid_ids[3] = get_voxel_id(grid, 1, 0, 1);
    solid_ids[4] = get_voxel_id(grid, 2, 0, 1);
    REQUIRE(failure_count == 0);
    for (index = 0; index < 5; ++index) {
        REQUIRE(ga_grid_set_state(
            grid,
            solid_ids[index],
            GA_GAS_STATE_SOLID) == GA_STATUS_SUCCESS);
    }

    REQUIRE(ga_classify_exterior(grid, &summary) == GA_STATUS_SUCCESS);
    REQUIRE(summary.solid_count == 5);
    REQUIRE(summary.outside_accessible_count == 3);
    REQUIRE(summary.closed_void_count == 1);

    target_id = get_voxel_id(grid, 1, 0, 0);
    source_id = get_voxel_id(grid, 1, 0, 2);
    REQUIRE(ga_grid_get_state(grid, target_id, &gas_state) == GA_STATUS_SUCCESS);
    REQUIRE(gas_state == GA_GAS_STATE_CLOSED_VOID);

    REQUIRE(ga_is_voxel_outside_accessible(
        grid,
        target_id,
        &is_accessible) == GA_STATUS_SUCCESS);
    REQUIRE(is_accessible == 0);
    REQUIRE(ga_is_voxel_outside_accessible(
        grid,
        source_id,
        &is_accessible) == GA_STATUS_SUCCESS);
    REQUIRE(is_accessible == 1);

    REQUIRE(ga_is_site_accessible(
        grid,
        (ga_point3){1.5, 0.5, 0.5},
        &is_accessible) == GA_STATUS_SUCCESS);
    REQUIRE(is_accessible == 0);
    REQUIRE(ga_is_site_accessible(
        grid,
        (ga_point3){1.5, 0.5, 2.5},
        &is_accessible) == GA_STATUS_SUCCESS);
    REQUIRE(is_accessible == 1);

    REQUIRE(ga_is_stencil_accessible(
        grid,
        &source_id,
        1,
        &is_accessible) == GA_STATUS_SUCCESS);
    REQUIRE(is_accessible == 1);

    ga_grid_destroy(grid);
}

static void test_c_atom_voxelization(void)
{
    ga_grid_spec grid_spec = make_grid_spec(3, 3, 3);
    ga_grid* grid = NULL;
    ga_atom atom;
    ga_voxel_id newly_solid_count = 0;
    ga_gas_state gas_state = GA_GAS_STATE_UNCLASSIFIED;

    atom.position.x = 1.5;
    atom.position.y = 1.5;
    atom.position.z = 1.5;
    atom.radius = 0.1;

    REQUIRE(ga_grid_create(&grid_spec, &grid) == GA_STATUS_SUCCESS);
    REQUIRE(ga_voxelize_atoms(
        grid,
        &atom,
        1,
        0.2,
        &newly_solid_count) == GA_STATUS_SUCCESS);
    REQUIRE(newly_solid_count == 1);
    REQUIRE(ga_grid_get_state_at(
        grid,
        (ga_voxel_coord){1, 1, 1},
        &gas_state) == GA_STATUS_SUCCESS);
    REQUIRE(gas_state == GA_GAS_STATE_SOLID);

    ga_grid_destroy(grid);
}

static void test_c_error_handling(void)
{
    ga_grid_spec invalid_spec = make_grid_spec(0, 1, 1);
    ga_grid* grid = NULL;
    ga_grid_spec valid_spec = make_grid_spec(1, 1, 1);
    ga_gas_state gas_state = GA_GAS_STATE_UNCLASSIFIED;
    uint8_t is_accessible = 0;

    REQUIRE(sizeof(ga_gas_state) == 1);
    REQUIRE(sizeof(ga_voxel_id) == 8);
    REQUIRE(GA_MAX_SITE_VOXEL_COUNT == 27);
    REQUIRE(ga_grid_create(&invalid_spec, &grid) == GA_STATUS_INVALID_ARGUMENT);
    REQUIRE(grid == NULL);
    REQUIRE(strlen(ga_last_error_message()) != 0);

    REQUIRE(ga_grid_create(&valid_spec, &grid) == GA_STATUS_SUCCESS);
    REQUIRE(strlen(ga_last_error_message()) == 0);
    REQUIRE(ga_grid_set_state(grid, 0, (ga_gas_state)99)
        == GA_STATUS_INVALID_ARGUMENT);
    REQUIRE(ga_grid_get_state(grid, 1, &gas_state) == GA_STATUS_OUT_OF_RANGE);
    REQUIRE(ga_is_stencil_accessible(
        grid,
        NULL,
        1,
        &is_accessible) == GA_STATUS_INVALID_ARGUMENT);
    REQUIRE(ga_is_stencil_accessible(
        grid,
        NULL,
        GA_MAX_SITE_VOXEL_COUNT + 1,
        &is_accessible) == GA_STATUS_INVALID_ARGUMENT);

    ga_grid_destroy(grid);
    ga_grid_destroy(NULL);
}

int main(void)
{
    test_c_grid_classification_and_queries();
    if (failure_count == 0) {
        printf("[PASS] C grid classification and queries\n");
    }

    test_c_atom_voxelization();
    if (failure_count == 0) {
        printf("[PASS] C atom voxelization\n");
    }

    test_c_error_handling();
    if (failure_count == 0) {
        printf("[PASS] C error handling\n");
    }

    if (failure_count != 0) {
        fprintf(stderr, "%d C API test failure(s)\n", failure_count);
        return 1;
    }

    printf("3 C API test groups passed\n");
    return 0;
}
