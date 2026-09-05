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
    grid_spec.spacing = (ga_grid_spacing){1.0, 1.0, 1.0};
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
    ga_voxel_id state_count = 0;
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
    REQUIRE(ga_grid_get_state_count(
        grid,
        GA_GAS_STATE_SOLID,
        &state_count) == GA_STATUS_SUCCESS);
    REQUIRE(state_count == summary.solid_count);
    REQUIRE(ga_grid_get_state_count(
        grid,
        GA_GAS_STATE_OUTSIDE_ACCESSIBLE,
        &state_count) == GA_STATUS_SUCCESS);
    REQUIRE(state_count == summary.outside_accessible_count);
    REQUIRE(ga_grid_get_state_count(
        grid,
        GA_GAS_STATE_CLOSED_VOID,
        &state_count) == GA_STATUS_SUCCESS);
    REQUIRE(state_count == summary.closed_void_count);

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

static void test_c_deposition_update(void)
{
    ga_grid_spec grid_spec = make_grid_spec(3, 1, 1);
    ga_voxel_coord explicit_source = {1, 0, 0};
    ga_grid* grid = NULL;
    ga_atom atom;
    ga_update_result* update_result = NULL;
    ga_deposition_update_summary update_summary;
    const ga_voxel_id* changed_voxel_ids = NULL;
    size_t changed_voxel_count = 0;

    grid_spec.periodic.x = 1;
    grid_spec.explicit_source_voxels = &explicit_source;
    grid_spec.explicit_source_count = 1;
    atom.position.x = 1.5;
    atom.position.y = 0.5;
    atom.position.z = 0.5;
    atom.radius = 0.0;

    REQUIRE(ga_grid_create(&grid_spec, &grid) == GA_STATUS_SUCCESS);
    REQUIRE(ga_classify_exterior(grid, &(ga_classification_summary){0})
        == GA_STATUS_SUCCESS);
    REQUIRE(ga_apply_deposition(
        grid,
        &atom,
        1,
        0.0,
        &update_result) == GA_STATUS_SUCCESS);
    REQUIRE(update_result != NULL);
    REQUIRE(ga_update_result_get_summary(update_result, &update_summary)
        == GA_STATUS_SUCCESS);
    REQUIRE(update_summary.newly_solid_count == 1);
    REQUIRE(update_summary.changed_voxel_count == 3);
    REQUIRE(update_summary.repair_visited_voxel_count == 2);
    REQUIRE(update_summary.repair_closed_voxel_count == 2);
    REQUIRE(update_summary.full_reclassification_performed == 0);
    REQUIRE(update_summary.affected_region_repair_performed == 1);
    REQUIRE(update_summary.classification.solid_count == 1);
    REQUIRE(update_summary.classification.outside_accessible_count == 0);
    REQUIRE(update_summary.classification.closed_void_count == 2);
    REQUIRE(ga_update_result_get_changed_voxels(
        update_result,
        &changed_voxel_ids,
        &changed_voxel_count) == GA_STATUS_SUCCESS);
    REQUIRE(changed_voxel_count == 3);
    REQUIRE(changed_voxel_ids != NULL);
    REQUIRE(changed_voxel_ids[0] < changed_voxel_ids[1]);
    REQUIRE(changed_voxel_ids[1] < changed_voxel_ids[2]);

    ga_update_result_destroy(update_result);
    update_result = NULL;
    ga_grid_destroy(grid);

    grid = NULL;
    REQUIRE(ga_grid_create(&grid_spec, &grid) == GA_STATUS_SUCCESS);
    REQUIRE(ga_classify_exterior(grid, &(ga_classification_summary){0})
        == GA_STATUS_SUCCESS);
    REQUIRE(ga_apply_deposition_with_mode(
        grid,
        &atom,
        1,
        0.0,
        GA_CONNECTIVITY_REPAIR_FULL_RECLASSIFICATION,
        &update_result) == GA_STATUS_SUCCESS);
    REQUIRE(ga_update_result_get_summary(update_result, &update_summary)
        == GA_STATUS_SUCCESS);
    REQUIRE(update_summary.full_reclassification_performed == 1);
    REQUIRE(update_summary.affected_region_repair_performed == 0);
    REQUIRE(update_summary.repair_visited_voxel_count == 0);
    REQUIRE(update_summary.repair_closed_voxel_count == 0);

    ga_update_result_destroy(update_result);
    ga_update_result_destroy(NULL);
    ga_grid_destroy(grid);
}

static void test_c_reversible_updates(void)
{
    ga_grid_spec grid_spec = make_grid_spec(5, 1, 1);
    ga_grid* grid = NULL;
    ga_atom initial_atom = {{2.5, 0.5, 0.5}, 0.0};
    ga_atom moved_atom = {{1.5, 0.5, 0.5}, 0.0};
    ga_atom cancellation_atom = {{4.5, 0.5, 0.5}, 0.0};
    ga_voxel_id newly_solid_count = 0;
    ga_classification_summary classification;
    ga_atom_change_batch changes;
    ga_atom_change_result* change_result = NULL;
    ga_atom_change_update_summary summary;
    const ga_voxel_id* changed_voxel_ids = NULL;
    size_t changed_voxel_count = 0;

    grid_spec.reservoir_faces.x_low = 1;
    REQUIRE(ga_grid_create(&grid_spec, &grid) == GA_STATUS_SUCCESS);
    REQUIRE(ga_voxelize_atoms(
        grid,
        &initial_atom,
        1,
        0.0,
        &newly_solid_count) == GA_STATUS_SUCCESS);
    REQUIRE(newly_solid_count == 1);
    REQUIRE(ga_classify_exterior(grid, &classification) == GA_STATUS_SUCCESS);
    REQUIRE(classification.solid_count == 1);
    REQUIRE(classification.outside_accessible_count == 2);
    REQUIRE(classification.closed_void_count == 2);

    changes.added_atoms = (ga_atom_view){&moved_atom, 1};
    changes.removed_atoms = (ga_atom_view){&initial_atom, 1};
    REQUIRE(ga_apply_atom_changes(
        grid,
        &changes,
        0.0,
        &change_result) == GA_STATUS_SUCCESS);
    REQUIRE(change_result != NULL);
    REQUIRE(ga_atom_change_result_get_summary(change_result, &summary)
        == GA_STATUS_SUCCESS);
    REQUIRE(summary.blocker_count_changed_voxel_count == 2);
    REQUIRE(summary.newly_solid_count == 1);
    REQUIRE(summary.newly_gas_count == 1);
    REQUIRE(summary.changed_voxel_count == 2);
    REQUIRE(summary.repair_kind == GA_ACCESSIBILITY_REPAIR_MIXED);
    REQUIRE(summary.full_reclassification_performed == 0);
    REQUIRE(summary.classification.solid_count == 1);
    REQUIRE(summary.classification.outside_accessible_count == 1);
    REQUIRE(summary.classification.closed_void_count == 3);
    REQUIRE(ga_atom_change_result_get_changed_voxels(
        change_result,
        &changed_voxel_ids,
        &changed_voxel_count) == GA_STATUS_SUCCESS);
    REQUIRE(changed_voxel_count == 2);
    REQUIRE(changed_voxel_ids[0] == get_voxel_id(grid, 1, 0, 0));
    REQUIRE(changed_voxel_ids[1] == get_voxel_id(grid, 2, 0, 0));
    ga_atom_change_result_destroy(change_result);
    change_result = NULL;

    REQUIRE(ga_apply_desorption(
        grid,
        &moved_atom,
        1,
        0.0,
        &change_result) == GA_STATUS_SUCCESS);
    REQUIRE(ga_atom_change_result_get_summary(change_result, &summary)
        == GA_STATUS_SUCCESS);
    REQUIRE(summary.newly_solid_count == 0);
    REQUIRE(summary.newly_gas_count == 1);
    REQUIRE(summary.changed_voxel_count == 4);
    REQUIRE(summary.repair_kind == GA_ACCESSIBILITY_REPAIR_OPENING);
    REQUIRE(summary.opening_repair_performed == 1);
    REQUIRE(summary.repair_opened_voxel_count == 4);
    REQUIRE(summary.classification.solid_count == 0);
    REQUIRE(summary.classification.outside_accessible_count == 5);
    REQUIRE(summary.classification.closed_void_count == 0);
    ga_atom_change_result_destroy(change_result);
    change_result = NULL;

    changes.added_atoms = (ga_atom_view){&cancellation_atom, 1};
    changes.removed_atoms = (ga_atom_view){&cancellation_atom, 1};
    REQUIRE(ga_apply_atom_changes_with_mode(
        grid,
        &changes,
        0.0,
        GA_ATOM_CHANGE_REPAIR_FULL_RECLASSIFICATION,
        &change_result) == GA_STATUS_SUCCESS);
    REQUIRE(ga_atom_change_result_get_summary(change_result, &summary)
        == GA_STATUS_SUCCESS);
    REQUIRE(summary.blocker_count_changed_voxel_count == 0);
    REQUIRE(summary.changed_voxel_count == 0);
    REQUIRE(summary.repair_kind == GA_ACCESSIBILITY_REPAIR_NONE);
    REQUIRE(summary.full_reclassification_performed == 0);
    ga_atom_change_result_destroy(change_result);
    change_result = NULL;

    changes.added_atoms = (ga_atom_view){&initial_atom, 1};
    changes.removed_atoms = (ga_atom_view){NULL, 0};
    REQUIRE(ga_apply_atom_changes_with_mode(
        grid,
        &changes,
        0.0,
        GA_ATOM_CHANGE_REPAIR_FULL_RECLASSIFICATION,
        &change_result) == GA_STATUS_SUCCESS);
    REQUIRE(ga_atom_change_result_get_summary(change_result, &summary)
        == GA_STATUS_SUCCESS);
    REQUIRE(summary.newly_solid_count == 1);
    REQUIRE(summary.newly_gas_count == 0);
    REQUIRE(summary.repair_kind
        == GA_ACCESSIBILITY_REPAIR_FULL_RECLASSIFICATION);
    REQUIRE(summary.full_reclassification_performed == 1);
    REQUIRE(summary.classification.solid_count == 1);
    ga_atom_change_result_destroy(change_result);
    ga_atom_change_result_destroy(NULL);
    ga_grid_destroy(grid);
}

static void test_c_error_handling(void)
{
    ga_grid_spec invalid_spec = make_grid_spec(0, 1, 1);
    ga_grid* grid = NULL;
    ga_grid_spec valid_spec = make_grid_spec(1, 1, 1);
    ga_gas_state gas_state = GA_GAS_STATE_UNCLASSIFIED;
    ga_voxel_id state_count = 0;
    uint8_t is_accessible = 0;
    ga_atom atom = {{0.5, 0.5, 0.5}, 0.0};
    ga_update_result* update_result = NULL;
    ga_atom_change_result* change_result = NULL;
    ga_atom_change_batch atom_changes = {{&atom, 1}, {NULL, 0}};

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
    REQUIRE(ga_grid_get_state_count(grid, (ga_gas_state)99, &state_count)
        == GA_STATUS_INVALID_ARGUMENT);
    REQUIRE(ga_grid_get_state_count(grid, GA_GAS_STATE_SOLID, NULL)
        == GA_STATUS_INVALID_ARGUMENT);
    REQUIRE(ga_grid_get_state(grid, 1, &gas_state) == GA_STATUS_OUT_OF_RANGE);
    REQUIRE(ga_apply_deposition(
        grid,
        &atom,
        1,
        0.0,
        &update_result) == GA_STATUS_INVALID_ARGUMENT);
    REQUIRE(update_result == NULL);
    REQUIRE(ga_apply_atom_changes_with_mode(
        grid,
        &atom_changes,
        0.0,
        (ga_atom_change_repair_mode)99,
        &change_result) == GA_STATUS_INVALID_ARGUMENT);
    REQUIRE(change_result == NULL);
    REQUIRE(ga_apply_atom_changes(
        grid,
        NULL,
        0.0,
        &change_result) == GA_STATUS_INVALID_ARGUMENT);
    REQUIRE(change_result == NULL);
    REQUIRE(ga_apply_desorption(
        grid,
        NULL,
        1,
        0.0,
        &change_result) == GA_STATUS_INVALID_ARGUMENT);
    REQUIRE(change_result == NULL);
    REQUIRE(ga_apply_deposition_with_mode(
        grid,
        &atom,
        1,
        0.0,
        (ga_connectivity_repair_mode)99,
        &update_result) == GA_STATUS_INVALID_ARGUMENT);
    REQUIRE(update_result == NULL);
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

static void test_c_non_cubic_spacing(void)
{
    ga_grid_spec grid_spec = make_grid_spec(2, 2, 2);
    ga_grid* grid = NULL;
    ga_voxel_id outside_id;
    uint8_t is_accessible = 0;

    grid_spec.spacing = (ga_grid_spacing){0.5, 1.0, 2.0};
    REQUIRE(ga_grid_create(&grid_spec, &grid) == GA_STATUS_SUCCESS);
    outside_id = get_voxel_id(grid, 1, 1, 1);
    REQUIRE(ga_grid_set_state(
        grid,
        outside_id,
        GA_GAS_STATE_OUTSIDE_ACCESSIBLE) == GA_STATUS_SUCCESS);
    REQUIRE(ga_is_site_accessible(
        grid,
        (ga_point3){0.75, 1.5, 3.0},
        &is_accessible) == GA_STATUS_SUCCESS);
    REQUIRE(is_accessible == 1);

    ga_grid_destroy(grid);
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

    test_c_deposition_update();
    if (failure_count == 0) {
        printf("[PASS] C deposition update\n");
    }

    test_c_reversible_updates();
    if (failure_count == 0) {
        printf("[PASS] C reversible updates\n");
    }

    test_c_non_cubic_spacing();
    if (failure_count == 0) {
        printf("[PASS] C non-cubic spacing\n");
    }

    if (failure_count != 0) {
        fprintf(stderr, "%d C API test failure(s)\n", failure_count);
        return 1;
    }

    printf("6 C API test groups passed\n");
    return 0;
}
