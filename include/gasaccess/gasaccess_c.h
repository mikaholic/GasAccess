#ifndef GASACCESS_GASACCESS_C_H
#define GASACCESS_GASACCESS_C_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ga_grid ga_grid;
typedef struct ga_update_result ga_update_result;

typedef uint64_t ga_voxel_id;
typedef uint8_t ga_gas_state;
typedef uint8_t ga_connectivity_repair_mode;
typedef int32_t ga_status;

#define GA_STATUS_SUCCESS ((ga_status)0)
#define GA_STATUS_INVALID_ARGUMENT ((ga_status)1)
#define GA_STATUS_OUT_OF_RANGE ((ga_status)2)
#define GA_STATUS_OVERFLOW ((ga_status)3)
#define GA_STATUS_ALLOCATION_FAILED ((ga_status)4)
#define GA_STATUS_INTERNAL_ERROR ((ga_status)5)

#define GA_GAS_STATE_UNCLASSIFIED ((ga_gas_state)0)
#define GA_GAS_STATE_SOLID ((ga_gas_state)1)
#define GA_GAS_STATE_OUTSIDE_ACCESSIBLE ((ga_gas_state)2)
#define GA_GAS_STATE_CLOSED_VOID ((ga_gas_state)3)

#define GA_CONNECTIVITY_REPAIR_AFFECTED_REGION \
    ((ga_connectivity_repair_mode)0)
#define GA_CONNECTIVITY_REPAIR_FULL_RECLASSIFICATION \
    ((ga_connectivity_repair_mode)1)

#define GA_MAX_SITE_VOXEL_COUNT ((size_t)27)

typedef struct ga_voxel_coord {
    int64_t x;
    int64_t y;
    int64_t z;
} ga_voxel_coord;

typedef struct ga_point3 {
    double x;
    double y;
    double z;
} ga_point3;

typedef struct ga_grid_dimensions {
    uint64_t x;
    uint64_t y;
    uint64_t z;
} ga_grid_dimensions;

typedef struct ga_periodic_axes {
    uint8_t x;
    uint8_t y;
    uint8_t z;
} ga_periodic_axes;

typedef struct ga_reservoir_faces {
    uint8_t x_low;
    uint8_t x_high;
    uint8_t y_low;
    uint8_t y_high;
    uint8_t z_low;
    uint8_t z_high;
} ga_reservoir_faces;

typedef struct ga_grid_spec {
    ga_point3 origin;
    double spacing;
    ga_grid_dimensions dimensions;
    ga_periodic_axes periodic;
    ga_reservoir_faces reservoir_faces;
    const ga_voxel_coord* explicit_source_voxels;
    size_t explicit_source_count;
} ga_grid_spec;

typedef struct ga_atom {
    ga_point3 position;
    double radius;
} ga_atom;

typedef struct ga_classification_summary {
    ga_voxel_id solid_count;
    ga_voxel_id outside_accessible_count;
    ga_voxel_id closed_void_count;
} ga_classification_summary;

typedef struct ga_deposition_update_summary {
    ga_voxel_id newly_solid_count;
    size_t changed_voxel_count;
    /* Repair traversal work and the number of gas voxels relabeled closed. */
    ga_voxel_id repair_visited_voxel_count;
    ga_voxel_id repair_closed_voxel_count;
    /* At most one of these method flags is nonzero for an update. */
    uint8_t full_reclassification_performed;
    uint8_t affected_region_repair_performed;
    ga_classification_summary classification;
} ga_deposition_update_summary;

/* The returned thread-local message remains valid until the next C API call
 * on the same thread. It is empty after a successful status-returning call. */
const char* ga_last_error_message(void);

/* ga_grid_create copies the grid specification and explicit source array. */
ga_status ga_grid_create(const ga_grid_spec* grid_spec, ga_grid** out_grid);
void ga_grid_destroy(ga_grid* grid);

ga_status ga_grid_get_voxel_count(const ga_grid* grid, ga_voxel_id* out_count);
ga_status ga_grid_get_voxel_id(
    const ga_grid* grid,
    ga_voxel_coord voxel_coord,
    ga_voxel_id* out_voxel_id);
ga_status ga_grid_get_state(
    const ga_grid* grid,
    ga_voxel_id voxel_id,
    ga_gas_state* out_gas_state);
ga_status ga_grid_get_state_at(
    const ga_grid* grid,
    ga_voxel_coord voxel_coord,
    ga_gas_state* out_gas_state);
ga_status ga_grid_get_state_count(
    const ga_grid* grid,
    ga_gas_state gas_state,
    ga_voxel_id* out_count);
ga_status ga_grid_set_state(
    ga_grid* grid,
    ga_voxel_id voxel_id,
    ga_gas_state gas_state);

ga_status ga_voxelize_atoms(
    ga_grid* grid,
    const ga_atom* atoms,
    size_t atom_count,
    double precursor_radius,
    ga_voxel_id* out_newly_solid_count);

ga_status ga_classify_exterior(
    ga_grid* grid,
    ga_classification_summary* out_summary);

/* The grid must already be fully classified. ga_apply_deposition copies the
 * small event batch and returns an owned result handle. */
ga_status ga_apply_deposition(
    ga_grid* grid,
    const ga_atom* deposited_atoms,
    size_t atom_count,
    double precursor_radius,
    ga_update_result** out_update_result);
/* Equivalent to ga_apply_deposition with an explicit incremental or full
 * reference repair mode. */
ga_status ga_apply_deposition_with_mode(
    ga_grid* grid,
    const ga_atom* deposited_atoms,
    size_t atom_count,
    double precursor_radius,
    ga_connectivity_repair_mode repair_mode,
    ga_update_result** out_update_result);
void ga_update_result_destroy(ga_update_result* update_result);
ga_status ga_update_result_get_summary(
    const ga_update_result* update_result,
    ga_deposition_update_summary* out_summary);
/* The returned array is sorted, unique, and valid until the result is
 * destroyed. It may be null when out_count is zero. */
ga_status ga_update_result_get_changed_voxels(
    const ga_update_result* update_result,
    const ga_voxel_id** out_voxel_ids,
    size_t* out_count);

ga_status ga_is_voxel_outside_accessible(
    const ga_grid* grid,
    ga_voxel_id voxel_id,
    uint8_t* out_is_accessible);
ga_status ga_is_site_accessible(
    const ga_grid* grid,
    ga_point3 site_position,
    uint8_t* out_is_accessible);
/* The explicit stencil is allocation-free and is capped at
 * GA_MAX_SITE_VOXEL_COUNT. */
ga_status ga_is_stencil_accessible(
    const ga_grid* grid,
    const ga_voxel_id* voxel_ids,
    size_t voxel_count,
    uint8_t* out_is_accessible);

#ifdef __cplusplus
}
#endif

#endif
