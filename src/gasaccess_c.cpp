#include "gasaccess/gasaccess_c.h"

#include "gasaccess/accessibility_query.hpp"
#include "gasaccess/atom_change_updater.hpp"
#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/adsorption_updater.hpp"
#include "gasaccess/exterior_classifier.hpp"
#include "gasaccess/gas_grid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <memory>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

struct ga_grid {
    explicit ga_grid(gasaccess::GridSpec grid_spec)
        : gas_grid(std::move(grid_spec))
    {
    }

    gasaccess::GasGrid gas_grid;
    std::unique_ptr<gasaccess::AdsorptionUpdater> adsorption_updater;
    std::unique_ptr<gasaccess::AtomChangeUpdater> atom_change_updater;
};

struct ga_update_result {
    gasaccess::AdsorptionUpdateResult update_result;
};

struct ga_atom_change_result {
    gasaccess::AtomChangeUpdateResult update_result;
};

namespace {

thread_local std::array<char, 512> last_error_message{};

void set_last_error(const char* message) noexcept
{
    if (message == nullptr) {
        last_error_message[0] = '\0';
        return;
    }
    static_cast<void>(std::snprintf(
        last_error_message.data(),
        last_error_message.size(),
        "%s",
        message));
}

template <typename Function>
ga_status protect_c_api(Function&& function) noexcept
{
    try {
        std::forward<Function>(function)();
        set_last_error(nullptr);
        return GA_STATUS_SUCCESS;
    } catch (const std::bad_alloc& error) {
        set_last_error(error.what());
        return GA_STATUS_ALLOCATION_FAILED;
    } catch (const std::overflow_error& error) {
        set_last_error(error.what());
        return GA_STATUS_OVERFLOW;
    } catch (const std::out_of_range& error) {
        set_last_error(error.what());
        return GA_STATUS_OUT_OF_RANGE;
    } catch (const std::invalid_argument& error) {
        set_last_error(error.what());
        return GA_STATUS_INVALID_ARGUMENT;
    } catch (const std::length_error& error) {
        set_last_error(error.what());
        return GA_STATUS_OVERFLOW;
    } catch (const std::exception& error) {
        set_last_error(error.what());
        return GA_STATUS_INTERNAL_ERROR;
    } catch (...) {
        set_last_error("unknown C++ exception");
        return GA_STATUS_INTERNAL_ERROR;
    }
}

void require_pointer(const void* pointer, const char* name)
{
    if (pointer == nullptr) {
        throw std::invalid_argument(name);
    }
}

gasaccess::VoxelCoord convert_voxel_coord(ga_voxel_coord voxel_coord)
{
    return {voxel_coord.x, voxel_coord.y, voxel_coord.z};
}

gasaccess::Point3 convert_point(ga_point3 point)
{
    return {point.x, point.y, point.z};
}

gasaccess::GridSpec convert_grid_spec(const ga_grid_spec& c_grid_spec)
{
    if (c_grid_spec.explicit_source_count != 0
        && c_grid_spec.explicit_source_voxels == nullptr) {
        throw std::invalid_argument(
            "explicit source pointer is null with nonzero count");
    }

    gasaccess::GridSpec grid_spec{};
    grid_spec.origin = convert_point(c_grid_spec.origin);
    grid_spec.spacing = {
        c_grid_spec.spacing.x,
        c_grid_spec.spacing.y,
        c_grid_spec.spacing.z
    };
    grid_spec.dimensions = {
        c_grid_spec.dimensions.x,
        c_grid_spec.dimensions.y,
        c_grid_spec.dimensions.z
    };
    grid_spec.periodic = {
        c_grid_spec.periodic.x != 0,
        c_grid_spec.periodic.y != 0,
        c_grid_spec.periodic.z != 0
    };
    grid_spec.reservoir_faces = {
        c_grid_spec.reservoir_faces.x_low != 0,
        c_grid_spec.reservoir_faces.x_high != 0,
        c_grid_spec.reservoir_faces.y_low != 0,
        c_grid_spec.reservoir_faces.y_high != 0,
        c_grid_spec.reservoir_faces.z_low != 0,
        c_grid_spec.reservoir_faces.z_high != 0
    };
    grid_spec.explicit_source_voxels.reserve(c_grid_spec.explicit_source_count);
    for (std::size_t index = 0;
         index < c_grid_spec.explicit_source_count;
         ++index) {
        grid_spec.explicit_source_voxels.push_back(
            convert_voxel_coord(c_grid_spec.explicit_source_voxels[index]));
    }
    return grid_spec;
}

gasaccess::GasState convert_gas_state(ga_gas_state gas_state)
{
    switch (gas_state) {
    case GA_GAS_STATE_UNCLASSIFIED:
        return gasaccess::GasState::Unclassified;
    case GA_GAS_STATE_SOLID:
        return gasaccess::GasState::Solid;
    case GA_GAS_STATE_OUTSIDE_ACCESSIBLE:
        return gasaccess::GasState::OutsideAccessible;
    case GA_GAS_STATE_CLOSED_VOID:
        return gasaccess::GasState::ClosedVoid;
    default:
        throw std::invalid_argument("invalid gas state value");
    }
}

gasaccess::ConnectivityRepairMode convert_repair_mode(
    ga_connectivity_repair_mode repair_mode)
{
    switch (repair_mode) {
    case GA_CONNECTIVITY_REPAIR_AFFECTED_REGION:
        return gasaccess::ConnectivityRepairMode::AffectedRegion;
    case GA_CONNECTIVITY_REPAIR_FULL_RECLASSIFICATION:
        return gasaccess::ConnectivityRepairMode::FullReclassification;
    default:
        throw std::invalid_argument("invalid connectivity repair mode");
    }
}

gasaccess::AtomChangeRepairMode convert_atom_change_repair_mode(
    ga_atom_change_repair_mode repair_mode)
{
    switch (repair_mode) {
    case GA_ATOM_CHANGE_REPAIR_INCREMENTAL:
        return gasaccess::AtomChangeRepairMode::Incremental;
    case GA_ATOM_CHANGE_REPAIR_FULL_RECLASSIFICATION:
        return gasaccess::AtomChangeRepairMode::FullReclassification;
    default:
        throw std::invalid_argument("invalid atom-change repair mode");
    }
}

ga_accessibility_repair_kind convert_repair_kind(
    gasaccess::AccessibilityRepairKind repair_kind)
{
    switch (repair_kind) {
    case gasaccess::AccessibilityRepairKind::None:
        return GA_ACCESSIBILITY_REPAIR_NONE;
    case gasaccess::AccessibilityRepairKind::Closing:
        return GA_ACCESSIBILITY_REPAIR_CLOSING;
    case gasaccess::AccessibilityRepairKind::Opening:
        return GA_ACCESSIBILITY_REPAIR_OPENING;
    case gasaccess::AccessibilityRepairKind::Mixed:
        return GA_ACCESSIBILITY_REPAIR_MIXED;
    case gasaccess::AccessibilityRepairKind::FullReclassification:
        return GA_ACCESSIBILITY_REPAIR_FULL_RECLASSIFICATION;
    }
    throw std::logic_error("invalid atom-change repair kind");
}

gasaccess::AdsorptionUpdater& adsorption_updater(
    ga_grid& grid,
    double precursor_radius,
    gasaccess::ConnectivityRepairMode repair_mode)
{
    if (!grid.adsorption_updater
        || grid.adsorption_updater->precursor_radius() != precursor_radius
        || grid.adsorption_updater->repair_mode() != repair_mode) {
        grid.adsorption_updater = std::make_unique<gasaccess::AdsorptionUpdater>(
            precursor_radius,
            repair_mode);
    }
    return *grid.adsorption_updater;
}

gasaccess::AtomChangeUpdater& atom_change_updater(
    ga_grid& grid,
    double precursor_radius,
    gasaccess::AtomChangeRepairMode repair_mode)
{
    if (!grid.atom_change_updater
        || grid.atom_change_updater->precursor_radius() != precursor_radius
        || grid.atom_change_updater->repair_mode() != repair_mode) {
        grid.atom_change_updater =
            std::make_unique<gasaccess::AtomChangeUpdater>(
                precursor_radius,
                repair_mode);
    }
    return *grid.atom_change_updater;
}

void validate_c_atoms(
    const ga_grid& grid,
    const ga_atom* atoms,
    std::size_t atom_count,
    double precursor_radius)
{
    if (atom_count != 0 && atoms == nullptr) {
        throw std::invalid_argument("atom pointer is null with nonzero count");
    }
    if (!std::isfinite(precursor_radius) || precursor_radius < 0.0) {
        throw std::invalid_argument("precursor radius must be finite and nonnegative");
    }

    for (std::size_t index = 0; index < atom_count; ++index) {
        const auto& atom = atoms[index];
        if (!std::isfinite(atom.position.x)
            || !std::isfinite(atom.position.y)
            || !std::isfinite(atom.position.z)) {
            throw std::invalid_argument("atom position must be finite");
        }
        if (!std::isfinite(atom.radius) || atom.radius < 0.0) {
            throw std::invalid_argument("atom radius must be finite and nonnegative");
        }
        const double excluded_radius = atom.radius + precursor_radius;
        if (!std::isfinite(excluded_radius)
            || !std::isfinite(excluded_radius * excluded_radius)) {
            throw std::invalid_argument("excluded radius is not finite");
        }
        if (!grid.gas_grid.locate_voxel(convert_point(atom.position))) {
            throw std::invalid_argument(
                "atom position is outside a non-periodic grid boundary");
        }
    }
}

std::vector<gasaccess::Atom> convert_c_atoms(
    const ga_atom* atoms,
    std::size_t atom_count)
{
    std::vector<gasaccess::Atom> converted_atoms;
    converted_atoms.reserve(atom_count);
    for (std::size_t index = 0; index < atom_count; ++index) {
        converted_atoms.push_back({
            convert_point(atoms[index].position),
            atoms[index].radius
        });
    }
    return converted_atoms;
}

}  // namespace

extern "C" {

const char* ga_last_error_message(void)
{
    return last_error_message.data();
}

ga_status ga_grid_create(const ga_grid_spec* grid_spec, ga_grid** out_grid)
{
    if (out_grid != nullptr) {
        *out_grid = nullptr;
    }
    return protect_c_api([&]() {
        require_pointer(grid_spec, "grid specification pointer is null");
        require_pointer(out_grid, "output grid pointer is null");
        auto new_grid = std::make_unique<ga_grid>(convert_grid_spec(*grid_spec));
        *out_grid = new_grid.release();
    });
}

void ga_grid_destroy(ga_grid* grid)
{
    delete grid;
}

ga_status ga_grid_get_voxel_count(const ga_grid* grid, ga_voxel_id* out_count)
{
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        require_pointer(out_count, "output voxel-count pointer is null");
        *out_count = grid->gas_grid.voxel_count();
    });
}

ga_status ga_grid_get_voxel_id(
    const ga_grid* grid,
    ga_voxel_coord voxel_coord,
    ga_voxel_id* out_voxel_id)
{
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        require_pointer(out_voxel_id, "output voxel-id pointer is null");
        *out_voxel_id = grid->gas_grid.voxel_id(convert_voxel_coord(voxel_coord));
    });
}

ga_status ga_grid_get_state(
    const ga_grid* grid,
    ga_voxel_id voxel_id,
    ga_gas_state* out_gas_state)
{
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        require_pointer(out_gas_state, "output gas-state pointer is null");
        *out_gas_state = static_cast<ga_gas_state>(grid->gas_grid.gas_state(voxel_id));
    });
}

ga_status ga_grid_get_state_at(
    const ga_grid* grid,
    ga_voxel_coord voxel_coord,
    ga_gas_state* out_gas_state)
{
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        require_pointer(out_gas_state, "output gas-state pointer is null");
        *out_gas_state = static_cast<ga_gas_state>(
            grid->gas_grid.gas_state(convert_voxel_coord(voxel_coord)));
    });
}

ga_status ga_grid_get_state_count(
    const ga_grid* grid,
    ga_gas_state gas_state,
    ga_voxel_id* out_count)
{
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        require_pointer(out_count, "output state-count pointer is null");
        *out_count = grid->gas_grid.gas_state_count(convert_gas_state(gas_state));
    });
}

ga_status ga_grid_set_state(
    ga_grid* grid,
    ga_voxel_id voxel_id,
    ga_gas_state gas_state)
{
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        grid->gas_grid.set_gas_state(voxel_id, convert_gas_state(gas_state));
    });
}

ga_status ga_voxelize_atoms(
    ga_grid* grid,
    const ga_atom* atoms,
    size_t atom_count,
    double precursor_radius,
    ga_voxel_id* out_newly_solid_count)
{
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        require_pointer(
            out_newly_solid_count,
            "output newly-solid-count pointer is null");
        validate_c_atoms(*grid, atoms, atom_count, precursor_radius);

        constexpr std::size_t chunk_capacity = 256;
        std::array<gasaccess::Atom, chunk_capacity> atom_chunk{};
        gasaccess::AtomVoxelizer voxelizer(precursor_radius);
        ga_voxel_id newly_solid_count = 0;

        std::size_t offset = 0;
        while (offset < atom_count) {
            const auto chunk_size = std::min(chunk_capacity, atom_count - offset);
            for (std::size_t index = 0; index < chunk_size; ++index) {
                const auto& atom = atoms[offset + index];
                atom_chunk[index] = {convert_point(atom.position), atom.radius};
            }
            newly_solid_count += voxelizer.voxelize(
                grid->gas_grid,
                {atom_chunk.data(), chunk_size});
            offset += chunk_size;
        }
        *out_newly_solid_count = newly_solid_count;
    });
}

ga_status ga_classify_exterior(
    ga_grid* grid,
    ga_classification_summary* out_summary)
{
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        require_pointer(out_summary, "output classification-summary pointer is null");
        const auto summary = gasaccess::ExteriorClassifier{}.classify(grid->gas_grid);
        out_summary->solid_count = summary.solid_count;
        out_summary->outside_accessible_count = summary.outside_accessible_count;
        out_summary->closed_void_count = summary.closed_void_count;
    });
}

ga_status ga_apply_adsorption(
    ga_grid* grid,
    const ga_atom* adsorbed_atoms,
    size_t atom_count,
    double precursor_radius,
    ga_update_result** out_update_result)
{
    return ga_apply_adsorption_with_mode(
        grid,
        adsorbed_atoms,
        atom_count,
        precursor_radius,
        GA_CONNECTIVITY_REPAIR_AFFECTED_REGION,
        out_update_result);
}

ga_status ga_apply_adsorption_with_mode(
    ga_grid* grid,
    const ga_atom* adsorbed_atoms,
    size_t atom_count,
    double precursor_radius,
    ga_connectivity_repair_mode repair_mode,
    ga_update_result** out_update_result)
{
    if (out_update_result != nullptr) {
        *out_update_result = nullptr;
    }
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        require_pointer(out_update_result, "output update-result pointer is null");
        const auto converted_repair_mode = convert_repair_mode(repair_mode);
        validate_c_atoms(*grid, adsorbed_atoms, atom_count, precursor_radius);
        const auto converted_atoms = convert_c_atoms(adsorbed_atoms, atom_count);
        auto update_result = std::make_unique<ga_update_result>();
        update_result->update_result = adsorption_updater(
            *grid,
            precursor_radius,
            converted_repair_mode).apply_adsorption(
                grid->gas_grid,
                {converted_atoms.data(), converted_atoms.size()});
        *out_update_result = update_result.release();
    });
}

void ga_update_result_destroy(ga_update_result* update_result)
{
    delete update_result;
}

ga_status ga_update_result_get_summary(
    const ga_update_result* update_result,
    ga_adsorption_update_summary* out_summary)
{
    return protect_c_api([&]() {
        require_pointer(update_result, "update-result pointer is null");
        require_pointer(out_summary, "output adsorption-summary pointer is null");
        const auto& result = update_result->update_result;
        out_summary->newly_solid_count = result.newly_solid_count;
        out_summary->changed_voxel_count = result.changed_voxel_ids.size();
        out_summary->repair_visited_voxel_count = result.repair_visited_voxel_count;
        out_summary->repair_closed_voxel_count = result.repair_closed_voxel_count;
        out_summary->full_reclassification_performed =
            result.used_full_reclassification() ? 1U : 0U;
        out_summary->affected_region_repair_performed =
            result.used_affected_region_repair() ? 1U : 0U;
        out_summary->classification.solid_count = result.classification.solid_count;
        out_summary->classification.outside_accessible_count =
            result.classification.outside_accessible_count;
        out_summary->classification.closed_void_count =
            result.classification.closed_void_count;
    });
}

ga_status ga_update_result_get_changed_voxels(
    const ga_update_result* update_result,
    const ga_voxel_id** out_voxel_ids,
    size_t* out_count)
{
    return protect_c_api([&]() {
        require_pointer(update_result, "update-result pointer is null");
        require_pointer(out_voxel_ids, "output changed-voxel pointer is null");
        require_pointer(out_count, "output changed-voxel-count pointer is null");
        const auto& changed_voxel_ids = update_result->update_result.changed_voxel_ids;
        *out_voxel_ids = changed_voxel_ids.data();
        *out_count = changed_voxel_ids.size();
    });
}

ga_status ga_apply_desorption(
    ga_grid* grid,
    const ga_atom* removed_atoms,
    size_t atom_count,
    double precursor_radius,
    ga_atom_change_result** out_change_result)
{
    return ga_apply_desorption_with_mode(
        grid,
        removed_atoms,
        atom_count,
        precursor_radius,
        GA_ATOM_CHANGE_REPAIR_INCREMENTAL,
        out_change_result);
}

ga_status ga_apply_desorption_with_mode(
    ga_grid* grid,
    const ga_atom* removed_atoms,
    size_t atom_count,
    double precursor_radius,
    ga_atom_change_repair_mode repair_mode,
    ga_atom_change_result** out_change_result)
{
    const ga_atom_change_batch atom_changes = {
        {nullptr, 0},
        {removed_atoms, atom_count}};
    return ga_apply_atom_changes_with_mode(
        grid,
        &atom_changes,
        precursor_radius,
        repair_mode,
        out_change_result);
}

ga_status ga_apply_atom_changes(
    ga_grid* grid,
    const ga_atom_change_batch* atom_changes,
    double precursor_radius,
    ga_atom_change_result** out_change_result)
{
    return ga_apply_atom_changes_with_mode(
        grid,
        atom_changes,
        precursor_radius,
        GA_ATOM_CHANGE_REPAIR_INCREMENTAL,
        out_change_result);
}

ga_status ga_apply_atom_changes_with_mode(
    ga_grid* grid,
    const ga_atom_change_batch* atom_changes,
    double precursor_radius,
    ga_atom_change_repair_mode repair_mode,
    ga_atom_change_result** out_change_result)
{
    if (out_change_result != nullptr) {
        *out_change_result = nullptr;
    }
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        require_pointer(atom_changes, "atom-change batch pointer is null");
        require_pointer(
            out_change_result,
            "output atom-change-result pointer is null");
        const auto converted_repair_mode =
            convert_atom_change_repair_mode(repair_mode);
        validate_c_atoms(
            *grid,
            atom_changes->added_atoms.atoms,
            atom_changes->added_atoms.count,
            precursor_radius);
        validate_c_atoms(
            *grid,
            atom_changes->removed_atoms.atoms,
            atom_changes->removed_atoms.count,
            precursor_radius);
        const auto added_atoms = convert_c_atoms(
            atom_changes->added_atoms.atoms,
            atom_changes->added_atoms.count);
        const auto removed_atoms = convert_c_atoms(
            atom_changes->removed_atoms.atoms,
            atom_changes->removed_atoms.count);

        auto change_result = std::make_unique<ga_atom_change_result>();
        change_result->update_result = atom_change_updater(
            *grid,
            precursor_radius,
            converted_repair_mode).apply_atom_changes(
                grid->gas_grid,
                {{added_atoms.data(), added_atoms.size()},
                 {removed_atoms.data(), removed_atoms.size()}});
        *out_change_result = change_result.release();
    });
}

void ga_atom_change_result_destroy(ga_atom_change_result* change_result)
{
    delete change_result;
}

ga_status ga_atom_change_result_get_summary(
    const ga_atom_change_result* change_result,
    ga_atom_change_update_summary* out_summary)
{
    return protect_c_api([&]() {
        require_pointer(change_result, "atom-change-result pointer is null");
        require_pointer(
            out_summary,
            "output atom-change-summary pointer is null");
        const auto& result = change_result->update_result;
        out_summary->blocker_count_changed_voxel_count =
            result.blocker_count_changed_voxel_count;
        out_summary->newly_solid_count = result.newly_solid_count;
        out_summary->newly_gas_count = result.newly_gas_count;
        out_summary->changed_voxel_count = result.changed_voxel_ids.size();
        out_summary->repair_kind = convert_repair_kind(result.repair_kind);
        out_summary->full_reclassification_performed =
            result.used_full_reclassification() ? 1U : 0U;
        out_summary->closing_repair_performed =
            result.used_closing_repair() ? 1U : 0U;
        out_summary->opening_repair_performed =
            result.used_opening_repair() ? 1U : 0U;
        out_summary->closing_visited_voxel_count =
            result.closing_visited_voxel_count;
        out_summary->opening_visited_voxel_count =
            result.opening_visited_voxel_count;
        out_summary->repair_closed_voxel_count =
            result.repair_closed_voxel_count;
        out_summary->repair_opened_voxel_count =
            result.repair_opened_voxel_count;
        out_summary->classification.solid_count =
            result.classification.solid_count;
        out_summary->classification.outside_accessible_count =
            result.classification.outside_accessible_count;
        out_summary->classification.closed_void_count =
            result.classification.closed_void_count;
    });
}

ga_status ga_atom_change_result_get_changed_voxels(
    const ga_atom_change_result* change_result,
    const ga_voxel_id** out_voxel_ids,
    size_t* out_count)
{
    return protect_c_api([&]() {
        require_pointer(change_result, "atom-change-result pointer is null");
        require_pointer(out_voxel_ids, "output changed-voxel pointer is null");
        require_pointer(out_count, "output changed-voxel-count pointer is null");
        const auto& changed_voxel_ids =
            change_result->update_result.changed_voxel_ids;
        *out_voxel_ids = changed_voxel_ids.data();
        *out_count = changed_voxel_ids.size();
    });
}

ga_status ga_is_voxel_outside_accessible(
    const ga_grid* grid,
    ga_voxel_id voxel_id,
    uint8_t* out_is_accessible)
{
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        require_pointer(out_is_accessible, "output accessibility pointer is null");
        const gasaccess::GasAccessibilityQuery query(grid->gas_grid);
        *out_is_accessible = query.is_voxel_outside_accessible(voxel_id) ? 1U : 0U;
    });
}

ga_status ga_is_site_accessible(
    const ga_grid* grid,
    ga_point3 site_position,
    uint8_t* out_is_accessible)
{
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        require_pointer(out_is_accessible, "output accessibility pointer is null");
        const gasaccess::GasAccessibilityQuery query(grid->gas_grid);
        *out_is_accessible = query.is_site_accessible(convert_point(site_position))
            ? 1U
            : 0U;
    });
}

ga_status ga_is_stencil_accessible(
    const ga_grid* grid,
    const ga_voxel_id* voxel_ids,
    size_t voxel_count,
    uint8_t* out_is_accessible)
{
    return protect_c_api([&]() {
        require_pointer(grid, "grid pointer is null");
        require_pointer(out_is_accessible, "output accessibility pointer is null");
        const gasaccess::GasAccessibilityQuery query(grid->gas_grid);
        *out_is_accessible = query.is_site_accessible({voxel_ids, voxel_count})
            ? 1U
            : 0U;
    });
}

}  // extern "C"
