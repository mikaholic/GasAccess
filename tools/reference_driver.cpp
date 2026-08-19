#include "gasaccess/accessibility_query.hpp"
#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/deposition_updater.hpp"
#include "gasaccess/exterior_classifier.hpp"
#include "gasaccess/gas_grid.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using gasaccess::Atom;
using gasaccess::AtomView;
using gasaccess::ClassificationSummary;
using gasaccess::GasGrid;
using gasaccess::GasState;
using gasaccess::GridSpacing;
using gasaccess::GridSpec;
using gasaccess::Point3;
using gasaccess::VoxelCoord;
using gasaccess::VoxelId;

constexpr std::uint64_t fnv_offset_basis = UINT64_C(14695981039346656037);
constexpr std::uint64_t fnv_prime = UINT64_C(1099511628211);

struct Options {
    std::string scenario = "open-trench";
    std::uint64_t nx = 64;
    std::uint64_t ny = 32;
    std::uint64_t nz = 64;
    GridSpacing spacing{1.0, 1.0, 1.0};
    double atom_radius = 0.2;
    double precursor_radius = 0.25;
    std::uint64_t atom_count = 10000;
    std::uint64_t query_count = 100000;
    std::uint64_t update_count = 5;
    std::uint64_t seed = UINT64_C(5489);
    bool periodic_x = false;
    bool periodic_y = false;
    bool periodic_z = false;
};

class SplitMix64 {
public:
    explicit SplitMix64(std::uint64_t seed) noexcept
        : state_(seed)
    {
    }

    std::uint64_t next() noexcept
    {
        state_ += UINT64_C(0x9e3779b97f4a7c15);
        std::uint64_t value = state_;
        value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
        value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
        return value ^ (value >> 31U);
    }

    std::uint64_t bounded(std::uint64_t upper_bound) noexcept
    {
        const std::uint64_t threshold = (UINT64_C(0) - upper_bound) % upper_bound;
        for (;;) {
            const std::uint64_t value = next();
            if (value >= threshold) {
                return value % upper_bound;
            }
        }
    }

private:
    std::uint64_t state_ = 0;
};

double elapsed_ms(
    std::chrono::steady_clock::time_point start,
    std::chrono::steady_clock::time_point end)
{
    return std::chrono::duration<double, std::milli>(end - start).count();
}

double elapsed_us(
    std::chrono::steady_clock::time_point start,
    std::chrono::steady_clock::time_point end)
{
    return std::chrono::duration<double, std::micro>(end - start).count();
}

std::uint64_t parse_uint64(const std::string& text, const std::string& option_name)
{
    std::size_t parsed_count = 0;
    unsigned long long parsed_value = 0;
    try {
        parsed_value = std::stoull(text, &parsed_count, 10);
    } catch (const std::exception&) {
        throw std::invalid_argument(option_name + " requires a nonnegative integer");
    }
    if (parsed_count != text.size() || (!text.empty() && text.front() == '-')) {
        throw std::invalid_argument(option_name + " requires a nonnegative integer");
    }
    if (parsed_value > std::numeric_limits<std::uint64_t>::max()) {
        throw std::out_of_range(option_name + " is too large");
    }
    return static_cast<std::uint64_t>(parsed_value);
}

double parse_double(const std::string& text, const std::string& option_name)
{
    std::size_t parsed_count = 0;
    double parsed_value = 0.0;
    try {
        parsed_value = std::stod(text, &parsed_count);
    } catch (const std::exception&) {
        throw std::invalid_argument(option_name + " requires a finite number");
    }
    if (parsed_count != text.size() || !std::isfinite(parsed_value)) {
        throw std::invalid_argument(option_name + " requires a finite number");
    }
    return parsed_value;
}

std::size_t checked_size(std::uint64_t value, const std::string& field_name)
{
    if (value > std::numeric_limits<std::size_t>::max()) {
        throw std::length_error(field_name + " exceeds the platform size limit");
    }
    return static_cast<std::size_t>(value);
}

void print_help(const char* executable_name)
{
    std::cout
        << "Usage: " << executable_name << " [options]\n"
        << "\n"
        << "Deterministic GasAccess correctness and timing driver.\n"
        << "\n"
        << "Options:\n"
        << "  --scenario NAME       open-trench, sealed-trench, pinch-off, or bulk\n"
        << "  --nx N --ny N --nz N  voxel dimensions (defaults: 64 32 64)\n"
        << "  --spacing VALUE       set all voxel spacings (default: 1.0)\n"
        << "  --spacing-x VALUE     override x voxel spacing\n"
        << "  --spacing-y VALUE     override y voxel spacing\n"
        << "  --spacing-z VALUE     override z voxel spacing\n"
        << "  --atom-radius VALUE   solid atom radius (default: 0.2)\n"
        << "  --precursor-radius V  spherical steric radius (default: 0.25)\n"
        << "  --atom-count N        generated atom records for bulk (default: 10000)\n"
        << "  --query-count N       cached site queries (default: 100000)\n"
        << "  --update-count N      one-atom deposition updates (default: 5)\n"
        << "  --seed N              deterministic generator seed (default: 5489)\n"
        << "  --periodic-x           enable periodic x boundary\n"
        << "  --periodic-y           enable periodic y boundary\n"
        << "  --periodic-z           enable periodic z boundary\n"
        << "  --help                 display this help\n";
}

Options parse_options(int argc, char* argv[])
{
    Options options{};
    for (int argument_index = 1; argument_index < argc; ++argument_index) {
        const std::string option = argv[argument_index];
        const auto require_value = [&]() -> std::string {
            if (argument_index + 1 >= argc) {
                throw std::invalid_argument(option + " requires a value");
            }
            ++argument_index;
            return argv[argument_index];
        };

        if (option == "--scenario") {
            options.scenario = require_value();
        } else if (option == "--nx") {
            options.nx = parse_uint64(require_value(), option);
        } else if (option == "--ny") {
            options.ny = parse_uint64(require_value(), option);
        } else if (option == "--nz") {
            options.nz = parse_uint64(require_value(), option);
        } else if (option == "--spacing") {
            const double spacing = parse_double(require_value(), option);
            options.spacing = {spacing, spacing, spacing};
        } else if (option == "--spacing-x") {
            options.spacing.x = parse_double(require_value(), option);
        } else if (option == "--spacing-y") {
            options.spacing.y = parse_double(require_value(), option);
        } else if (option == "--spacing-z") {
            options.spacing.z = parse_double(require_value(), option);
        } else if (option == "--atom-radius") {
            options.atom_radius = parse_double(require_value(), option);
        } else if (option == "--precursor-radius") {
            options.precursor_radius = parse_double(require_value(), option);
        } else if (option == "--atom-count") {
            options.atom_count = parse_uint64(require_value(), option);
        } else if (option == "--query-count") {
            options.query_count = parse_uint64(require_value(), option);
        } else if (option == "--update-count") {
            options.update_count = parse_uint64(require_value(), option);
        } else if (option == "--seed") {
            options.seed = parse_uint64(require_value(), option);
        } else if (option == "--periodic-x") {
            options.periodic_x = true;
        } else if (option == "--periodic-y") {
            options.periodic_y = true;
        } else if (option == "--periodic-z") {
            options.periodic_z = true;
        } else if (option == "--help") {
            print_help(argv[0]);
            std::exit(EXIT_SUCCESS);
        } else {
            throw std::invalid_argument("unknown option: " + option);
        }
    }
    return options;
}

void validate_options(const Options& options)
{
    if (options.scenario != "open-trench"
        && options.scenario != "sealed-trench"
        && options.scenario != "pinch-off"
        && options.scenario != "bulk") {
        throw std::invalid_argument(
            "scenario must be open-trench, sealed-trench, pinch-off, or bulk");
    }
    if (options.nx == 0 || options.ny == 0 || options.nz == 0) {
        throw std::invalid_argument("all voxel dimensions must be positive");
    }
    if (!(options.spacing.x > 0.0)
        || !(options.spacing.y > 0.0)
        || !(options.spacing.z > 0.0)) {
        throw std::invalid_argument("all spacing components must be positive");
    }
    if (options.atom_radius < 0.0 || options.precursor_radius < 0.0) {
        throw std::invalid_argument("atom and precursor radii must be nonnegative");
    }
    if (options.scenario != "bulk" && (options.nx < 5 || options.nz < 4)) {
        throw std::invalid_argument("trench scenarios require nx >= 5 and nz >= 4");
    }
}

GridSpec make_grid_spec(const Options& options)
{
    GridSpec grid_spec{};
    grid_spec.spacing = options.spacing;
    grid_spec.dimensions = {options.nx, options.ny, options.nz};
    grid_spec.periodic = {
        options.periodic_x,
        options.periodic_y,
        options.periodic_z
    };

    const bool is_trench = options.scenario == "open-trench"
        || options.scenario == "sealed-trench"
        || options.scenario == "pinch-off";
    if (!options.periodic_z) {
        grid_spec.reservoir_faces.z_high = true;
    } else if (is_trench) {
        grid_spec.explicit_source_voxels.push_back({
            0,
            0,
            static_cast<std::int64_t>(options.nz - 1)
        });
    } else if (!options.periodic_y) {
        grid_spec.reservoir_faces.y_high = true;
    } else if (!options.periodic_x) {
        grid_spec.reservoir_faces.x_high = true;
    } else {
        grid_spec.explicit_source_voxels.push_back({0, 0, 0});
    }
    return grid_spec;
}

Atom atom_at(const GasGrid& gas_grid, const VoxelCoord& voxel_coord, double radius)
{
    return {gas_grid.voxel_center(voxel_coord), radius};
}

std::vector<Atom> make_trench_atoms(
    const GasGrid& gas_grid,
    const Options& options,
    bool sealed)
{
    std::vector<Atom> atoms;
    const std::uint64_t left_wall = options.nx / 3;
    const std::uint64_t right_wall = options.nx - left_wall - 1;
    const std::uint64_t roof_z = options.nz - 2;

    for (std::uint64_t z = 0; z <= roof_z; ++z) {
        for (std::uint64_t y = 0; y < options.ny; ++y) {
            atoms.push_back(atom_at(
                gas_grid,
                {static_cast<std::int64_t>(left_wall),
                 static_cast<std::int64_t>(y),
                 static_cast<std::int64_t>(z)},
                options.atom_radius));
            atoms.push_back(atom_at(
                gas_grid,
                {static_cast<std::int64_t>(right_wall),
                 static_cast<std::int64_t>(y),
                 static_cast<std::int64_t>(z)},
                options.atom_radius));
        }
    }

    if (options.periodic_z) {
        for (std::uint64_t y = 0; y < options.ny; ++y) {
            for (std::uint64_t x = left_wall + 1; x < right_wall; ++x) {
                atoms.push_back(atom_at(
                    gas_grid,
                    {static_cast<std::int64_t>(x),
                     static_cast<std::int64_t>(y),
                     0},
                    options.atom_radius));
            }
        }
    }

    if (sealed) {
        for (std::uint64_t y = 0; y < options.ny; ++y) {
            for (std::uint64_t x = left_wall + 1; x < right_wall; ++x) {
                atoms.push_back(atom_at(
                    gas_grid,
                    {static_cast<std::int64_t>(x),
                     static_cast<std::int64_t>(y),
                     static_cast<std::int64_t>(roof_z)},
                    options.atom_radius));
            }
        }
    }
    return atoms;
}

std::vector<Atom> make_bulk_atoms(const GasGrid& gas_grid, const Options& options)
{
    std::vector<Atom> atoms;
    atoms.reserve(checked_size(options.atom_count, "atom count"));
    SplitMix64 generator(options.seed);
    for (std::uint64_t atom_index = 0;
         atom_index < options.atom_count;
         ++atom_index) {
        const VoxelCoord voxel_coord{
            static_cast<std::int64_t>(generator.bounded(options.nx)),
            static_cast<std::int64_t>(generator.bounded(options.ny)),
            static_cast<std::int64_t>(generator.bounded(options.nz))
        };
        atoms.push_back(atom_at(gas_grid, voxel_coord, options.atom_radius));
    }
    return atoms;
}

std::vector<Atom> make_atoms(const GasGrid& gas_grid, const Options& options)
{
    if (options.scenario == "open-trench"
        || options.scenario == "pinch-off") {
        return make_trench_atoms(gas_grid, options, false);
    }
    if (options.scenario == "sealed-trench") {
        return make_trench_atoms(gas_grid, options, true);
    }
    return make_bulk_atoms(gas_grid, options);
}

std::vector<Atom> make_update_atoms(
    const GasGrid& gas_grid,
    const Options& options)
{
    const auto update_capacity = options.update_count < gas_grid.voxel_count()
        ? options.update_count
        : gas_grid.voxel_count();
    const auto update_capacity_size = checked_size(update_capacity, "update count");
    std::vector<Atom> update_atoms;
    update_atoms.reserve(update_capacity_size);

    if (options.scenario == "pinch-off") {
        const std::uint64_t left_wall = options.nx / 3;
        const std::uint64_t right_wall = options.nx - left_wall - 1;
        const std::uint64_t roof_z = options.nz - 2;
        for (std::uint64_t y = 0;
             y < options.ny && update_atoms.size() < update_capacity_size;
             ++y) {
            for (std::uint64_t x = left_wall + 1;
                 x < right_wall && update_atoms.size() < update_capacity_size;
                 ++x) {
                update_atoms.push_back(atom_at(
                    gas_grid,
                    {static_cast<std::int64_t>(x),
                     static_cast<std::int64_t>(y),
                     static_cast<std::int64_t>(roof_z)},
                    options.atom_radius));
            }
        }
        return update_atoms;
    }

    for (VoxelId voxel_id = 0;
         voxel_id < gas_grid.voxel_count()
             && update_atoms.size() < update_capacity_size;
         ++voxel_id) {
        if (gas_grid.gas_state(voxel_id) != GasState::Solid) {
            update_atoms.push_back(
                {gas_grid.voxel_center(voxel_id), options.atom_radius});
        }
    }
    return update_atoms;
}

double percentile_us(
    const std::vector<double>& sorted_latencies_us,
    double percentile)
{
    if (sorted_latencies_us.empty()) {
        return 0.0;
    }
    const double rank = std::ceil(
        percentile * static_cast<double>(sorted_latencies_us.size()));
    const auto index = static_cast<std::size_t>(rank) - 1;
    return sorted_latencies_us[index];
}

void print_latency_summary(
    const std::string& prefix,
    std::vector<double> latencies_us)
{
    std::sort(latencies_us.begin(), latencies_us.end());
    std::cout << prefix << "_sample_count=" << latencies_us.size() << '\n';
    std::cout << prefix << "_median_us="
              << percentile_us(latencies_us, 0.50) << '\n';
    std::cout << prefix << "_p95_us="
              << percentile_us(latencies_us, 0.95) << '\n';
    std::cout << prefix << "_max_us="
              << (latencies_us.empty() ? 0.0 : latencies_us.back()) << '\n';
}

std::uint64_t state_checksum(const GasGrid& gas_grid)
{
    std::uint64_t checksum = fnv_offset_basis;
    for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
        checksum ^= static_cast<std::uint8_t>(gas_grid.gas_state(voxel_id));
        checksum *= fnv_prime;
    }
    return checksum;
}

void validate_summary(const GasGrid& gas_grid, const ClassificationSummary& summary)
{
    if (summary.solid_count > gas_grid.voxel_count()
        || summary.outside_accessible_count
            > gas_grid.voxel_count() - summary.solid_count
        || summary.closed_void_count
            != gas_grid.voxel_count() - summary.solid_count
                - summary.outside_accessible_count
        || gas_grid.gas_state_count(GasState::Unclassified) != 0
        || gas_grid.gas_state_count(GasState::Solid) != summary.solid_count
        || gas_grid.gas_state_count(GasState::OutsideAccessible)
            != summary.outside_accessible_count
        || gas_grid.gas_state_count(GasState::ClosedVoid)
            != summary.closed_void_count) {
        throw std::logic_error("classification summary does not match voxel count");
    }
}

std::optional<std::uint64_t> read_peak_rss_bytes()
{
#if defined(__linux__)
    std::ifstream status_file("/proc/self/status");
    std::string line;
    while (std::getline(status_file, line)) {
        if (line.compare(0, 6, "VmHWM:") != 0) {
            continue;
        }
        std::istringstream line_stream(line.substr(6));
        std::uint64_t value_kib = 0;
        std::string unit;
        if (line_stream >> value_kib >> unit && unit == "kB"
            && value_kib <= std::numeric_limits<std::uint64_t>::max() / 1024U) {
            return value_kib * 1024U;
        }
    }
#endif
    return std::nullopt;
}

void print_checksum(const std::string& name, std::uint64_t checksum)
{
    std::cout << name << "=0x" << std::hex << std::setw(16) << std::setfill('0')
              << checksum << std::dec << std::setfill(' ') << '\n';
}

int run(const Options& options)
{
    validate_options(options);

    const auto grid_start = std::chrono::steady_clock::now();
    GasGrid gas_grid(make_grid_spec(options));
    const auto grid_end = std::chrono::steady_clock::now();

    const auto atom_start = std::chrono::steady_clock::now();
    const auto atoms = make_atoms(gas_grid, options);
    const auto atom_end = std::chrono::steady_clock::now();

    const auto voxelization_start = std::chrono::steady_clock::now();
    const VoxelId newly_solid_count = gasaccess::AtomVoxelizer(
        options.precursor_radius).voxelize(
            gas_grid,
            {atoms.data(), atoms.size()});
    const auto voxelization_end = std::chrono::steady_clock::now();

    const auto classification_start = std::chrono::steady_clock::now();
    const auto initial_summary = gasaccess::ExteriorClassifier{}.classify(gas_grid);
    const auto classification_end = std::chrono::steady_clock::now();
    validate_summary(gas_grid, initial_summary);
    const auto initial_checksum = state_checksum(gas_grid);

    SplitMix64 query_generator(options.seed ^ UINT64_C(0xd1b54a32d192ed03));
    gasaccess::GasAccessibilityQuery accessibility_query(gas_grid);
    VoxelId accessible_query_count = 0;
    const auto query_start = std::chrono::steady_clock::now();
    for (std::uint64_t query_index = 0;
         query_index < options.query_count;
         ++query_index) {
        const auto voxel_id = query_generator.bounded(gas_grid.voxel_count());
        if (accessibility_query.is_site_accessible(gas_grid.voxel_center(voxel_id))) {
            ++accessible_query_count;
        }
    }
    const auto query_end = std::chrono::steady_clock::now();

    const auto update_atoms = make_update_atoms(gas_grid, options);

    VoxelId geometry_changing_update_count = 0;
    VoxelId locally_safe_update_count = 0;
    VoxelId full_reclassification_update_count = 0;
    VoxelId affected_region_repair_update_count = 0;
    VoxelId repair_visited_voxel_count = 0;
    VoxelId repair_closed_voxel_count = 0;
    VoxelId update_newly_solid_count = 0;
    VoxelId changed_state_count = 0;
    std::vector<double> update_latencies_us;
    std::vector<double> safe_update_latencies_us;
    std::vector<double> repair_update_latencies_us;
    std::vector<double> full_update_latencies_us;
    update_latencies_us.reserve(update_atoms.size());
    safe_update_latencies_us.reserve(update_atoms.size());
    repair_update_latencies_us.reserve(update_atoms.size());
    full_update_latencies_us.reserve(update_atoms.size());
    ClassificationSummary final_summary = initial_summary;
    const gasaccess::DepositionUpdater deposition_updater(options.precursor_radius);
    const auto update_start = std::chrono::steady_clock::now();
    for (const auto& update_atom : update_atoms) {
        const auto single_update_start = std::chrono::steady_clock::now();
        const auto result = deposition_updater.apply_deposition(
            gas_grid,
            AtomView{&update_atom, 1});
        const auto single_update_end = std::chrono::steady_clock::now();
        const double update_latency_us = elapsed_us(
            single_update_start,
            single_update_end);
        update_latencies_us.push_back(update_latency_us);
        if (result.geometry_changed()) {
            ++geometry_changing_update_count;
        }
        if (result.used_full_reclassification()) {
            ++full_reclassification_update_count;
        }
        if (result.used_affected_region_repair()) {
            ++affected_region_repair_update_count;
            repair_update_latencies_us.push_back(update_latency_us);
        } else if (result.used_full_reclassification()) {
            full_update_latencies_us.push_back(update_latency_us);
        } else if (result.geometry_changed()) {
            ++locally_safe_update_count;
            safe_update_latencies_us.push_back(update_latency_us);
        }
        repair_visited_voxel_count += result.repair_visited_voxel_count;
        repair_closed_voxel_count += result.repair_closed_voxel_count;
        update_newly_solid_count += result.newly_solid_count;
        changed_state_count += static_cast<VoxelId>(result.changed_voxel_ids.size());
        final_summary = result.classification;
    }
    const auto update_end = std::chrono::steady_clock::now();
    validate_summary(gas_grid, final_summary);
    const auto final_checksum = state_checksum(gas_grid);

    const double query_ms = elapsed_ms(query_start, query_end);
    const double query_throughput = query_ms > 0.0
        ? static_cast<double>(options.query_count) * 1000.0 / query_ms
        : 0.0;
    const std::uint64_t persistent_state_bytes = gas_grid.voxel_count()
        * static_cast<std::uint64_t>(sizeof(GasState));
    const std::uint64_t traversal_payload_upper_bound_bytes = gas_grid.voxel_count()
        <= std::numeric_limits<std::uint64_t>::max() / sizeof(VoxelId)
        ? gas_grid.voxel_count() * static_cast<std::uint64_t>(sizeof(VoxelId))
        : std::numeric_limits<std::uint64_t>::max();
    const auto peak_rss_bytes = read_peak_rss_bytes();

    std::cout << std::boolalpha << std::fixed << std::setprecision(6);
    std::cout << "driver_version=3\n";
    std::cout << "scenario=" << options.scenario << '\n';
    std::cout << "nx=" << options.nx << '\n';
    std::cout << "ny=" << options.ny << '\n';
    std::cout << "nz=" << options.nz << '\n';
    std::cout << "spacing_x=" << options.spacing.x << '\n';
    std::cout << "spacing_y=" << options.spacing.y << '\n';
    std::cout << "spacing_z=" << options.spacing.z << '\n';
    std::cout << "atom_radius=" << options.atom_radius << '\n';
    std::cout << "precursor_radius=" << options.precursor_radius << '\n';
    std::cout << "periodic_x=" << options.periodic_x << '\n';
    std::cout << "periodic_y=" << options.periodic_y << '\n';
    std::cout << "periodic_z=" << options.periodic_z << '\n';
    std::cout << "seed=" << options.seed << '\n';
    std::cout << "atom_count=" << atoms.size() << '\n';
    std::cout << "voxel_count=" << gas_grid.voxel_count() << '\n';
    std::cout << "voxelization_newly_solid_count=" << newly_solid_count << '\n';
    std::cout << "initial_solid_count=" << initial_summary.solid_count << '\n';
    std::cout << "initial_outside_accessible_count="
              << initial_summary.outside_accessible_count << '\n';
    std::cout << "initial_closed_void_count=" << initial_summary.closed_void_count << '\n';
    print_checksum("initial_state_checksum", initial_checksum);
    std::cout << "query_count=" << options.query_count << '\n';
    std::cout << "accessible_query_count=" << accessible_query_count << '\n';
    std::cout << "update_count_requested=" << options.update_count << '\n';
    std::cout << "update_count_performed=" << update_atoms.size() << '\n';
    std::cout << "geometry_changing_update_count="
              << geometry_changing_update_count << '\n';
    std::cout << "locally_safe_update_count="
              << locally_safe_update_count << '\n';
    std::cout << "full_reclassification_update_count="
              << full_reclassification_update_count << '\n';
    std::cout << "affected_region_repair_update_count="
              << affected_region_repair_update_count << '\n';
    std::cout << "repair_visited_voxel_count="
              << repair_visited_voxel_count << '\n';
    std::cout << "repair_closed_voxel_count="
              << repair_closed_voxel_count << '\n';
    std::cout << "update_newly_solid_count=" << update_newly_solid_count << '\n';
    std::cout << "changed_state_count=" << changed_state_count << '\n';
    std::cout << "final_solid_count=" << final_summary.solid_count << '\n';
    std::cout << "final_outside_accessible_count="
              << final_summary.outside_accessible_count << '\n';
    std::cout << "final_closed_void_count=" << final_summary.closed_void_count << '\n';
    print_checksum("final_state_checksum", final_checksum);
    std::cout << "grid_construction_ms=" << elapsed_ms(grid_start, grid_end) << '\n';
    std::cout << "atom_generation_ms=" << elapsed_ms(atom_start, atom_end) << '\n';
    std::cout << "voxelization_ms="
              << elapsed_ms(voxelization_start, voxelization_end) << '\n';
    std::cout << "classification_ms="
              << elapsed_ms(classification_start, classification_end) << '\n';
    std::cout << "query_ms=" << query_ms << '\n';
    std::cout << "query_throughput_per_second=" << query_throughput << '\n';
    std::cout << "deposition_update_ms=" << elapsed_ms(update_start, update_end) << '\n';
    print_latency_summary("update_latency", update_latencies_us);
    print_latency_summary("safe_update_latency", safe_update_latencies_us);
    print_latency_summary("repair_update_latency", repair_update_latencies_us);
    print_latency_summary("full_update_latency", full_update_latencies_us);
    std::cout << "persistent_state_storage_bytes=" << persistent_state_bytes << '\n';
    std::cout << "persistent_state_bytes_per_voxel=" << sizeof(GasState) << '\n';
    std::cout << "traversal_frontier_payload_upper_bound_bytes="
              << traversal_payload_upper_bound_bytes << '\n';
    if (peak_rss_bytes) {
        std::cout << "peak_rss_bytes=" << *peak_rss_bytes << '\n';
        std::cout << "peak_rss_bytes_per_voxel="
                  << static_cast<double>(*peak_rss_bytes)
                      / static_cast<double>(gas_grid.voxel_count())
                  << '\n';
    } else {
        std::cout << "peak_rss_bytes=unavailable\n";
        std::cout << "peak_rss_bytes_per_voxel=unavailable\n";
    }
    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char* argv[])
{
    try {
        return run(parse_options(argc, argv));
    } catch (const std::exception& exception) {
        std::cerr << "error=" << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
