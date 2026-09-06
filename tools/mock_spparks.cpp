#include "mock_spparks.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace gasaccess::testing {
namespace {

struct AtomWire {
    std::int64_t id = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double radius = 0.0;
};

void check_mpi(int error_code, const char* operation)
{
    if (error_code == MPI_SUCCESS) {
        return;
    }
    throw std::runtime_error(std::string(operation) + " failed");
}

int checked_int(std::size_t value, const char* description)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(std::string(description) + " exceeds MPI int range");
    }
    return static_cast<int>(value);
}

bool is_finite_point(const Point3& point) noexcept
{
    return std::isfinite(point.x) && std::isfinite(point.y)
        && std::isfinite(point.z);
}

std::int64_t normalized_location(
    std::int64_t location,
    int process_count,
    bool periodic,
    bool& valid)
{
    if (location >= 0 && location < process_count) {
        return location;
    }
    if (!periodic) {
        valid = false;
        return 0;
    }
    if (location < 0) {
        return static_cast<std::int64_t>(process_count - 1);
    }
    return 0;
}

std::vector<int> axis_offsets(
    double coordinate,
    double lower,
    double upper,
    double ghost_distance)
{
    std::vector<int> offsets{0};
    if (coordinate - lower <= ghost_distance) {
        offsets.push_back(-1);
    }
    if (upper - coordinate <= ghost_distance) {
        offsets.push_back(1);
    }
    return offsets;
}

}  // namespace

MockSpparksDomain::MockSpparksDomain(
    MPI_Comm communicator,
    const Point3& global_lower,
    const Point3& global_upper,
    const GridDimensions& process_grid,
    const PeriodicAxes& periodic)
    : world(communicator)
{
    if (communicator == MPI_COMM_NULL) {
        throw std::invalid_argument("mock SPPARKS communicator is null");
    }
    if (!is_finite_point(global_lower) || !is_finite_point(global_upper)
        || global_upper.x <= global_lower.x
        || global_upper.y <= global_lower.y
        || global_upper.z <= global_lower.z) {
        throw std::invalid_argument("mock SPPARKS global bounds are invalid");
    }
    const std::array<std::uint64_t, 3> process_counts{
        process_grid.x, process_grid.y, process_grid.z};
    for (const auto process_count : process_counts) {
        if (process_count == 0
            || process_count
                > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
            throw std::invalid_argument("mock SPPARKS process grid is invalid");
        }
    }

    check_mpi(MPI_Comm_rank(world, &me), "MPI_Comm_rank");
    check_mpi(MPI_Comm_size(world, &nprocs), "MPI_Comm_size");
    const auto product = process_grid.x * process_grid.y * process_grid.z;
    if (product != static_cast<std::uint64_t>(nprocs)) {
        throw std::invalid_argument(
            "mock SPPARKS process-grid product differs from communicator size");
    }

    procgrid[0] = static_cast<int>(process_grid.x);
    procgrid[1] = static_cast<int>(process_grid.y);
    procgrid[2] = static_cast<int>(process_grid.z);
    std::copy(std::begin(procgrid), std::end(procgrid), std::begin(user_procgrid));

    myloc[0] = me % procgrid[0];
    myloc[1] = (me / procgrid[0]) % procgrid[1];
    myloc[2] = me / (procgrid[0] * procgrid[1]);

    periodicity[0] = periodic.x ? 1 : 0;
    periodicity[1] = periodic.y ? 1 : 0;
    periodicity[2] = periodic.z ? 1 : 0;
    xperiodic = periodicity[0];
    yperiodic = periodicity[1];
    zperiodic = periodicity[2];
    nonperiodic = (periodic.x && periodic.y && periodic.z) ? 0 : 1;

    boxxlo = global_lower.x;
    boxxhi = global_upper.x;
    boxylo = global_lower.y;
    boxyhi = global_upper.y;
    boxzlo = global_lower.z;
    boxzhi = global_upper.z;
    xprd = boxxhi - boxxlo;
    yprd = boxyhi - boxylo;
    zprd = boxzhi - boxzlo;

    const std::array<double, 3> lower{boxxlo, boxylo, boxzlo};
    const std::array<double, 3> period{xprd, yprd, zprd};
    std::array<double, 3> sublo{};
    std::array<double, 3> subhi{};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const auto width = period[axis] / static_cast<double>(procgrid[axis]);
        sublo[axis] = lower[axis] + static_cast<double>(myloc[axis]) * width;
        subhi[axis] = sublo[axis] + width;

        for (int side = 0; side < 2; ++side) {
            VoxelCoord neighbor_location{
                static_cast<std::int64_t>(myloc[0]),
                static_cast<std::int64_t>(myloc[1]),
                static_cast<std::int64_t>(myloc[2])};
            auto* component = axis == 0
                ? &neighbor_location.x
                : (axis == 1 ? &neighbor_location.y : &neighbor_location.z);
            *component += side == 0 ? -1 : 1;
            bool valid = true;
            *component = normalized_location(
                *component,
                procgrid[axis],
                periodicity[axis] != 0,
                valid);
            procneigh[axis][side] = valid
                ? rank_for_location(neighbor_location)
                : MPI_PROC_NULL;
        }
    }
    subxlo = sublo[0];
    subxhi = subhi[0];
    subylo = sublo[1];
    subyhi = subhi[1];
    subzlo = sublo[2];
    subzhi = subhi[2];
}

bool MockSpparksDomain::owns(const Point3& position) const noexcept
{
    return position.x >= subxlo && position.x < subxhi
        && position.y >= subylo && position.y < subyhi
        && position.z >= subzlo && position.z < subzhi;
}

int MockSpparksDomain::rank_for_location(
    const VoxelCoord& process_location) const
{
    if (process_location.x < 0 || process_location.x >= procgrid[0]
        || process_location.y < 0 || process_location.y >= procgrid[1]
        || process_location.z < 0 || process_location.z >= procgrid[2]) {
        throw std::out_of_range("mock SPPARKS process location is out of range");
    }
    const auto rank = process_location.x
        + static_cast<std::int64_t>(procgrid[0])
            * (process_location.y
               + static_cast<std::int64_t>(procgrid[1]) * process_location.z);
    return static_cast<int>(rank);
}

MockSpparksApp::MockSpparksApp(const MockSpparksDomain& domain)
    : domain_(domain)
{
}

void MockSpparksApp::set_owned_atoms(
    std::vector<MockSpparksAtom> owned_atoms)
{
    if (owned_atoms.size()
        > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("mock SPPARKS owned atom count exceeds int");
    }
    for (const auto& atom : owned_atoms) {
        if (!is_finite_point(atom.position) || !std::isfinite(atom.radius)
            || atom.radius < 0.0) {
            throw std::invalid_argument("mock SPPARKS atom is invalid");
        }
        if (!domain_.owns(atom.position)) {
            throw std::invalid_argument(
                "mock SPPARKS owned atom lies outside its subdomain");
        }
    }
    std::sort(
        owned_atoms.begin(),
        owned_atoms.end(),
        [](const MockSpparksAtom& left, const MockSpparksAtom& right) {
            return left.id < right.id;
        });
    atoms_ = std::move(owned_atoms);
    nlocal = checked_int(atoms_.size(), "mock SPPARKS owned atom count");
    nghost = 0;
    refresh_public_arrays();
}

void MockSpparksApp::synchronize_ghost_atoms(double ghost_distance)
{
    if (!std::isfinite(ghost_distance) || ghost_distance < 0.0) {
        throw std::invalid_argument(
            "mock SPPARKS ghost distance must be finite and nonnegative");
    }
    const std::array<double, 3> subdomain_widths{
        domain_.subxhi - domain_.subxlo,
        domain_.subyhi - domain_.subylo,
        domain_.subzhi - domain_.subzlo};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (domain_.procgrid[axis] > 1
            && ghost_distance > subdomain_widths[axis]) {
            throw std::invalid_argument(
                "mock SPPARKS ghost distance spans more than one subdomain");
        }
    }
    atoms_.resize(static_cast<std::size_t>(nlocal));

    std::vector<std::vector<AtomWire>> send_by_rank(
        static_cast<std::size_t>(domain_.nprocs));
    for (const auto& atom : atoms_) {
        const auto x_offsets = axis_offsets(
            atom.position.x, domain_.subxlo, domain_.subxhi, ghost_distance);
        const auto y_offsets = axis_offsets(
            atom.position.y, domain_.subylo, domain_.subyhi, ghost_distance);
        const auto z_offsets = axis_offsets(
            atom.position.z, domain_.subzlo, domain_.subzhi, ghost_distance);
        std::set<int> target_ranks;
        for (const auto z_offset : z_offsets) {
            for (const auto y_offset : y_offsets) {
                for (const auto x_offset : x_offsets) {
                    if (x_offset == 0 && y_offset == 0 && z_offset == 0) {
                        continue;
                    }
                    bool valid = true;
                    VoxelCoord location{
                        normalized_location(
                            static_cast<std::int64_t>(domain_.myloc[0]) + x_offset,
                            domain_.procgrid[0],
                            domain_.periodicity[0] != 0,
                            valid),
                        normalized_location(
                            static_cast<std::int64_t>(domain_.myloc[1]) + y_offset,
                            domain_.procgrid[1],
                            domain_.periodicity[1] != 0,
                            valid),
                        normalized_location(
                            static_cast<std::int64_t>(domain_.myloc[2]) + z_offset,
                            domain_.procgrid[2],
                            domain_.periodicity[2] != 0,
                            valid)};
                    if (valid) {
                        const auto target_rank = domain_.rank_for_location(location);
                        if (target_rank != domain_.me) {
                            target_ranks.insert(target_rank);
                        }
                    }
                }
            }
        }
        const AtomWire wire{
            atom.id,
            atom.position.x,
            atom.position.y,
            atom.position.z,
            atom.radius};
        for (const auto target_rank : target_ranks) {
            send_by_rank[static_cast<std::size_t>(target_rank)].push_back(wire);
        }
    }

    std::vector<int> send_counts(static_cast<std::size_t>(domain_.nprocs), 0);
    std::vector<int> receive_counts(static_cast<std::size_t>(domain_.nprocs), 0);
    for (int rank = 0; rank < domain_.nprocs; ++rank) {
        send_counts[static_cast<std::size_t>(rank)] = checked_int(
            send_by_rank[static_cast<std::size_t>(rank)].size(),
            "mock SPPARKS ghost send count");
    }
    check_mpi(
        MPI_Alltoall(
            send_counts.data(), 1, MPI_INT,
            receive_counts.data(), 1, MPI_INT,
            domain_.world),
        "MPI_Alltoall");

    std::vector<int> send_displacements(send_counts.size(), 0);
    std::vector<int> receive_displacements(receive_counts.size(), 0);
    int send_total = 0;
    int receive_total = 0;
    for (int rank = 0; rank < domain_.nprocs; ++rank) {
        const auto index = static_cast<std::size_t>(rank);
        send_displacements[index] = send_total;
        receive_displacements[index] = receive_total;
        if (send_counts[index] > std::numeric_limits<int>::max() - send_total
            || receive_counts[index]
                > std::numeric_limits<int>::max() - receive_total) {
            throw std::overflow_error("mock SPPARKS ghost exchange exceeds int");
        }
        send_total += send_counts[index];
        receive_total += receive_counts[index];
    }

    std::vector<AtomWire> send_buffer;
    send_buffer.reserve(static_cast<std::size_t>(send_total));
    for (const auto& rank_buffer : send_by_rank) {
        send_buffer.insert(
            send_buffer.end(), rank_buffer.begin(), rank_buffer.end());
    }
    std::vector<AtomWire> receive_buffer(
        static_cast<std::size_t>(receive_total));

    MPI_Datatype atom_wire_type = MPI_DATATYPE_NULL;
    check_mpi(
        MPI_Type_contiguous(
            checked_int(sizeof(AtomWire), "mock SPPARKS atom wire size"),
            MPI_BYTE,
            &atom_wire_type),
        "MPI_Type_contiguous");
    check_mpi(MPI_Type_commit(&atom_wire_type), "MPI_Type_commit");
    AtomWire dummy{};
    const auto exchange_error = MPI_Alltoallv(
        send_buffer.empty() ? &dummy : send_buffer.data(),
        send_counts.data(), send_displacements.data(), atom_wire_type,
        receive_buffer.empty() ? &dummy : receive_buffer.data(),
        receive_counts.data(), receive_displacements.data(), atom_wire_type,
        domain_.world);
    const auto free_error = MPI_Type_free(&atom_wire_type);
    check_mpi(exchange_error, "MPI_Alltoallv");
    check_mpi(free_error, "MPI_Type_free");

    std::sort(
        receive_buffer.begin(),
        receive_buffer.end(),
        [](const AtomWire& left, const AtomWire& right) {
            return left.id < right.id;
        });
    for (const auto& wire : receive_buffer) {
        atoms_.push_back({wire.id, {wire.x, wire.y, wire.z}, wire.radius});
    }
    nghost = checked_int(receive_buffer.size(), "mock SPPARKS ghost atom count");
    refresh_public_arrays();
}

void MockSpparksApp::refresh_public_arrays()
{
    id_storage_.resize(atoms_.size());
    xyz_storage_.resize(atoms_.size());
    xyz_rows_.resize(atoms_.size());
    radius_storage_.resize(atoms_.size());
    for (std::size_t atom_index = 0; atom_index < atoms_.size(); ++atom_index) {
        const auto& atom = atoms_[atom_index];
        id_storage_[atom_index] = atom.id;
        xyz_storage_[atom_index] = {
            atom.position.x, atom.position.y, atom.position.z};
        xyz_rows_[atom_index] = xyz_storage_[atom_index].data();
        radius_storage_[atom_index] = atom.radius;
    }
    id = id_storage_.empty() ? nullptr : id_storage_.data();
    xyz = xyz_rows_.empty() ? nullptr : xyz_rows_.data();
    radius = radius_storage_.empty() ? nullptr : radius_storage_.data();
}

}  // namespace gasaccess::testing
