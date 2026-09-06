#ifndef GASACCESS_TOOLS_MOCK_SPPARKS_HPP
#define GASACCESS_TOOLS_MOCK_SPPARKS_HPP

#include "gasaccess/gas_grid.hpp"

#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX 1
#endif
#include <mpi.h>

#include <array>
#include <cstdint>
#include <vector>

namespace gasaccess::testing {

struct MockSpparksAtom {
    std::int64_t id = 0;
    Point3 position{};
    double radius = 0.0;
};

class MockSpparksDomain {
public:
    MockSpparksDomain(
        MPI_Comm communicator,
        const Point3& global_lower,
        const Point3& global_upper,
        const GridDimensions& process_grid,
        const PeriodicAxes& periodic);

    bool owns(const Point3& position) const noexcept;
    int rank_for_location(const VoxelCoord& process_location) const;

    MPI_Comm world = MPI_COMM_NULL;
    int me = 0;
    int nprocs = 0;
    int procgrid[3]{};
    int user_procgrid[3]{};
    int myloc[3]{};
    int procneigh[3][2]{};

    int box_exist = 1;
    int dimension = 3;
    int nonperiodic = 0;
    int xperiodic = 1;
    int yperiodic = 1;
    int zperiodic = 1;
    int periodicity[3]{1, 1, 1};

    double xprd = 0.0;
    double yprd = 0.0;
    double zprd = 0.0;
    double boxxlo = 0.0;
    double boxxhi = 0.0;
    double boxylo = 0.0;
    double boxyhi = 0.0;
    double boxzlo = 0.0;
    double boxzhi = 0.0;
    double subxlo = 0.0;
    double subxhi = 0.0;
    double subylo = 0.0;
    double subyhi = 0.0;
    double subzlo = 0.0;
    double subzhi = 0.0;
};

class MockSpparksApp {
public:
    explicit MockSpparksApp(const MockSpparksDomain& domain);

    MockSpparksApp(const MockSpparksApp&) = delete;
    MockSpparksApp& operator=(const MockSpparksApp&) = delete;

    void set_owned_atoms(std::vector<MockSpparksAtom> owned_atoms);
    void synchronize_ghost_atoms(double ghost_distance);

    int nlocal = 0;
    int nghost = 0;
    std::int64_t* id = nullptr;
    double** xyz = nullptr;
    double* radius = nullptr;

private:
    void refresh_public_arrays();

    const MockSpparksDomain& domain_;
    std::vector<MockSpparksAtom> atoms_{};
    std::vector<std::int64_t> id_storage_{};
    std::vector<std::array<double, 3>> xyz_storage_{};
    std::vector<double*> xyz_rows_{};
    std::vector<double> radius_storage_{};
};

}  // namespace gasaccess::testing

#endif
