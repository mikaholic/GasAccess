#include "gasaccess/spparks_adapter.hpp"

#include "app.h"
#include "domain.h"

#include <cstddef>

namespace {

// Compilation of this function proves that the adapter consumes the public
// field contract of the actual SPPARKS Domain and App base classes. It is not
// executed because those objects require a complete SPPARKS runtime.
void compile_real_spparks_contract(
    const SPPARKS_NS::Domain& domain,
    const SPPARKS_NS::App& app,
    MPI_Comm communicator)
{
    gasaccess::SpparksAdapterConfig config{};
    config.requested_spacing = {1.0, 1.0, 1.0};
    const auto grid_spec = gasaccess::make_spparks_grid_spec(domain, config);
    const auto decomposition_spec = gasaccess::make_spparks_decomposition_spec(
        domain, communicator, config);
    gasaccess::SpparksAtomBuffer atom_buffer;
    atom_buffer.assign(
        app,
        [](std::size_t) {
            return 0.0;
        });
    (void)grid_spec;
    (void)decomposition_spec;
    (void)atom_buffer;
}

}  // namespace

int main()
{
    (void)&compile_real_spparks_contract;
    return 0;
}
