#ifndef GASACCESS_ATOM_CHANGE_EVENT_BUFFER_HPP
#define GASACCESS_ATOM_CHANGE_EVENT_BUFFER_HPP

#include "gasaccess/atom_voxelizer.hpp"

#include <cstddef>
#include <vector>

namespace gasaccess {

// Owns synchronized KMC event records until apply_atom_changes() returns.
// Removed atoms are retained at their old positions and radii. Views returned
// by atom_changes() remain valid until the buffer is modified or destroyed.
class AtomChangeEventBuffer {
public:
    void clear() noexcept;
    bool empty() const noexcept;
    std::size_t added_atom_count() const noexcept;
    std::size_t removed_atom_count() const noexcept;

    void reserve(
        std::size_t added_atom_capacity,
        std::size_t removed_atom_capacity);
    void record_adsorption(const Atom& added_atom);
    void record_desorption(const Atom& removed_atom);
    void record_move(const Atom& old_atom, const Atom& new_atom);
    void append(const AtomChangeBatch& atom_changes);

    AtomView added_atoms() const noexcept;
    AtomView removed_atoms() const noexcept;
    AtomChangeBatch atom_changes() const noexcept;

private:
    std::vector<Atom> added_atoms_{};
    std::vector<Atom> removed_atoms_{};
};

}  // namespace gasaccess

#endif
