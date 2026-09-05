#include "gasaccess/atom_change_event_buffer.hpp"

#include <limits>
#include <stdexcept>

namespace gasaccess {
namespace {

void validate_view(AtomView atom_view)
{
    if (atom_view.count != 0 && atom_view.atoms == nullptr) {
        throw std::invalid_argument(
            "atom-change event pointer is null with nonzero count");
    }
}

std::size_t checked_size_sum(std::size_t lhs, std::size_t rhs)
{
    if (rhs > std::numeric_limits<std::size_t>::max() - lhs) {
        throw std::length_error("atom-change event count exceeds size_t");
    }
    return lhs + rhs;
}

}  // namespace

void AtomChangeEventBuffer::clear() noexcept
{
    added_atoms_.clear();
    removed_atoms_.clear();
}

bool AtomChangeEventBuffer::empty() const noexcept
{
    return added_atoms_.empty() && removed_atoms_.empty();
}

std::size_t AtomChangeEventBuffer::added_atom_count() const noexcept
{
    return added_atoms_.size();
}

std::size_t AtomChangeEventBuffer::removed_atom_count() const noexcept
{
    return removed_atoms_.size();
}

void AtomChangeEventBuffer::reserve(
    std::size_t added_atom_capacity,
    std::size_t removed_atom_capacity)
{
    added_atoms_.reserve(added_atom_capacity);
    removed_atoms_.reserve(removed_atom_capacity);
}

void AtomChangeEventBuffer::record_deposition(const Atom& added_atom)
{
    added_atoms_.push_back(added_atom);
}

void AtomChangeEventBuffer::record_desorption(const Atom& removed_atom)
{
    removed_atoms_.push_back(removed_atom);
}

void AtomChangeEventBuffer::record_move(
    const Atom& old_atom,
    const Atom& new_atom)
{
    const auto next_added_size = checked_size_sum(added_atoms_.size(), 1);
    const auto next_removed_size = checked_size_sum(removed_atoms_.size(), 1);
    added_atoms_.reserve(next_added_size);
    removed_atoms_.reserve(next_removed_size);
    removed_atoms_.push_back(old_atom);
    added_atoms_.push_back(new_atom);
}

void AtomChangeEventBuffer::append(const AtomChangeBatch& atom_changes)
{
    validate_view(atom_changes.added_atoms);
    validate_view(atom_changes.removed_atoms);
    const auto next_added_size = checked_size_sum(
        added_atoms_.size(),
        atom_changes.added_atoms.count);
    const auto next_removed_size = checked_size_sum(
        removed_atoms_.size(),
        atom_changes.removed_atoms.count);
    added_atoms_.reserve(next_added_size);
    removed_atoms_.reserve(next_removed_size);
    if (atom_changes.added_atoms.count != 0) {
        added_atoms_.insert(
            added_atoms_.end(),
            atom_changes.added_atoms.atoms,
            atom_changes.added_atoms.atoms
                + atom_changes.added_atoms.count);
    }
    if (atom_changes.removed_atoms.count != 0) {
        removed_atoms_.insert(
            removed_atoms_.end(),
            atom_changes.removed_atoms.atoms,
            atom_changes.removed_atoms.atoms
                + atom_changes.removed_atoms.count);
    }
}

AtomView AtomChangeEventBuffer::added_atoms() const noexcept
{
    return {added_atoms_.data(), added_atoms_.size()};
}

AtomView AtomChangeEventBuffer::removed_atoms() const noexcept
{
    return {removed_atoms_.data(), removed_atoms_.size()};
}

AtomChangeBatch AtomChangeEventBuffer::atom_changes() const noexcept
{
    return {added_atoms(), removed_atoms()};
}

}  // namespace gasaccess
