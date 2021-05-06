#include "node_entry.h"
#include "node_entry_view.h"
#include <iostream>
#include <cmath>

using namespace xpas;

node_entry::node_entry(std::string _id, vector_type columns)
    : _branch_label{ std::move(_id) }
    , _columns{std::move(columns) }
{}

void node_entry::add_column(column_type column)
{
    _columns.push_back(column);
}

size_t node_entry::get_alignment_size() const
{
    return _columns.size();
}

std::string node_entry::get_label() const
{
    return _branch_label;
}

const proba_pair& node_entry::at(size_t column, size_t row) const
{
    return _columns[column][row];
}

double entropy(const xpas::column_type& array)
{
    double entropy = 0.0;
    for (const auto& proba_pair : array)
    {
        entropy -= proba_pair.score * std::log2(proba_pair.score);
    }
    return entropy;
}

void node_entry::calculate_entropies()
{
    for (const auto& log_column : _columns)
    {
        auto column = log_column;
        for (size_t i = 0; i < log_column.size(); ++i)
        {
            column[i].score = static_cast<xpas::phylo_kmer::score_type>(std::pow(10, log_column[i].score));
        }

        _entropies.push_back(entropy(column));
    }
}

double node_entry::get_column_entropy(size_t column_idx) const
{
    return _entropies[column_idx];
}

bool operator==(const node_entry& lhs, const node_entry& rhs)
{
    return lhs.get_label() == rhs.get_label();
}

using namespace impl;

node_entry_view make_empty_view()
{
    return node_entry_view{ nullptr, 0.0, 0, 0 };
}

chain_window_iterator::chain_window_iterator(node_entry_view view,
                                             size_t kmer_size, phylo_kmer::score_type threshold) noexcept
    : _view{ std::move(view) }
    , _kmer_size{ kmer_size }
    , _threshold{ threshold }
    , _first_view_pos{ 0 }
{
}

chain_window_iterator& chain_window_iterator::operator++()
{
    const auto entry = _view.get_entry();

    /// Prefix size is the suffix size from the previous window,
    /// see dac_kmer_iterator::_next_phylokmer
    const auto prefix_size = _view.get_prefix_size();
    /// Suffix size for the next window
    const auto suffix_size = _kmer_size - prefix_size;

    /// in the branch-and-bound, there are no saved suffixes, and we
    /// chain windows just iteratively
    if (prefix_size == 0)
    {
        if (size_t(_view.get_end_pos() + 1) < entry->get_alignment_size())
        {
            _view.set_start_pos(_view.get_start_pos() + 1);
            _view.set_end_pos(_view.get_end_pos() + 1);
        }
        else
        {
            _view = make_empty_view();
        }
    }
    /// in the divide-and-conquer, we chain windows in a special way
    else
    {
        /// The start position of the last chain
        const auto last_chain_pos = (_kmer_size % 2) ? _kmer_size / 2 : _kmer_size / 2 - 1;

        /// continue the chain if possible
        if (size_t(_view.get_end_pos() + suffix_size) < entry->get_alignment_size())
        {
            _view.set_start_pos(_view.get_start_pos() + suffix_size);
            _view.set_end_pos(_view.get_end_pos() + suffix_size);
        }
        /// if the chain is over, start the next one if possible
        else if (_first_view_pos + 1 <= last_chain_pos)
        {
            ++_first_view_pos;

            _view.set_start_pos(_first_view_pos);
            _view.set_end_pos(_first_view_pos + _kmer_size - 1);
            _view.set_prefixes({});
            _view.set_prefix_size(_kmer_size / 2);
        }
        /// otherwise, the iterator is over
        else
        {
            _view = make_empty_view();
        }
    }

    return *this;
}

bool chain_window_iterator::operator==(const chain_window_iterator& rhs) const noexcept
{
    return _view == rhs._view;
}

bool chain_window_iterator::operator!=(const chain_window_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

chain_window_iterator::reference chain_window_iterator::operator*() noexcept
{
    return _view;
}

chain_windows::chain_windows(const node_entry& entry, size_t kmer_size, xpas::phylo_kmer::score_type threshold)
    : _entry{ entry }
    , _kmer_size{ kmer_size }
    , _start_pos{ 0 }
    , _threshold{ threshold }
{}

chain_windows::const_iterator chain_windows::begin() const
{
    auto first_window_view = node_entry_view(&_entry, _threshold, 0, _kmer_size - 1);
    return { std::move(first_window_view), _kmer_size, _threshold };
}

chain_windows::const_iterator chain_windows::end() const noexcept
{
    return { make_empty_view(), 0, 0.0 };
}