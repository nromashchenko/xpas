#ifndef RAPPAS_CPP_NODE_ENTRY_H
#define RAPPAS_CPP_NODE_ENTRY_H

#include <vector>
#include "node_entry_view.h"

namespace xpas
{
    class view_iterator;

    /// \brief A submatrix of posterior probabilities matrix (fixed branch, all the positions of input alignment)
    class node_entry final
    {
    public:
        using const_iterator = view_iterator;
        using vector_type = std::vector<xpas::row_type>;

        explicit node_entry() noexcept = default;
        node_entry(std::string _id, vector_type&& rows);
        node_entry(const node_entry&) = delete;
        node_entry(node_entry&&) = default;
        node_entry& operator=(const node_entry&) = delete;
        node_entry& operator=(node_entry&&) = default;
        ~node_entry() noexcept = default;

        const_iterator begin(size_t kmer_size, xpas::phylo_kmer::score_type threshold) const;
        const_iterator end() const;

        void push_back(xpas::row_type&& row);

        size_t get_alignment_size() const;
        std::string get_label() const;

        const xpas::proba_pair& at(size_t position, size_t variant) const;

    private:
        std::string _branch_label;
        vector_type _rows;
    };

    class view_iterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using const_reference = const node_entry_view&;
        using const_pointer = const node_entry_view*;

        view_iterator(node_entry_view view) noexcept;
        view_iterator(const view_iterator& view) = delete;
        view_iterator(view_iterator&& view) = delete;
        view_iterator& operator=(const view_iterator&) = delete;
        view_iterator& operator=(view_iterator&&) = delete;
        ~view_iterator() = default;

        view_iterator& operator++();

        bool operator==(const view_iterator& rhs) const noexcept;
        bool operator!=(const view_iterator& rhs) const noexcept;
        const_reference operator*() const noexcept;
        const_pointer operator->() const noexcept;
    private:
        node_entry_view _view;
    };

    bool operator==(const node_entry& lhs, const node_entry& rhs);
}


#endif //RAPPAS_CPP_NODE_ENTRY_H
