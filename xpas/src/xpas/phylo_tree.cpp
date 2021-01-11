#include <xpas/phylo_kmer.h>
#include <xpas/phylo_tree.h>
#include <xpas/newick.h>
#include <boost/filesystem.hpp>
#include <algorithm>

namespace fs = boost::filesystem;
using namespace xpas;
using namespace xpas::impl;
using std::vector;
using std::string;
using std::move;
using std::begin, std::end;


phylo_tree::phylo_tree(phylo_node* root)
    : _root{ root }, _node_count{ 0 }
{
    index();
}

phylo_tree::phylo_tree(phylo_tree&& other) noexcept
{
    _root = other._root;
    other._root = nullptr;

    _node_count = other._node_count;
    other._node_count = 0;

    _preorder_id_to_node = std::move(other._preorder_id_to_node);
    _postorder_id_node_mapping = std::move(other._postorder_id_node_mapping);
    _label_to_node = std::move(other._label_to_node);
}

phylo_tree::~phylo_tree() noexcept
{
    delete _root;
}

phylo_tree::const_iterator xpas::phylo_tree::begin() const noexcept
{
    return visit_subtree(_root).begin();
}

phylo_tree::const_iterator xpas::phylo_tree::end() const noexcept
{
    return visit_subtree(_root).end();
}

phylo_tree::iterator xpas::phylo_tree::begin() noexcept
{
    return visit_subtree<postorder_tree_iterator<false>>(_root).begin();
}

phylo_tree::iterator xpas::phylo_tree::end() noexcept
{
    return visit_subtree<postorder_tree_iterator<false>>(_root).end();
}

size_t phylo_tree::get_node_count() const noexcept
{
    return _node_count;
}

phylo_tree::value_pointer phylo_tree::get_root() const noexcept
{
    return _root;
}

void phylo_tree::set_root(value_pointer root)
{
    _root = root;
}

bool phylo_tree::is_rooted() const noexcept
{
    return _root && _root->get_children().size() < 3;
}

void phylo_tree::index()
{
    if (_root->get_parent())
    {
        throw std::invalid_argument{ "Can not create a tree from non-root node: "
                                     "the parent of the root must be nullptr."};
    }

    /// Recreate all search maps
    _index_preorder_id();
    _index_postorder_id();
    _index_labels();

    _count_nodes();
}

optional<const phylo_node*> phylo_tree::get_by_preorder_id(phylo_node::id_type preorder_id) const noexcept
{
    if (const auto it = _preorder_id_to_node.find(preorder_id); it != _postorder_id_node_mapping.end())
    {
        return { it->second };
    }
    else
    {
        return { nullopt };
    }
}

optional<const phylo_node*> phylo_tree::get_by_postorder_id(phylo_node::id_type postorder_id) const noexcept
{
    if (const auto it = _postorder_id_node_mapping.find(postorder_id); it != _postorder_id_node_mapping.end())
    {
        return { it->second };
    }
    else
    {
        return { nullopt };
    }
}

optional<const xpas::phylo_node*> phylo_tree::get_by_label(const std::string& label) const noexcept
{
    if (const auto it = _label_to_node.find(label); it != _label_to_node.end())
    {
        return { it->second };
    }
    else
    {
        return { nullopt };
    }
}

void phylo_tree::_index_preorder_id()
{
    _preorder_id_to_node.clear();

    phylo_node::id_type preorder_id = 0;
    for (auto& node : xpas::visit_subtree<preorder_tree_iterator<false>>(_root))
    {
        node._preorder_id = preorder_id;
        _preorder_id_to_node[preorder_id] = &node;
        ++preorder_id;
    }
}

void phylo_tree::_index_postorder_id()
{
    _postorder_id_node_mapping.clear();

    phylo_node::id_type postorder_id = 0;
    for (auto& node : xpas::visit_subtree<postorder_tree_iterator<false>>(_root))
    {
        node._postorder_id = postorder_id;
        _postorder_id_node_mapping[postorder_id] = &node;
        ++postorder_id;
    }
}

void phylo_tree::_index_labels()
{
    _label_to_node.clear();

    for (auto& node : xpas::visit_subtree(_root))
    {
        _label_to_node[node.get_label()] = &node;
    }
}

void phylo_tree::_count_nodes()
{
    _node_count = 0;

    for (auto& node : xpas::visit_subtree(_root))
    {
        ++_node_count;
    }
}

namespace xpas
{
    void save_tree(const phylo_tree& tree, const std::string& filename)
    {
        std::ofstream out(filename);
        out << xpas::io::to_newick(tree);
    }
}