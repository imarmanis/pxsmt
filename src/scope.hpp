#ifndef _SCOPE_HPP
#define _SCOPE_HPP

#include "graph.hpp"

class Scope
{
public:
    unsigned count;
    /* Tracks inserted (adjacency lists) graph edges */
    std::vector<Edge *> edges;
    /* Tracks TC pairs */
    std::vector<upair> tc_pairs;
    /* Tracks propagated A->F pair (persistency) */
    std::vector<npair> propagated_paths;
    /* Tracks fixed identifiers */
    std::vector<z3::expr> fixed;
    /* Tracks porf TC pairs */
    std::vector<upair> porf_tc_pairs;
    /* Tracks nodes that are porf-before a Rec */
    std::vector<Node*> porf_rec;
    /* Track crashed assume nodes */
    std::vector<Node*> crashed;
    /* Track po_cl pairs */
    std::vector<std::pair<Node*, Node*>> pocl;
    Scope();
    bool dead();
    bool empty();
    void clear();
    void incr_count();
    unsigned decr_count(unsigned);
};

#endif
