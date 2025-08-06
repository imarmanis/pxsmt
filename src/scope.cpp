#include "scope.hpp"

Scope::Scope() : count(1) {};

bool Scope::dead() { return (count == 0); }

bool Scope::empty() {
    return edges.empty() && tc_pairs.empty() && propagated_paths.empty()
        && fixed.empty() && porf_tc_pairs.empty() && porf_rec.empty()
        && crashed.empty();// && pocl.empty();
};

void Scope::incr_count()
{
    SASSERT(empty());
    count++;
};

/* Returns how many scopes were removed */
unsigned Scope::decr_count(unsigned d)
{
    SASSERT(empty());
    if (d > count) {
        unsigned decr = count;
        count = 0;
        return decr;
    }
    count -= d;
    return d;
}

void Scope::clear()
{
    edges.clear();
    tc_pairs.clear();
    propagated_paths.clear();
    fixed.clear();
    porf_tc_pairs.clear();
    porf_rec.clear();
    crashed.clear();
    //pocl.clear();
}
