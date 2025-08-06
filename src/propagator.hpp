#ifndef _PROPAGATOR_HPP
#define _PROPAGATOR_HPP

#include "stack"
#include "vector"
#include "tuple"
#include "unordered_set"
#include "unordered_map"
#include "z3++.h"
#include "graph.hpp"
#include "scope.hpp"
#include "debug.hpp"
#include "expr_hasher.hpp"
#include "llvm/ADT/STLExtras.h"

using namespace z3;

class Propagator : public user_propagator_base
{
private:
    Graph g;
    std::list<Scope> scopes;
    bool pop_pending = false;
    bool after_final = false;
    std::unordered_set<z3::expr> fixed_exprs;
    //std::map<unsigned, z3::expr> expr_from_id;
    /* Tracked edge between two nodes to edge id */
    /* FIXME : multimap, not map */
    std::unordered_map<std::pair<Node*, Node*>, z3::expr, llvm::pair_hash<Node *, Node*>> tracked;
    /* Tracks nodes that are porf before a Rec */
    std::set<Node*, np_cmp> porf_rec;
    /* Crashed nodes */
    std::unordered_map<z3::expr, Node*> cnode_from_id;
    std::map<Node*, z3::expr, np_cmp> crashed;
    /* track crash propagations : id of node -> [[id of reasons]] */
    std::multimap<Node *, std::list<z3::expr>, np_cmp> c_propagations;
    /* track dom(rf;[REC]) and rf */
    std::map<Node*, Edge*, np_cmp> recovered_writes;
    /* track rng([REC];fr) and fr */
    //std::map<Node*, Edge*, np_cmp> lost_writes;
    /* Node to pending in edges */
    std::unordered_multimap<Node *, Edge *> pendings_in;
    /* Node to pending out edges */
    std::unordered_multimap<Node *, Edge *> pendings_out;
    #ifdef DEBUG
    bool is_fixed(const z3::expr &);
    bool is_fixed(Edge *);
    bool verifyStack();
    #endif
    unsigned insertions = 0;
    unsigned path_conflicts = 0;
    unsigned consequences = 0;
    unsigned tracked_variables = 0;

    /*
     * Maps write nodes to incoming coe/fre edges.
     * Needed for coh_r{r,w}, and atomic axioms.
     */
    std::multimap<Node *, Edge *, np_cmp> inc_obsw;
    std::unordered_set<std::pair<Node *, Node*>, llvm::pair_hash<Node *, Node *>> rmws;
    /* po : (static) and transitive */
    bool ** po = nullptr;
    /* sppo : static and transitive */
    bool ** sppo = nullptr;
    /* (po U rf)+ : ITC */
    bool ** porf = nullptr;
    std::vector<std::pair<Node *, Node*>> poi;
    std::vector<std::pair<Node *, Node*>> sppoi;
    bool new_edge(Edge *);
    std::vector<Edge*> derive_frs(Edge *e);
    void propagate_conflict(std::list<Edge*>*, std::map<Node*, z3::expr> &);
    void propagate_consequence(std::list<Edge*>*, z3::expr);
    bool future_read(Edge *);
    bool coh_rr(Edge *);
    bool coh_rw(Edge *);
    bool coh_wr(Edge *);
    bool coh_ww(Edge *);
    bool atomicity(Edge *);
    Node *new_node(llvm::User *, z3::expr, z3::expr, z3::expr, NODE_KIND);
    /* Add rf edge to porf relation, return set of new FL / FO_helper node */
    std::list<Node*> porf_rf_add(Edge *);
    bool porf_reachable(Node *, Node*);
    void porf_meld(std::list<Node*> &, Node *, Node *);
    void process_porf_rec_paths(Node *, reasons_cb);
    void process_porf_path(Node *x, Node *y, reasons_cb cb);
    void porf_paths(std::list<z3::expr> &, bool, Node *, Node *, reasons_cb);
    bool redundant_cprop(Node *, const reasons&);

public:
    std::vector<z3::expr> *same_loc;
    std::vector<z3::expr> *same_cl;
    std::vector<z3::expr> *ob_lt;
    void emitStats();
    void conflict_wrapper(const z3::expr_vector &a, std::string);
    Propagator(solver *s);
    std::vector<Node *> *get_nodes();
    solver *s;
    std::set<Node *, np_cmp> reads;
    std::set<Node *, np_cmp> writes;
    std::set<Node *, np_cmp> recs;
    /* All FL and FO node */
    std::list<Node *> flushes;
    /* Maps FO node to it's po-later FO_helpers */
    std::multimap<Node *, Node *, np_cmp> fo_helpers;
    std::set<Node *, np_cmp> fo_helps;
    void push() override;
    void pop(unsigned int) override;
    void fixed(expr const &, expr const &) override;
    bool process_edge(Edge *);
    void decide(expr &, unsigned &, Z3_lbool &) override;
    user_propagator_base *fresh(z3::context &ctx) override { assert(false); }
    Node *read_node(llvm::User *, z3::expr, z3::expr, z3::expr);
    Node *read_ex_node(llvm::User *, z3::expr, z3::expr, z3::expr);
    Node *write_node(llvm::User *, z3::expr, z3::expr, z3::expr);
    Node *write_ex_node(llvm::User *, z3::expr, z3::expr, z3::expr);
    Node *flush_node(llvm::User *, z3::expr, z3::expr);
    Node *flush_opt_node(llvm::User *, z3::expr, z3::expr);
    Node *rec_node(llvm::User *, z3::expr, z3::expr, z3::expr);
    Node *mfence_node(llvm::User *, z3::expr);
    Node *sfence_node(llvm::User *, z3::expr);
    Node *assume_node(llvm::User *, z3::expr, z3::expr);
    void track_rf(Node *, Node *, const expr &);
    void track_co(Node *, Node *, const expr &);
    void track_fb(Node *, Node *, const expr &);
    void track_dppoi(Node *, Node *, const expr &);
    void track_crashed(Node *, const expr &);
    void track(const expr &);
    void add_sppoi(Node *, Node *);
    void add_poi(Node *, Node *);
    void add_rmw(Node *, Node *);
    bool in_po(Node *, Node *);
    bool external(Node *, Node *);
    bool is_rmw_pair(Node *, Node *);
    /* Compute the TC of static relations, as well as initial porf (= po) */
    void static_tc();
    void finalize_nodes();

};

#endif
