#ifndef _GRAPH_HPP
#define _GRAPH_HPP

#include "list"
#include "unordered_map"
#include "unordered_set"
#include "map"
#include "set"
#include "vector"
#include "llvm/IR/User.h"
#include "z3++.h"
#include "debug.hpp"
#include "expr_hasher.hpp"
#include "llvm/ADT/STLExtras.h"

enum BW_S_RES {
    BW_OK,
    BW_OOF,
    BW_CYCLE
};

enum NODE_KIND {
    P_READ,
    P_WRITE,
    EX_READ,
    EX_WRITE,
    SFENCE,
    MFENCE,
    FLUSH,
    FLUSH_OPT,
    REC_READ,
    ASSUME
};

class Propagator;
class Graph;
class Node;

enum EdgeKind {
    SPPO = 0,
    DPPO = 1,
    RF = 2,
    CO = 3,
    FR = 4,
    /* alive FL/FO and same cache-line as the Rec read */
    A_CL = 5,
    FB = 6,
    PO = 7
};

#define MAX_EDGE_COST 1
#define EDGE_K_TEMP SPPO

std::string to_string(EdgeKind);

class Edge
{
public:
    const EdgeKind kind;

private:
    Node *x;
    Node *y;
    const std::vector<z3::expr> reasons;

public:
    Edge(EdgeKind, Node *, Node *);
    Edge(EdgeKind, Node *, Node *, const z3::expr &);
    Edge(EdgeKind, Node *, Node *, const z3::expr &, const z3::expr &);
    Edge(EdgeKind, Node *, Node *, Edge *);
    unsigned getCost();
    friend std::ostream &operator<<(std::ostream &, const Edge &);
    friend class Graph;
    friend class Propagator;
};

extern std::pair<Edge*, unsigned> not_desc;

typedef std::list<Edge*> path;
typedef std::list<z3::expr> reasons;
typedef std::function<void(path&)> path_cb;
typedef std::function<void(reasons&)> reasons_cb;

typedef std::function<void(path&, const z3::expr &)> extra_path_cb;

typedef std::tuple<unsigned, unsigned> upair;
typedef std::pair<Node*, Node*> npair;
typedef std::vector<upair> graph_scope;
typedef std::vector<npair> path_scope;

struct np_cmp;

class Node
{
private:
    /* ICD specific */
    unsigned level = 1;
    /* All PPO U OBS outgoing edges */
    std::list<Edge *> out;
    /* All PPO U OBS incoming edges */
    std::list<Edge *> in;
    /* Incoming edges (in) that have same-level endpoints */
    std::list<Edge *> icd_in;
    unsigned dist;
    /* Outgoing rfi and coi edges, needed for fr derivation */
    std::list<Edge *> out_i;
    /* Outgoing po U rf edges */
    std::list<Edge *> porf_out;
    Node(z3::expr porf_clock, z3::expr ob_block, unsigned nid, llvm::User *user,
         z3::expr guard, z3::expr loc, z3::expr val, NODE_KIND kind,
         z3::expr crash, z3::expr rf_w_id);
    /* Temporary, used to find porf path */
    Edge * parent;

public:
    bool is_plain_read();
    bool is_plain_write();
    bool is_mfence();
    bool is_sfence();
    bool is_flush();
    bool is_flush_opt();
    bool is_rec();
    bool is_rex();
    bool is_wex();
    bool is_fo_helper();
    bool is_assume();
    bool has_location();
    const unsigned nid;
    llvm::User *user;
    z3::expr guard;
    /* For assume nodes, loc holds the np expr */
    z3::expr loc;
    z3::expr val;
    z3::expr porf_clock;
    z3::expr ob_clock;
    NODE_KIND kind;
    z3::expr crash;
    /* Used for rec nodes */
    z3::expr rf_w_id;
    #ifdef DEBUG
    bool hasOut(Edge *);
    bool hasOut_i(Edge *);
    #endif
    void log_edges();
    friend std::ostream &operator<<(std::ostream &, const Node &);
    friend class Graph;
    friend class Propagator;
    friend struct np_cmp;
    void print(z3::model m);
    void print();
};

struct np_cmp {
    bool operator() (Node *x, Node *y) const {
        return x->nid < y->nid;
    }
};


class Graph
{
private:
    unsigned edges;
    void register_edge(EdgeKind, Node *, Node *, const z3::expr &);

public:
    std::unordered_set<z3::expr> *fixed;
    std::unordered_map<std::pair<Node*, Node*>, z3::expr, llvm::pair_hash<Node *, Node*>> *tracked;
    std::pair<Edge *, unsigned> **desc = nullptr;
    bool reachable(Node *, Node *);
    void meld(std::list<std::tuple<Node*, Node*, z3::expr>> &, Node *, Node *, Edge *, graph_scope *, extra_path_cb,
        unsigned);
    void process_paths(Node *, Node *, path_cb, bool **);
    void find_short_paths(unsigned &, bool, Node *, Node *, bool **, path *, path_cb);
    #ifdef DEBUG
    bool verifyEdges();
    bool find_node(Node *target, Node *curr, std::set<Node*> &v);
    bool visit(Node *, std::map<Node*, bool> *);
    bool acyclic();
    bool contains(Edge *);
    #endif

    std::unordered_map<z3::expr, Edge*> edge_from_id;
    /* debug : to print path found */
    std::map<Node*, Edge*> prev;
    Graph();
    std::vector<Node *> nodes;
    unsigned num_nodes();
    unsigned num_edges();
    void add_node(Node *);
    void register_dppoi_edge(Node *, Node*, const z3::expr &);
    void register_rf_edge(Node *, Node*, const z3::expr &);
    void register_co_edge(Node *, Node*, const z3::expr &);
    void register_fb_edge(Node *, Node*, const z3::expr &);
    void register_cl_edge(Node *, Node*, const z3::expr &);
    void insert_edge(Edge *, graph_scope *, extra_path_cb);
    void remove_edge(Edge *);
    Edge* get_tracked_edge(const z3::expr &);
    BW_S_RES backward_search(std::unordered_set<Node*> &, Edge *);
    bool forward_search(std::unordered_set<Node*> &, Edge *);
    unsigned delta();
};

bool compare_nodep (Node*&, Node*&);

#endif
