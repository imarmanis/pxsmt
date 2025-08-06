#include "string"
#include "algorithm"
#include "unordered_set"
#include "set"
#include "map"
#include "tuple"
#include "stack"
#include "graph.hpp"
#include "config.hpp"
#include "EdgePriorityQueue.hpp"
#include "llvm/Support/Casting.h"
#include "llvm/IR/GlobalVariable.h"
#include "llvm/IR/Instructions.h"
#include <llvm/Support/raw_ostream.h>

bool compare_nodep (Node*& first, Node*& second)
{
    return first->nid < second->nid;
}

void Graph::add_node(Node *n) {
    nodes.push_back(n);
    std::string s;
    llvm::raw_string_ostream ss(s);
    ss << *n->user;
    LOG("node " << *n << " : " << ss.str() << "\n");
}

unsigned Graph::num_nodes() { return nodes.size(); }

unsigned Graph::num_edges() { return edges; }

void Graph::register_dppoi_edge(Node *x, Node*y, const z3::expr &reason)
{
    register_edge(DPPO, x, y , reason);
}

void Graph::register_rf_edge(Node *x, Node *y, const z3::expr &reason)
{
    register_edge(RF, x, y, reason);
}

void Graph::register_fb_edge(Node *x, Node *y, const z3::expr &reason)
{
    register_edge(FB, x, y, reason);
}

void Graph::register_co_edge(Node *x, Node *y, const z3::expr &reason)
{
    register_edge(CO, x, y, reason);
}

void Graph::register_cl_edge(Node *x, Node *y, const z3::expr &reason)
{
    assert(x->is_flush() || x->is_flush_opt());
    assert(y->is_rec());
    register_edge(A_CL, x, y, reason);
}

void Graph::register_edge(EdgeKind kind, Node *x, Node *y, const z3::expr &reason)
{
    Edge *e = new Edge(kind, x, y, reason);
    LOG("edge : " << reason << " = " << *e << "\n");
    edge_from_id.insert({reason,  e});
}

/* Just adds the edge, MUST be preceded by a call to reachable if ICD is used */
void Graph::insert_edge(Edge *e, graph_scope *scope, extra_path_cb ecb)
{
    Node *s = e->x;
    Node *t = e->y;

    if (algk == ITC) {
        /* Might not hold because of cycles containing recovery nodes */
        assert(!reachable(t, s));

        std::list<std::tuple<Node*, Node*, z3::expr>> pairs;
        if (!reachable(s, t))
            for(Node *v : nodes)
                if (reachable(v, s) &&
                !reachable(v, t))
                        meld(pairs, v, t, e, scope, ecb, max_pairs);

        DEBUG_CODE(LOG("n_tracked : " << pairs.size() << "\n"););
        for (auto p : pairs) {
            Node *x, *y;
            // FIXME
            std::tie(x, y, std::ignore) = p;
            z3::expr id = std::get<2>(p);
            Edge *not_fixed = edge_from_id.at(id);
            DEBUG_CODE(
                LOG("   while adding : " << *x << " -> " << *y << " to tc,\n");
                LOG("       found tracked but not fixed : " << id << " = " << *not_fixed << ",\n");
                LOG("       ");
            );
            path pa;
            Node *curr = y;
            do {
                Edge *ce = desc[x->nid][curr->nid].first;
                DEBUG_CODE(LOG(*ce<< ", "););
                pa.push_back(ce);
                curr = ce->x;
            } while (curr != x);
            DEBUG_CODE(LOG("\n"););
            assert(not_fixed->reasons.size() == 1);
            ecb(pa, not_fixed->reasons[0]);
        }
    } else {
        /* ICD */
        if (s->level == t->level) {
            assert(e->x->level == e->y->level);
            t->icd_in.push_back(e);
        }
    }
    s->out.push_back(e);
    t->in.push_back(e);
    edges++;
}

/* Add a -> b to tc, using edge h : a ->* c ->(e) b */
void Graph::meld(std::list<std::tuple<Node *, Node*, z3::expr>> &pairs, Node *a, Node *b, Edge *h, graph_scope *scope,
                     extra_path_cb ecbl, unsigned n_pairs)
{
    std::stack<std::tuple<Node*, Node*, Edge *>> W;
    W.push(std::make_tuple(a, b, h));

    while (!W.empty()) {
        Node *v, *t;
        Edge *par;
        auto tup = W.top();
        W.pop();
        std::tie(v, t, par) = tup;
        unsigned vs_dist = desc[v->nid][par->x->nid].second;
        assert(vs_dist != UINT_MAX);
        unsigned vt_dist = vs_dist + par->getCost();
        desc[v->nid][t->nid] = std::make_pair(par, vt_dist);

        if (n_pairs > 0) {
            auto it = tracked->find(std::make_pair(t, v));
            /* Check that there exists a t->v tracked edge (e.g. RF) */
            if (it != tracked->end()) {
                auto tv_tid = it->second;
                if ((vt_dist == 1) && (edge_from_id.at(tv_tid)->kind == RF)) {
                    n_pairs--;
                    pairs.push_back(std::make_tuple(v, t, tv_tid));
                }
            }
        }
        scope->push_back(std::make_pair(v->nid, t->nid));
        for (Edge *e : t->out)
            if (!reachable(v, e->y))
                W.push(std::make_tuple(v, e->y, e));
    }

    return;
}

void Graph::remove_edge(Edge *e)
{
    e->x->out.remove(e);
    e->y->in.remove(e);
    if ((algk == ICD) && (e->x->level == e->y->level)) {
        auto &l = e->y->icd_in;
        auto I = std::find(l.begin(), l.end(), e);
        assert(I != l.end());
        l.erase(I);
    }
    edges--;
}

Graph::Graph() : edges(0) {}

Edge* Graph::get_tracked_edge(const z3::expr &reason)
{
    return edge_from_id.at(reason);
}

Node::Node(z3::expr porf_clock, z3::expr ob_clock, unsigned nid, llvm::User* user,
    z3::expr guard, z3::expr loc, z3::expr val, NODE_KIND kind,
    z3::expr crash, z3::expr rf_w_id)
    : nid(nid), user(user), guard(guard), loc(loc), val(val),
    porf_clock(porf_clock), ob_clock(ob_clock), kind(kind),
    crash(crash), rf_w_id(rf_w_id) {}

Edge::Edge(EdgeKind kind, Node *x, Node *y)
    : kind(kind), x(x), y(y) {}

Edge::Edge(EdgeKind kind, Node *x, Node *y, const z3::expr &reason)
    : kind(kind), x(x), y(y), reasons({reason}) {}

Edge::Edge(EdgeKind kind, Node *x, Node *y, const z3::expr &reason1, const z3::expr &reason2)
    : kind(kind), x(x), y(y), reasons({reason1, reason2}) {}

/* Change MAX_EDGE_COST accordingly */
unsigned Edge::getCost()
{
    switch (kind) {
    case SPPO:
        return 0;

    default:
        return 1;
    }
}
std::ostream &operator<<(std::ostream &out, const Node &n)
{
    out << "#" << n.nid;
    return out;
}

std::ostream &operator<<(std::ostream &out, const Edge &e)
{
    out << to_string(e.kind) << " : [" << *e.x << " ~> " << *e.y << "]";
    return out;
}

std::string to_string(EdgeKind kind) {
    switch (kind) {
    case SPPO:
        return "SPPO";

    case DPPO:
        return "DPPO";

    case RF:
        return "RF";

    case CO:
        return "CO";

    case FR:
        return "FR";

    case A_CL:
        return "CL";

    case FB:
        return "FB";

    case PO:
        return "PO";
    default:
        assert(false);
    }
}

void Node::log_edges()
{
    TRACE_CODE(
        LOG("Node " << *this << ", edges:\n");
        if (!out.empty())
            LOG("out\n");
        for (auto *e : out)
            LOG("   " << *e << "\n");
        if (!out_i.empty())
            LOG("out_i\n");
        for (auto *e : out_i)
            LOG("   " << *e << "\n");
    );
}

void Graph::process_paths(Node *s, Node *t, path_cb cb, bool **sppo)
{
    assert(s != t);
    if (algk == ITC)
        assert(reachable(s, t));
    else
        /* ICD, avoid unnecessary reachable calls */
        SASSERT(reachable(s, t));

    for (auto *n : nodes)
        n->dist = UINT_MAX;
    s->dist = 0;

    EdgePriorityQueue W;
    std::unordered_set<Node*> V;
    for (Edge *e : s->out)
        W.push_back(e);

    while (!W.empty()) {
        Edge *curr_e = W.back();
        Node *curr = curr_e->y;
        auto dist = W.getCost();
        W.pop_back();

        /* Optimize path finding if ITC is used (fast queries) */
        if ((algk == ITC) && !reachable(curr, t))
            continue;

        if (V.find(curr) != V.end())
            continue;
        V.insert(curr);

        curr->dist = dist;
        if (curr == t)
            continue;

        for (Edge *out_e : curr->out)
            W.push_back(out_e);
    }

    assert(t->dist != UINT_MAX);
    path p;
    unsigned c = 0;
    find_short_paths(c, true, t, s, sppo, &p, cb);
}

void Graph::find_short_paths(unsigned &c, bool use_zero_cost, Node *curr, Node *target, bool **sppo, path *p, path_cb cb)
{
    if (c > max_paths)
        return;

    if (curr == target) {
        c++;
        cb(*p);
        return;
    }

    for (auto *e : curr->in) {
        /* Check if it's static ppo (zero cost) */
        if (e->kind == SPPO)
            continue;

        if (curr->dist == e->x->dist + e->getCost()) {
            p->push_back(e);
            find_short_paths(c, true, e->x, target, sppo, p, cb);
            assert(p->back() == e);
            p->pop_back();
        }
    }

    if (use_zero_cost) {
        for (auto *k : nodes) {
            if (curr == k)
                continue;
            if (sppo[k->nid][curr->nid] && (curr->dist == k->dist))
                find_short_paths(c, false, k, target, sppo, p, cb);
        }
    }
}


#ifdef DEBUG

bool Graph::contains(Edge *q)
{
    for (auto n : nodes) {
        for (auto *e : n->out)
            if (e == q)
                return true;
        for (auto *e : n->out_i)
            if (e == q)
                return true;
    }

    return false;
}

bool Node::hasOut(Edge *e)
{
    for (Edge *outg : out)
        if (outg == e)
            return true;
    return false;
}

bool Node::hasOut_i(Edge *e)
{
    for (Edge *outg : out_i)
        if (outg == e)
            return true;
    return false;
}

bool Graph::acyclic()
{

    std::set<Node*> unmarked;
    std::map<Node*, bool> marked;
    for (auto n : nodes)
        unmarked.insert(n);

    while (!unmarked.empty()) {
        auto n = *unmarked.begin();
        if (visit(n, &marked))
            return false;
        assert(marked.find(n) != marked.end());
        unmarked.erase(n);
    }

    return true;
}

/* return true if not DAG */
bool Graph::visit(Node *n, std::map<Node*, bool> *marked)
{
    auto m = marked->find(n);
    /* false : temporary, true : permanent */
    if (m != marked->end()) {
        if (!m->second)
            return true;
        return false;
    }

    marked->insert({n, false});

    for (auto *outg : n->out) {
        if (visit(outg->y, marked))
            return true;
    }

    marked->erase(n);
    marked->insert({n, true});
    return false;
}

bool Graph::verifyEdges()
{
    unsigned c = 0;
    for (auto *n : nodes) {
        for (auto *e : n->out) {
            if (e->x != n) return false;
            c++;
        }

        for (auto *e : n->out_i)
            if (e->x != n) return false;
    }

    return (c == edges);
}
#endif

std::pair<Edge*, unsigned> not_desc = std::make_pair(nullptr, UINT_MAX);

bool Node::is_plain_read()
{
    if (llvm::dyn_cast<llvm::GlobalVariable>(user))
        return false;
    auto *i = llvm::cast<llvm::Instruction>(user);
    return i->getOpcode() == llvm::Instruction::Load;
}

bool Node::is_plain_write()
{
    if (llvm::dyn_cast<llvm::GlobalVariable>(user))
        return true;
    auto *i = llvm::cast<llvm::Instruction>(user);
    return i->getOpcode() == llvm::Instruction::Store;
}

bool Node::is_mfence()
{
    return kind == MFENCE;
}

bool Node::is_sfence()
{
    return kind == SFENCE;
}

bool Node::is_flush()
{
    return kind == FLUSH;
}

bool Node::is_flush_opt()
{
    return kind == FLUSH_OPT;
}

/* Update accordingly when adding new event types */
bool Node::has_location()
{
    return !is_mfence() && !is_sfence() && !is_assume();
}

BW_S_RES Graph::backward_search(std::unordered_set<Node*> &B, Edge *e)
{
    Node *u = e->x;
    Node *target = e->y;
    /* First iteration doesn't count */
    unsigned fuel = 1 + (use_delta ? delta() : u->level);

    EdgePriorityQueue W;
    W.push_back(e);

    while(!W.empty() && (fuel > 0)) {
        fuel--;

        auto *curr_e = W.back();
        W.pop_back();
        Node *curr = curr_e->x;

        if (curr == target)
            return BW_CYCLE;

        /* Already visited */
        if (B.find(curr) != B.end())
            continue;

        B.insert(curr);

        for (Edge *inc_e : curr->icd_in) {
            SASSERT(inc_e->x->level == inc_e->y->level);
            W.push_back(inc_e);
        }
    }

    /* Not interrupted */
    if (W.empty())
        return BW_OK;

    /* Out of fuel */
    return BW_OOF;
}

/* Returns whether there exists a path towards a node of B */
bool Graph::forward_search(std::unordered_set<Node*> &B, Edge *e)
{
    Node *w = e->y;
    unsigned w_level = w->level;

    /* Instead of node u, use the edge e s.t. e->y = u */
    EdgePriorityQueue F;
    for (Edge *outg : w->out)
        F.push_back(outg);

    bool ret = false;
    while (!F.empty()) {
        Edge *curr_e = F.back();
        F.pop_back();
        Node *curr = curr_e->y;

        if (B.find(curr) != B.end())
            ret = true;

        if (curr->level == w_level) {
            curr->icd_in.push_back(curr_e);
            SASSERT(curr_e->x->level == curr_e->y->level);
        } else if (curr->level < w_level) {
            curr->level = w_level;
            curr->icd_in.clear();
            curr->icd_in.push_back(curr_e);
            SASSERT(curr_e->x->level == curr_e->y->level);
            for (Edge *outg : curr->out)
                F.push_back(outg);
        }
    }

    return ret;
}

bool Graph::reachable(Node *x, Node *y)
{
    if (algk == ITC)
        return desc[x->nid][y->nid].second != UINT_MAX;

    /* Determine whether a x -> y path exists by checking acyclicity of G+<y,x> */
    Edge __e(EDGE_K_TEMP, y, x);
    Edge *e = &__e;
    Node *u = e->x;
    Node *w = e->y;
    std::unordered_set<Node *> B;

    if (u->level < w->level)
        return false;

    auto bw_res = backward_search(B, e);

    if (bw_res == BW_CYCLE)
        return true;
    else if (bw_res == BW_OOF) {
        w->level = u->level + 1;
        w->icd_in.clear();
    } else if ((bw_res == BW_OK) && (u->level > w->level)) {
        w->level = u->level;
        w->icd_in.clear();
    } else
        return false;

    return forward_search(B, e);
}

unsigned Graph::delta()
{
    return std::min({ sqrt(edges), pow(nodes.size(), 2.0 / 3) });
}

bool Node::is_rec()
{
    return kind == REC_READ;
}

bool Node::is_fo_helper()
{
    return is_rex() || is_mfence() || is_sfence();
}

bool Node::is_rex()
{
    return kind == EX_READ;
}

bool Node::is_wex()
{
    return kind == EX_WRITE;
}

bool Node::is_assume()
{
    return kind == ASSUME;
}

void Node::print(z3::model m)
{
    std::cout << "     " << "N " << nid << ", loc : " << m.eval(loc).as_uint64()
        << ", guard : " << m.eval(guard)
        << ", crash : " << m.eval(crash)
        << ", val : " << m.eval(val).as_uint64() << "\n";
}

void Node::print()
{
    std::cout << "N " << *this << " : " << guard.simplify() << "\n";
    std::string istr;
    llvm::raw_string_ostream ss(istr);
    user->print(ss);
    std::cout << istr << "\n";
    for (auto *e : porf_out)
        if (e->kind == RF)
            std::cout << "  " << *e << "\n";
}
