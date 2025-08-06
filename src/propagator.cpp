#include "vector"
#include "tuple"
#include "z3++.h"
#include "propagator.hpp"
#include "config.hpp"
#include "debug.hpp"
#include "graph.hpp"

using namespace z3;

bool prefix(const reasons &a, const reasons &b) {
    unsigned n = a.size();
    if (n > b.size())
        return false;
    for (auto ita = a.begin(), itb = b.begin(); ita != a.end(); ++ita, ++itb)
        if (ita->id() != itb->id())
            return false;
    return true;
}

void Propagator::decide(expr &e, unsigned &n, Z3_lbool &b)
{
}

void Propagator::fixed(const expr &n, const expr &expr)
{
    SASSERT(verifyStack());
    SASSERT(g.verifyEdges());
    SASSERT(!is_fixed(n));
    assert(!after_final);
    scopes.back().fixed.push_back(n);
    fixed_exprs.insert(n);

    if (pop_pending) {
        LOG("while a pop is pending : set " << n << " to " << expr << "\n");
        return;
    }

    if (!expr.is_true()) {
        return;
    }

    auto cit = cnode_from_id.find(n);
    if (cit != cnode_from_id.end()) {
        assert(crash_encoding == PARTIAL);
        assert(!porf_idl);
        auto *cnode = cit->second;
        crashed.emplace(cnode, n);
        scopes.back().crashed.push_back(cnode);
        if (porf_rec.find(cnode) != porf_rec.end()) {
            auto cb = [this, n] (reasons& rs) {
                reasons rl = rs;
                rl.sort([](const z3::expr &x, const z3::expr &y) {
                        return x.id() < y.id();
                });
                LOG("porf+crash conflict\n");
                z3::expr_vector confls(s->ctx());
                for (auto i : rl)
                    confls.push_back(i);
                confls.push_back(n);
                conflict_wrapper(confls, "crash");
            };
            process_porf_rec_paths(cnode, cb);
        }
        return;
    }

    /* An edge is set */
    Edge *e = g.get_tracked_edge(n);
    LOG("set : " << n << " :: " << *e << "\n");
    process_edge(e);
}

bool Propagator::process_edge(Edge *e)
{
    /* If node is not porf-rec reachable, then add edge to pending and stop */
    if (crash_encoding == PARTIAL) {
        auto alive_x = porf_rec.find(e->x) != porf_rec.end();
        auto alive_y = porf_rec.find(e->y) != porf_rec.end();
        if ((!alive_x || !alive_y) && !((e->kind == RF) && alive_y)) {
            pendings_out.insert({e->x, e});
            pendings_in.insert({e->y, e});
            return true;
        }
    }

    /* Check coherence */
    if (!(future_read(e) && coh_ww(e) && coh_rw(e) && coh_rr(e) && atomicity(e))) {
        return false;
    }
    auto frs = derive_frs(e);
    for (Edge *fr : frs)
        if (!(coh_wr(fr) && coh_rr(fr) && atomicity(fr))) {
            return false;
        }

    if (!new_edge(e)) {
        return false;
    }

    /* Add derived fr edges */
    for (Edge *fr : frs) {
        if (!new_edge(fr)) {
            return false;
        }
    }

    auto n_porf_recs = porf_rf_add(e);
    for (auto p : n_porf_recs)
            scopes.back().porf_rec.push_back(p);

    if ((e->kind == RF) && (crash_encoding == PARTIAL) && check_pers && !porf_idl) {
        for (auto p : n_porf_recs) {
        /* In case of partial crash encoding, only FL, FOH, and assumes are needed */
            if (!p->is_flush() && !p->is_fo_helper() && !p->is_assume())
                    continue;
            auto cb = [this, p](reasons& rs) {
                    reasons rl = rs;
                    // FIXME
                    rl.sort([](const z3::expr &x, const z3::expr &y) {
                                    return x.id() < y.id();
                                    });
                    z3::expr_vector porf_crash_prop(s->ctx());
                    for (auto i : rl)
                            porf_crash_prop.push_back(i);
                    auto it = crashed.find(p);
                    if (it != crashed.end()) {
                            LOG("crash+porf conflict\n");
                            porf_crash_prop.push_back(it->second);
                            conflict_wrapper(porf_crash_prop, "crash");
                    } else if (!redundant_cprop(p, rl)) {
                            LOG("porf+crash consequence :" << *p << " (" <<
                                            p->kind << ")\n");
                            consequences++;
                            TRACE_CODE(
                                            for (auto x : porf_crash_prop)
                                            LOG(x << ", ");
                                            LOG("\n");
                                      );
                            propagate(porf_crash_prop, !p->crash);
                            c_propagations.emplace(p, rl);
                    }
            };
            process_porf_rec_paths(p, cb);
        }
    }
    if (pop_pending)
            return false;


    for (auto p : n_porf_recs) {
            auto itp_in = pendings_in.equal_range(p);
            for (auto it = itp_in.first; it != itp_in.second; it++) {
                auto *e  = it->second;
                if (porf_rec.find(e->x) == porf_rec.end())
                        continue;
                if (!process_edge(e))
                        return false;
            }

            auto itp_out = pendings_out.equal_range(p);
            for (auto it = itp_out.first; it != itp_out.second; it++) {
                auto *e  = it->second;
                if (porf_rec.find(e->y) == porf_rec.end())
                        continue;
                if (!process_edge(e))
                        return false;
            }
    }
    return true;
}

Node *Propagator::new_node(llvm::User *user, z3::expr guard, z3::expr loc,
    z3::expr val, NODE_KIND kind)
{
    std::string str1 = "porfcl" + std::to_string(g.nodes.size());
    std::string str2 = "obcl" + std::to_string(g.nodes.size());
    auto porf_clock = s->ctx().int_const(str1.c_str());
    auto ob_clock = s->ctx().int_const(str2.c_str());
    auto nid = g.nodes.size();
    auto ncr = s->ctx().bool_const((("crash_") + std::to_string(nid)).c_str());
    auto rf_w_id = s->ctx().bv_const((("rf_w_id") + std::to_string(nid)).c_str(), 64);
    Node *n = new Node(porf_clock, ob_clock, nid, user, guard, loc,
                        val, kind, ncr, rf_w_id);
    /* FL and assume always have their crash tracked, but only some fo_helpers */
    if ((crash_encoding == PARTIAL) && check_pers && (n->is_flush() || n->is_assume())
        && !porf_idl)
        track_crashed(n, ncr);

    LOG("node " << *n);
    if (idl)
        LOG(", ob_clock " << ob_clock);
    if (porf_idl)
        LOG(", porf_clock " << porf_clock);
    LOG("\n");
    g.add_node(n);
    if (n->is_rec())
            porf_rec.insert(n);
    return n;
}

Node *Propagator::read_node(llvm::User *user, z3::expr guard, z3::expr loc,
    z3::expr val)
{
    auto n = new_node(user, guard, loc, val, P_READ);
    reads.insert(n);
    return n;
}

Node *Propagator::read_ex_node(llvm::User *user, z3::expr guard, z3::expr loc,
    z3::expr val)
{
    auto n = new_node(user, guard, loc, val, EX_READ);
    reads.insert(n);
    return n;
}

Node *Propagator::write_node(llvm::User *user, z3::expr guard, z3::expr loc,
    z3::expr val)
{
    auto n = new_node(user, guard, loc, val, P_WRITE);
    writes.insert(n);
    return n;
}

Node *Propagator::write_ex_node(llvm::User *user, z3::expr guard, z3::expr loc,
    z3::expr val)
{
    auto n = new_node(user, guard, loc, val, EX_WRITE);
    writes.insert(n);
    return n;
}

Node *Propagator::flush_node(llvm::User *user, z3::expr guard, z3::expr loc)
{
    auto n = new_node(user, guard, loc, z3::expr(s->ctx()), FLUSH);
    flushes.push_back(n);
    return n;
}

Node *Propagator::flush_opt_node(llvm::User *user, z3::expr guard, z3::expr loc)
{
    auto n = new_node(user, guard, loc, z3::expr(s->ctx()), FLUSH_OPT);
    flushes.push_back(n);
    return n;
}

Node *Propagator::mfence_node(llvm::User *user, z3::expr guard)
{
    return new_node(user, guard, z3::expr(s->ctx()), z3::expr(s->ctx()), MFENCE);
}

Node *Propagator::sfence_node(llvm::User *user, z3::expr guard)
{
    return new_node(user, guard, z3::expr(s->ctx()), z3::expr(s->ctx()), SFENCE);
}

Node *Propagator::rec_node(llvm::User *user, z3::expr guard, z3::expr loc,
    z3::expr val)
{
    auto n = new_node(user, guard, loc, val, REC_READ);
    reads.insert(n);
    recs.insert(n);
    return n;
}

Node *Propagator::assume_node(llvm::User *user, z3::expr guard, z3::expr val)
{
    return new_node(user, guard, z3::expr(s->ctx()), val, ASSUME);
}

void Propagator::track_crashed(Node *cn, const z3::expr &np)
{
    assert(cn->is_flush() || cn->is_fo_helper() || cn->is_assume());
    LOG("crash track : " << *cn);
    track(np);
    cnode_from_id.emplace(np, cn);
}

void Propagator::track(const z3::expr &e)
{
        tracked_variables++;
        add(e);
}

void Propagator::track_dppoi(Node *x, Node *y, const expr &e)
{
    LOG("dppoi : " << *x << " ~> " << *y << "\n");
    if (idl) {
        /* Can it be the case? If so, implication shouldn't be added */
        assert(!is_rmw_pair(x, y));
        s->add(z3::implies(e, ob_lt[x->nid][y->nid]));
        return;
    }
    std::tuple<Node*, Node*, expr> t(x, y, e);
    track(e);
    LOG("tracked " << e << "\n");
    tracked.insert({std::make_pair(x, y), e});
    g.register_dppoi_edge(x, y, e);
}

void Propagator::track_rf(Node *x, Node *y, const expr &e)
{
    std::tuple<Node*, Node*, expr> t(x, y, e);
    track(e);
    LOG("tracked " << e << "\n");
    tracked.insert({std::make_pair(x, y), e});
    g.register_rf_edge(x, y, e);
}

void Propagator::track_fb(Node *x, Node *y, const expr &e)
{
    std::tuple<Node*, Node*, expr> t(x, y, e);
    track(e);
    LOG("tracked " << e << "\n");
    tracked.insert({std::make_pair(x, y), e});
    g.register_fb_edge(x, y, e);
}

void Propagator::track_co(Node *x, Node *y, const expr &e)
{
    if (idl) {
        s->add(z3::implies(e, ob_lt[x->nid][y->nid]));
        return;
    }
    std::tuple<Node*, Node*, expr> t(x, y, e);
    track(e);
    LOG("tracked " << e << "\n");
    tracked.insert({std::make_pair(x, y), e});
    g.register_co_edge(x, y, e);
}

/* Keeps track of the edge and inserts it into the graph */
void Propagator::add_sppoi(Node *x, Node *y)
{
    sppoi.push_back(std::make_pair(x, y));
    LOG("SPPOI : " << *x << " ~> " << * y << "\n");
    if (idl) {
        /* For <R,W> in RMW pairs, their ob_clock are equal */
        if (!is_rmw_pair(x, y))
            s->add(ob_lt[x->nid][y->nid]);
        return;
    }
    Edge *po = new Edge(SPPO, x, y);
    /* In case ICD is used, reachable MUST be called before insert_edge */
    assert(!g.reachable(y, x));
    auto ecb = [](path&, unsigned){ return; };
    g.insert_edge(po, &scopes.back().tc_pairs, ecb);
}

/* Adds to poi, sppoi, and keeps track of the rmw */
void Propagator::add_rmw(Node *x, Node *y)
{
    add_poi(x, y);
    add_sppoi(x, y);
    rmws.insert(std::make_pair(x, y));
    /* RMW is atomic wrt crash */
    if ((crash_encoding == FULL) && check_pers)
        s->add(z3::implies(y->crash, x->crash));
    if (idl)
        s->add(y->ob_clock == x->ob_clock);
}

void Propagator::static_tc() {
    unsigned n = g.nodes.size();
    Node *x, *y;
    for (auto p : poi) {
        std::tie(x, y) = p;
        po[x->nid][y->nid] = true;
    }
    for (auto p : sppoi) {
        std::tie(x, y) = p;
        sppo[x->nid][y->nid] = true;
    }

    for (unsigned k = 0; k < n; k++)
        for (unsigned i = 0; i < n; i++)
            for (unsigned j = 0; j < n; j++)
                po[i][j] = po[i][j] || (po[i][k] && po[k][j]);

    for (unsigned k = 0; k < n; k++)
        for (unsigned i = 0; i < n; i++)
            for (unsigned j = 0; j < n; j++)
                sppo[i][j] = sppo[i][j] || (sppo[i][k] && sppo[k][j]);

    for (unsigned i = 0; i < n; i++)
        for (unsigned j = 0; j < n; j++)
            /* Initially porf = po, since no rf is set */
            porf[i][j] = po[i][j];

    for (unsigned i = 0; i < n; i++)
        porf[i][i] = true;

    for (unsigned i = 0; i < n; i++)
        SASSERT(!po[i][i] && !sppo[i][i]);
    PLOG("nodes : " << n << "\n");
}

bool Propagator::new_edge(Edge *e)
{
    /* 
     * Don't add rfi and coi to (ppo U obs)-graph
     * SPPO edges are not added via this function
     */
    SASSERT(e->kind != SPPO);
    auto k = e->kind;
    bool obs = !in_po(e->x, e->y) && ((k == RF) || (k == CO) || (k == FR));
    if ((k == DPPO) || (k == FB) || obs) {
        if (g.reachable(e->y, e->x)) {
            std::map<Node*, z3::expr> cache;
            auto cb = [this, e, &cache](path &p) {
                p.push_back(e);
                propagate_conflict(&p, cache);
                p.pop_back();
            };
            g.process_paths(e->y, e->x, cb, sppo);
            /* If a cycle is found, don't bother updating additional structures */
            return false;
        }
        auto ecb = [this](path &p, const z3::expr &id) {
            propagate_consequence(&p, id);
        };
        g.insert_edge(e, &scopes.back().tc_pairs, ecb);
        insertions++;
        LOG("new graph edge : " << *e << "\n");
    }

    scopes.back().edges.push_back(e);

    if (!in_po(e->x, e->y) && ((k == CO) || (k == FR)))
        inc_obsw.insert({e->y, e});

    if (in_po(e->x, e->y) && ((k == RF) || (k == CO)))
        e->x->out_i.push_back(e);

    return true;
}

std::vector<Edge*> Propagator::derive_frs(Edge *e)
{
    /* fr_reasons = { rf_reason, co_reason } */

    std::vector<Edge *> frs;
    switch (e->kind) {
    case RF:
        for (std::list<Edge*> *l : {&e->x->out, &e->x->out_i})
            for (Edge *ep : *l)
                if (ep->kind == CO)
                    frs.push_back(new Edge(FR, e->y, ep->y, e->reasons[0],
                                                           ep->reasons[0]));
        break;

    case CO:
        for (std::list<Edge*> *l : {&e->x->out, &e->x->out_i})
            for (Edge *ep : *l)
                if (ep->kind == RF)
                    frs.push_back(new Edge(FR, ep->y, e->y, ep->reasons[0],
                                                           e->reasons[0]));
        break;

    default:
        break;
    }

    return frs;
}

void Propagator::propagate_consequence(std::list<Edge *> *path, z3::expr conseq)
{
    consequences++;
    DEBUG_CODE(LOG("path consequence\n"););
    z3::expr_vector causes(s->ctx());
    for (Edge *e : *path) {
        DEBUG_CODE(
            LOG(*e << ", ");
            LOG("reason : ");
        );
        for (auto reason : e->reasons) {
            causes.push_back(reason);
            DEBUG_CODE(LOG(reason << ", "););
        }
        DEBUG_CODE(LOG("\n"););
    }

    propagate(causes, conseq);

}

/* Propagate a path conflict, cache reason for which a Fl node is alive */
void Propagator::propagate_conflict(std::list<Edge *> *path, std::map<Node*, z3::expr> &cache)
{
    /* TODO : cache results, usefull when propagating more than one path conflict */
    path_conflicts++;
    DEBUG_CODE(LOG("path conflict\n"););
    z3::expr_vector confls(s->ctx());
    for (Edge *e : *path) {
        DEBUG_CODE(
            LOG(*e << ", ");
            LOG("reason : ");
        );
        for (auto reason : e->reasons) {
            confls.push_back(reason);
            DEBUG_CODE(LOG(reason << ", "););
        }
        DEBUG_CODE(LOG("\n"););
    }

    conflict_wrapper(confls, "path");
}

/* irreflexive(po;rf) */
bool Propagator::future_read(Edge *e)
{
    if ((e->kind == RF) && in_po(e->y, e->x)) {
        assert(false);
        return false;
    }
    return true;
}

/* irreflexive(po;co) */
bool Propagator::coh_ww(Edge *e)
{
    if ((e->kind == CO) && in_po(e->y, e->x)) {
        assert(false);
        return false;
    }
    return true;
}

/* irreflexive(po;fr) */
bool Propagator::coh_wr(Edge *e)
{
    if ((e->kind == FR) && in_po(e->y, e->x)) {
        z3::expr_vector confls(s->ctx());
        for (auto e : e->reasons)
                confls.push_back(e);
        conflict_wrapper(confls, "coh_wr");
        return false;
    }
    return true;
}

/* irreflexive(po;coe;rfe) */
bool Propagator::coh_rw(Edge *e)
{
    if (in_po(e->x, e->y))
        return true;

    switch (e->kind) {
    case CO: {
        for (Edge *rf : e->y->out) {
            if (in_po(rf->y, e->x)) {
                z3::expr_vector confls(s->ctx());
                confls.push_back(e->reasons[0]);
                confls.push_back(rf->reasons[0]);
                conflict_wrapper(confls, "coh_rw-co");
                return false;
            }
        }
        break;
    }

    case RF: {
        bool ret = true;
        auto itp = inc_obsw.equal_range(e->x);
        for (auto it = itp.first; it != itp.second; it++) {
            if ((it->second->kind == CO) && (in_po(e->y, it->second->x))) {
                z3::expr_vector confls(s->ctx());
                confls.push_back(it->second->reasons[0]);
                confls.push_back(e->reasons[0]);
                conflict_wrapper(confls, "coh_rw-rf");
                ret = false;
            }
        }
        return ret;
    }

    default:
        break;
    }

    return true;
}

/* irreflexive(po;fre;rfe) */
bool Propagator::coh_rr(Edge *e)
{
    if (in_po(e->x, e->y))
        return true;

    switch (e->kind) {
    case FR: {
        for (Edge *rf : e->y->out) {
            if (in_po(rf->y, e->x)) {
                z3::expr_vector confls(s->ctx());
                confls.push_back(e->reasons[0]);
                confls.push_back(e->reasons[1]);
                confls.push_back(rf->reasons[0]);
                conflict_wrapper(confls, "coh_rr-fr");
                return false;
            }
        }
        break;
    }

    case RF: {
        bool ret = true;
        auto itp = inc_obsw.equal_range(e->x);
        for (auto it = itp.first; it != itp.second; it++) {
            if ((it->second->kind == FR) && (in_po(e->y, it->second->x))) {
                z3::expr_vector confls(s->ctx());
                confls.push_back(it->second->reasons[0]);
                confls.push_back(it->second->reasons[1]);
                confls.push_back(e->reasons[0]);
                conflict_wrapper(confls, "coh_rr-rf");
                ret = false;
            }
        }
        return ret;
    }

    default:
        break;
    }

    return true;
}

/* irreflexive(rmw-1;fre;coe) */
bool Propagator::atomicity(Edge *e)
{
    if (in_po(e->x, e->y))
        return true;

    switch (e->kind) {
    case FR: {
        for (Edge *co : e->y->out) {
            if ((co->kind == CO) &&
                (rmws.find(std::make_pair(e->x, co->y)) != rmws.end())) {
                z3::expr_vector confls(s->ctx());
                confls.push_back(e->reasons[0]);
                confls.push_back(e->reasons[1]);
                confls.push_back(co->reasons[0]);
                conflict_wrapper(confls, "atom-fr");
                return false;
            }
        }
        break;
    }

    case CO: {
        auto itp = inc_obsw.equal_range(e->x);
        for (auto it = itp.first; it != itp.second; it++) {
            if ((it->second->kind == FR) &&
                (rmws.find(std::make_pair(it->second->x, e->y)) != rmws.end())) {
                z3::expr_vector confls(s->ctx());
                confls.push_back(it->second->reasons[0]);
                confls.push_back(it->second->reasons[1]);
                confls.push_back(e->reasons[0]);
                conflict_wrapper(confls, "atomc-co");
                return false;
            }
        }
        break;
    }

    default:
        break;
    }

    return true;
}

/* Keeps track of the poi edge */
void Propagator::add_poi(Node *x, Node *y)
{
    if ((crash_encoding == FULL) && check_pers)
        s->add(z3::implies(x->guard && y->guard && x->crash, y->crash));
    if (porf_idl)
        s->add(x->porf_clock < y->porf_clock);
    poi.push_back(std::make_pair(x, y));
    /* PO is added only once in porf_out */
    x->porf_out.push_back(new Edge(PO, x, y));
}

bool Propagator::in_po(Node *x, Node *y)
{
    SASSERT(po);
    return po[x->nid][y->nid];
}

bool Propagator::external(Node *x, Node *y)
{
    return !(in_po(x,y) || in_po(y, x));
}

std::vector<Node *> *Propagator::get_nodes()
{
    return &g.nodes;
}

void Propagator::push()
{
    LOG("scope : push\n");
    assert(!after_final);

    auto &tscope = scopes.back();
    if (tscope.empty())
        tscope.incr_count();
    else
        scopes.emplace_back();
    SASSERT(verifyStack());
}

void Propagator::pop(unsigned d)
{
    LOG("scope : pop " << d << "\n");

    after_final = false;
    pop_pending = false;
    while(d > 0) {
        auto &tscope = scopes.back();
        for(Edge *e : tscope.edges) {
            LOG("unset : " << *e << "\n";);

            auto k = e->kind;
            SASSERT(k != SPPO);
            bool obs = !in_po(e->x, e->y) && ((k == RF) || (k == CO) || (k == FR));
            if ((k == DPPO) || (k == FB) || obs) {
                SASSERT(g.contains(e));
                g.remove_edge(e);
            }

            if (external(e->x, e->y) && ((k == CO) || (k == FR)))
            {
                auto itp = inc_obsw.equal_range(e->y);
                bool found = false;
                for (auto it = itp.first; it != itp.second; it++) {
                    if (it->second == e) {
                        inc_obsw.erase(it);
                        found = true;
                        break;
                    }
                }
                SASSERT(found);
                std::ignore = found;
                DEBUG_CODE(
                    auto itp = inc_obsw.equal_range(e->y);
                    for (auto it = itp.first; it != itp.second; it++)
                        assert(it->second != e);
                );
            }

            if (in_po(e->x, e->y) && ((k == RF) || (k == CO))) {
                SASSERT(e->x->hasOut_i(e));
                e->x->out_i.remove(e);
            }

            if (k == RF)
                e->x->porf_out.remove(e);

            if ((k == RF) && e->y->is_rec())
                recovered_writes.erase(e->x);
            /*
            else if ((k == FR) && e->x->is_rec())
                lost_writes.erase(e->y);
            */

            SASSERT(!g.contains(e));

            if (k == FR)
                delete e;
        }
        for(auto t : tscope.tc_pairs) {
            unsigned x, y;
            std::tie(x, y) = t;
            g.desc[x][y] = not_desc;
        };
        for (auto fs : tscope.fixed) {
            auto cit = cnode_from_id.find(fs);
            if (cit != cnode_from_id.end()) {
                    // TODO: Replace crashed (scope)
            } else if (crash_encoding == PARTIAL) {
                Edge *e = g.get_tracked_edge(fs);
                auto itpx = pendings_out.equal_range(e->x);
                for (auto it = itpx.first; it != itpx.second; it++) {
                    if (it->second == e) {
                        pendings_out.erase(it);
                        break;
                    }
                }
                auto itpy = pendings_in.equal_range(e->y);
                for (auto it = itpy.first; it != itpy.second; it++) {
                    if (it->second == e) {
                        pendings_in.erase(it);
                        break;
                    }
                }
            }
            fixed_exprs.erase(fs);
        }
        for (auto t : tscope.porf_tc_pairs) {
            unsigned x, y;
            std::tie(x, y) = t;
            porf[x][y] = false;
        }
        for (auto p : tscope.porf_rec)
            porf_rec.erase(p);
        for (auto cn : tscope.crashed)
            crashed.erase(cn);
        SASSERT(!scopes.empty());

        tscope.clear();
        auto decr = tscope.decr_count(d);
        d -= decr;
        if (tscope.dead())
            scopes.pop_back();
    }

    SASSERT(verifyStack());
}

Propagator::Propagator(solver *s) : user_propagator_base(s), s(s)
{
    scopes.emplace_back();
    g.fixed = &fixed_exprs;
    g.tracked = &tracked;
    register_fixed();
    //register_decide();
}

void Propagator::conflict_wrapper(const z3::expr_vector &a, std::string str)
{
    pop_pending = true;
    LOG("conflict : " << str << " : ");
    for (const z3::expr &exp : a) {
        TRACE_CODE(
            auto cit = cnode_from_id.find(exp);
            if (cit != cnode_from_id.end()) {
                SASSERT(is_fixed(exp));
                LOG(" " << exp << " (" << *cit->second << "), ");
            } else {
                SASSERT(is_fixed(exp));
                LOG(" " << exp << " (" << *g.edge_from_id.at(exp) << "), ");
            }
        );
    }
    LOG("\n");

    conflict(a);
}

#ifdef DEBUG
bool Propagator::is_fixed(const z3::expr &n)
{
    for (auto &scope : scopes) {
        for (auto id : scope.fixed) {
            if (id == n) {
                assert(fixed_exprs.find(n) != fixed_exprs.end());
                return true;
            }
        }
    }
    assert(fixed_exprs.find(n) == fixed_exprs.end());
    return false;
}

bool Propagator::is_fixed(Edge *e)
{
    for (auto &scope : scopes)
        for (auto ep : scope.edges)
            if (ep == e)
                return true;

    return false;
}

bool Propagator::verifyStack()
{
    for (auto n : g.nodes) {
        for (auto *e : n->out) {
            for (auto r : e->reasons)
                if (!is_fixed(r))
                    return false;
            if ((e->kind != SPPO) && !is_fixed(e))
                return false;
        }
    }

    return true;
}

#endif /* DEBUG */

void Propagator::emitStats()
{
    std:: cout << "stats : " << insertions << " insertions, "\
    << path_conflicts << " conflicts, "\
    << consequences << " consequences, "
    << get_nodes()->size() << " nodes\n";
}

void Propagator::finalize_nodes()
{
    unsigned n = g.num_nodes();
    g.desc = new std::pair<Edge*, unsigned>*[n];
    po = new bool*[n];
    sppo = new bool*[n];
    porf = new bool*[n];
    for (unsigned i = 0; i < n; i++) {
        g.desc[i] = new std::pair<Edge*, unsigned>[n];
        po[i] = new bool[n];
        sppo[i] = new bool[n];
        porf[i] = new bool[n];
    }
    for (unsigned i = 0; i < n; i++)
        for (unsigned j = 0; j < n; j++) {
            g.desc[i][j] = not_desc;
            po[i][j] = false;
            sppo[i][j] = false;
            porf[i][j] = false;
        }

    for (unsigned i = 0; i < n; i++) {
        g.desc[i][i] = std::make_pair(nullptr, 0);
        porf[i][i] = true;
    }

    same_loc = new std::vector<z3::expr>[n];
    same_cl = new std::vector<z3::expr>[n];
    ob_lt = new std::vector<z3::expr>[n];
    for (unsigned i = 0; i < n; i++) {
        auto *v = &same_loc[i];
        v->reserve(n);
        /* Bypass copy constructor of z3::expr problem */
        for (unsigned j = 0; j < n; j++)
            v->push_back(z3::expr(s->ctx()));
    }
    for (unsigned i = 0; i < n; i++) {
        auto *v = &same_cl[i];
        v->reserve(n);
        /* Bypass copy constructor of z3::expr problem */
        for (unsigned j = 0; j < n; j++)
            v->push_back(z3::expr(s->ctx()));
    }

    for (auto n1 : *get_nodes()) {
        auto *v2 = &ob_lt[n1->nid];
        for (auto n2 : *get_nodes()) {
            z3::expr ltporf(s->ctx());
            auto ltob = n1->ob_clock < n2->ob_clock;
            v2->push_back(ltob);
        }
    }
}

std::list<Node*> Propagator::porf_rf_add(Edge *rf)
{
    assert(!porf_idl);
    Node * s = rf->x, *t = rf->y;
    std::list<Node*> n_porf_rec;
    if (rf->kind != RF)
            return n_porf_rec;

    /* PORF cycle, it will be detected later, don't do anything */
    //assert(!porf_reachable(t,s));
    if (porf_reachable(t, s))
        return n_porf_rec;
    s->porf_out.push_back(rf);
    assert(porf_reachable(s, s));
    if (!porf_reachable(s, t))
        for (Node *v : g.nodes)
            if (porf_reachable(v, s) &&
            !porf_reachable(v, t))
                porf_meld(n_porf_rec, v, t);
    assert(porf_reachable(s, t));
    return n_porf_rec;
}

void Propagator::porf_meld(std::list<Node*> &n_porf_rec, Node *a, Node *b)
{
    std::stack<std::pair<Node*, Node*>> W;
    W.push(std::make_pair(a, b));
    auto &s = scopes.back();

    while (!W.empty()) {
        Node *v, *t;
        auto p = W.top();
        W.pop();
        std::tie(v, t) = p;
        if (porf_reachable(v, t))
            continue;
        porf[v->nid][t->nid] = true;
        s.porf_tc_pairs.push_back(std::make_pair(v->nid, t->nid));
        if (t->is_rec() && (porf_rec.find(v) == porf_rec.end())) {
                porf_rec.insert(v);
                n_porf_rec.push_back(v);
        }
        for (Edge *e : t->porf_out)
            if (!porf_reachable(v, e->y))
                W.push(std::make_pair(v, e->y));
    }
}

bool Propagator::porf_reachable(Node *x, Node *y)
{
    return porf[x->nid][y->nid];
}

void Propagator::process_porf_path(Node *x, Node *y, reasons_cb cb)
{
    std::list<std::pair<Edge*, unsigned>> W;
    std::unordered_set<Node *> V;
    for (auto e : x->porf_out)
        if (e->kind == PO)
            W.push_back(std::make_pair(e, 0));
        else if (e->kind == RF)
            W.push_back(std::make_pair(e, 1));
        else
            assert(false);
    for (auto *k : *get_nodes())
        k->dist = UINT_MAX;
    x->dist = 0;
    V.insert(x);

    x->parent = nullptr;
    while (!W.empty()) {
        auto ewt = W.front();
        auto e = ewt.first;
        auto d = ewt.second;
        W.pop_front();
        auto n = e->y;

        if (V.find(n) != V.end())
            continue;

        V.insert(n);
        n->dist = d;
        n->parent = e;

        if (n == y)
                break;

        for (auto ep : n->porf_out)
            if (ep->kind == PO)
                W.push_front(std::make_pair(ep, d + 0));
            else if (ep->kind == RF)
                W.push_back(std::make_pair(ep, d + 1));
            else assert(false);
    }

    // Might break if rf can contradict po, shouldn't happen tought
    assert(y->dist != 0);
        Node *curr = y;
        std::list<z3::expr> reasons;
        while (curr->parent != nullptr) {
            auto e = curr->parent;
            if (e->kind == RF)
                reasons.push_back(e->reasons[0]);
            curr = e->x;
        }
        cb(reasons);
        return;

}

void Propagator::process_porf_rec_paths(Node *n, reasons_cb cb)
{
    assert(n->is_flush() || n->is_assume() || n->is_fo_helper());
    std::list<std::pair<Edge*, unsigned>> W;
    std::unordered_set<Node *> V;
    for (auto e : n->porf_out)
        if (e->kind == PO)
            W.push_back(std::make_pair(e, 0));
        else if (e->kind == RF)
            W.push_back(std::make_pair(e, 1));
        else
            assert(false);
    for (auto *k : *get_nodes())
        k->dist = UINT_MAX;
    n->dist = 0;
    V.insert(n);

    std::list<Node *> recs;
    n->parent = nullptr;
    while (!W.empty()) {
        auto ewt = W.front();
        auto e = ewt.first;
        auto d = ewt.second;
        W.pop_front();
        auto n = e->y;

        if (V.find(n) != V.end())
            continue;

        V.insert(n);
        n->dist = d;
        n->parent = e;

        if (n->kind == REC_READ) {
            recs.push_back(n);
            continue;
        }

        for (auto ep : n->porf_out)
            if (ep->kind == PO)
                W.push_front(std::make_pair(ep, d + 0));
            else if (ep->kind == RF)
                W.push_back(std::make_pair(ep, d + 1));
            else assert(false);
    }

    assert(!recs.empty());
    if (!all_porf_paths) {
        Node *curr = recs.front();
        std::list<z3::expr> reasons;
        while (curr->parent != nullptr) {
            auto e = curr->parent;
            if (e->kind == RF)
                reasons.push_back(e->reasons[0]);
            curr = e->x;
        }
        cb(reasons);
        return;
    }

    std::list<z3::expr> path;
    for (auto rn : recs)
        porf_paths(path, true, rn, n, cb);
}

void Propagator::porf_paths(std::list<z3::expr> &p, bool use_po, Node * curr,
    Node *target, reasons_cb cb)
{
    if (curr == target) {
        cb(p);
        return;
    }
    if (curr->parent->kind == RF) {
        p.push_back(curr->parent->reasons[0]);
        porf_paths(p, true, curr->parent->x, target, cb);
        p.pop_back();
    }
    if (use_po) {
        for (auto *k : *get_nodes()) {
            if (curr == k)
                continue;
            if (po[k->nid][curr->nid] && (curr->dist == k->dist))
                porf_paths(p, false, k, target, cb);
        }
    }
}

bool Propagator::is_rmw_pair(Node *x, Node *y)
{
    return x->is_rex() && y->is_wex() && (x->user == y->user);
}

bool Propagator::redundant_cprop(Node *n, const reasons &l) {
    SASSERT(std::is_sorted(l.begin(), l.end()));
    auto itp = c_propagations.equal_range(n);
    for (auto it = itp.first; it != itp.second; it++)
        if (prefix(l, it->second))
            return true;
    return false;
}
