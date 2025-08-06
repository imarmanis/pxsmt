#include "llvm/IR/CFG.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Verifier.h"
#include "llvm/IR/User.h"
#include "llvm/IR/InstrTypes.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/Operator.h"
#include "llvm/Analysis/InlineCost.h"
#include "llvm/Transforms/Utils/Cloning.h"
#include "llvm/AsmParser/Parser.h"
#include "llvm/ADT/SCCIterator.h"
#include "llvm/ADT/APInt.h"
#include "llvm/Support/SourceMgr.h"
#include "string"
#include "deque"
#include "ctime"
#include "execution.hpp"
#include "propagator.hpp"
#include "debug.hpp"
#include "config.hpp"
#include "verifier.hpp"
#include "z3++.h"

using namespace llvm;

Verifier::Verifier(std::unique_ptr<Module> _mod, AnalysisInfo &analysis)
    : s(z3::solver(c, (struct z3::solver::simple) {})), prop(&s),
    assumptions(c), assertions(c), no_diverge_assertions(c),
    data_layout(_mod.get()), alloc_counter(64)
{
    if (porf_idl) {
        z3::params p(c);
        p.set("auto_config", false);
        p.set("smt.arith.solver", (unsigned) 1);
        s.set(p);
    }
    std::srand(std::time(nullptr));
    unsigned seed = std::rand();
    z3::set_param("smt.random_seed", (int) seed);
    z3::set_param("sat.random_seed", (int) seed);
    s.set("random_seed", (unsigned) seed);
    s.set("smt.random_seed", (unsigned) seed);
    s.set("sat.random_seed", (unsigned) seed);
    mod = std::move(_mod);
    if (llvm::verifyModule(*mod))
        exit_msg("Broken module");

    Function *main = nullptr, *rec = nullptr;
    for (auto &F : *mod) {
        auto fname = F.getName();
        if (!(fname.compare("main")))
            main = &F;
        else if (!(fname.compare("__VERIFIER_recovery")))
            rec = &F;
    }

    check_pers = rec != nullptr;
    assert(main);

    node_list inits;
    for (auto &g : mod->globals()) {

        auto *initializer = g.getInitializer();
        unsigned alloc_size = data_layout.getTypeAllocSize(initializer->getType());
        auto *ptr_type = cast<PointerType>(g.getType());
        unsigned ptr_width = pointer_size(ptr_type) * BITS_IN_BYTE;

        auto align = g.getAlignment();
        auto pval = c.bv_val(getNewAllocPtr(alloc_size, align), ptr_width);
        scontext.emplace(&g, new SymValue(pval));
    }

    z3::expr_vector empty(c);
    assert(main->arg_empty());
    auto main_exec = newFExecution(main, c.bool_val(true), false);
    main_exec->process(empty, inits);

    FExecution *rec_exec = nullptr;
    if (rec) {
        assert(rec->arg_empty());
        rec_exec = newFExecution(rec, c.bool_val(true), true);
        rec_exec->process(empty, inits);
    }

    prop.finalize_nodes();
    for (auto t : pois) {
        Node *x, *y;
        std::tie(x, y) = t;
        prop.add_poi(x, y);
        if (mmk == SC)
            prop.add_sppoi(x, y);
    }
    for (auto t : static_rmws) {
        Node *x, *y;
        std::tie(x, y) = t;
        prop.add_rmw(x, y);
    }
    for (auto *x : *prop.get_nodes()) {
        if (!x->has_location()) continue;
        for (auto *y : *prop.get_nodes()) {
            if (!y->has_location()) continue;
            assert(!y->is_mfence());
            auto aar = analysis.getAliasInfo(std::make_pair(x->user, y->user));
            if ((aar == MustAlias) || (x == y))
                prop.same_loc[x->nid][y->nid] = c.bool_val(true);
            else if (aar == NoAlias)
                prop.same_loc[x->nid][y->nid] = c.bool_val(false);
            else
                prop.same_loc[x->nid][y->nid] = (x->loc == y->loc).simplify();
        }
    }
    for (auto *x : *prop.get_nodes()) {
        if (!x->has_location()) continue;
        for (auto *y : *prop.get_nodes()) {
            if (!y->has_location()) continue;
            if (prop.same_loc[x->nid][y->nid].is_true())
                prop.same_cl[x->nid][y->nid] = c.bool_val(true);
            else
                prop.same_cl[x->nid][y->nid] =
                    ((x->loc / 64) == (y->loc / 64)).simplify();
        }
    }

    Predecessors init_preds;
    if (mmk == TSO) {
        main_exec->process_dppo(init_preds);
        for (auto t : ppos) {
            Node *x, *y;
            z3::expr cond(c);
            std::tie(x, y, cond) = t;
            if (y->is_plain_read()) {
                auto *yp = cast<Instruction>(y->user)->getParent();
                auto *xi = dyn_cast<Instruction>(x->user);
                if ((xi && (xi->getParent() == yp))
                        || analysis.dominates(std::make_pair(xi, yp))) {
                    prop.add_sppoi(x, y);
                } else {
                    auto dppoi_e = c.bool_const(getUID("dppo").c_str());
                    if (crash_encoding == FULL)
                        cond = cond && !x->crash && !y->crash;
                    s.add(dppoi_e == cond);
                    prop.track_dppoi(x, y, dppoi_e);
                }
            } else {
                prop.add_sppoi(x, y);
            }
        }
    }

    prop.static_tc();

    /* FO helpers */
    for (auto *fo : prop.flushes) {
        if (!fo->is_flush_opt())
            continue;
        std::list<Node*> helper_candidates;
        for (auto *h : *prop.get_nodes()) {
            if (!(h->is_fo_helper() && prop.in_po(fo, h)))
                continue;
            helper_candidates.push_back(h);
        }
        for (auto *hc : helper_candidates) {
            /* helper can only be an instruction, not an init write */
            auto *hci = cast<Instruction>(hc->user);
            bool skip = false;
            /* Check if hc is dominated by any other hc' */
            for (auto *hcp : helper_candidates) {
                if (hcp == hc)
                    continue;
                auto *hcip = cast<Instruction>(hcp->user);
                auto same_block_po = (hci->getParent() == hcip->getParent())
                        && prop.in_po(hcp, hc);
                if (same_block_po ||
                    analysis.dominates(std::make_pair(hcip, hci->getParent()))) {
                    skip = true;
                    break;
                }
            }
            if (!skip) {
                prop.fo_helpers.emplace(fo, hc);
                prop.fo_helps.emplace(hc);
            }
        }
    }
    if ((crash_encoding == PARTIAL) && check_pers && !porf_idl) {
        for (auto *h : prop.fo_helps)
            prop.track_crashed(h, h->crash);
    }

    std::set<Node *, np_cmp> relevant_fohs;
    if (check_pers) {
        for (auto *f : prop.flushes) {
            if (f->is_flush()) {
                for (auto *r : prop.recs) {
                    auto sl = prop.same_cl[f->nid][r->nid];
                    if (sl.is_false())
                        continue;
                    z3::expr fb = c.bool_const(getUID("fb").c_str());
                    auto cond = ((crash_encoding == GUARD) ?
                                    f->guard : f->guard && !f->crash)
                                    && sl;
                    s.add(fb == cond);
                    prop.track_fb(f, r, fb);
                }
            } else if (f->is_flush_opt()) {
                z3::expr_vector vec(c);
                z3::expr foh = c.bool_const(getUID("foh").c_str());
                auto itp = prop.fo_helpers.equal_range(f);
                for (auto it = itp.first; it != itp.second; it++) {
                    auto h = it->second;
                    relevant_fohs.emplace(h);
                    auto cond = ((crash_encoding == GUARD) ? h->guard :
                                    h->guard && !h->crash);
                    vec.push_back(cond);
                }
                s.add(foh == z3::mk_or(vec));
                for (auto *r : prop.recs) {
                    auto sl = prop.same_cl[f->nid][r->nid];
                    if (sl.is_false())
                        continue;
                    z3::expr fb = c.bool_const(getUID("fb").c_str());
                    /* !f->crash not necessary because it is implied by helpers */
                    s.add(fb == (f->guard && foh && sl));
                    prop.track_fb(f, r, fb);
                }
            } else
                assert(false);
        }
    }

    if ((crash_encoding == PARTIAL) && check_pers && porf_idl) {
        z3::expr crash_sink_clock = c.int_const(getUID("crash_sink_clock").c_str());
        for (auto *r : prop.recs)
            s.add(r->porf_clock < crash_sink_clock);
        for (auto *x : *prop.get_nodes()) {
            if (x->is_flush() || x->is_assume())
                s.add(z3::implies(x->porf_clock < crash_sink_clock, !x->crash));
        }
        for (auto *h : relevant_fohs)
                s.add(z3::implies(h->porf_clock < crash_sink_clock, !h->crash));
    }

    SASSERT(s.check() == z3::sat);

    /* Maps read node to its rf candidates (node + rf expr) */
    std::multimap<Node*, std::pair<z3::expr, Node*>, np_cmp> rf_candidates;
    /* rf */
    for (auto r : prop.reads) {
        std::vector<Node*> po_before;
        std::vector<Node*> candidates;

        for (auto w : prop.writes) {
            /*
             * TODO : Remove silent assumption that all pointers have
             * the same width, i.e., locs have the same sort
             */
            z3::expr sym_sl = prop.same_loc[r->nid][w->nid];
            if (sym_sl.is_false())
                continue;

            /* Assumption : no mized-size accesses */
            if (r->val.get_sort().id() != w->val.get_sort().id())
                continue;

            /* Avoid future_read violations */
            if (prop.in_po(r, w))
                continue;

            if (prop.in_po(w, r)) po_before.push_back(w);
            else candidates.push_back(w);
        }

        std::list<Node*> min_po_before;
        for (auto pb : po_before) {
            bool skip = false;
            for (auto pbp : po_before) {
                /* Don't add those that are po-overwritten */
                if (prop.in_po(pb, pbp) && z3::implies(r->guard, z3::implies(pb->guard, pbp->guard
                    && (prop.same_loc[pb->nid][pbp->nid]))).simplify().is_true()) {
                    LOG(*pb << " is overwritten by " << *pbp << " wrt" << *r << "\n");
                    skip = true;
                    break;
                }
            }
            if (!skip)
                candidates.push_back(pb);
        }

        z3::expr_vector rfs(c);

        if (candidates.empty())
            exit_msg("All values must be initialized\n");

        for (auto w : candidates) {
            z3::expr rf = c.bool_const(getUID("rf").c_str());

            auto sv = r->val == w->val;
            auto sl = prop.same_loc[r->nid][w->nid];
            auto rf_rel = r->guard && w->guard && sv && sl;
            if (idl)
                rf_rel = rf_rel && (w->ob_clock < r->ob_clock);
            if (porf_idl)
                rf_rel = rf_rel && (w->porf_clock < r->porf_clock);
            if ((crash_encoding == FULL) && check_pers)
                rf_rel = rf_rel && (r->is_rec() ? !w->crash : !w->crash && !r->crash);
            if (r->is_rec())
                rf_rel = rf_rel && (r->rf_w_id == c.bv_val(w->nid, 64));

            /* Skip nodes that can't possibly refer to the same location */
            if (sl.is_false())
                continue;
            s.add(z3::implies(rf, rf_rel));
            rfs.push_back(rf);
            if (!idl)
                prop.track_rf(w, r, rf);
            rf_candidates.emplace(r, std::make_pair(rf, w));
        }

        auto rguard = (r->is_rec() || !check_pers || (crash_encoding != FULL))
                        ? r->guard : r->guard && !r->crash;
        s.add(z3::implies(rguard, z3::mk_or(rfs)));
    }

    /* Maps write node to its co-after candidates (node + rf expr) */
    std::multimap<Node*, std::pair<z3::expr, Node*>, np_cmp> co_candidates;
    /* co */
    for (auto w = prop.writes.begin(), we = prop.writes.end(); w != we; ++w) {
        auto wp = w;
        for (++wp; wp != we; ++wp) {
            auto w1 = *w, w2 = *wp;
            z3::expr co1 = c.bool_const(getUID("co").c_str());
            z3::expr co2 = c.bool_const(getUID("co").c_str());

            auto sl = prop.same_loc[w1->nid][w2->nid];
            /* Skip nodes that can't possibly refer to the same location */
            if (sl.is_false())
                continue;

            auto co_rel = w1->guard && w2->guard && sl;
            if ((crash_encoding == FULL) && check_pers)
                co_rel = co_rel && !w1->crash && !w2->crash;
            /* Avoid coh_ww violation */
            if (prop.in_po(w1, w2)) {
                s.add(co1 == co_rel);
                prop.track_co(w1, w2, co1);
                co_candidates.emplace(w1, std::make_pair(co1, w2));
            } else if (prop.in_po(w2, w1)) {
                s.add(co2 == co_rel);
                prop.track_co(w2, w1, co2);
                co_candidates.emplace(w2, std::make_pair(co2, w1));
            } else {
                s.add((co1 || co2) == co_rel);
                prop.track_co(w1, w2, co1);
                prop.track_co(w2, w1, co2);
                co_candidates.emplace(w1, std::make_pair(co1, w2));
                co_candidates.emplace(w2, std::make_pair(co2, w1));
            }
        }
    }

    for (auto w = prop.writes.begin(), we = prop.writes.end(); w != we; ++w) {
        auto wp = w;
        for (++wp; wp != we; ++wp) {
            auto w1 = *w, w2 = *wp;
            auto scl = prop.same_cl[w1->nid][w2->nid];
            if (!prop.in_po(w1, w2))
                continue;

            if (scl.is_false())
                continue;

            /*
            if (scl.is_true())
                prop.pocl.emplace(std::make_pair(w1, w2), -1);
            else
                prop.track_pocl(w1, w2);
                */
        }
    }

    /* fr */
    if (idl)
        for (auto r : prop.reads) {
            auto itp = rf_candidates.equal_range(r);
            for (auto it = itp.first; it != itp.second; it++) {
                z3::expr rf_w(c);
                Node *w;
                std::tie(rf_w, w) = it->second;
                auto itp2 = co_candidates.equal_range(w);
                for (auto it2 = itp2.first; it2 != itp2.second; it2++) {
                    z3::expr co_w(c);
                    Node *wp;
                    std::tie(co_w, wp) = it2->second;
                    /* r ->(fr) wp */
                    if (prop.is_rmw_pair(r, wp))
                        /* r <_ob wp shouldn't be added since it is r =_ob wp */
                        continue;
                    if (prop.in_po(wp, r))
                        s.add(!(rf_w && co_w));
                    else
                        s.add(z3::implies(rf_w && co_w, prop.ob_lt[r->nid][wp->nid]));
                }
            }
        }

    /* rf_w_id (for recovery nodes) */
    for (auto r = prop.recs.begin(), re = prop.recs.end(); r != re; ++r) {
        auto rp = r;
        for (++rp; rp != re; ++rp) {
            auto r1 = *r, r2 = *rp;
            auto sl = prop.same_loc[r1->nid][r2->nid];
            if (sl.is_false())
                continue;
            s.add(z3::implies(sl, r1->rf_w_id == r2->rf_w_id));
        }
    }

    auto assume_clause = z3::mk_and(assumptions);
    s.add(assume_clause);

    if (sanity_check)
        assert(s.check() == z3::sat);
    else
        SASSERT(s.check() == z3::sat);

    auto assert_no_diverge = check_unroll ?
        !z3::mk_and(no_diverge_assertions) : c.bool_val(false);
    auto assert_clause = !z3::mk_and(assertions) || assert_no_diverge;
    s.add(assert_clause);

    LOG(s);

    clock_t begin = std::clock();
    auto res = s.check();
    clock_t end = std::clock();
    switch (res) {
    case z3::unsat:
        PLOG("SAFE\n");
        std::cout << "SAFE\n";
        break;
    case z3::sat: {
        PLOG("UNSAFE\n");
        auto m = s.get_model();
        if (m.eval(assert_no_diverge).is_true())
            exit_msg("Insufficient unrolling\n");
        std::cout << "UNSAFE\n";
        for (auto n : *prop.get_nodes())
            n->log_edges();
        LOG(m);
        main_exec->print_ce(m);
        if (rec_exec)
            rec_exec->print_ce(m);
        print_graph(m);
        break;
    }
    case z3::unknown:
        exit_msg(("UNKNOWN : " + s.reason_unknown() + "\n").c_str());
        break;
    }

    if (emit_stats) {
        auto stats = s.statistics();
        std::cout << stats << "\n";
        std::ignore = double(end - begin) / CLOCKS_PER_SEC;
        prop.emitStats();
    }
}

/*
 * nl : list with po|imm-before events, each event appears at most once.
 * snl : list of fence events, along with the condition
 * upon which they are before the current instruction.
 */

std::string Verifier::getUID(const char *pref)
{
    std::string spref(pref);
    return (spref + std::to_string(counter++));
}

z3::expr Verifier::symVal(Value *val)
{
    auto it = scontext.find(val);
    if (it != scontext.end())
        return it->second->val;
    if (auto *cd = dyn_cast<ConstantInt>(val)) {
        unsigned width = cd->getBitWidth();

        return c.bv_val(cd->getSExtValue(), width);
    }

    if (isa<UndefValue>(val)) {
        LOG("Using undefined value\n");
        auto *ut = cast<IntegerType>(val->getType());
        unsigned width = ut->getBitWidth();

        return c.bv_const(getUID("uv").c_str(), width);
    }

    if (isa<ConstantPointerNull>(val)) {
        auto width = pointer_size(cast<PointerType>(val->getType())) * BITS_IN_BYTE;
        return c.bv_val(0, width);
    }

    if (auto *bc = dyn_cast<BitCastOperator>(val)) {
        if (bc->getSrcTy()->isPointerTy())
            return symVal(bc->getOperand(0));
        /* TODO : other types */
        assert(false);
    }

    /* Should be checked after checking context, since GEPInstr isa<GEPOperator> */
    if (isa<GEPOperator>(val))
        return gep_sval(val);

    /* For ConstExpr IntToPtr casts */
    if (auto *op = dyn_cast<Operator>(val)) {
        CastInst *ci = nullptr;
        if (auto *x = dyn_cast<CastInst>(op))
            ci = x;
        else if (auto *ce = dyn_cast<ConstantExpr>(op))
            ci = cast<CastInst>(ce->getAsInstruction());
        else
            assert(false);
        if (ci->isNoopCast(data_layout))
            return symVal(ci->getOperand(0));
        assert(false);
    }

    assert(false);
}

z3::expr Verifier::compare(CmpInst::Predicate p,  const z3::expr &e1, const z3::expr &e2)
{
    switch (p) {
    case CmpInst::ICMP_EQ :
        return (e1 == e2);

    case CmpInst::ICMP_NE :
        return (e1 != e2);

    case CmpInst::ICMP_SGE :
        return (e1 >= e2);

    case CmpInst::ICMP_ULE :
        return (ule(e1, e2));

    case CmpInst::ICMP_UGE :
        return (uge(e1, e2));

    case CmpInst::ICMP_ULT :
        return (ult(e1, e2));

    case CmpInst::ICMP_UGT :
        return (ugt(e1, e2));

    case CmpInst::ICMP_SGT :
        return (e1 > e2);

    case CmpInst::ICMP_SLE :
        return (e1 <= e2);

    case CmpInst::ICMP_SLT :
        return (e1 < e2);

    default:
        assert(false);
        break;
    }
}

/* Does not include alive conditions */
z3::expr Verifier::getEdgeCondition(BasicBlock *bh, BasicBlock *bt)
{
    auto *br = cast<BranchInst>(bh->getTerminator());

    if (br->isUnconditional()) {
        SASSERT(br->getSuccessor(0) == bt);
        return c.bool_val(true);
    }

    auto val = symVal(br->getOperand(0));
    auto w = val.get_sort().bv_size();
    auto scond = val != c.bv_val(0, w);

    z3::expr ret(c);
    if (br->getSuccessor(0) == bt)
        ret = scond;
    else if (br->getSuccessor(1) == bt)
        ret = !scond;
    else
        assert(false);
    return ret;
}

unsigned Verifier::getNewAllocPtr(unsigned d, unsigned alignment)
{
    assert(alignment);
    unsigned rem = alloc_counter % alignment;
    /* Round up if needed to satisfy alignment */
    if (rem != 0)
        alloc_counter += alignment - rem;
    auto ret = alloc_counter;
    alloc_counter += d;
    return ret;
}

void Verifier::print_graph(z3::model m)
{
    for (auto n : *prop.get_nodes())
        n->print();
}

/* Size of pointer, in bytes */
unsigned Verifier::pointer_size(PointerType *ptr_t)
{
    auto size = data_layout.getPointerSize(ptr_t->getAddressSpace());
    assert(size == 8);
    return size;
}

z3::expr Verifier::get_uninterpreted(Type *t)
{
    if (t->isVoidTy())
        return z3::expr(c);

    if (auto it = dyn_cast<IntegerType>(t)) {
        auto width = it->getBitWidth();
        auto val = c.bv_const(getUID("ui").c_str(), width);
        return val;
    } else if (auto pt = dyn_cast<PointerType>(t)) {
        auto width = pointer_size(pt) * BITS_IN_BYTE;
        auto val = c.bv_const(getUID("up").c_str(), width);
        return val;
    }

    assert(false);
}

z3::expr Verifier::gep_sval(Value *val)
{
    auto *gep = cast<GEPOperator>(val);
    auto *t = gep->getPointerOperandType();
    auto ptr_width = pointer_size(cast<PointerType>(gep->getType())) * BITS_IN_BYTE;
    z3::expr offset = c.bv_val(0, ptr_width);
    for (auto I = gep->idx_begin(), IE = gep->idx_end(); I != IE; ++I) {
        auto sindx = symVal(*I);
        auto indw = sindx.get_sort().bv_size();
        if (indw < ptr_width)
            sindx = sext(sindx, ptr_width - indw);

        if (auto *st = dyn_cast<StructType>(t)) {
            /* struct offset shouldn't be symbolic */
            auto indx = sindx.simplify();
            assert(indx.is_const());
            auto *sl = data_layout.getStructLayout(st);
            auto ti = indx.get_numeral_uint();
            auto el_off = sl->getElementOffset(ti);
            offset = offset + c.bv_val(el_off, ptr_width);
            t = st->getTypeAtIndex(ti);
            continue;
        }

        if (auto *pt = dyn_cast<PointerType>(t))
            t = pt->getElementType();
        else if (auto *at = dyn_cast<ArrayType>(t))
            t = at->getElementType();
        else
            assert(false);

        auto add = sindx * c.bv_val(data_layout.getTypeAllocSize(t), ptr_width);
        offset = offset + add;
    }
    auto base = symVal(gep->getPointerOperand());
    return base + offset;
}

std::pair<z3::func_decl, z3::func_decl_vector> Verifier::getCustom(StructType * st)
{
    auto res = custom_types.find(st);
    if (res != custom_types.end())
        return res->second;

    unsigned n = st->getNumContainedTypes();
    assert(n == 2);
    std::vector<z3::sort> sorts;
    char *names[2];
    unsigned i = 0;
    for (auto *I = st->element_begin(), *IE = st->element_end(); I != IE; I++) {
        auto s = to_sort(*I);
        sorts.push_back(s);
        auto str = s.name().str();
        auto l = str.length();
        names[i] = new char[l + 1];
        str.copy(names[i], l);
        names[i][l] = '\0';
        i++;
    }
    assert(i == n);
    func_decl_vector v(c);
    func_decl fd = c.tuple_sort(getUID("ts").c_str(), n, names, &sorts[0], v);
    for (unsigned i = 0; i < n; i++)
        delete[] names[i];
    auto p = std::make_pair(fd, v);
    custom_types.emplace(st, p);
    return p;
}

z3::sort Verifier::to_sort(Type *t)
{
    if (auto *it = dyn_cast<IntegerType>(t)) {
        return c.bv_sort(it->getBitWidth());
    } else if (auto *pt = dyn_cast<PointerType>(t)) {
        return c.bv_sort(pointer_size(pt) * BITS_IN_BYTE);
    }
    assert(false);
}

FExecution * Verifier::newFExecution(Function *f, z3::expr alive, bool rec)
{
    auto fexec = new FExecution(f, this, alive, rec);
    auto empl = to_fexec.emplace(f, fexec);
    assert(empl.second);
    return fexec;
}

SymValue * Verifier::newSymValue(Value *val, z3::expr sval)
{
    auto ret = new SymValue(sval);
    auto empl = scontext.emplace(val, ret);
    assert(empl.second);
    return ret;
}

bool Verifier::in_rec(Node *n)
{
    auto *instr = cast<Instruction>(n->user);
    if (!instr)
        return false;
    return instr->getParent()->getName() == "__VERIFIER_recovery";
}
