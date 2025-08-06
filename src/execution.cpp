#include "z3++.h"
#include "execution.hpp"
#include "verifier.hpp"
#include "llvm/IR/BasicBlock.h"
#include "llvm/IR/CFG.h"
#include "llvm/ADT/SCCIterator.h"
#include "llvm/Transforms/Utils/Cloning.h"
#include "deque"
#include "config.hpp"

using namespace llvm;

void append_unique(std::list<Node*> &ps, Node *n)
{
    if (std::find(ps.begin(), ps.end(), n) == ps.end())
        ps.push_back(n);
}

/* Add conditional event, merge conditions (||) if event exists already */
void add_cond(snode_list &ps, Node *n, z3::expr ncond)
{
    snode_list::iterator I = ps.begin(), IE = ps.end();
    while ((I != IE) && (I->first != n)) I++;

    if (I != IE)
        I->second = I->second || ncond;
    else
        ps.push_back(std::make_pair(n, ncond));

    return;
}

void FExecution::process(z3::expr_vector params, node_list preds)
{
    for (auto I = llvm::scc_begin(func), IE = llvm::scc_end(func); I != IE; ++I) {
        auto bvec = *I;
        assert(bvec.size() == 1);
        auto *bbexec = newBBExecution(bvec[0]);
        bb_executions.push_front(bbexec);
    }

    auto *eb = &func->getEntryBlock();
    auto *eb_exec = to_bbexec[eb];

    eb_exec->preds = preds;
    ver->s.add(alive == eb_exec->alive_in);

    unsigned i = 0;
    for (auto &I : func->args())
        std::ignore = ver->newSymValue(&I, params[i++]);

    for (auto bb_ex : bb_executions)
        bb_ex->process();

    for (auto bb_ex : bb_executions) {
        auto bb = bb_ex->bb;
        auto balive = bb_ex->alive_in;
        z3::expr_vector inc_edges(ver->c);
        for (auto *pred : predecessors(bb))
            /* Exclude diverge block */
            if (pred != bb) {
                auto pred_exec = to_bbexec.at(pred);
                inc_edges.push_back(pred_exec->alive_out && getEdgeCondition(pred_exec, bb_ex));
            }
        if (bb != &func->getEntryBlock()) {
            auto cond = balive == mk_or(inc_edges);
            ver->s.add(cond);
        }
    }

}

BBExecution * FExecution::newBBExecution(BasicBlock *bb)
{
    auto bbexec = new BBExecution(bb, this);
    auto empl = to_bbexec.emplace(bb, bbexec);
    assert(empl.second);
    return bbexec;
}

void BBExecution::process()
{
    auto ipreds = preds;
    auto reached = alive_in;
    for (auto &i : *bb) {
        auto t = process_instruction(&i, ipreds, reached);
        instr_executions.push_back(t.first);
        reached = t.second;
    }
    alive_out = reached;
}

FExecution::FExecution(Function *func, Verifier *ver, z3::expr alive, bool rec)
    : func(func), ver(ver), alive(alive),
    ret_val(ver->get_uninterpreted(func->getReturnType())), rec(rec),
    terminated(ver->c)

{ }

BBExecution::BBExecution(BasicBlock *bb, FExecution *parent)
    : bb(bb), parent(parent), ver(parent->ver),
    alive_in(ver->c.bool_const(ver->getUID("b").c_str())),
    alive_out(alive_in) { }

std::pair<InstrExecution *, z3::expr> BBExecution::process_instruction(
    Instruction *instr, node_list &preds, z3::expr reached)
{
    auto reachedp = reached;
    auto &i = *instr;
    if (i.isBinaryOp()) {
        z3::expr (*op)(z3::expr const &, z3:: expr const &);
        switch (i.getOpcode()) {
        case Instruction::Add : {
            op = &z3::operator+;
            break;
        }

        case Instruction::Sub : {
            op = &z3::operator-;
            break;
        }

        case Instruction::And : {
            op = &z3::operator&;
            break;
        }

        case Instruction::Or : {
            op = &z3::operator|;
            break;
        }

        case Instruction::Xor : {
            op = &z3::operator^;
            break;
        }

        case Instruction::Mul : {
            op = &z3::operator*;
            break;
        }

        case Instruction::SDiv : {
            op = &z3::operator/;
            break;
        }

        case Instruction::UDiv : {
            op = &z3::udiv;
            break;
        }

        case Instruction::SRem : {
            op = &z3::srem;
            break;
        }

        case Instruction::URem : {
            op = &z3::urem;
            break;
        }

        case Instruction::Shl : {
            op = &z3::shl;
            break;
        }

        case Instruction::LShr : {
            op = &z3::lshr;
            break;
        }

        case Instruction::AShr : {
            op = &z3::ashr;
            break;
        }

        default:
            assert(false);

        }

        auto op1 = i.getOperand(0), op2 = i.getOperand(1);
        auto sop1 = ver->symVal(op1), sop2 = ver->symVal(op2);
        auto res = op(sop1, sop2);
        auto iexec = newIExecution(instr);
        iexec->sval = ver->newSymValue(&i, res);
        return std::make_pair(iexec, reachedp);
    }

    InstrExecution *iexec = nullptr;
    switch (i.getOpcode()) {

    case Instruction::ICmp : {
        auto sval1 = ver->symVal(i.getOperand(0)), sval2 = ver->symVal(i.getOperand(1));
        auto *icmp_i = cast<ICmpInst>(&i);
        auto cmp = ver->compare(icmp_i->getPredicate(), sval1, sval2);
        auto sval = z3::ite(cmp, ver->c.bv_val(1, 1), ver->c.bv_val(0, 1));
        iexec = newIExecution(instr);
        iexec->sval = ver->newSymValue(&i, sval);
        break;
    }

    case Instruction::Alloca: {
        auto *ai = cast<AllocaInst>(&i);
        auto a_size_opt = ai->getAllocationSizeInBits(ver->data_layout);
        assert(a_size_opt.hasValue());
        auto a_size = a_size_opt.getValue();
        assert(!(a_size % BITS_IN_BYTE));
        auto ptr_size = ver->pointer_size(ai->getType());
        /* Always align to cache-line size to avoid false cache-line sharing */
        auto align = 64; //ai->getAlignment();
        assert(ptr_size == 8);
        auto pval = ver->c.bv_val(ver->getNewAllocPtr(a_size / BITS_IN_BYTE,
            align), ptr_size * BITS_IN_BYTE);
        iexec = newIExecution(instr);
        iexec->sval = ver->newSymValue(&i, pval);
        break;
    }

    case Instruction::Load : {

        auto *li = cast<LoadInst>(&i);
        unsigned width;
        if (auto *t = dyn_cast<IntegerType>(li->getType()))
            width = t->getBitWidth();
        else if (auto *t = dyn_cast<PointerType>(li->getType()))
            width = ver->pointer_size(t) * BITS_IN_BYTE;
        else
            assert(false);

        assert(!(width % BITS_IN_BYTE));
        auto val = ver->c.bv_const(ver->getUID("l").c_str(), width);

        Node *rnode;
        if (!parent->rec)
            reachedp = mk_guard_crash(reached);
        if (parent->rec)
            rnode = ver->prop.rec_node(li, reachedp,
                ver->symVal(li->getPointerOperand()), val);
        else
            rnode = ver->prop.read_node(li, reachedp,
                ver->symVal(li->getPointerOperand()), val);

        iexec = newIExecution(instr, rnode);
        iexec->sval = ver->newSymValue(&i, val);

        if (parent->rec)
            break;

        for (auto x : preds)
            ver->pois.push_back(std::make_tuple(x, rnode));

        preds = {rnode};

        break;
    }

    case Instruction::Store : {

        assert(!parent->rec);

        auto *si = cast<StoreInst>(&i);
        auto sloc = ver->symVal(si->getPointerOperand());
        auto sval = ver->symVal(si->getValueOperand());

        /* no writes in the recovery routine */
        reachedp = mk_guard_crash(reached);
        auto wnode = ver->prop.write_node(si, reachedp, sloc, sval);
        iexec = newIExecution(instr, wnode);

        for (auto x : preds)
            ver->pois.push_back(std::make_tuple(x, wnode));
        preds = {wnode};
        break;
    }

    case Instruction::BitCast :
    case Instruction::PtrToInt :
    case Instruction::IntToPtr :
    {
        auto *ci = cast<CastInst>(&i);
        if (ci->isNoopCast(ver->data_layout)) {
            iexec = newIExecution(instr);
            iexec->sval = ver->newSymValue(&i, ver->symVal(ci->getOperand(0)));
            break;
        }
        assert(false);
    }

    case Instruction::ZExt : {
        auto *ze_i = cast<ZExtInst>(&i);
        SASSERT(ze_i->isIntegerCast());
        auto srct = cast<IntegerType>(ze_i->getSrcTy()), destt = cast<IntegerType>(ze_i->getDestTy());
        auto val = z3::zext(ver->symVal(ze_i->getOperand(0)), destt->getBitWidth() - srct->getBitWidth());
        iexec = newIExecution(instr);
        iexec->sval = ver->newSymValue(&i, val);
        break;
    }

    case Instruction::SExt : {
        auto *sext_i = cast<SExtInst>(&i);
        auto destw = cast<IntegerType>(sext_i->getDestTy())->getBitWidth();
        auto srctw = cast<IntegerType>(sext_i->getSrcTy())->getBitWidth();
        auto val = sext(ver->symVal(sext_i->getOperand(0)), destw - srctw);
        iexec = newIExecution(instr);
        iexec->sval = ver->newSymValue(&i, val);
        break;
    }

    case Instruction::Trunc : {
        auto *tr_i = cast<TruncInst>(&i);
        auto destt = cast<IntegerType>(tr_i->getDestTy());
        auto destw = destt->getBitWidth();
        auto val = ver->symVal(tr_i->getOperand(0)).extract(destw - 1, 0);
        iexec = newIExecution(instr);
        iexec->sval = ver->newSymValue(&i, val);
        break;
    }

    case Instruction::PHI : {
        auto *phi_i = cast<PHINode>(&i);
        auto t = phi_i->getType();
        unsigned w;
        if (auto *pt = dyn_cast<PointerType>(t))
            w = ver->pointer_size(pt) * BITS_IN_BYTE;
        else if (auto *it = dyn_cast<IntegerType>(t))
            w = it->getBitWidth();
        else
            assert(false);

        auto phi_val = ver->c.bv_const(ver->getUID("p").c_str(), w);

        iexec = newIExecution(instr);
        iexec->sval = ver->newSymValue(&i, phi_val);

        for (unsigned i = 0, ie = phi_i->getNumIncomingValues(); i < ie ; ++i) {
            auto *b = phi_i->getIncomingBlock(i);
            auto b_ex = parent->to_bbexec[b];
            auto *v = phi_i->getIncomingValue(i);
            auto c = b_ex->alive_out && ver->getEdgeCondition(b, phi_i->getParent());
            ver->s.add(z3::implies(c, phi_val == ver->symVal(v)));
        }
        break;
    }

    case Instruction::AtomicRMW : {
        assert(!parent->rec);

        auto *ai = cast<AtomicRMWInst>(&i);
        assert(ai->getOperation() == AtomicRMWInst::BinOp::Add);
        auto *at = cast<IntegerType>(ai->getType());
        unsigned width = at->getBitWidth();
        assert(width == 32);

        auto rval = ver->c.bv_const(ver->getUID("l").c_str(), width);

        auto loc = ver->symVal(ai->getPointerOperand());
        auto add = ver->symVal(ai->getValOperand());
        reachedp = mk_guard_crash(reached);
        auto rnode = ver->prop.read_ex_node(ai, reachedp, loc, rval);
        auto wnode = ver->prop.write_ex_node(ai, reachedp, loc, rval + add);
        ver->static_rmws.push_back(std::make_tuple(rnode, wnode));

        for (auto x : preds)
            ver->pois.push_back(std::make_tuple(x, rnode));

        iexec = newIExecution(&i, rnode, wnode);
        iexec->sval = ver->newSymValue(&i, rval);
        preds = {wnode};
        break;
    }

    case Instruction::AtomicCmpXchg : {
        assert(!parent->rec);

        auto *casi = cast<AtomicCmpXchgInst>(&i);
        auto *p_op = casi->getPointerOperand();
        auto *t = cast<PointerType>(p_op->getType())->getElementType();
        unsigned width;
        if (auto *it = dyn_cast<IntegerType>(t)) {
            width = it->getBitWidth();
        } else if (auto *pt = dyn_cast<PointerType>(t)) {
            width = ver->pointer_size(pt) * BITS_IN_BYTE;
        } else
            assert(false);

        auto rval = ver->c.bv_const(ver->getUID("l").c_str(), width);

        auto loc = ver->symVal(p_op);
        auto test_val = ver->symVal(casi->getCompareOperand());
        auto new_val = ver->symVal(casi->getNewValOperand());
        auto success = (rval == test_val);
        auto success_bv = z3::ite(success, ver->c.bv_val(1, 1), ver->c.bv_val(0, 1));
        reachedp = mk_guard_crash(reached);
        auto rnode = ver->prop.read_ex_node(&i, reachedp, loc, rval);
        auto wnode = ver->prop.write_ex_node(&i, reachedp && success, loc, new_val);
        /* Note : seems OK to always add the static RMW edge */
        ver->static_rmws.push_back(std::make_tuple(rnode, wnode));

        auto st = cast<StructType>(casi->getType());
        auto ct = ver->getCustom(st);

        for (auto x : preds)
            ver->pois.push_back(std::make_tuple(x, rnode));

        iexec = newIExecution(&i, rnode, wnode);
        iexec->sval = ver->newSymValue(&i, ct.first(rval, success_bv));
        /* Include both since the write might not always be alive */
        preds = {rnode, wnode};

        break;
    }

    case Instruction::ExtractValue : {
        auto *ev = cast<ExtractValueInst>(&i);
        auto *op = ev->getAggregateOperand();
        assert(ev->getNumIndices() == 1);
        auto stp = ver->symVal(op);
        auto ind = ev->getIndices()[0];
        assert(stp.is_app());
        auto val = stp.arg(ind);
        iexec = newIExecution(&i);
        iexec->sval = ver->newSymValue(&i, val);

        /* Avoid using projections since we know the index
        auto *ot = cast<StructType>(op->getType());
        auto ct = getCustom(ot);
        auto val = ct.second[ind](stp);
        */
        break;
    }

    case Instruction::Call : {
        auto *ci = cast<CallInst>(&i);
        auto fname = ci->getCalledFunction()->getName().str();
        if (!fname.compare("__VERIFIER_assert")) {
            auto cond = ci->getOperand(0);

            auto scond = ver->symVal(cond) != ver->c.bv_val(0,1);

            reachedp = mk_guard_crash(reached);
            auto assertion = z3::implies(reachedp, scond);
            ver->assertions.push_back(assertion);
            iexec = newIExecution(&i);
            /* assert are only allowed inside the recovery routine,
             * in case we are checking persistency */
            assert(!check_pers || parent->rec);
            break;
        }
        if (!fname.compare("__VERIFIER_assume")) {
            auto cond = ci->getOperand(0);

            auto scond = ver->symVal(cond) != ver->c.bv_val(0,1);

            if (!parent->rec)
                reachedp = mk_guard_crash(reached);
            auto guard = reachedp;
            auto anode = ver->prop.assume_node(&i, guard, scond);
            auto hyp =
                (crash_encoding == GUARD) ? guard : (
                    (!parent->rec && check_pers) ?
                    !anode->crash && reached :
                    reached
                );
            auto assumption = z3::implies(hyp, scond);
            ver->assumptions.push_back(assumption);
            iexec = newIExecution(&i);
            for (auto x : preds)
                ver->pois.push_back(std::make_tuple(x, anode));
            preds = {anode};
            break;
        }

        if (!fname.compare("__VERIFIER_end_loop")) {
            iexec = newIExecution(&i);
            break;
        }

        if (!fname.compare("__VERIFIER_parallel")) {

            assert(!parent->rec);
            ValueToValueMapTy empty;
            std::list<Function *> par_fs;

            for (auto I = ci->arg_begin(), IE = ci->arg_end(); I != IE; ++I) {
                auto fval = cast<Value>(I);
                auto func = cast<Function>(fval->stripPointerCasts());
                par_fs.push_back(func);
            }

            node_list npreds;
            auto icexec = newCallIExecution(&i);
            z3::expr_vector pterm(ver->c);
            for (auto *f : par_fs) {
                node_list f_term;
                z3::expr_vector empty(ver->c);
                auto f_exec = ver->newFExecution(f, reached, parent->rec);
                icexec->called_execs.push_back(f_exec);
                f_exec->process(empty, preds);
                for (auto pr : f_exec->terminators)
                    append_unique(npreds, pr);
                pterm.push_back(z3::mk_or(f_exec->terminated));
            }
            if (crash_encoding != PARTIAL) {
                reachedp = ver->c.bool_const(ver->getUID("pt").c_str());
                ver->s.add(reachedp == mk_and(pterm));
            }
            preds = npreds;
            iexec = icexec;
            break;
        }

        if (!fname.compare("llvm.x86.sse2.clflush")) {

            assert(!parent->rec);
            auto loc = ver->symVal(ci->getArgOperand(0));
            reachedp = mk_guard_crash(reached);
            Node *fnode = ver->prop.flush_node(ci, reachedp, loc);


            for (auto x : preds)
                ver->pois.push_back(std::make_tuple(x, fnode));
            preds = {fnode};
            iexec = newIExecution(&i, fnode);

            break;
        }

        if (!fname.compare("llvm.x86.clflushopt")) {

            assert(!parent->rec);
            auto loc = ver->symVal(ci->getArgOperand(0));
            reachedp = mk_guard_crash(reached);
            Node *fnode = ver->prop.flush_opt_node(ci, reachedp, loc);

            for (auto x : preds)
                ver->pois.push_back(std::make_tuple(x, fnode));

            preds = {fnode};
            iexec = newIExecution(&i, fnode);

            break;
        }

        if (!fname.compare("__VERIFIER_nondet_int")) {

            auto *t = cast<IntegerType>(ci->getType());
            unsigned width = t->getBitWidth();

            auto val = ver->c.bv_const(ver->getUID("nd").c_str(), width);
            iexec = newIExecution(&i);
            iexec->sval = ver->newSymValue(&i, val);

            break;
        }

        if (!fname.rfind("__VERIFIER_", 0)) {
            std::cerr << "Unhandled " << fname << "\n";
            exit(1);
        }

        if (!fname.compare("llvm.x86.sse2.mfence")) {
            assert(!parent->rec);

            /* Maybe we don't need that */
            reachedp = mk_guard_crash(reached);
            auto node = ver->prop.mfence_node(&i, reachedp);
            for (auto x : preds)
                ver->pois.push_back(std::make_tuple(x, node));

            preds = {node};
            iexec = newIExecution(&i, node);
            break;
        }

        if (!fname.compare("llvm.x86.sse.sfence")) {
            assert(!parent->rec);

            reachedp = mk_guard_crash(reached);
            auto node = ver->prop.sfence_node(&i, reachedp);
            for (auto x : preds)
                ver->pois.push_back(std::make_tuple(x, node));

            preds = {node};
            iexec = newIExecution(&i, node);
            break;
        }

        /* Skip the rest of llvm intrinsic functions */
        if (!fname.rfind("llvm", 0)) {
            iexec = newIExecution(&i);
            break;
        }

        /* Regular function, dead code for now since they are inlined */
        assert(false);
        std::ignore = ci->getCalledFunction();
        auto n = ci->getNumArgOperands();
        z3::expr_vector params(ver->c);
        for (unsigned i = 0; i < n; i++)
            params.push_back(ver->symVal(ci->getArgOperand(i)));
        /* TODO */

        break;
    }

    case Instruction::Br : {
        auto *bri = cast<BranchInst>(&i);
        auto succs = bri->successors();
        iexec = newIExecution(&i);

        /* Diverge block */
        auto head = succs.begin();
        if (*head == i.getParent()) {
            /* If full encoding is used, we are missing
             * the crash-po axiom since asserts do not
             * corresspond to nodes (implementation detail).
             * Note that reached does not include any crash
             * variable. */
            if ((crash_encoding == FULL) && check_pers) {
                z3::expr_vector v(ver->c);
                for (auto p : preds)
                    v.push_back(p->guard && !p->crash);
                ver->no_diverge_assertions.push_back(
                    z3::implies(z3::mk_or(v), !reached)
                );
            } else
                ver->no_diverge_assertions.push_back(!reached);
            break;
        }

        for (auto s : succs) {
            auto b_exec = parent->to_bbexec[s];
            b_exec->add_predecessors(preds);
        }

        break;
    }

    case Instruction::Select : {
        auto *si = cast<SelectInst>(&i);
        auto scond = ver->symVal(si->getCondition()) == ver->c.bv_val(1, 1);
        auto val = z3::ite(scond, ver->symVal(si->getTrueValue()),
            ver->symVal(si->getFalseValue()));

        iexec = newIExecution(&i);
        iexec->sval = ver->newSymValue(&i, val);
        break;
    }

    case Instruction::GetElementPtr : {
        auto gepval = ver->gep_sval(&i);
        iexec = newIExecution(&i);
        iexec->sval = ver->newSymValue(&i, gepval);
        break;
    }

    case Instruction::Ret : {
        parent->terminated.push_back(reached);
        for (auto p : preds)
            append_unique(parent->terminators, p);
        auto *reti = cast<ReturnInst>(&i);
        auto rval = reti->getReturnValue();
        auto sret = parent->ret_val;
        if (sret) {
            assert(rval);
            auto ret_cond = z3::implies(reached, sret == ver->symVal(rval));
            ver->s.add(ret_cond);
        } else
            assert(!rval);
        iexec = newIExecution(&i);
        break;

    }

    case Instruction::Unreachable : {
        assert(false);
    }

    default:
        outs() << "Unhandled instruction : " << i << "\n";
        assert(false);
        break;
    }
    assert(iexec);
    return std::make_pair(iexec, reachedp);
}

SymValue::SymValue(z3::expr val) : val(val) {};

InstrExecution::InstrExecution(Instruction *instr, BBExecution *parent)
    : InstrExecution(K_INSTR, instr, parent) {};

InstrExecution::InstrExecution(Kind k, Instruction *instr, BBExecution *parent)
    : kind(k), instr(instr), parent(parent), ver(parent->ver) {};

InstrExecution * BBExecution::newIExecution(Instruction *instr)
{
    auto iexec = new InstrExecution(instr, this);
    auto empl = to_iexec.emplace(instr, iexec);
    assert(empl.second);
    return iexec;
}

InstrExecution * BBExecution::newIExecution(Instruction *instr, Node *n)
{
    auto iexec = new MemoryInstrExecution(instr, this, n);
    auto empl = to_iexec.emplace(instr, iexec);
    assert(empl.second);
    return iexec;
}

InstrExecution * BBExecution::newIExecution(Instruction *instr, Node *n1, Node *n2)
{
    auto iexec = new MemoryInstrExecution(instr, this, n1, n2);
    auto empl = to_iexec.emplace(instr, iexec);
    assert(empl.second);
    return iexec;
}

CallInstrExecution * BBExecution::newCallIExecution(Instruction *instr)
{
    auto iexec = new CallInstrExecution(instr, this);
    auto empl = to_iexec.emplace(instr, iexec);
    assert(empl.second);
    return iexec;
}

MemoryInstrExecution::MemoryInstrExecution(Instruction *instr, BBExecution *parent,
    Node *n) : InstrExecution(K_MEM, instr, parent), nodes{n} {};
MemoryInstrExecution::MemoryInstrExecution(Instruction *instr, BBExecution *parent,
    Node *n1, Node *n2): InstrExecution(K_MEM, instr, parent), nodes{n1, n2} {};

CallInstrExecution::CallInstrExecution(Instruction *instr, BBExecution *parent)
    : InstrExecution(K_INSTR, instr, parent) {};


z3::expr FExecution::getEdgeCondition(BBExecution *bbh, BBExecution *bbt)
{
    return ver->getEdgeCondition(bbh->bb, bbt->bb);
}

void BBExecution::add_predecessors(std::list<Node*> prds)
{
    for (auto prd : prds)
        append_unique(preds, prd);
}

void BBExecution::add_c_predecessors(const Predecessors &ncpreds)
{ c_preds.merge_with(ncpreds); }

void Predecessors::merge_with(const Predecessors &nps)
{
    for (auto np : nps.all)
        add_cond(all, np.first, np.second);
}

void Predecessors::merge_with(const Predecessors &nps, z3::expr cond)
{
    for (auto np : nps.all)
        add_cond(all, np.first, np.second && cond);
}

void BBExecution::add_c_predecessors(const Predecessors &ncpreds, z3::expr cond)
{ c_preds.merge_with(ncpreds, cond); }

void FExecution::add_c_terminators(const Predecessors &ncpred)
{ c_terminators.merge_with(ncpred); }

void FExecution::process_dppo(Predecessors c_preds)
{
    auto *eb = &func->getEntryBlock();
    auto *eb_exec = to_bbexec[eb];
    eb_exec->c_preds = c_preds;
    for (auto bb_ex : bb_executions)
        bb_ex->process_dppo();
}

void BBExecution::process_dppo()
{
    Predecessors ic_preds = c_preds;
    for (auto iexec : instr_executions)
        iexec->process_dppo(ic_preds);
}

void InstrExecution::process_dppo(Predecessors &c_preds)
{
    auto &i = *instr;
    switch (i.getOpcode()) {
    case Instruction::Load : {
        auto *mi = dyn_cast<MemoryInstrExecution>(this);
        Node *rnode = mi->nodes[0];
        snode_list nc_preds;
        for (auto x : c_preds.all) {
            Node *pred;
            z3::expr cond(ver->c);
            std::tie(pred, cond) = x;
            auto t = std::make_tuple(pred, rnode, cond);
            switch (pred->kind) {
            /* These are ordered before the read, and everything is ordered before a read
             * so we don't need to keep them */
            case P_READ :
            case EX_READ :
            case EX_WRITE :
            case MFENCE :
                ver->ppos.push_back(t);
                break;

            /* These are not ordered before the read, we keep them for following (non-read) events */
            case P_WRITE :
            case SFENCE :
            case FLUSH :
            case FLUSH_OPT :
                nc_preds.push_back(x);
                break;

            default:
                assert(false);
            }
        }
        nc_preds.push_back({std::make_pair(rnode, rnode->guard)});
        c_preds.all = nc_preds;
        break;
    }

    case Instruction::Store : {
        auto *mi = dyn_cast<MemoryInstrExecution>(this);
        Node *wnode = mi->nodes[0];
        snode_list nc_preds;
        for (auto x : c_preds.all) {
            Node *pred;
            z3::expr cond(ver->c);
            std::tie(pred, cond) = x;
            auto t = std::make_tuple(pred, wnode, cond);
            /* Everything is ordered before a write */
            ver->ppos.push_back(t);
            /* Following events may not be ordered after the write
             * but they might need to be ordered before the preds */
            nc_preds.push_back(x);
        }
        nc_preds.push_back({std::make_pair(wnode, wnode->guard)});
        c_preds.all = nc_preds;
        break;
    }

    case Instruction::AtomicRMW : {
        auto *mi = dyn_cast<MemoryInstrExecution>(this);
        Node *rnode = mi->nodes[0], *wnode = mi->nodes[1];
        for (auto x : c_preds.all) {
            Node *pred;
            z3::expr cond(ver->c);
            std::tie(pred, cond) = x;
            auto t = std::make_tuple(pred, rnode, cond);
            ver->ppos.push_back(t);
        }
        c_preds.all = {std::make_pair(wnode, wnode->guard)};
        break;
    }

    case Instruction::AtomicCmpXchg : {
        auto *mi = dyn_cast<MemoryInstrExecution>(this);
        Node *rnode = mi->nodes[0], *wnode = mi->nodes[1];
        for (auto x : c_preds.all) {
            Node *pred;
            z3::expr cond(ver->c);
            std::tie(pred, cond) = x;
            auto t = std::make_tuple(pred, rnode, cond);
            ver->ppos.push_back(t);
        }
        c_preds.all = {std::make_pair(rnode, rnode->guard),
                      std::make_pair(wnode, wnode->guard)};
        break;
    }

    case Instruction::Call : {
        snode_list nc_preds;
        if (auto *cie = dyn_cast<CallInstrExecution>(this)) {
            assert(!cie->called_execs.empty());
            for (auto f_exec : cie->called_execs) {
                f_exec->process_dppo(c_preds);
                for (auto spr : f_exec->c_terminators.all)
                    add_cond(nc_preds, spr.first, spr.second);
            }
            c_preds.all = nc_preds;
        } else if (auto *mie = dyn_cast<MemoryInstrExecution>(this)) {
            auto *node = mie->nodes[0];
            switch (node->kind) {
            case FLUSH_OPT :
                for (auto x : c_preds.all) {
                    Node *pred;
                    z3::expr cond(ver->c);
                    std::tie(pred, cond) = x;
                    auto t = std::make_tuple(pred, node, cond);
                    switch (pred->kind) {
                    case P_WRITE :
                    case FLUSH :
                    case FLUSH_OPT : {
                        auto sl = ver->prop.same_cl[pred->nid][node->nid];
                        ver->ppos.push_back(std::make_tuple(pred, node, cond && sl));
                        /* If there is a certain W/FL/FL|x -> FO|x edge, don't keep the events
                         * since the current FO|x suffices */
                        if (!sl.is_true())
                            nc_preds.push_back(x);
                        break;
                    }

                    case P_READ :
                    case SFENCE :
                    case MFENCE :
                    case EX_READ :
                    case EX_WRITE :
                        /* Keep it for a following FO */
                        ver->ppos.push_back(std::make_tuple(pred, node, cond));
                        nc_preds.push_back(x);
                        break;

                    default:
                        assert(false);
                    }
                }
                nc_preds.push_back(std::make_pair(node, node->guard));
                c_preds.all = nc_preds;
                break;

            case FLUSH :
            case SFENCE :
                for (auto x : c_preds.all) {
                    Node *pred;
                    z3::expr cond(ver->c);
                    std::tie(pred, cond) = x;
                    auto t = std::make_tuple(pred, node, cond);
                    ver->ppos.push_back(t);
                    nc_preds.push_back(x);
                }
                nc_preds.push_back(std::make_pair(node, node->guard));
                c_preds.all = nc_preds;
                break;

            case MFENCE :
                for (auto x : c_preds.all) {
                    Node *pred;
                    z3::expr cond(ver->c);
                    std::tie(pred, cond) = x;
                    auto t = std::make_tuple(pred, node, cond);
                    ver->ppos.push_back(t);
                }
                c_preds.all = {std::make_pair(node, node->guard)};
                break;

            default:
                assert(false);
            }
        } else {
            /* TODO : regular function */;
        }
        break;
    }

    case Instruction::Br : {
        auto *bri = cast<BranchInst>(instr);
        if (bri->isUnconditional()) {
            auto *tb = bri->getSuccessor(0);
            auto tb_exec = parent->parent->to_bbexec[tb];
            tb_exec->add_c_predecessors(c_preds);
        }else {
            for (unsigned i = 0; i < bri->getNumSuccessors(); i++) {
                auto *tb = bri->getSuccessor(i);
                auto bcond = ver->getEdgeCondition(bri->getParent(), tb);
                auto tb_exec = parent->parent->to_bbexec[tb];
                tb_exec->add_c_predecessors(c_preds, bcond);
            }
        }
        break;
    }

    case Instruction::Ret : {
        parent->parent->add_c_terminators(c_preds);
        break;
    }

    default:
        break;
    }
    return;
}

void FExecution::print_ce(z3::model m)
{
    std::cout << "Start function " << func->getName().str() << "\n";
    auto *be = bb_executions.front();
    while (true) {
        auto bb = be->print_ce(m);
        if (!bb)
            break;
        be = to_bbexec[bb];
    }
    std::cout << "End function " << func->getName().str() << "\n";
}

BasicBlock *BBExecution::print_ce(z3::model m)
{
    if (crash_encoding != GUARD)
        assert(m.eval(alive_in).is_true());
    for (auto *ie : instr_executions)
        ie->print(m);
    auto *tb = bb->getTerminator();
    if (auto *bi = dyn_cast<BranchInst>(tb)) {
        if (!bi->isUnconditional()) {
            auto *ti = cast<Instruction>(bi->getOperand(0));
            auto *tie = to_iexec[ti];
            auto bt = m.eval(tie->sval->val).as_uint64();
            /* true is the first, must flip */
            return bi->getSuccessor(1 - bt);
        } else {
            auto ss = bi->getSuccessor(0);
            if (ss != bb)
                return ss;
            else
                return nullptr;
        }
            return bi->getSuccessor(0);
    } else if (isa<ReturnInst>(tb)) {
        std::cout << "ret\n";
        return nullptr;
    }

    assert(false);
}

void InstrExecution::print(z3::model m)
{
    std::string istr;
    raw_string_ostream ss(istr);
    instr->print(ss);
    std::cout << istr << "\n";
    if (auto *cie = dyn_cast<CallInstrExecution>(this)) {
        for(auto c : cie->called_execs)
            c->print_ce(m);
        return;
    }

    if (auto *mie = dyn_cast<MemoryInstrExecution>(this)) {
        for (auto *n : mie->nodes)
            n->print(m);
    }

    if (sval)
        std::cout << "     " << m.eval(sval->val).as_uint64() << "\n";
}

z3::expr BBExecution::mk_guard_crash(z3::expr reached)
{
    if ((!check_pers) || (crash_encoding != GUARD))
        return reached;
    auto crash = ver->c.bool_const(ver->getUID("g_crash").c_str());
    auto guard = ver->c.bool_const(ver->getUID("guard").c_str());
    ver->s.add(guard == (reached && !crash));
    return guard;
}
