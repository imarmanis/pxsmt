#include "InlineClonePass.hpp"
#include "llvm/ADT/SCCIterator.h"
#include <llvm/Transforms/Utils/Cloning.h>
#include "llvm/IR/CFG.h"
#include <llvm/IR/Module.h>
#include "llvm/IR/Instructions.h"
#include <list>
#include <algorithm>

using namespace llvm;

void InlineClonePass::handleFunction(Function *F)
{
    std::list<CallInst*> calls;
    for (auto &b : *F)
        for (auto &i : b)
            if (auto *ci = dyn_cast<CallInst>(&i)) {
                auto func = ci->getCalledFunction();
                auto fname = func->getName().str();
                if (!fname.compare("__VERIFIER_parallel")) {
                    for (auto I = ci->arg_begin(), IE = ci->arg_end(); I != IE; ++I) {
                        auto func = cast<Function>(cast<Value>(I)->stripPointerCasts());
                        handleFunction(func);
                    }
                    calls.push_back(ci);
                } else if (fname.rfind("__VERIFIER_", 0)) {
                    handleFunction(func);
                    calls.push_back(ci);
                }
            }
    for  (auto ci : calls) {
        InlineFunctionInfo ifi;
        auto func = ci->getCalledFunction();
        auto fname = func->getName().str();
        if (!fname.compare("__VERIFIER_parallel")) {
            for (unsigned i = 0, ie = ci->getNumArgOperands(); i < ie; i++) {
                auto arg = ci->getArgOperand(i);
                auto func = cast<Function>(arg->stripPointerCasts());
                ValueToValueMapTy empty;
                auto clone = CloneFunction(func, empty);
                ci->setArgOperand(i, clone);
            }
        } else
            assert(InlineFunction(*ci, ifi).isSuccess() || ci->getName().rfind("llvm.x86"));
    }

}

bool InlineClonePass::hasLoop (const std::vector<CallGraphNode *> &nv)
{
    if (nv.size() > 1)
        return true;

    auto n = nv[0];
    for (auto x : *n)
        if (n == x.second)
            return true;

    return false;
}

bool InlineClonePass::runOnModule(Module &M)
{
    Function *main = nullptr, *rec = nullptr;
    CallGraph CG(M);
    for (auto I = llvm::scc_begin(&CG), IE = llvm::scc_end(&CG); I != IE; ++I) {
        auto fvec = *I;
        assert(!hasLoop(fvec));

        auto fnode = fvec[0];
        auto func = fnode->getFunction();
        if (!func)
            continue;

        auto fname = func->getName();
        if (!(fname.compare("main")))
            main = func;
        else if (!(fname.compare("__VERIFIER_recovery")))
            rec = func;
    }
    assert(main);

    handleFunction(main);
    if (rec)
        handleFunction(rec);

    return true;
}

Pass *createInlineClonePass()
{
	auto *p = new InlineClonePass();
	return p;
}

char InlineClonePass::ID = 42;
