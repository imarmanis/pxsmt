#ifndef __INLINE_CLONE_PASS_H__
#define __INLINE_CLONE_PASS_H__

#include <llvm/Pass.h>
#include "llvm/Analysis/CallGraph.h"

using namespace llvm;

class InlineClonePass : public llvm::ModulePass {

public:
	static char ID;

    void cloneParFunctions(Function *F);
    bool hasLoop (const std::vector<CallGraphNode *> &nv);
    void handleFunction(Function *);

	InlineClonePass() : llvm::ModulePass(ID) {}

	bool runOnModule(llvm::Module &M) override;

private:

};

#endif
