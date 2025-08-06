#ifndef __MANALYSIS_PASS_H__
#define __MANALYSIS_PASS_H__

#include <llvm/Pass.h>
#include "AnalysisInfo.hpp"

class MAnalysisPass : public llvm::ModulePass {

public:
	static char ID;

	MAnalysisPass(AnalysisInfo &analysis) :
		llvm::ModulePass(ID), analysis(analysis) {}

	bool runOnModule(llvm::Module &M) override;

	virtual void getAnalysisUsage(llvm::AnalysisUsage &au) const override;

private:
	AnalysisInfo &analysis;

};

#endif
