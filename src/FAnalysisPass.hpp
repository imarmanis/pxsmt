#ifndef __F_ANALYSIS_PASS_H__
#define __F_ANALYSIS_PASS_H__

#include <llvm/Pass.h>
#include "AnalysisInfo.hpp"

class FAnalysisPass : public llvm::FunctionPass {

public:
	static char ID;

	FAnalysisPass(AnalysisInfo &analysis) :
		llvm::FunctionPass(ID), analysis(analysis) {}

	bool runOnFunction(llvm::Function &F) override;

	virtual void getAnalysisUsage(llvm::AnalysisUsage &au) const override;

private:
	AnalysisInfo &analysis;

};

#endif
