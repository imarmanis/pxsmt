#include "FAnalysisPass.hpp"
#include "llvm/Analysis/AliasAnalysis.h"
#include <llvm/Analysis/PostDominators.h>
#include <llvm/Analysis/BasicAliasAnalysis.h>
#include <llvm/Analysis/MemoryLocation.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/Function.h>
#include "AnalysisInfo.hpp"

using namespace llvm;

void FAnalysisPass::getAnalysisUsage(AnalysisUsage &AU) const
{
	FunctionPass::getAnalysisUsage(AU);
	AU.addRequiredTransitive<BasicAAWrapperPass>();
	AU.addRequiredTransitive<DominatorTreeWrapperPass>();
	AU.setPreservesAll();
}

bool FAnalysisPass::runOnFunction(Function &F)
{
	auto &AAB = getAnalysis<BasicAAWrapperPass>().getResult();
	std::vector<std::pair<Instruction *, MemoryLocation>> pairs;
	for (auto &b : F)
		for (auto &i : b) {
			auto opt = MemoryLocation::getOrNone(&i);
			if (opt.hasValue())
				pairs.push_back(std::make_pair(&i, opt.getValue()));
		}
	llvm::AAQueryInfo aaq;
	for (unsigned i = 0; i < pairs.size(); i++) {
		MemoryLocation a_ml;
		Instruction *a_i;
		std::tie(a_i, a_ml) = pairs[i];
		for (unsigned j = i + 1; j < pairs.size(); j++) {
			MemoryLocation b_ml;
			Instruction *b_i;
			std::tie(b_i, b_ml) = pairs[j];
			auto pair = std::make_pair(a_i, b_i);

			if (analysis.getAliasInfo(pair) != MayAlias)
				continue;

			if (a_i == b_i) {
				analysis.registerMustAlias(pair);
				continue;
			}

			auto aar_b = AAB.alias(a_ml, b_ml, aaq);
			if (aar_b == NoAlias)
				analysis.registerNoAlias(pair);
			else if (aar_b == MustAlias)
				analysis.registerMustAlias(pair);
		}
	}

	auto &DT = getAnalysis<DominatorTreeWrapperPass>().getDomTree();
	for (auto &b : F)
		for (auto &i : b)
			for (auto &bp : F)
				if (DT.dominates(&i, &bp))
					analysis.registerDomination(std::make_pair(&i, &bp));
	return false;
}

Pass *createFAnalysisPass(AnalysisInfo &analysis)
{
	auto *p = new FAnalysisPass(analysis);
	return p;
}

char FAnalysisPass::ID = 42;
