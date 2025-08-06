#include "MAnalysisPass.hpp"
#include "llvm/Analysis/AliasAnalysis.h"
#include <llvm/Analysis/PostDominators.h>
#include <llvm/Analysis/CFLAndersAliasAnalysis.h>
#include <llvm/Analysis/CFLSteensAliasAnalysis.h>
#include <llvm/Analysis/TypeBasedAliasAnalysis.h>
#include <llvm/Analysis/GlobalsModRef.h>
#include "llvm/Analysis/MemoryLocation.h"
#include <llvm/IR/Module.h>

using namespace llvm;

void MAnalysisPass::getAnalysisUsage(AnalysisUsage &AU) const
{
	ModulePass::getAnalysisUsage(AU);
	AU.addRequiredTransitive<CFLAndersAAWrapperPass>();
	AU.addRequiredTransitive<CFLSteensAAWrapperPass>();
	AU.addRequiredTransitive<TypeBasedAAWrapperPass>();
	AU.addRequiredTransitive<GlobalsAAWrapperPass>();
	AU.setPreservesAll();
}

bool MAnalysisPass::runOnModule(Module &M)
{
	auto &AAA = getAnalysis<CFLAndersAAWrapperPass>().getResult();
	auto &AAS = getAnalysis<CFLSteensAAWrapperPass>().getResult();
	auto &AAT = getAnalysis<TypeBasedAAWrapperPass>().getResult();
	auto &AAG = getAnalysis<GlobalsAAWrapperPass>().getResult();
	std::vector<std::pair<Instruction *, MemoryLocation>> pairs;
	for (auto &f : M)
		for (auto &b : f)
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

			auto aar = AAA.alias(a_ml, b_ml, aaq);
			if (aar == MayAlias)
				aar = AAS.alias(a_ml, b_ml, aaq);
			if (aar == MayAlias)
				aar = AAT.alias(a_ml, b_ml, aaq);
			if (aar == MayAlias)
				aar = AAG.alias(a_ml, b_ml, aaq);

			if (aar == NoAlias)
				analysis.registerNoAlias(pair);
			else if (aar == MustAlias)
				analysis.registerMustAlias(pair);
		}
	}

	return false;
}

Pass *createMAnalysisPass(AnalysisInfo &analysis)
{
	auto *p = new MAnalysisPass(analysis);
	return p;
}

char MAnalysisPass::ID = 42;
