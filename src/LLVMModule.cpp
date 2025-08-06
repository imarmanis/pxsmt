/*
 * GenMC -- Generic Model Checking.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, you can access it online at
 * http://www.gnu.org/licenses/gpl-3.0.html.
 *
 * Author: Michalis Kokologiannakis <michalis@mpi-sws.org>
 */

#include "Passes.hpp"
#include "LLVMModule.hpp"
#include <iostream>
#include <llvm/InitializePasses.h>
#include <llvm/IR/PassManager.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/IR/LLVMContext.h>
#include <llvm/IR/IRPrintingPasses.h>
#include <llvm/IR/Verifier.h>
#include <llvm/IRReader/IRReader.h>
#include <llvm/Support/MemoryBuffer.h>
#include <llvm/Support/SourceMgr.h>
#include <llvm/Support/Debug.h>
#include <llvm/Support/FileSystem.h>
#include <llvm/Support/raw_ostream.h>
#include <llvm/Transforms/IPO.h>
#include <llvm/Transforms/Scalar.h>
#include <llvm/Transforms/Utils.h>
#include <llvm/Analysis/AliasAnalysis.h>
#include "AnalysisInfo.hpp"

# define PassManager llvm::legacy::PassManager

namespace LLVMModule {

    llvm::LLVMContext *globalContext = nullptr;

	llvm::LLVMContext &getLLVMContext(void)
	{
		if (!globalContext)
			globalContext = new llvm::LLVMContext();
		return *globalContext;
	}


	std::unique_ptr<llvm::Module> getLLVMModule(std::string &src)
	{
		llvm::MemoryBuffer *buf;
		llvm::SMDiagnostic err;

		buf = llvm::MemoryBuffer::getMemBuffer(src, "", false).release();
		auto mod = llvm::parseIR(buf->getMemBufferRef(), err, getLLVMContext()).release();
		if (!mod) {
            std::cerr << "Could not parse LLVM IR!\n";
            err.print("err", llvm::dbgs());
            exit(1);
		}
		return std::unique_ptr<llvm::Module>(mod);
	}

	bool transformLLVMModule(llvm::Module &mod, int bound, AnalysisInfo &analysis)
	{
		llvm::PassRegistry &Registry = *llvm::PassRegistry::getPassRegistry();
		PassManager OptPM, BndPM, AnalysisPM;
		bool modified;

		llvm::initializeCore(Registry);
		llvm::initializeScalarOpts(Registry);
		llvm::initializeObjCARCOpts(Registry);
		llvm::initializeVectorization(Registry);
		llvm::initializeIPO(Registry);
		llvm::initializeAnalysis(Registry);
		llvm::initializeTransformUtils(Registry);
		llvm::initializeInstCombine(Registry);
		llvm::initializeInstrumentation(Registry);
		llvm::initializeTarget(Registry);

		OptPM.add(createDeclareInternalsPass());
		OptPM.add(llvm::createPromoteMemoryToRegisterPass());
		OptPM.add(llvm::createDeadArgEliminationPass());
        OptPM.add(createInlineClonePass());
		OptPM.add(llvm::createCFGSimplificationPass());

        if (bound >= 0)
            BndPM.add(createBLoopUnrollPass(bound));

		AnalysisPM.add(createMAnalysisPass(analysis));
		AnalysisPM.add(createFAnalysisPass(analysis));

		modified = BndPM.run(mod);

		modified |= OptPM.run(mod);

		modified |= AnalysisPM.run(mod);

		assert(!llvm::verifyModule(mod, &llvm::dbgs()));
		return modified;
	}

	void printLLVMModule(llvm::Module &mod, const std::string &out)
	{
		PassManager PM;
		std::error_code errs;

		llvm::raw_ostream *os = new llvm::raw_fd_ostream(out.c_str(), errs,
								 llvm::sys::fs::F_None);

		if (errs) {
			delete os;
            std::cerr << "Failed to write transformed module to file "
			     + out + ": " + errs.message();
			return;
		}
		PM.add(llvm::createPrintModulePass(*os));
		PM.run(mod);
		return;
	}

}
