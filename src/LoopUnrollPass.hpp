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

#ifndef __LOOP_UNROLL_PASS_HPP__
#define __LOOP_UNROLL_PASS_HPP__

#include <llvm/Pass.h>

#include <llvm/ADT/StringRef.h>
#include <llvm/Analysis/LoopPass.h>
#include <llvm/Transforms/Utils/Cloning.h>

class LoopUnrollPass : public llvm::LoopPass {

protected:
	int unrollDepth;

	llvm::BasicBlock *makeDivergeBlock(llvm::Loop *l);
	void redirectBranch(int bodyIdx, int blockIdx, int unrollDepth,
			    llvm::BasicBlock *divergeBlock,
			    std::map<llvm::BasicBlock const *, int> &loopBlockIdx,
			    std::vector<std::vector<llvm::BasicBlock *> > &loopBodies);
	void redirectPHIOrValue(int bodyIdx, int blockIdx,
				std::vector<llvm::ValueToValueMapTy> &VMaps,
				std::map<llvm::BasicBlock const *, int> &loopBlockIdx,
				std::vector<std::vector<llvm::BasicBlock *> > &loopBodies);

public:
	static char ID;

	LoopUnrollPass(int depth) : llvm::LoopPass(ID), unrollDepth(depth) {
		if (unrollDepth < 0)
			unrollDepth = 0;
	};

	virtual llvm::StringRef getPassName() const { return "LoopUnrollPass"; } ;
	virtual void getAnalysisUsage(llvm::AnalysisUsage &au) const;
	virtual bool runOnLoop(llvm::Loop *l, llvm::LPPassManager &LPM);
};

#endif /* __LOOP_UNROLL_PASS_HPP__ */
