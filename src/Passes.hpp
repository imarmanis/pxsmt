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

#ifndef __PASSES_HPP__
#define __PASSES_HPP__

#include <llvm/Pass.h>
#include "AnalysisInfo.hpp"

llvm::ModulePass *createDeclareInternalsPass();

llvm::Pass *createBLoopUnrollPass(int depth);

llvm::ModulePass *create();

llvm::ModulePass *createMAnalysisPass(AnalysisInfo &analysis);

llvm::FunctionPass *createFAnalysisPass(AnalysisInfo &analysis);

llvm::ModulePass *createInlineClonePass();

#endif /* __PASSES_HPP__ */
