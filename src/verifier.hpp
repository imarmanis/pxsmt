#ifndef _VERIFIER_HPP
#define _VERIFIER_HPP

#include "llvm/IR/Module.h"
#include "llvm/Analysis/CallGraph.h"
#include "map"
#include "deque"
#include "z3++.h"
#include "propagator.hpp"
#include "execution.hpp"
#include "graph.hpp"
#include "AnalysisInfo.hpp"

using namespace llvm;

#define BITS_IN_BYTE 8

class Verifier
{
private:
    z3::context c;
    z3::solver s;
    unsigned counter = 0;
    std::unique_ptr<Module> mod;
    std::map<Value*, SymValue*> scontext;
    std::list<std::tuple<Node*, Node*>> pois;
    std::list<std::tuple<Node*, Node*, z3::expr>> ppos;
    std::list<std::tuple<Node*, Node*>> static_rmws;
    Propagator prop;
    z3::expr_vector assumptions;
    z3::expr_vector assertions;
    z3::expr_vector no_diverge_assertions;
    DataLayout data_layout;
     /* Counter for allocas, in bytes. */
    unsigned alloc_counter;
    z3::expr symVal(Value *);
    z3::expr compare(CmpInst::Predicate,  const z3::expr &, const z3::expr &);
    z3::expr getEdgeCondition(BasicBlock *, BasicBlock *);
    unsigned getNewAllocPtr(unsigned size, unsigned alignment);
    std::string getUID(const char*);
    unsigned pointer_size(PointerType *);
    z3::expr gep_sval(Value *);
    std::map<StructType*,
        std::pair<z3::func_decl, z3::func_decl_vector>> custom_types;
    z3::sort to_sort(Type*);
    FExecution * newFExecution(Function *f, z3::expr alive, bool rec);
    /* Check whether node belongs to the recovery routine */
    bool in_rec(Node *);

public:
    SymValue * newSymValue(Value *val, z3::expr sval);
    std::map<Function *, FExecution*> to_fexec;
    z3::expr get_uninterpreted(Type *);
    z3::expr get_uninterpreted(z3::sort);
    z3::expr get_uninterpreted_bool();
    std::pair<z3::func_decl, z3::func_decl_vector> getCustom(StructType *);
    Verifier(std::unique_ptr<Module>, AnalysisInfo&);
    void print_graph(z3::model);
    friend class FExecution;
    friend class BBExecution;
    friend class InstrExecution;
};

#endif
