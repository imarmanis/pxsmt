#ifndef __EXECUTION_H__
#define __EXECUTION_H__

#include "z3++.h"
#include "graph.hpp"
#include "llvm/IR/Module.h"
#include "llvm/IR/Value.h"
#include "llvm/IR/Type.h"
#include "llvm/IR/GlobalVariable.h"

using namespace llvm;

typedef std::list<Node*> node_list;
typedef std::list<std::pair<Node*, z3::expr>> snode_list;

class Predecessors
{
public :
    snode_list all;
    /* used just for po+cl check */
    node_list writes;
    void merge_with(const Predecessors&);
    void merge_with(const Predecessors&, z3::expr);
};

class Verifier;
class FExecution;
class BBExecution;

class SymValue
{
public:
    z3::expr val;
    SymValue(z3::expr val);
};

class InstrExecution
{
public:
    enum Kind { K_INSTR, K_MEM, K_CALL };
    Kind kind;
    Instruction *instr;
    BBExecution *parent;
    Verifier *ver;
    SymValue *sval = nullptr;
    InstrExecution(Instruction *instr, BBExecution *parent);
    InstrExecution(Kind k, Instruction *instr, BBExecution *parent);
    InstrExecution() = delete;
    void process_dppo(Predecessors &c_preds);
    virtual ~InstrExecution() = default;
    void print(z3::model m);
    static bool classof(const InstrExecution *b) { return b->kind == K_INSTR; }
};

class MemoryInstrExecution : public InstrExecution
{
public:
    std::vector<Node*> nodes;
    MemoryInstrExecution(Instruction *instr, BBExecution *parent, Node *n);
    MemoryInstrExecution(Instruction *instr, BBExecution *parent, Node *n1, Node *n2);
    static bool classof(const InstrExecution *b) { return b->kind == K_MEM; }
};

class CallInstrExecution : public InstrExecution
{
public:
    std::vector<FExecution*> called_execs;
    CallInstrExecution(Instruction *instr, BBExecution *parent);
    static bool classof(const InstrExecution *b) { return b->kind == K_CALL; }
};

class BBExecution
{
public:
    BasicBlock *bb;
    FExecution *parent;
    Verifier *ver;
    z3::expr alive_in;
    z3::expr alive_out;
    BBExecution(BasicBlock *bb, FExecution *fexec);
    std::map<Instruction *, InstrExecution*> to_iexec;
    InstrExecution * newIExecution(Instruction *instr);
    InstrExecution * newIExecution(Instruction *instr, Node *n);
    InstrExecution * newIExecution(Instruction *instr, Node *n1, Node *n2);
    CallInstrExecution * newCallIExecution(Instruction *instr);
    void add_predecessors(node_list npreds);
    void add_c_predecessors(const Predecessors &cpreds);
    void add_c_predecessors(const Predecessors &cpreds, z3::expr);
    std::list<InstrExecution*> instr_executions;
    node_list preds;
    Predecessors c_preds;
    void process();
    std::pair<InstrExecution *, z3::expr> process_instruction(Instruction *i,
        node_list &preds, z3::expr reached);
    void process_dppo();
    BasicBlock *print_ce(z3::model m);
    z3::expr mk_guard_crash(z3::expr);
};

class FExecution
{
public:
    Function *func;
    Verifier *ver;
    z3::expr alive;
    z3::expr ret_val;
    bool rec = false;
    BBExecution * newBBExecution(BasicBlock *bb);
    std::map<BasicBlock *, BBExecution*> to_bbexec;
    /* Topollogically ordered */
    std::list<BBExecution*> bb_executions;
    node_list terminators;
    z3::expr_vector terminated;
    Predecessors c_terminators;
    void add_c_terminators(const Predecessors&);
    FExecution(Function *func, Verifier *ver, z3::expr alive, bool rec);
    void process(z3::expr_vector params, node_list preds);
    void process_dppo(Predecessors c_preds);
    void print_ce(z3::model m);
    z3::expr getEdgeCondition(BBExecution *bbh, BBExecution *bbt);
};

#endif
