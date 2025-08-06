#ifndef __ANALYSIS_INFO_H__
#define __ANALYSIS_INFO_H__

#include <llvm/Analysis/MemoryLocation.h>
#include "llvm/Analysis/AliasAnalysis.h"
#include "set"

using namespace llvm;

typedef std::pair<User*, User*> user_pair;
typedef std::set<user_pair> alias_set;
typedef std::pair<Instruction*, BasicBlock*> dpair;

class AnalysisInfo
{
public:
    AliasResult getAliasInfo(user_pair);
    void registerMustAlias(user_pair);
    void registerNoAlias(user_pair);
    bool dominates(dpair);
    void registerDomination(dpair);
private:
    alias_set must_alias;
    alias_set no_alias;
    std::set<dpair> dominator;
    user_pair ordered(user_pair);

};

#endif
