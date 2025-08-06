#include "llvm/ADT/APInt.h"
#include "AnalysisInfo.hpp"

AliasResult AnalysisInfo::getAliasInfo(user_pair pair)
{
    auto opair = ordered(pair);
    auto must = must_alias.find(opair);
    if (must != must_alias.end())
        return MustAlias;
    auto no = no_alias.find(opair);
    if (no != no_alias.end())
        return NoAlias;
    return MayAlias;
}

bool AnalysisInfo::dominates(dpair pair)
{
    return dominator.find(pair) != dominator.end();
}

void AnalysisInfo::registerMustAlias(user_pair pair)
{
    must_alias.insert(ordered(pair));
}

void AnalysisInfo::registerNoAlias(user_pair pair)
{
    no_alias.insert(ordered(pair));
}

void AnalysisInfo::registerDomination(dpair pair)
{
    dominator.insert(pair);
}

user_pair AnalysisInfo::ordered(user_pair pair)
{
    return std::minmax(pair.first, pair.second);
}
