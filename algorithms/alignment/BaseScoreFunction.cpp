#include "BaseScoreFunction.hpp"

BaseScoreFunction::BaseScoreFunction() {
    ins = del = substitutionPrior = globalDeletionPrior = affineExtend = 0;
}

BaseScoreFunction::BaseScoreFunction(int insP, int delP, int subPriorP, 
        int delPriorP, int affineExtensionP)  {
    ins = insP;
    del = delP;
    substitutionPrior = subPriorP;
    globalDeletionPrior = delPriorP;
    affineExtend = affineExtensionP;
}
