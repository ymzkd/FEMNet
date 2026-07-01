#include "FEAnalysis.h"

Displacement FEModeOperator::GetBeamDeform(int eid, int mode_id, double p)
{
    BeamElement* elm = model->GetBeamElement(eid);
    std::vector<std::vector<Displacement>> mv = ModeVectors();
    Displacement disp = elm->DisplaceAt(
        mv[mode_id][elm->Nodes[0]->id], mv[mode_id][elm->Nodes[1]->id], p);

    return disp;
}
