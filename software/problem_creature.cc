#include "problem_creature.h"

ProblemCreature::ProblemCreature(ProblemDescription* pd, int massId, double t,
                                 double tStep) {
    this->problemDescription = pd;
    this->t = t;
    this->tStep = tStep;
    this->massId = massId;
    this->hasGetCostCache = false;

    // dna = [k1, k2, ..., kn, c1, c2, ..., cn]
    this->dna = std::vector<Bounded>(pd->NumberOfSpringsAndDampers());
    for (int i = 0; i < int(this->dna.size()); i++) {
        (this->dna.at(i)).Set(Random());
    }
}

double ProblemCreature::GetCost() {
    if (this->hasGetCostCache == true) {
        return this->getCostCache;
    }
    // TODO: A function that edits an existing Problem by setting a dna
    // could be slightly better, since BuildFromDNA creates the masses again,
    // for example, which we don't alter at the dna.
    auto p = this->problemDescription->BuildFromDNA(dna).val;
    auto e0 = p->Integrate(0.0, this->t, this->tStep);
    if (e0.isError) {
        throw("Unexpected Integrate error in GetCost");
    }
    auto e1 = p->GetMassMaxAbsAccel(this->massId);
    if (e1.isError) {
        throw("Unexpected GetMassMaxAbsAccel error in GetCost");
    }
    this->getCostCache = e1.val;
    this->hasGetCostCache = true;
    return e1.val;
}