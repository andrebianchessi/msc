#include "problem_creature.h"

ProblemCreature::ProblemCreature(ProblemDescription* pd, int massId, double t,
                                 double tStep) {
    this->problemDescription = pd;
    this->t = t;
    this->tStep = tStep;
    this->massId = massId;

    // dna = [k1, k2, ..., kn, c1, c2, ..., cn]
    this->dna = std::vector<Bounded>(pd->NumberOfSpringsAndDampers());
    for (int i = 0; i < this->dna.size(); i++) {
        (this->dna.at(i)).Set(Random());
    }
    auto p = pd->BuildFromDNA(dna).val;
}

double ProblemCreature::GetCost() {
    // TODO: A function that edits an existing Problem by setting a dna
    // could be slightly better, since BuildFromDNA creates the masses again,
    // for example, which we don't alter at the dna.
    auto p = this->problemDescription->BuildFromDNA(dna).val;
    p->Integrate(0.0, this->t, this->tStep);
    return (p->GetMassMaxAbsAccel(this->massId)).val;
}