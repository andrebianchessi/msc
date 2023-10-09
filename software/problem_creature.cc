#include "problem_creature.h"

ProblemCreature::ProblemCreature(ProblemDescription* pd, int massId, double T) {
    this->problemDescription = pd;
    this->T = T;
    this->massId = massId;
    this->hasGetCostCache = false;

    this->pimodels = nullptr;
    this->timesToCheck = -1;

    // dna = [k1, k2, ..., kn, c1, c2, ..., cn]
    this->dna = std::vector<Bounded>(pd->NumberOfSpringsAndDampers());
    for (int i = 0; i < int(this->dna.size()); i++) {
        (this->dna.at(i)).Set(Random());
    }
}

ProblemCreature::ProblemCreature(ProblemDescription* pd, int massId,
                                 Pimodels* pimodels, int timesToCheck) {
    this->problemDescription = pd;
    this->T = T;
    this->massId = massId;
    this->hasGetCostCache = false;

    this->pimodels = pimodels;
    this->timesToCheck = timesToCheck;

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
    if (this->pimodels == nullptr) {
        // TODO: A function that edits an existing Problem by setting a dna
        // could be slightly better, since BuildFromDNA creates the masses
        // again, for example, which we don't alter at the dna.
        auto p = this->problemDescription->BuildFromDNA(dna);
        if (p.isError) {
            throw std::runtime_error(
                "Unexpected BuildFromDNA error in GetCost");
        }
        auto ei = p.val.Integrate(this->T);
        if (ei.isError) {
            throw std::runtime_error("Unexpected Integrate error in GetCost");
        }
        auto e1 = p.val.GetMassMaxAbsAccel(this->massId);
        if (e1.isError) {
            throw std::runtime_error(
                "Unexpected GetMassMaxAbsAccel error in GetCost");
        }
        this->getCostCache = e1.val;
    } else {
        this->getCostCache =
            this->pimodels->GetMaxAbsAccel(massId, timesToCheck, dna);
    }
    this->hasGetCostCache = true;
    return this->getCostCache;
}