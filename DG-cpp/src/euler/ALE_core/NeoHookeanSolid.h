#ifndef NEO_HOOKEAN_SOLID_H
#define NEO_HOOKEAN_SOLID_H

#include "ElasticSolid.h"

namespace euler_ale {

struct NeoHookeanMaterial {
    double density = 70.0;
    double thickness = 1.0;
    double young = 2.8e2;
    double poisson = 0.34;
    double damping = 0.55;

    static NeoHookeanMaterial fromSolidMaterial(const SolidMaterial& material);
    double mu() const;
    double lambda() const;
};

class NeoHookeanSolidModel {
public:
    explicit NeoHookeanSolidModel(const NeoHookeanMaterial& material);

    const NeoHookeanMaterial& material() const;
    MatrixXd internalForces(const ElasticSolid2D& solid,
                            double* strainEnergy = nullptr) const;
    double strainEnergy(const ElasticSolid2D& solid) const;
    double kineticEnergy(const ElasticSolid2D& solid) const;
    double stableTimeStep(const ElasticSolid2D& solid, double cfl) const;
    void advanceExplicit(ElasticSolid2D& solid, double dt) const;

private:
    NeoHookeanMaterial material_;
};

} // namespace euler_ale

#endif
