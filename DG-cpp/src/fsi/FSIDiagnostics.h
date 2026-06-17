#ifndef FSI_DIAGNOSTICS_H
#define FSI_DIAGNOSTICS_H

#include "EulerALE1D.h"

#include <string>

namespace euler_fsi {

struct FSIEnergySnapshot {
    double gasWallWork = 0.0;
    double dampingLoss = 0.0;
    double gasWorkDrift = 0.0;
    double coupledEnergyDrift = 0.0;
};

class FSIEnergyBudget {
public:
    void reset(double gasEnergy, const SpringPiston& piston);
    FSIEnergySnapshot advance(double dt, double gasEnergy, const SpringPiston& piston,
                              double pWall, double xOld, double vOld);

    double gasWallWork() const { return gasWallWork_; }
    double dampingLoss() const { return dampingLoss_; }
    double maxGasWorkDrift() const { return maxGasWorkDrift_; }
    double maxCoupledEnergyDrift() const { return maxCoupledEnergyDrift_; }

    static const char* csvHeader();
    std::string csvValues(const FSIEnergySnapshot& s) const;

private:
    double initialGasEnergy_ = 0.0;
    double initialCoupledEnergy_ = 0.0;
    double gasWallWork_ = 0.0;
    double dampingLoss_ = 0.0;
    double maxGasWorkDrift_ = 0.0;
    double maxCoupledEnergyDrift_ = 0.0;
};

} // namespace euler_fsi

#endif
