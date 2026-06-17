#include "FSIDiagnostics.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace euler_fsi {

void FSIEnergyBudget::reset(double gasEnergy, const SpringPiston& piston) {
    initialGasEnergy_ = gasEnergy;
    initialCoupledEnergy_ = gasEnergy + piston.kineticEnergy() + piston.springEnergy()
                          + piston.externalPressurePotential();
    gasWallWork_ = 0.0;
    dampingLoss_ = 0.0;
    maxGasWorkDrift_ = 0.0;
    maxCoupledEnergyDrift_ = 0.0;
}

FSIEnergySnapshot FSIEnergyBudget::advance(double dt, double gasEnergy,
                                           const SpringPiston& piston,
                                           double pWall, double xOld, double vOld) {
    gasWallWork_ += piston.area * pWall * (piston.x - xOld);
    dampingLoss_ += 0.5 * piston.damping * (vOld * vOld + piston.v * piston.v) * dt;

    FSIEnergySnapshot s;
    s.gasWallWork = gasWallWork_;
    s.dampingLoss = dampingLoss_;
    s.gasWorkDrift = (gasEnergy + gasWallWork_) - initialGasEnergy_;
    double coupledEnergy = gasEnergy + piston.kineticEnergy() + piston.springEnergy()
                         + piston.externalPressurePotential() + dampingLoss_;
    s.coupledEnergyDrift = coupledEnergy - initialCoupledEnergy_;
    maxGasWorkDrift_ = std::max(maxGasWorkDrift_, std::abs(s.gasWorkDrift));
    maxCoupledEnergyDrift_ = std::max(maxCoupledEnergyDrift_, std::abs(s.coupledEnergyDrift));
    return s;
}

const char* FSIEnergyBudget::csvHeader() {
    return "gas_wall_work,damping_loss,gas_work_drift,coupled_energy_drift";
}

std::string FSIEnergyBudget::csvValues(const FSIEnergySnapshot& s) const {
    std::ostringstream out;
    out << std::setprecision(16)
        << s.gasWallWork << "," << s.dampingLoss << ","
        << s.gasWorkDrift << "," << s.coupledEnergyDrift;
    return out.str();
}

} // namespace euler_fsi
