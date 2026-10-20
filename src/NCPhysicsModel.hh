////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef NCPlugin_PhysicsModel_hh
#define NCPlugin_PhysicsModel_hh
#include "NCrystal/NCPluginBoilerplate.hh"
#include <memory>
#include <vector>
namespace NCPluginNamespace {
  class PhysicsModel final {
  public:
    struct EnergyCache {
      struct QInterval { double lower,upper,mass,inverseMass=-1; };
      double energyEV=-1, sigma=0;
      std::vector<QInterval> intervals;
      std::vector<double> cumulative;
    };
    struct ScatEvent { double ekin_final,mu; };
    static bool isApplicable(const NC::Info&);
    static PhysicsModel createFromInfo(const NC::Info&);
    PhysicsModel(PhysicsModel&&) noexcept;
    PhysicsModel& operator=(PhysicsModel&&) noexcept;
    ~PhysicsModel();
    double maxEnergyEV() const noexcept;
    double calcCrossSection(double,EnergyCache&) const;
    ScatEvent sampleScatteringEvent(NC::RNG&,double,EnergyCache&) const;
    std::string descriptionJSON() const;
  private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
    explicit PhysicsModel(std::unique_ptr<Impl>);
  };
}
#endif
