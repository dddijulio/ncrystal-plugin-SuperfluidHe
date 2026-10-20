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

#include "NCPluginFactory.hh"
#include "NCPhysicsModel.hh"
#include "NCrystal/internal/freegas/NCFreeGas.hh"

namespace NCPluginNamespace {
  //The helium model up to 1 eV; above it NCrystal's free gas for He-4 at the same
  //temperature (its impulse limit: within 0.11% of the model at 1 eV, 0.1-1.3 K).
  class PluginScatter final : public NC::ProcImpl::ScatterIsotropicMat {
    struct Cache final : public NC::CacheBase {
      PhysicsModel::EnergyCache energy;
      NC::CachePtr freeGas;//the free gas keeps its own cache
      void invalidateCache() override { energy.energyEV=-1;freeGas.reset(); }
    };
    PhysicsModel m_model;
    NC::FreeGas m_freeGas;
  public:
    PluginScatter(PhysicsModel&& model,NC::Temperature temperature,const NC::AtomData& atom)
      :m_model(std::move(model)),m_freeGas(temperature,atom) {}
    const char* name() const noexcept override { return "SuperfluidHeFiniteWidth"; }
    NC::CrossSect crossSectionIsotropic(NC::CachePtr& ptr,NC::NeutronEnergy e) const override {
      auto& cache=accessCache<Cache>(ptr);
      if (e.dbl()>m_model.maxEnergyEV()) return m_freeGas.crossSectionIsotropic(cache.freeGas,e);
      return NC::CrossSect{m_model.calcCrossSection(e.dbl(),cache.energy)};
    }
    NC::ScatterOutcomeIsotropic sampleScatterIsotropic(NC::CachePtr& ptr,NC::RNG& rng,
                                                       NC::NeutronEnergy e) const override {
      auto& cache=accessCache<Cache>(ptr);
      if (e.dbl()>m_model.maxEnergyEV()) return m_freeGas.sampleScatterIsotropic(cache.freeGas,rng,e);
      auto out=m_model.sampleScatteringEvent(rng,e.dbl(),cache.energy);
      return {NC::NeutronEnergy{out.ekin_final},NC::CosineScatAngle{out.mu}};
    }
  protected:
    NC::Optional<std::string> specificJSONDescription() const override {
      return m_model.descriptionJSON();
    }
  };
}
const char* NCP::PluginFactory::name() const noexcept {return NCPLUGIN_NAME_CSTR "Factory";}
NC::Priority NCP::PluginFactory::query(const NC::FactImpl::ScatterRequest& cfg) const {
  if (!PhysicsModel::isApplicable(cfg.info())) return NC::Priority::Unable;
  return NC::Priority{999};
}
NC::ProcImpl::ProcPtr NCP::PluginFactory::produce(const NC::FactImpl::ScatterRequest& cfg) const {
  //Complete helium model (it has no elastic part). Do not append an elastic channel or the
  //standard SAB/free-gas machinery: that would double-count the supplied model. With the
  //inelastic component disabled (e.g. comp=elas) nothing scatters: NCrystal's own elastic
  //scattering from the DOS (that of a solid) does not apply to superfluid helium.
  if (cfg.get_inelas()=="0") return NC::ProcImpl::getGlobalNullScatter();
  const auto& info=cfg.info();
  auto model=PhysicsModel::createFromInfo(info);//checks the material before it is used below
  return NC::makeSO<PluginScatter>(std::move(model),info.getTemperature(),
                                   info.getComposition().front().atom.data());
}
