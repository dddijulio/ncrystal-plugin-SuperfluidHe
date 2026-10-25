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

#include "NCTestPlugin.hh"
#include "NCrystal/NCPluginBoilerplate.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include <cmath>
#include <string>

// Self-test, run by "ncrystal-pluginmanager --test SuperfluidHe" (which afterwards also loads
// every shipped file). Two shipped materials span the line widths: 0.1 K (built-in width, the
// narrowest line) and 1.3 K (input width, the broadest). The stored values come from the
// plugin itself (2026-10-19) and must be updated whenever the physics changes.
namespace {
  struct Case {
    const char* file;
    double sigma[13];//barn at 1e-11, 1e-10, ..., 1 eV and at 2 eV (NCrystal free gas)
    double meanEkin[3],meanMu[3];//at 1e-7, 1e-3 and 1e-2 eV: means of 2e6 sampled events
  };
  const Case cases[]={
    {"He4_Superfluid_T0.1K.ncmat",
#if NCRYSTAL_VERSION >= 4004007
     //Finer G1 binning in NCrystal 4.4.7 (VDOSGn) changes sigma below 1e-3 eV by up to 0.7%:
     {2.5629871349497987e-10,8.105111580804839e-11,2.5634147414061317e-11,8.118494494693076e-12,
      2.60625163260337e-12,9.544209614072062e-13,1.0171410867945624e-12,2.4971493121146246e-10,
      0.02845230394890344,0.7615001948595991,0.8428925198080393,0.8510808991367165,
      0.8519911827778801},
#else
     {2.549195322947928e-10,8.061268908028805e-11,2.5495965859337426e-11,8.074758785875311e-12,
      2.5921426373711426e-12,9.48981770763614e-13,1.0097289649386999e-12,2.478963964791386e-10,
      0.028452280614790503,0.7615001373435621,0.8428925408213456,0.8510809026377465,0.8519911827778801},
#endif
     {4.976978514195556e-05,0.00023550177447217034,0.006782611584951787},
     {-0.07287935240217178,0.7456612311438169,0.08855729811243404}},
    {"He4_Superfluid_T1.3K.ncmat",
     {0.8427145343466804,0.26649024300528185,0.08427319334958068,0.026654884036153563,
      0.008447022484334998,0.0026877817985224107,0.0008896052970701619,0.0005440303000820469,
      0.028857715038930306,0.7614228111513015,0.8428921560639818,0.8510809269593408,0.8519967332922316},
     {0.0004195318166732618,0.00025135110533765133,0.00678558088974789},
     {0.024125667322514924,0.7076375807949362,0.08863863701655006}}
  };
  const double sampledEnergies[3]={1e-7,1e-3,1e-2};
  const unsigned events=20000;

  // True when loading the material fails with an error message containing the given text.
  bool refused(std::string text,const char* expected) {
    try {
      NC::createScatter(NC::MatCfg::createFromRawData(std::move(text)));
    } catch (const NC::Error::BadInput& e) {
      return std::string(e.what()).find(expected)!=std::string::npos;
    }
    return false;
  }
}

void NCP::customPluginTest() {
  for (const auto& c:cases) {
    auto scatter=NC::createScatter(std::string("plugins::SuperfluidHe/")+c.file);
    //Cross sections within 1e-4, the accuracy sigma is computed to at these energies (fixed
    //table energies, or below the table).
    for (int k=0;k<13;++k) {
      double e=(k<12 ? std::pow(10.0,k-11) : 2.0);
      double xs=scatter.crossSectionIsotropic(NC::NeutronEnergy{e}).dbl();
      if (!(std::abs(xs/c.sigma[k]-1)<=1e-4))
        NCRYSTAL_THROW2(CalcError,"SuperfluidHe self-test: "<<c.file<<" sigma at "<<e<<" eV is "
                        <<xs<<" barn, stored "<<c.sigma[k]);
    }
    //The hand-over to NCrystal's free gas just above 1 eV.
    double below=scatter.crossSectionIsotropic(NC::NeutronEnergy{1.0}).dbl();
    double above=scatter.crossSectionIsotropic(NC::NeutronEnergy{1.000001}).dbl();
    if (!(std::abs(above/below-1)<=2e-3))
      NCRYSTAL_THROW2(CalcError,"SuperfluidHe self-test: "<<c.file<<" sigma steps from "<<below
                      <<" to "<<above<<" barn at 1 eV");
    //Sampling: every event physical; the means of E' and mu within 5 standard errors of the
    //stored ones (the stored means carry a 10 times smaller error).
    for (int j=0;j<3;++j) {
      auto e=NC::NeutronEnergy{sampledEnergies[j]};
      double sum[2]={0,0},squares[2]={0,0};
      for (unsigned i=0;i<events;++i) {
        auto event=scatter.sampleScatterIsotropic(e);
        double v[2]={event.ekin.dbl(),event.mu.dbl()};
        if (!(std::isfinite(v[0]) && v[0]>0 && std::isfinite(v[1]) && std::abs(v[1])<=1))
          NCRYSTAL_THROW2(CalcError,"SuperfluidHe self-test: "<<c.file<<" unphysical event at "
                          <<e.dbl()<<" eV: E'="<<v[0]<<" eV, mu="<<v[1]);
        for (int q=0;q<2;++q) {sum[q]+=v[q];squares[q]+=v[q]*v[q];}
      }
      const double stored[2]={c.meanEkin[j],c.meanMu[j]};
      for (int q=0;q<2;++q) {
        double mean=sum[q]/events,error=std::sqrt(std::max(0.0,squares[q]/events-mean*mean)/events);
        if (!(std::abs(mean-stored[q])<=5*error))
          NCRYSTAL_THROW2(CalcError,"SuperfluidHe self-test: "<<c.file<<" mean "<<(q ? "mu" : "E' (eV)")
                          <<" at "<<e.dbl()<<" eV is "<<mean<<" +- "<<error<<", stored "<<stored[q]);
      }
    }
  }
  //Temperatures outside 0.1-1.3 K are refused, also the 293.15 K NCrystal assumes for a file
  //without @TEMPERATURE (with the hint to set it in the cfg-string).
  auto text=NC::FactImpl::createTextData(NC::TextDataPath("plugins::SuperfluidHe/He4_Superfluid_T1.3K.ncmat"))
    ->rawDataCopy();
  const std::string temperature="@TEMPERATURE\n 1.3\n";
  auto at=text.find(temperature);
  if (at==std::string::npos)
    NCRYSTAL_THROW(CalcError,"SuperfluidHe self-test: no @TEMPERATURE 1.3 in the 1.3 K file");
  if (!refused(std::string(text).replace(at,temperature.size(),"@TEMPERATURE\n 1.4\n"),"0.1-1.3 K"))
    NCRYSTAL_THROW(CalcError,"SuperfluidHe self-test: a material at 1.4 K was not refused");
  if (!refused(std::string(text).erase(at,temperature.size()),";temp="))
    NCRYSTAL_THROW(CalcError,"SuperfluidHe self-test: a material without @TEMPERATURE was not refused");
  NCPLUGIN_MSG("Cross sections, free-gas hand-over, sampling and temperature checks passed");
}
