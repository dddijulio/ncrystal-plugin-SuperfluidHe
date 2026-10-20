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

// SuperfluidHe: superfluid helium-4 model of Granada et al., NIM A 1053 (2023)
// 168284: a finite-width Lorentzian single-phonon line on the measured dispersion plus a
// DOS-derived multiphonon part (orders n>=2). Inside this file energies are meV and Q is
// 1/Aa; the public interface uses eV. Material data are immutable after setup (apart from
// a once-only warning flag); the incident-energy caches of the exact method belong to the
// NCrystal streams.
// Energy transfer follows the convention of Granada et al., NIM A 1053 (2023) 168284:
// hbar*omega = E_in-E_out (positive = neutron energy loss; the line sits at
// +hbar*omega_p(Q), their Eq. 2, and the second sum rule is their Eq. 12). NCrystal's G_n
// use E_out-E_in; they are converted where they enter (section 2) and checked at load.
//
// Sections:
//   1. Helpers
//   2. Multiphonon spectra G_n
//   3. Lorentzian line pieces
//   4. Model data and numerical settings
//   5. The kernel at one Q: u(Q) and Poisson weights
//   6. Masses in an energy-transfer window
//   7. Method A (exact): adaptive Q integration at one incident energy
//   8. Method B (transport above 1e-10 eV): energy table built at load
//   9. Energy-transfer sampling at a chosen Q
//  10. Setup and reading the input
//  11. Public interface: cross section and scattering event
//  12. Summary
#include "NCPhysicsModel.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/vdos/NCVDOSGn.hh"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <sstream>

namespace NCPluginNamespace {
namespace {
// ============================================================================ 1. Helpers
constexpr double cn=NC::ksq2ekin(1.0)*1000;//neutron E=cn*k^2 (meV, 1/Aa)
using Curve=std::vector<std::pair<double,double>>;
// Linear interpolation, constant beyond the end points.
double interp(const Curve& c,double x) {
  if (x<=c.front().first) return c.front().second;
  if (x>=c.back().first) return c.back().second;
  auto hi=std::upper_bound(c.begin(),c.end(),x,
    [](double a,const std::pair<double,double>& b){return a<b.first;});
  auto lo=hi-1;
  return lo->second+(hi->second-lo->second)*(x-lo->first)/(hi->first-lo->first);
}
double logistic(double x) {
  if (x>=0) return 1/(1+std::exp(-x));
  double e=std::exp(x); return e/(1+e);
}
double randOpen(NC::RNG& r) {
  double eps=std::numeric_limits<double>::epsilon();
  return NC::ncclamp(r.generate(),eps,1-eps);
}
// Gauss-Legendre points and weights on [-1,1].
std::pair<std::vector<double>,std::vector<double>> gauss(unsigned n) {
  std::vector<double> x(n),w(n);
  for (unsigned i=0;i<(n+1)/2;++i) {
    double z=std::cos(NC::kPi*(i+.75)/(n+.5)),dp=0;
    for (unsigned iter=0;iter<100;++iter) {
      double p0=1,p1=z;
      for (unsigned j=2;j<=n;++j) {
        double p=((2.0*j-1)*z*p1-(j-1.0)*p0)/j; p0=p1;p1=p;
      }
      dp=n*(z*p1-p0)/(z*z-1);
      double delta=p1/dp;z-=delta;
      if (std::abs(delta)<2e-15) break;
      if (iter==99) NCRYSTAL_THROW(CalcError,"Gauss quadrature failed");
    }
    x[i]=-z;x[n-1-i]=z;w[i]=w[n-1-i]=2/((1-z*z)*dp*dp);
  }
  return {std::move(x),std::move(w)};
}
// Mean order of a Poisson(w) distribution restricted to n>=2.
double meanN(double w) {
  if (w<1e-3) return 2+w/3+w*w/18+w*w*w/270;
  double a=-std::expm1(-w);return w*a/(a-w*std::exp(-w));
}
// Stable real roots of a*x^2+b*x+c (a!=0) strictly inside (lo,hi).
void quadraticCuts(std::vector<double>& cuts,double a,double b,double c,double lo,double hi) {
  auto add=[&](double x){if (std::isfinite(x) && x>lo && x<hi) cuts.push_back(x);};
  double disc=b*b-4*a*c;if (disc<0) return;
  double v=-.5*(b+std::copysign(std::sqrt(disc),b));
  if (v==0) {add(-b/(2*a));return;}
  add(v/a);add(c/v);
}
double abeMorishimaFWHM(double temperature) {
  // Abe and Morishima, NIM A 459 (2001), Eqs. A.6--A.7. T, Gamma_R and Delta are in
  // kelvin; A.1 identifies kB*Gamma_R as the small-Q energy FWHM (no factor of two).
  if (!(std::isfinite(temperature) && temperature>0))
    NCRYSTAL_THROW(BadInput,"Abe--Morishima linewidth requires a positive finite temperature");
  double prefactor=(1+.0603*std::sqrt(temperature))*std::sqrt(temperature);
  // s=(Delta(0)-Delta(T))/T solves s*exp(-s)=r. The branch continuous from T=0 has
  // 0<=s<=1 and exists only for r<=1/e.
  double logR=std::log(24.72)+std::log(prefactor)-std::log(temperature)-8.625/temperature;
  if (logR>-1)
    NCRYSTAL_THROW(BadInput,"No real low-temperature Abe--Morishima gap; supply fwhm_meV explicitly");
  double lower=0,upper=1;
  for (unsigned i=0;i<128;++i) {
    double mid=.5*(lower+upper);
    if (std::log(mid)-mid<logR) lower=mid;else upper=mid;
  }
  double width=NC::constant_boltzmann*1000*83.2*prefactor*std::exp(-8.625/temperature+.5*(lower+upper));
  if (!(std::isfinite(width) && width>0))
    NCRYSTAL_THROW(BadInput,"Abe--Morishima FWHM is outside the positive representable range");
  return width;
}

// ============================================================================ 2. Multiphonon spectra G_n
// NCrystal's n-phonon spectra G_n (A. Sjolander, Ark. Fys. 14 (1958) 315; the tau_n of
// Granada et al., NIM A 1053 (2023) 168284, Eq. 7), flipped to the sign convention of
// Granada et al. and stored as normalised piecewise-linear densities, with the cumulative
// and inverse used for the multiphonon masses and sampling.
struct Spectrum {
  double emin,step,mean=0;
  std::vector<double> density,cdf;
  Spectrum(const NC::VDOSGn& gn,unsigned n,double kT) {
    const auto& raw=gn.getRawSpectrum(n);auto range=gn.eRange(n);
    step=gn.binWidth(n)*1000;
    //VDOSGn's grid weights have area step*sum(Gn). Zero-padded linear interpolation keeps
    //exactly that area and first moment (trapezoids without padding would lose half of
    //each endpoint weight and break the scaling with order).
    emin=-range.second*1000-step;
    density.assign(raw.size()+2,0.0);
    std::reverse_copy(raw.begin(),raw.end(),density.begin()+1);
    double peak=*std::max_element(density.begin(),density.end());
    if (*std::min_element(density.begin(),density.end()) < -1e-12*peak)
      NCRYSTAL_THROW(CalcError,"Significant negative Gn spectrum");
    for (auto& d:density) d=std::max(0.0,d);
    // FFT roundoff can swamp the exponentially small energy-gain side: below the
    // precision floor, rebuild it from the loss side by detailed balance at the grid points.
    double precisionFloor=64*std::numeric_limits<double>::epsilon()*peak;
    for (std::size_t i=0;i<density.size();++i) {
      double e=emin+i*step;
      if (e>=0) break;
      double x=(-e-emin)/step,value=0;
      if (x>=0 && x<density.size()-1) {
        auto j=static_cast<std::size_t>(x);double fraction=x-j;
        value=(1-fraction)*density[j]+fraction*density[j+1];
      }
      value*=std::exp(e/kT);
      if (value<precisionFloor) density[i]=value;
    }
    cdf.resize(density.size(),0);
    for (std::size_t i=1;i<density.size();++i)
      cdf[i]=cdf[i-1]+.5*step*(density[i-1]+density[i]);
    double norm=cdf.back();
    if (!(norm>0)) NCRYSTAL_THROW(CalcError,"Empty Gn spectrum");
    for (auto& d:density) d/=norm;
    for (auto& v:cdf) v/=norm;
    cdf.back()=1;
    for (std::size_t i=1;i<density.size();++i) {
      double e0=emin+(i-1)*step,d0=density[i-1],d1=density[i];
      //Exact first moment of the piecewise-linear density.
      mean+=step*(.5*e0*(d0+d1)+step*(d0/6+d1/3));
    }
  }
  double cumulative(double e) const {
    if (e<=emin) return 0;
    if (e>=emin+step*(density.size()-1)) return 1;
    auto i=std::min(density.size()-2,static_cast<std::size_t>((e-emin)/step));
    double d=e-(emin+i*step);
    return cdf[i]+density[i]*d+.5*(density[i+1]-density[i])*d*d/step;
  }
  double inverse(double p) const {
    p=std::max(0.0,std::min(1.0,p));
    if (p==1) return emin+step*(density.size()-1);
    auto it=std::upper_bound(cdf.begin(),cdf.end(),p);
    auto i=std::min(density.size()-2,static_cast<std::size_t>(it-cdf.begin()-1));
    double area=p-cdf[i],slope=(density[i+1]-density[i])/step;
    double denom=density[i]+std::sqrt(std::max(0.0,density[i]*density[i]+2*slope*area));
    double delta=denom>0 ? 2*area/denom : 0;
    return emin+i*step+std::max(0.0,std::min(step,delta));
  }
};

// ============================================================================ 3. Lorentzian line pieces
// The single-phonon Lorentzian of Granada et al., NIM A 1053 (2023) 168284, Eq. 2 (centre
// ep, HWHM gamma), handled as uniform in the angle atan((E-ep)/gamma): an interval
// [lower,upper] is stored as its start angle and angular span; mass is filled in by the caller.
struct LineInterval {
  double a=0,span=0,lower=0,upper=0,mass=0;
  int side=0;//0: interval crosses the peak; +/-1: above/below it (complementary tail angle)
};
struct LineWindow {
  std::vector<LineInterval> intervals;
  double mass=0;
};
LineInterval lorentzWindow(double ep,double gamma,double lower,double upper) {
  LineInterval r;r.lower=lower;r.upper=upper;
  if (!(upper>lower)) return r;
  //atan((upper-ep)/gamma)-atan((lower-ep)/gamma) without subtracting nearly equal
  //+/-pi/2 values; atan2 also gives the correct branch across the peak.
  r.span=std::atan2(gamma*(upper-lower),gamma*gamma+(lower-ep)*(upper-ep));
  if (lower>=ep) {r.side=+1;r.a=std::atan2(gamma,upper-ep);}
  else if (upper<=ep) {r.side=-1;r.a=std::atan2(gamma,ep-lower);}
  else r.a=std::atan2(lower-ep,gamma);
  return r;
}
// Energy at the fraction 0..1 of the interval's angular span.
double lorentzEnergy(double ep,double gamma,const LineInterval& win,double fraction) {
  double angle=win.a+win.span*fraction;
  //In a far tail the angle is O(gamma/distance), so gamma/tan(angle) keeps the energy
  //even where atan(distance/gamma) rounds to +/-pi/2.
  double e=win.side ? ep+win.side*gamma/std::tan(angle) : ep+gamma*std::tan(angle);
  return std::max(win.lower,std::min(win.upper,e));
}

// Kernel quantities at one Q (section 5) and Q-density values (sections 6-8).
struct Weights { unsigned first=2;std::vector<double> values; };//orders first, first+1, ...
struct QPoint { double q,z,h,ep,w;Weights weights;double lineNorm; };//w: Poisson mean
struct QMass {
  double single=0,multi=0;
  double total() const {return single+multi;}
};
}

// ============================================================================ 4. Model data and numerical settings
// Everything a loaded material holds: the input, what setup derives from it, and the fixed
// numerical settings. The energy table and the warning flag are declared with their code in
// section 8.
struct PhysicsModel::Impl {
  // Material input and fixed model settings.
  double temperature=0,kT=0,gamma=0;//gamma: HWHM of the line (meV), Q independent
  double anchor=.7,emax=1000;//u(Q) anchor (1/Aa); upper end of the model (meV), also where the line ends
  bool widthFromInput=false;
  unsigned qOrder=16;
  double qTolerance=1e-8;//relative tolerance of method A in transport
  Curve z,h,phonons;
  // Derived at setup.
  double msd=0,recoil=0,sigmaBound=0,staticMax=0;
  double seedMean=0,seedLaplaceUpper=0,anchorU=0,uSlope=0;
  unsigned maxOrder=2;
  std::vector<Spectrum> spectra;//order n at index n-1
  std::vector<double> inputQ;//input Q points: all curves are linear between neighbours
  std::vector<double> uGridQ,uGridW;//Poisson mean w(Q) from the anchor to qmax (u(Q) table)
  std::pair<std::vector<double>,std::vector<double>> qFine,qCoarse;
  std::vector<QPoint> points;//Q grid: input Q points plus qOrder Gauss points between them
  std::vector<QPoint> coarsePoints;//the coarse rule's points on the same intervals
  // Numerical settings (accuracy choices; loop-count safety limits stay inline).
  static constexpr unsigned lineRuleOrder=8;//Gauss points per panel of the line mass
  static constexpr double tailRatio=8;//ratio of successive Lorentzian-tail cuts
  static constexpr double tailTolerance=1e-6;//relative cost of the unresolved inner tail
  static constexpr double uTableStep=5e-4;//u(Q) table spacing between input Q points (1/Aa)
  static constexpr double gainSkipTolerance=1e-8;//see components()
  static constexpr double tableMinEV=1e-10;//table from here to emax, method A below
  static constexpr double tablePerDecade=10;//initial table energies per decade
  static constexpr unsigned tableDepth=8;//bisection levels at most
  static constexpr double tableTolerance=1e-3;//bisect while log-log sigma at a midpoint is off by more
  static constexpr double tableSigmaTolerance=1e-4;//method A tolerance of sigma at each table energy
  static constexpr double tableBound=4;//rejection bound: Q density <= tableBound*proposal
  static constexpr double cutNodeNear=1e-7;//Q points at q_cut*(1-+cutNodeNear*4^k),
  static constexpr unsigned cutNodeLevels=9;//k<cutNodeLevels (out to 0.66%), see gridNode

  // ========================================================================== 5. The kernel at one Q
  // point(Q) collects what S(Q,E) needs at one Q: Z, H, the line centre and normalisation,
  // and the multiphonon Poisson mean w(Q) with its order weights. u(Q) (Granada et al., NIM A
  // 1053 (2023) 168284, Eq. 11) is set from their second sum rule (Eq. 12) at and above the
  // anchor and continued below it (see u()).

  // First moment Z*<E*tanh(E/2kT)> of the normalised line (0..emax), on panels as
  // in lineMass (section 6) except that the tail cuts run on to the line's end: there the
  // weight grows like E, so the far tails add about (gamma/pi)*ln(emax/Ep), which
  // a single Gauss rule over the whole line would miss.
  double lineMoment(double q) const {
    double strength=interp(z,q);if (strength==0) return 0;
    double ep=interp(phonons,q),top=emax;
    static const auto rule=gauss(lineRuleOrder);
    double cuts[128];int n=0;
    cuts[n++]=0;cuts[n++]=top;
    auto add=[&](double e){if (e>0 && e<top && n<128) cuts[n++]=e;};
    for (double f=.25;f<=64;f*=2) add(f*kT);
    double reach=std::max(ep,top-ep),logStart=tailTolerance*kT/gamma;//ln(d/gamma), as in lineMass
    if (logStart<std::log(reach/gamma))
      for (double d=gamma*std::exp(logStart);d<reach;d*=tailRatio) {add(ep-d);add(ep+d);}
    std::sort(cuts,cuts+n);
    double sum=0;
    for (int i=1;i<n;++i) {
      if (!(cuts[i]>cuts[i-1])) continue;
      auto win=lorentzWindow(ep,gamma,cuts[i-1],cuts[i]);
      if (!(win.span>0)) continue;
      double s=0;
      for (std::size_t j=0;j<rule.first.size();++j) {
        double e=lorentzEnergy(ep,gamma,win,.5*(1+rule.first[j]));
        s+=rule.second[j]*e*std::tanh(e/(2*kT));
      }
      sum+=.5*win.span*s;
    }
    return strength*sum/lorentzWindow(ep,gamma,0.0,top).span;
  }
  // Poisson mean w(Q) of the multiphonon orders from the second sum rule at Q (H(Q)>0).
  double solveW(double q) const {
    double target=(recoil*q*q-lineMoment(q))/interp(h,q)/seedMean;
    if (!(target>2)) NCRYSTAL_THROW2(BadInput,"No u(Q) solution at Q="<<q
      <<" 1/Aa (second sum rule); check the input curves");
    if (!(target<50)) return target;
    double lo=0,hi=target;
    for (unsigned j=0;j<60;++j) {
      double mid=.5*(lo+hi);
      if (meanN(mid)<target) lo=mid;else hi=mid;
    }
    return .5*(lo+hi);
  }
  // Why an anchor: the second sum rule fixes the mean phonon order of the multiphonon part,
  // (recoil*Q^2-line moment)/(H*<G_1>), but the expansion starts at n=2. Below Q~0.69/Aa
  // (shipped data) the rule asks for less than 2, so no u satisfies it. At and above the
  // anchor u=msd*Q^2/(H*w), with w(Q) interpolated from the table built at load (w is smooth
  // even where u is steep just above the anchor) or solved beyond it. Below the anchor u
  // continues exponentially with the anchor's value and slope: u becomes very large (1e13-1e15
  // at Q=0.15) and w tiny, so the multiphonon part is in effect the pure two-phonon G_2, as
  // close to the rule as the expansion gets (first moment up to 37% high near Q=0.3).
  double u(double q) const {
    if (q<anchor) return anchorU*std::exp(uSlope/anchorU*(q-anchor));
    double hv=interp(h,q);if (hv<=0) return 1;
    double w;
    if (uGridQ.size()<2 || q>uGridQ.back()) w=solveW(q);
    else {
      auto it=std::upper_bound(uGridQ.begin(),uGridQ.end(),q);
      std::size_t k=std::min(uGridQ.size()-2,static_cast<std::size_t>(it-uGridQ.begin()-1));
      w=uGridW[k]+(uGridW[k+1]-uGridW[k])*(q-uGridQ[k])/(uGridQ[k+1]-uGridQ[k]);
    }
    return msd*q*q/(hv*w);
  }
  // Poisson(w) weights of the orders that matter, conditioned on n>=2.
  Weights weights(double w) const {
    Weights r;if (w==0) {r.values={1};return r;}
    r.first=static_cast<unsigned>(std::max(2.0,std::floor(w-10*std::sqrt(w+1)-12)));
    unsigned last=static_cast<unsigned>(std::ceil(w+10*std::sqrt(w+1)+12));
    if (last>maxOrder) NCRYSTAL_THROW(CalcError,"Insufficient phonon orders");
    double v=std::exp(r.first*std::log(w)-w-std::lgamma(r.first+1.0)),sum=0;
    for (unsigned n=r.first;n<=last;++n) {r.values.push_back(v);sum+=v;v*=w/(n+1.0);}
    if (!(sum>0)) NCRYSTAL_THROW(CalcError,"Poisson weights underflowed");
    for (auto& p:r.values) p/=sum;
    return r;
  }
  QPoint point(double q) const {
    double hv=interp(h,q),w=(q>0 && hv>0) ? msd*q*q/(u(q)*hv) : 0;
    double ep=interp(phonons,q),normalisation=lorentzWindow(ep,gamma,0.0,emax).span;
    //Safety net for other inputs: as H->0, w=msd*Q^2/(u*H) can need more orders than were
    //built (not reached with the shipped inputs). Such many-phonon events are inaccessible
    //below emax: their whole accessible probability is bounded by the
    //compound-Poisson Laplace transform before they are omitted.
    if (std::ceil(w+10*std::sqrt(w+1)+12)>maxOrder) {
      double upper=2*std::sqrt(cn*emax)*q-cn*q*q;
      double logBound=upper/(2*kT)+w*(seedLaplaceUpper-1);
      double denominator=1-(1+w)*std::exp(-w);
      if (denominator>0 && logBound-std::log(denominator)<std::log(1e-14))
        return {q,interp(z,q),hv,ep,w,Weights{},normalisation};
      NCRYSTAL_THROW(CalcError,"Uncovered phonon orders with non-negligible accessible weight");
    }
    return {q,interp(z,q),hv,ep,w,weights(w),normalisation};
  }

  // ========================================================================== 6. Masses in an energy-transfer window
  // How much of S(Q,E) lies in an energy-transfer window at one Q. Single phonon: the
  // Lorentzian of Granada et al., NIM A 1053 (2023) 168284, Eq. 2, truncated to excitation
  // energies 0..emax (the most a neutron in the model can lose) and normalised,
  // times their detailed-balance factor 1/(1+exp(-hbar*omega/kT)): logistic(+E/kT) for loss
  // (transfer +E), logistic(-E/kT) for gain (transfer -E). Multiphonon: the Poisson mixture
  // of G_n. Masses are unscaled (Z and H multiply them).

  // Line mass of one sign for excitation energies in [lower,upper], one method for every
  // width: the Gauss rule in the Lorentzian angle on panels cut where the integrand changes.
  // The thermal factor changes on the kT scale from the window's lower edge (cuts at
  // lower+kT/4..64kT, doubling). The Lorentzian tails change on the scale of the distance to
  // the peak (cuts at Ep+/-d*tailRatio^j out to the window's far end, from the distance to
  // the window when the peak is outside it). For a peak inside the window, tails closer than
  // d are left unresolved, costing about (gamma/kT)*ln(d/gamma) relative: d is where that
  // reaches tailTolerance, and a line that narrow needs no tail cuts at all.
  double lineMass(const QPoint& p,double lower,double upper,int sign) const {
    lower=std::max(lower,0.0);upper=std::min(upper,emax);
    if (!(upper>lower) || p.z==0) return 0;
    static const auto rule=gauss(lineRuleOrder);
    double cuts[128];int n=0;
    cuts[n++]=lower;cuts[n++]=upper;
    auto add=[&](double e){if (e>lower && e<upper && n<128) cuts[n++]=e;};
    for (double f=.25;f<=64;f*=2) add(lower+f*kT);
    double distance=p.ep<lower ? lower-p.ep : (p.ep>upper ? p.ep-upper : 0);
    double reach=std::max(p.ep-lower,upper-p.ep);
    if (distance>0) {
      for (double d=distance*tailRatio;d<reach;d*=tailRatio) {add(p.ep-d);add(p.ep+d);}
    } else {
      double logStart=tailTolerance*kT/gamma;//ln(d/gamma)
      if (logStart<std::log(64*kT/gamma))
        for (double d=gamma*std::exp(logStart);d<reach && d<=64*kT;d*=tailRatio) {add(p.ep-d);add(p.ep+d);}
    }
    std::sort(cuts,cuts+n);
    double sum=0;
    for (int i=1;i<n;++i) {
      if (!(cuts[i]>cuts[i-1])) continue;
      auto win=lorentzWindow(p.ep,gamma,cuts[i-1],cuts[i]);
      if (!(win.span>0)) continue;
      double s=0;
      for (std::size_t j=0;j<rule.first.size();++j)
        s+=rule.second[j]*logistic(sign*lorentzEnergy(p.ep,gamma,win,.5*(1+rule.first[j]))/kT);
      sum+=.5*win.span*s;
    }
    return sum/p.lineNorm;
  }
  // Upper bound on the energy-gain line mass over [lower,upper]: logistic(-E/kT) decreases,
  // so its value at lower bounds it up to lower+20kT and its value there bounds the rest.
  double gainBound(const QPoint& p,double lower,double upper) const {
    lower=std::max(lower,0.0);upper=std::min(upper,emax);
    if (!(upper>lower) || p.z==0) return 0;
    double split=std::min(upper,lower+20*kT);
    double bound=logistic(-lower/kT)*lorentzWindow(p.ep,gamma,lower,split).span;
    if (upper>split) bound+=logistic(-split/kT)*lorentzWindow(p.ep,gamma,split,upper).span;
    return bound/p.lineNorm;
  }
  // Component masses for transfers in [lower,upper]. The cross section and the event
  // sampler both use this one function. With logistic(E/kT)=1-logistic(-E/kT) the energy
  // loss over excitation energies [0,upper] is the line's span there minus the gain integral
  // over the same energies, which is also the first part of the energy-gain window
  // [max(0,-upper),-lower]. The rest of that window is integrated only where its bound
  // exceeds gainSkipTolerance of the other mass at this Q (always if there is none);
  // gainTop is where the kept gain window ends, so sampling uses the same support.
  struct Components {
    double positive=0,negative=0,multi=0;//line loss, line gain, multiphonon
    double gainTop=0;//upper excitation energy of the gain window kept
    std::vector<double> orderCDF;//cumulative multiphonon mass by order
  };
  Components components(const QPoint& p,double lower,double upper,bool orders=false) const {
    Components c;
    if (p.h>0) for (std::size_t j=0;j<p.weights.values.size();++j) {
      const auto& s=spectra[p.weights.first+j-1];
      c.multi+=p.weights.values[j]*std::max(0.0,s.cumulative(upper)-s.cumulative(lower));
      if (orders) c.orderCDF.push_back(c.multi);
    }
    if (p.z==0) return c;
    double shared=0;
    if (upper>0) {
      double top=std::min(upper,emax);
      shared=lineMass(p,0.0,upper,-1);
      c.positive=(top>0 ? lorentzWindow(p.ep,gamma,0.0,top).span/p.lineNorm : 0)-shared;
    }
    double from=std::abs(upper),others=p.z*c.positive+p.h*c.multi,rest=0;//rest starts at |upper|
    c.gainTop=-lower;
    if (-lower>from) {
      if (!(others>0) || p.z*gainBound(p,from,-lower)>gainSkipTolerance*others) rest=lineMass(p,from,-lower,-1);
      else c.gainTop=from;
    }
    c.negative=shared+rest;
    return c;
  }
  // Q density at incident wave number ki: Q times the masses over the kinematic window, so
  // sigma = sigma_b/(2 k_i^2) * integral of qDensity dQ, and Q is sampled from it.
  QMass qDensity(const QPoint& p,double ki) const {
    //These forms avoid subtracting E from nearly equal kinetic energies.
    double q=p.q;
    auto c=components(p,-cn*q*(q+2*ki),cn*q*(2*ki-q));
    return {q*(p.z*(c.positive+c.negative)),q*(p.h*c.multi)};
  }

  // ========================================================================== 7. Method A (exact)
  // Method A: sigma(E) = sigma_b/(2 k_i^2) * integral over Q of Q*(masses in the kinematic
  // window), by adaptive Gauss integration on Q intervals split at every kinematic cut of the
  // line; the intervals are kept for sampling Q by inverting the integral. Used to build the
  // energy table and below its first energy (1e-10 eV).
  QMass qDensity(double q,double ki) const {
    // Reuse only exactly identical Q grid points, never an interpolated point.
    for (const auto* prepared:{&points,&coarsePoints}) {
      auto it=std::lower_bound(prepared->begin(),prepared->end(),q,
        [](const QPoint& p,double value){return p.q<value;});
      if (it!=prepared->end() && it->q==q) return qDensity(*it,ki);
    }
    return qDensity(point(q),ki);
  }
  // Break points: Q values at this incident energy where qDensity can jump: the input Q
  // points, 2k_i, and the kinematic cuts, where the line centre (loss +eps or gain -eps)
  // meets a kinematic limit. The dispersion is linear between neighbouring input Q points, so
  // the cuts are quadratic roots; a narrow line then switches on or off at an interval
  // boundary, never inside one.
  std::vector<double> qBreakpoints(double ki) const {
    std::vector<double> cuts=inputQ;
    if (2*ki<inputQ.back()) cuts.push_back(2*ki);
    for (std::size_t i=1;i<inputQ.size();++i) {
      double lo=inputQ[i-1],hi=inputQ[i];
      if (interp(z,lo)==0 && interp(z,hi)==0) continue;
      double el=interp(phonons,lo),er=interp(phonons,hi);
      double slope=(er-el)/(hi-lo),intercept=el-slope*lo;
      for (int sign:{-1,1}) for (int direction:{-1,1})
        quadraticCuts(cuts,-cn,direction*2*cn*ki-sign*slope,
                      -sign*intercept,lo,hi);
    }
    std::sort(cuts.begin(),cuts.end());
    cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    return cuts;
  }
  // Rigorous bound on the integral of Q*mass(Q) over [q0,inf) for q0>=2ki: only energy
  // gains of at least cn*Q*(Q-2ki) remain, with thermal factor at most exp(-gain/kT);
  // integrated Z+H is at most staticMax, and a factor 10 covers the approximate detailed
  // balance of the numerical G_n. With x=Q-ki>=ki, Q<=2x bounds the integral by
  // exp(-a*q0*(q0-2ki))/a, a=cn/kT.
  double gainTailBound(double q0,double ki) const {
    double a=cn/kT;
    return 10*staticMax*std::exp(-a*q0*(q0-2*ki))/a;
  }
  QMass qQuadrature(double lo,double hi,double ki,
      const std::pair<std::vector<double>,std::vector<double>>& rule) const {
    double half=.5*(hi-lo),mid=lo+half;QMass result;
    for (std::size_t j=0;j<rule.first.size();++j) {
      auto value=qDensity(mid+half*rule.first[j],ki);
      result.single+=rule.second[j]*value.single;
      result.multi+=rule.second[j]*value.multi;
    }
    result.single*=half;result.multi*=half;return result;
  }
  // Adaptive bisection of [lo,hi]: fine against coarse rule per cell. The sum of the cell
  // errors is controlled (not a share per cell), so a linewidth-sized boundary layer of
  // negligible total error needs no unattainable local accuracy.
  // The fine-coarse difference is a reliable error estimate at the table's tolerance of 1e-4
  // (the actual error is at most 0.42 of it at 0.1-1.3 K); at 1e-8 the actual error can
  // reach ~6e-7.
  QMass integrateQ(double lo,double hi,double ki,double tolerance,QMass absolute,
                   unsigned depth=0,PhysicsModel::EnergyCache* cache=nullptr,
                   const QMass* initialFine=nullptr) const {
    struct Cell {double lo,hi;QMass value,error;unsigned depth;};
    auto evaluate=[&](double a,double b,unsigned level,bool initial) {
      auto fine=initial && initialFine ? *initialFine : qQuadrature(a,b,ki,qFine);
      auto coarse=qQuadrature(a,b,ki,qCoarse);
      if (!std::isfinite(fine.total()) || !std::isfinite(coarse.total()))
        NCRYSTAL_THROW(CalcError,"Nonfinite Q integrand");
      return Cell{a,b,fine,{std::abs(fine.single-coarse.single),
                           std::abs(fine.multi-coarse.multi)},level};
    };
    std::vector<Cell> cells{evaluate(lo,hi,depth,true)};
    for (;;) {
      QMass total,error;
      for (const auto& cell:cells) {
        total.single+=cell.value.single;total.multi+=cell.value.multi;
        error.single+=cell.error.single;error.multi+=cell.error.multi;
      }
      double singleLimit=absolute.single+tolerance*total.single;
      double multiLimit=absolute.multi+tolerance*total.multi;
      if (error.single<=singleLimit && error.multi<=multiLimit) {
        if (cache) {
          std::sort(cells.begin(),cells.end(),[](const Cell& a,const Cell& b){return a.lo<b.lo;});
          for (const auto& cell:cells) if (cell.value.total()>0) {
            cache->intervals.push_back({cell.lo,cell.hi,cell.value.total()});
            cache->cumulative.push_back(cache->cumulative.back()+cell.value.total());
          }
        }
        return total;
      }
      std::size_t worst=cells.size();double score=-1;
      for (std::size_t i=0;i<cells.size();++i) {
        const auto& cell=cells[i];
        if (cell.depth>=48 || cell.hi-cell.lo<=32*std::numeric_limits<double>::epsilon()
            *std::max(1.0,std::max(std::abs(cell.lo),std::abs(cell.hi)))) continue;
        double candidate=std::max(cell.error.single/std::max(1e-300,singleLimit),
                                  cell.error.multi/std::max(1e-300,multiLimit));
        if (candidate>score) {score=candidate;worst=i;}
      }
      if (worst==cells.size() || cells.size()>=16384)
        NCRYSTAL_THROW(CalcError,"Adaptive Q integral could not meet its summed error tolerance");
      auto cell=cells[worst];double mid=cell.lo+.5*(cell.hi-cell.lo);
      cells[worst]=evaluate(cell.lo,mid,cell.depth+1,false);
      cells.push_back(evaluate(mid,cell.hi,cell.depth+1,false));
    }
  }
  // Sigma and the Q intervals at one incident energy; tolerance is relative to the total;
  // sigma is 0 outside 0 < E <= emax (above it the factory uses the free gas).
  void fillCache(double energyEV,PhysicsModel::EnergyCache& c,double tolerance) const {
    if (c.energyEV==energyEV) return;
    c.energyEV=-1;c.sigma=0;
    c.intervals.clear();c.cumulative.assign(1,0);
    if (energyEV>0 && energyEV<=.001*emax) {
      double ki=std::sqrt(energyEV*1000/cn);auto breaks=qBreakpoints(ki);
      QMass estimate,total;std::vector<QMass> initial;
      initial.reserve(breaks.size()-1);
      for (std::size_t i=1;i<breaks.size();++i) {
        auto value=qQuadrature(breaks[i-1],breaks[i],ki,qFine);
        initial.push_back(value);
        estimate.single+=value.single;estimate.multi+=value.multi;
        //Stop once the whole remaining Q range is rigorously below 1e-3 of the requested
        //error: larger Q only reach suppressed energy gains.
        if (breaks[i]>=2*ki && gainTailBound(breaks[i],ki)<=1e-3*tolerance*estimate.total()) {
          breaks.resize(i+1);break;
        }
      }
      //A small absolute budget, shared relative to the total, stops roundoff-level
      //intervals from forcing unlimited refinement.
      double budget=std::max(1e-300,.1*tolerance*estimate.total()/(breaks.size()-1));
      QMass absolute{budget,budget};
      for (std::size_t i=1;i<breaks.size();++i) {
        auto value=integrateQ(breaks[i-1],breaks[i],ki,tolerance,absolute,0,&c,&initial[i-1]);
        total.single+=value.single;total.multi+=value.multi;
      }
      double factor=sigmaBound/(2*ki*ki);
      c.sigma=factor*total.single+factor*total.multi;
    }
    c.energyEV=energyEV;
  }
  // Q from the cached intervals: pick an interval by mass, then a safeguarded Newton
  // inversion of the integrated Q density inside it (same quadrature as the cross section).
  double sampleQ(NC::RNG& rng,double ki,PhysicsModel::EnergyCache& cache) const {
    double target=randOpen(rng)*cache.cumulative.back();
    auto it=std::upper_bound(cache.cumulative.begin(),cache.cumulative.end(),target);
    auto i=std::min(cache.intervals.size()-1,
      static_cast<std::size_t>(it-cache.cumulative.begin()-1));
    auto& interval=cache.intervals[i];
    //Quantile error 1e-8 of the interval mass: far below any sampling statistics.
    double tolerance=std::max(1e-300,interval.mass*1e-8);
    double integration=std::min(qTolerance,1e-9);
    QMass absolute{.025*tolerance,.025*tolerance};
    if (interval.inverseMass<0)
      interval.inverseMass=integrateQ(interval.lower,interval.upper,ki,
                                    integration,absolute).total();
    if (!(interval.inverseMass>0))
      NCRYSTAL_THROW(CalcError,"Empty Q interval in inverse cumulative distribution");
    double fraction=(target-cache.cumulative[i])/interval.mass;
    double area=fraction*interval.inverseMass;
    double lo=interval.lower,hi=interval.upper;
    double x=lo+(hi-lo)*fraction;
    for (unsigned iteration=0;iteration<60;++iteration) {
      auto mass=integrateQ(interval.lower,x,ki,integration,absolute);
      double difference=mass.total()-area;
      if (std::abs(difference)<=tolerance) return x;
      if (difference<0) lo=x;else hi=x;
      if (hi-lo<=8*std::numeric_limits<double>::epsilon()*std::max(1.0,std::abs(x)))
        return lo+.5*(hi-lo);
      double density=qDensity(x,ki).total();
      double next=density>0 ? x-difference/density : lo+.5*(hi-lo);
      x=(next>lo && next<hi) ? next : lo+.5*(hi-lo);
    }
    NCRYSTAL_THROW(CalcError,"Q inverse cumulative distribution did not converge");
  }

  // ========================================================================== 8. Method B (energy table)
  // Transport above 1e-10 eV: a table over incident energy, built once at load. Each
  // table energy holds sigma from method A and, as the proposal for sampling Q, the Q density
  // on the Q grid (piecewise linear). Sigma is interpolated log-log between table energies;
  // Q is drawn from the two neighbouring table energies and corrected by rejection (section
  // 11).
  struct TableNode {double e=0,sigma=0;std::vector<double> q,f,cumulative;};
  std::vector<TableNode> table;//in increasing energy
  mutable std::atomic<bool> boundWarned{false};//tableBound exceeded while sampling (warn once)
  // One table energy. Around each kinematic cut of the line (other than an input Q point),
  // Q points at q_cut*(1-+cutNodeNear*4^k) on both sides follow the line's step at any
  // width.
  TableNode gridNode(double energyEV) const {
    TableNode n;n.e=energyEV;n.cumulative.assign(1,0);
    double e=energyEV*1000,ki=std::sqrt(e/cn);
    std::vector<double> cuts;
    for (double x:qBreakpoints(ki)) if (!std::binary_search(inputQ.begin(),inputQ.end(),x)) cuts.push_back(x);
    std::size_t next=0;double single=0,multi=0,prevS=0,prevM=0;
    auto add=[&](double q,const QPoint& p) {
      if (!n.q.empty() && !(q>n.q.back())) return;//keep Q points strictly increasing
      auto parts=components(p,-cn*q*(q+2*ki),cn*q*(2*ki-q));
      double s=q*p.z*(parts.positive+parts.negative),mu=q*p.h*parts.multi;
      if (!n.q.empty()) {
        double w=.5*(q-n.q.back());
        single+=w*(prevS+s);multi+=w*(prevM+mu);
        n.cumulative.push_back(single+multi);
      }
      n.q.push_back(q);n.f.push_back(s+mu);prevS=s;prevM=mu;
    };
    for (const auto& p:points) {
      while (next<cuts.size() && cuts[next]<p.q) {
        double x=cuts[next++],limit=next<cuts.size() ? std::min(p.q,cuts[next]) : p.q;
        for (unsigned k=cutNodeLevels;k-->0;) {double y=x*(1-cutNodeNear*std::pow(4.,k));add(y,point(y));}
        for (unsigned k=0;k<cutNodeLevels;++k) {
          double y=x*(1+cutNodeNear*std::pow(4.,k));
          if (y<limit) add(y,point(y));
        }
        //A piece between two cuts with no Q grid point inside gets 7 points.
        if (next<cuts.size() && cuts[next]<p.q)
          for (unsigned k=1;k<8;++k) {double y=x+(cuts[next]-x)*k/8;add(y,point(y));}
      }
      add(p.q,p);
      if (p.q>=2*ki && single+multi>0 && gainTailBound(p.q,ki)<=1e-12*(single+multi)) break;
    }
    //The Q density stays finite as Q->0 (Lorentzian tail with the Bose factor of soft
    //phonons), while this proposal would vanish at the point Q=0: start it flat at the
    //first interior point instead. Proposal only; sigma is below.
    if (n.q.size()>1 && n.q.front()==0) {
      n.f.front()=n.f[1];
      for (std::size_t i=1;i<n.q.size();++i)
        n.cumulative[i]=n.cumulative[i-1]+.5*(n.q[i]-n.q[i-1])*(n.f[i]+n.f[i-1]);
    }
    PhysicsModel::EnergyCache exact;fillCache(energyEV,exact,tableSigmaTolerance);
    n.sigma=exact.sigma;
    return n;
  }
  // Table energies from tableMinEV to emax, tablePerDecade per decade; each interval is
  // bisected at its geometric midpoint until log-log interpolation reproduces sigma there
  // within tableTolerance, at most tableDepth levels deep.
  void buildTable() {
    double emaxEV=.001*emax,lo=std::min(tableMinEV,.5*emaxEV);
    unsigned n=std::max(1u,static_cast<unsigned>(std::ceil(tablePerDecade*std::log10(emaxEV/lo))));
    std::vector<TableNode> nodes;
    for (unsigned i=0;i<=n;++i) nodes.push_back(gridNode(i==n ? emaxEV : lo*std::pow(emaxEV/lo,double(i)/n)));
    std::function<void(double,double,double,double,unsigned)> refine=
      [&](double ea,double sa,double eb,double sb,unsigned depth) {
        double em=std::sqrt(ea*eb);
        nodes.push_back(gridNode(em));double sm=nodes.back().sigma;
        double guess=(sa>0 && sb>0) ? std::sqrt(sa*sb) : .5*(sa+sb);
        if (depth<tableDepth && std::abs(sm-guess)>tableTolerance*std::max(sm,guess)) {
          refine(ea,sa,em,sm,depth+1);refine(em,sm,eb,sb,depth+1);
        }
      };
    for (unsigned i=0;i<n;++i) {
      double ea=nodes[i].e,sa=nodes[i].sigma,eb=nodes[i+1].e,sb=nodes[i+1].sigma;
      refine(ea,sa,eb,sb,0);
    }
    std::sort(nodes.begin(),nodes.end(),[](const TableNode& a,const TableNode& b){return a.e<b.e;});
    table=std::move(nodes);
  }
  // Log-log interpolation; i is the lower table energy and w the weight of table energy i+1.
  double tableSigma(double energyEV,std::size_t& i,double& w) const {
    i=0;w=0;
    if (!(energyEV>=table.front().e && energyEV<=.001*emax)) return 0;
    auto it=std::upper_bound(table.begin(),table.end(),energyEV,
      [](double v,const TableNode& n){return v<n.e;});
    i=std::min<std::size_t>(table.size()-2,static_cast<std::size_t>(it-table.begin()-1));
    const auto& a=table[i];const auto& b=table[i+1];
    w=std::log(energyEV/a.e)/std::log(b.e/a.e);
    if (a.sigma>0 && b.sigma>0) return a.sigma*std::exp(w*std::log(b.sigma/a.sigma));
    return a.sigma+(b.sigma-a.sigma)*w;
  }
  // The normalised proposal density of one table energy at q, and a draw from it.
  static double nodeDensity(const TableNode& n,double q) {
    if (q<n.q.front() || q>n.q.back()) return 0;
    auto it=std::upper_bound(n.q.begin(),n.q.end(),q);
    std::size_t k=std::min(n.q.size()-2,static_cast<std::size_t>(it-n.q.begin()-1));
    double t=(q-n.q[k])/(n.q[k+1]-n.q[k]);
    return (n.f[k]+(n.f[k+1]-n.f[k])*t)/n.cumulative.back();
  }
  static double sampleNodeQ(NC::RNG& rng,const TableNode& n) {
    double target=randOpen(rng)*n.cumulative.back();
    auto it=std::upper_bound(n.cumulative.begin(),n.cumulative.end(),target);
    auto i=std::min(n.q.size()-2,static_cast<std::size_t>(it-n.cumulative.begin()-1));
    double width=n.q[i+1]-n.q[i],area=target-n.cumulative[i],f0=n.f[i];
    double slope=(n.f[i+1]-f0)/width;
    double denom=f0+std::sqrt(std::max(0.0,f0*f0+2*slope*area));
    return n.q[i]+std::max(0.0,std::min(width,denom>0 ? 2*area/denom : 0));
  }

  // ========================================================================== 9. Energy-transfer sampling
  // Given Q, the energy transfer is drawn exactly from the kernel inside the kinematic
  // window: single or multiphonon by their masses; for the line, loss or gain, an interval
  // by mass, then a Lorentzian energy accepted with the detailed-balance factor; for
  // multiphonon, an order by its Poisson weight, then an energy from that G_n.

  // The line in [lower,upper] split into intervals with their masses: one for energy
  // loss; bands of at most 2 kT for energy gain, so the thermal tail is resolved wherever
  // the window starts and the in-band rejection costs at most exp(2). The last band keeps
  // the rest of the window (beyond 746 kT the Bose factor is zero in double precision).
  LineWindow lineWindow(const QPoint& p,double lower,double upper,int sign) const {
    LineWindow r;
    lower=std::max(lower,0.0);upper=std::min(upper,emax);
    if (!(upper>lower) || p.z==0) return r;
    auto integrate=[&](double lo,double hi) {
      auto interval=lorentzWindow(p.ep,gamma,lo,hi);
      interval.mass=lineMass(p,lo,hi,sign);
      r.mass+=interval.mass;
      if (interval.mass>0) r.intervals.push_back(interval);
    };
    if (sign>0) integrate(lower,upper);
    else {
      double left=lower,thermalEnd=std::min(upper,746*kT);
      while (left<thermalEnd) {
        double right=std::min(thermalEnd,left+2*kT);
        if (!(right>left)) break;
        integrate(left,right);left=right;
      }
      if (upper>left) integrate(left,upper);
    }
    return r;
  }
  // Excitation energy from the line: interval by mass, then Lorentzian energy accepted
  // with the detailed-balance factor against the interval's largest one.
  double sampleLine(NC::RNG& rng,const QPoint& p,const LineWindow& win,int sign) const {
    if (win.intervals.empty()) NCRYSTAL_THROW(CalcError,"Empty single-phonon interval");
    double target=randOpen(rng)*win.mass,cumulative=0;
    const LineInterval* interval=&win.intervals.back();
    for (const auto& candidate:win.intervals) {
      cumulative+=candidate.mass;
      if (target<cumulative) {interval=&candidate;break;}
    }
    double envelope=logistic(sign*(sign>0 ? interval->upper : interval->lower)/kT);
    if (!(envelope>0)) NCRYSTAL_THROW(CalcError,"Empty single-phonon envelope");
    for (unsigned attempt=0;attempt<10000;++attempt) {
      double e=lorentzEnergy(p.ep,gamma,*interval,randOpen(rng));
      if (randOpen(rng)*envelope<=logistic(sign*e/kT)) return sign*e;
    }
    NCRYSTAL_THROW(CalcError,"Single-phonon rejection sampler did not converge");
  }
  // Energy transfer at Q from the kernel's conditional distribution in [lower,upper]:
  // single phonon (loss or gain) or multiphonon (order, then G_n inside the window).
  double sampleTransfer(NC::RNG& rng,const QPoint& p,const Components& parts,double lower,double upper) const {
    double single=p.z*(parts.positive+parts.negative),multi=p.h*parts.multi;
    if (randOpen(rng)*(single+multi)<single) {
      bool plus=randOpen(rng)*(parts.positive+parts.negative)<parts.positive;
      auto selected=lineWindow(p,plus ? lower : -upper,plus ? upper : parts.gainTop,plus ? +1 : -1);
      return sampleLine(rng,p,selected,plus ? +1 : -1);
    }
    const auto& orderCDF=parts.orderCDF;
    double value=randOpen(rng)*orderCDF.back();
    auto j=std::min(orderCDF.size()-1,static_cast<std::size_t>(
      std::upper_bound(orderCDF.begin(),orderCDF.end(),value)-orderCDF.begin()));
    const auto& s=spectra[p.weights.first+j-1];
    double a=s.cumulative(lower),b=s.cumulative(upper);
    return s.inverse(a+(b-a)*randOpen(rng));
  }

  // ========================================================================== 10. Setup
  // Everything derived from the input, in order: checks and quadrature rules; the DOS from
  // NCrystal's VDOS input (MSD, G_1); u(Q) at the anchor; the input Q points and the u(Q)
  // table; the Q grid, the phonon orders needed and their G_n; the kernel at the Q grid
  // points; the energy table.
  void setup(const NC::AtomData& atom,const NC::DI_VDOS& vdosInput) {
    sigmaBound=atom.scatteringXS().dbl();recoil=cn*NC::const_neutron_mass_amu/atom.averageMassAMU().dbl();
    if (z.front().first!=0 || h.front().first!=0 || phonons.front().first!=0)
      NCRYSTAL_THROW(BadInput,"Q curves must start at zero");
    for (const auto& v:z) if (v.first>phonons.back().first && v.second!=0)
      NCRYSTAL_THROW(BadInput,"Nonzero Z beyond dispersion coverage");
    //Largest integrated strength Z+H: both are piecewise linear, constant beyond.
    for (const auto* c:{&z,&h}) for (const auto& v:*c)
      staticMax=std::max(staticMax,interp(z,v.first)+interp(h,v.first));
    qFine=gauss(std::max(4u,qOrder));
    qCoarse=gauss(std::max(2u,qOrder/2));
    //The DOS as NCrystal read and regularised it (@DYNINFO type vdos).
    const auto& vdos=vdosInput.vdosData();
    NC::VDOSEval eval(vdos);msd=eval.getMSD();
    NC::VDOSGn gn(eval);spectra.emplace_back(gn,1,kT);seedMean=spectra.front().mean;
    if (!(seedMean>0)) NCRYSTAL_THROW(CalcError,"Wrong Gn transfer convention");
    const auto& seed=spectra.front();
    for (std::size_t i=1;i<seed.density.size();++i) {
      double area=.5*seed.step*(seed.density[i-1]+seed.density[i]);
      if (area>0) seedLaplaceUpper+=area*std::exp(-(seed.emin+(i-1)*seed.step)/(2*kT));
    }
    if (!(seedLaplaceUpper<1))
      NCRYSTAL_THROW(CalcError,"Cannot bound inaccessible high-order tail");
    anchorU=u(anchor);uSlope=(u(anchor+.001)-anchorU)/.001;
    //Input Q points: the input curves' Q values, 0, the anchor, and 0.5/Aa steps beyond the
    //data up to the largest Q reachable at emax (with 30 kT of energy gain).
    double dataQ=std::min(z.back().first,h.back().first),ki=std::sqrt(emax/cn);
    if (!(anchor<dataQ))
      NCRYSTAL_THROW2(BadInput,"Z(Q) and H(Q) must extend beyond Q="<<anchor<<" 1/Aa");
    double qmax=std::max(dataQ,ki+std::sqrt((emax+30*kT)/cn));
    if (qmax>dataQ && (z.back().second!=0 || std::abs(h.back().second-1)>1e-6))
      NCRYSTAL_THROW(BadInput,"High-Q extension requires terminal H=1 and Z=0");
    inputQ={0,anchor,qmax};
    for (const auto* c:{&z,&h,&phonons})
      for (const auto& v:*c) if (v.first<=qmax) inputQ.push_back(v.first);
    unsigned next=static_cast<unsigned>(std::ceil((qmax-dataQ)/.5));
    for (unsigned i=0;i<=next;++i) inputQ.push_back(next ? dataQ+(qmax-dataQ)*i/next : dataQ);
    std::sort(inputQ.begin(),inputQ.end());inputQ.erase(std::unique(inputQ.begin(),inputQ.end()),inputQ.end());
    //w(Q) once, from the anchor to qmax: each interval between input Q points subdivided to uTableStep.
    for (std::size_t i=1;i<inputQ.size();++i) {
      if (inputQ[i]<=anchor) continue;
      double a=std::max(anchor,inputQ[i-1]),b=inputQ[i];
      unsigned m=std::max(1u,static_cast<unsigned>(std::ceil((b-a)/uTableStep)));
      for (unsigned j=uGridQ.empty() ? 0 : 1;j<=m;++j) {
        double q=a+(b-a)*j/m;
        uGridQ.push_back(q);uGridW.push_back(interp(h,q)>0 ? solveW(q) : 0);
      }
    }
    //The Q grid uses the integration's fine rule, so method A reuses these kernel values.
    std::vector<double> grid;
    for (std::size_t i=1;i<inputQ.size();++i) {
      grid.push_back(inputQ[i-1]);double half=.5*(inputQ[i]-inputQ[i-1]),mid=inputQ[i-1]+half;
      for (double x:qFine.first) grid.push_back(mid+half*x);
    }
    grid.push_back(qmax);
    //Phonon orders: enough for the largest Poisson mean on the grid; G_n generated once.
    double wmax=1;
    for (double q:grid) {
      double hv=interp(h,q);
      if (q>0 && hv>0) wmax=std::max(wmax,msd*q*q/(u(q)*hv));
    }
    double nmax=std::ceil(wmax+12*std::sqrt(wmax+1)+24);
    if (nmax>=99999) NCRYSTAL_THROW(BadInput,"Too many phonon orders requested");
    maxOrder=static_cast<unsigned>(nmax);
    gn.growMaxOrder(maxOrder);spectra.reserve(maxOrder);
    for (unsigned n=2;n<=maxOrder;++n) spectra.emplace_back(gn,n,kT);
    //Each G_n must have n times the mean of G_1 (checks the numerical convolutions).
    for (unsigned n=1;n<=maxOrder;++n) {
      double error=std::abs(spectra[n-1].mean/(n*seedMean)-1);
      if (!(error<=1e-7))
        NCRYSTAL_THROW2(CalcError,"Gn moments do not scale with order: relative error "
          <<error<<" at n="<<n<<" (limit 1e-7)");
    }
    for (double q:grid) points.push_back(point(q));
    coarsePoints.reserve((inputQ.size()-1)*qCoarse.first.size());
    for (std::size_t i=1;i<inputQ.size();++i) {
      double half=.5*(inputQ[i]-inputQ[i-1]),mid=inputQ[i-1]+half;
      for (double x:qCoarse.first) coarsePoints.push_back(point(mid+half*x));
    }
    buildTable();
  }
};

// ============================================================================ 10. Setup (continued): reading the input
bool PhysicsModel::isApplicable(const NC::Info& info) {
  return info.countCustomSections(pluginNameUpperCase())>0;
}
PhysicsModel PhysicsModel::createFromInfo(const NC::Info& info) {
  if (info.countCustomSections(pluginNameUpperCase())!=1)
    NCRYSTAL_THROW(BadInput,"Need exactly one @CUSTOM_SUPERFLUIDHE section");
  if (!info.hasTemperature() || info.getComposition().size()!=1)
    NCRYSTAL_THROW(BadInput,"Need temperature and pure He4 composition");
  const auto& atom=info.getComposition().front().atom.data();
  if (!atom.isSingleIsotope() || atom.Z()!=2 || atom.A()!=4)
    NCRYSTAL_THROW(BadInput,"SuperfluidHe requires pure He4");
  //Lines of two words: a curve name and its point count followed by that many rows of
  //two values, or a scalar option and its value.
  static const std::set<std::string> curveNames{"z_q","h_q","phonons_meV"};
  std::set<std::string> options{"fwhm_meV"};
  auto data=info.getCustomSection(pluginNameUpperCase());
  std::map<std::string,double> scalars;std::map<std::string,Curve> curves;
  for (std::size_t i=0;i<data.size();++i) {
    if (data[i].size()!=2) NCRYSTAL_THROW(BadInput,"Expected two words in custom section");
    const auto& name=data[i][0];double v;
    if (!NC::safe_str2dbl(data[i][1],v) || !std::isfinite(v))
      NCRYSTAL_THROW(BadInput,"Invalid model value");
    if (curveNames.count(name)) {
      if (curves.count(name) || v<2 || v>100000 || std::floor(v)!=v)
        NCRYSTAL_THROW(BadInput,"Duplicate curve or invalid point count");
      auto count=static_cast<std::size_t>(v);
      if (count>data.size()-i-1) NCRYSTAL_THROW(BadInput,"Incomplete model curve");
      auto& c=curves[name];
      for (std::size_t j=0;j<count;++j) {
        const auto& row=data[++i];double x,y;
        if (row.size()!=2 || !NC::safe_str2dbl(row[0],x) || !NC::safe_str2dbl(row[1],y)
            || !std::isfinite(x) || !std::isfinite(y) || x<0
            || (!c.empty() && x<=c.back().first) || y<0)
          NCRYSTAL_THROW(BadInput,"Invalid or unordered curve");
        c.emplace_back(x,y);
      }
    } else {
      if (!options.count(name)) NCRYSTAL_THROW2(BadInput,"Unknown helium field: "<<name);
      if (!scalars.emplace(name,v).second) NCRYSTAL_THROW(BadInput,"Duplicate parameter");
    }
  }
  for (const auto& n:curveNames)
    if (!curves.count(n)) NCRYSTAL_THROW2(BadInput,"Missing curve "<<n);
  //The density of states comes the standard NCrystal way, @DYNINFO type vdos.
  const NC::DI_VDOS* vdosInput=nullptr;
  for (const auto& di:info.getDynamicInfoList())
    if ((vdosInput=dynamic_cast<const NC::DI_VDOS*>(di.get()))) break;
  if (!vdosInput)
    NCRYSTAL_THROW(BadInput,"SuperfluidHe needs the density of states in @DYNINFO (type vdos)");
  auto m=std::unique_ptr<Impl>(new Impl);
  m->temperature=info.getTemperature().dbl();m->kT=info.getTemperature().kT()*1000;
  //The temperature range of the shipped inputs, where the model was checked (also catches a
  //missing @TEMPERATURE, for which NCrystal assumes 293.15 K).
  if (!(m->temperature>=0.1 && m->temperature<=1.3))
    NCRYSTAL_THROW2(BadInput,"SuperfluidHe supports temperatures 0.1-1.3 K (material: "
                    <<m->temperature<<" K). A file without @TEMPERATURE takes its temperature "
                    "from the cfg-string, e.g. \";temp=0.6K\".");
  m->widthFromInput=scalars.count("fwhm_meV")!=0;
  m->gamma=.5*(m->widthFromInput ? scalars.at("fwhm_meV") : abeMorishimaFWHM(m->temperature));
  if (!(m->gamma>0)) NCRYSTAL_THROW(BadInput,"fwhm_meV must be positive");
  m->z=std::move(curves["z_q"]);
  m->h=std::move(curves["h_q"]);m->phonons=std::move(curves["phonons_meV"]);
  m->setup(atom,*vdosInput);
  return PhysicsModel(std::move(m));
}

// ============================================================================ 11. Public interface
// What NCrystal calls: the cross section and one scattering event per call, from the energy
// table between 1e-10 eV and emax (1 eV) and from method A below it. Above emax the
// factory (NCPluginFactory.cc) hands over to NCrystal's free gas.
PhysicsModel::PhysicsModel(std::unique_ptr<Impl> p):m_impl(std::move(p)) {}
PhysicsModel::PhysicsModel(PhysicsModel&&) noexcept=default;
PhysicsModel& PhysicsModel::operator=(PhysicsModel&&) noexcept=default;
PhysicsModel::~PhysicsModel()=default;
double PhysicsModel::maxEnergyEV() const noexcept {return .001*m_impl->emax;}
double PhysicsModel::calcCrossSection(double e,EnergyCache& c) const {
  const auto& m=*m_impl;
  //Below the first table energy the Q density still changes shape: exact path.
  if (!m.table.empty() && e>=m.table.front().e) {std::size_t i;double w;return m.tableSigma(e,i,w);}
  m.fillCache(e,c,m.qTolerance);return c.sigma;
}
// Q, then the energy transfer at that Q, then the outgoing energy and angle. Table: Q is
// proposed from an equal mixture of the Q densities of the two neighbouring table energies
// and accepted with probability density/(tableBound*proposal), where the density is the
// exact Q density at this energy (normalised by sigma). Where the density exceeds
// tableBound*proposal such Q are under-sampled: scans at 0.1-1.3 K found at most 3.22, and
// an excess is reported once per material. Below the table: Q by inverting the integrated Q
// density (no rejection).
PhysicsModel::ScatEvent PhysicsModel::sampleScatteringEvent(NC::RNG& rng,double energyEV,EnergyCache& c) const {
  const auto& m=*m_impl;
  bool exact=m.table.empty() || energyEV<m.table.front().e;
  std::size_t i=0;double w=0,sigma;
  if (exact) {m.fillCache(energyEV,c,m.qTolerance);sigma=c.sigma;}
  else sigma=m.tableSigma(energyEV,i,w);
  if (!(sigma>0)) return {energyEV,1};
  double e=energyEV*1000,ki=std::sqrt(e/cn);
  const Impl::TableNode* lowNode=exact ? nullptr : &m.table[i];
  const Impl::TableNode* highNode=(exact || !(w>0)) ? lowNode : &m.table[i+1];
  double mix=(lowNode==highNode) ? 0 : .5,norm=sigma*2*ki*ki/m.sigmaBound;
  for (unsigned retry=0;retry<100000;++retry) {
    double q,proposal=0;
    if (exact) q=m.sampleQ(rng,ki,c);
    else {
      q=Impl::sampleNodeQ(rng,randOpen(rng)<mix ? *highNode : *lowNode);
      proposal=(1-mix)*Impl::nodeDensity(*lowNode,q)+mix*Impl::nodeDensity(*highNode,q);
    }
    auto p=m.point(q);double lower=-cn*q*(q+2*ki),upper=cn*q*(2*ki-q);
    auto parts=m.components(p,lower,upper,true);
    double single=p.z*(parts.positive+parts.negative),multi=p.h*parts.multi;
    if (!(single+multi>0)) continue;//possible only at a floating-point boundary
    double density=q*(single+multi);
    if (!exact && density>Impl::tableBound*proposal*norm && !m.boundWarned.exchange(true))
      NCPLUGIN_WARN(m.temperature<<" K: Q density "<<density/(proposal*norm)<<" x the sampling "
                    "proposal (bound "<<Impl::tableBound<<") at E="<<energyEV<<" eV, Q="<<q
                    <<" 1/Aa; such Q are under-sampled. Reported once per material.");
    if (!exact && randOpen(rng)*Impl::tableBound*proposal*norm>density) continue;
    double ef=e-m.sampleTransfer(rng,p,parts,lower,upper);if (!(ef>0)) continue;
    double kf=std::sqrt(ef/cn),mu=(ki*ki+kf*kf-q*q)/(2*ki*kf);
    if (!std::isfinite(mu) || mu < -1-1e-10 || mu > 1+1e-10)
      NCRYSTAL_THROW(CalcError,"Sampled energy and Q violate neutron kinematics");
    return {ef*.001,std::max(-1.0,std::min(1.0,mu))};
  }
  NCRYSTAL_THROW(CalcError,"No accessible interval for sampled Q");
}

// ============================================================================ 12. Summary
// The summary NCrystal asks every process for (scatter.getSummary()["specific"]).
std::string PhysicsModel::descriptionJSON() const {
  const auto& m=*m_impl;
  std::ostringstream s;
  NC::streamJSONDictEntry(s,"summarystr","SuperfluidHe model of Granada et al., NIM A 1053 "
    "(2023) 168284: finite-width single phonon + n>=2 multiphonon",NC::JSONDictPos::FIRST);
  NC::streamJSONDictEntry(s,"temperature_K",m.temperature);
  NC::streamJSONDictEntry(s,"fwhm_meV",2*m.gamma);
  NC::streamJSONDictEntry(s,"fwhm_source",m.widthFromInput ? "input" : "abe_morishima_2001");
  s<<'}';
  return s.str();
}
}
