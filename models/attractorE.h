#ifndef ATTRACTORE_H // Usual macro guard to prevent multiple inclusion
#define ATTRACTORE_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// Generalized and deformed alpha-attractor E-model (alpha-Starobinsky), following
// Ellis, Garcia, Olive, Verner, arXiv:2510.18656, Sec. V (generalized and deformed attractors):
//
//   V(phi) = (3/4) lambda MPl^4 [ kappa (1 - cosh(b x)) + sinh(b x) ]^k ,   x = phi/MPl, b = sqrt(2/(3 alphaAtt))
//
// kappa = 1 gives the generalized E-model V = (3/4) lambda MPl^4 (1 - e^{-b x})^k,
// and additionally k = 2, alphaAtt = 1 gives the Starobinsky model.
// Near the minimum V ~ phi^k for any kappa. See notes/inflation_models.tex.
//
// A daughter field chi (scalar 1) is coupled through (1/2) g^2 phi^2 chi^2.

#include "CosmoInterface/cosmointerface.h"

namespace TempLat
{
  struct ModelPars : public TempLat::DefaultModelPars {
    static constexpr size_t NScalars = 2;
    static constexpr size_t NPotTerms = 2;

    using FloatType = double;
  };

#define MODELNAME attractorE

  template <class R> using Model = MakeModel(R, ModelPars);

  class MODELNAME : public Model<MODELNAME>
  {
  private:
    FloatType lambda, alphaAtt, k, kappa, g;
    // b = sqrt(2/(3 alphaAtt)); pref = (3/4) lambda MPl^2/omegaStar^2 converts V/V0 into program units;
    // qInt = g^2 MPl^2/omegaStar^2 is the program-units interaction coupling.
    FloatType b, pref, qInt;

  public:
    MODELNAME(ParameterParser &parser, RunParameters<FloatType> &runPar, auto toolBox)
        : Model<MODELNAME>(parser, runPar.getLatParams(), toolBox, runPar.dt, STRINGIFY(MODELLABEL))
    {
      /////////
      // Independent parameters of the model (read from parameters file)
      /////////

      lambda = parser.get<FloatType>("lambda");          // CMB normalization, ~ 24 alphaAtt pi^2 A_s / N_*^2
      alphaAtt = parser.get<FloatType>("alphaAtt", 1.);   // attractor parameter (Kahler curvature 2/(3 alphaAtt))
      k = parser.get<FloatType>("k", 2.);                 // V ~ phi^k near the minimum, even integer
      kappa = parser.get<FloatType>("kappa", 1.);         // plateau deformation, 1 = undeformed
      g = parser.get<FloatType>("g", 0.);                 // inflaton-daughter coupling (1/2) g^2 phi^2 chi^2

      if (k < 2 || !AlmostEqual(k / 2, std::round(k / 2), runParameterTolerance))
        throw(RunParametersInconsistent("attractorE: k must be an even integer >= 2."));
      if (alphaAtt <= 0) throw(RunParametersInconsistent("attractorE: alphaAtt must be positive."));

      fldS0 = parser.get<FloatType, 2>("initial_amplitudes"); // in GeV
      piS0 = parser.get<FloatType, 2>("initial_momenta", {0, 0}); // in GeV^2

      if (AlmostEqual(fldS0[0], 0.0, runParameterTolerance))
        throw(RunParametersInconsistent("attractorE: the initial inflaton amplitude must be non-zero."));

      b = sqrt(FloatType(2) / (3 * alphaAtt));

      /////////
      // Rescaling for program variables
      /////////

      // Near the minimum V ~ (lambda_k/k) MPl^(4-k) phi^k, with lambda_k = k (3/4) lambda b^k.
      // omegaStar is the oscillation frequency of a condensate of amplitude phi_i, and
      // alpha = 3(k-2)/(k+2) keeps that frequency constant in program time
      // (alpha = 0 for k = 2, 1 for k = 4).
      const FloatType xi = std::abs(fldS0[0]) / MPl;
      const FloatType lambdak = k * FloatType(0.75) * lambda * pow(b, k);

      alpha = 3 * (k - 2) / (k + 2);
      fStar = MPl;
      omegaStar = sqrt(lambdak) * MPl * pow(xi, (k - 2) / 2);

      pref = FloatType(0.75) * lambda * pow<2>(MPl / omegaStar);
      qInt = pow<2>(g * MPl / omegaStar);

      setInitialPotentialAndMassesFromPotential();
    }

    /////////
    // Program potential. The program field is x = phi/MPl.
    /////////

    // G = kappa (1 - cosh(b x)) + sinh(b x), written to avoid the cancellation in 1 - cosh near x = 0.
    auto G() const { return sinh(b * fldS(0_c)) - 2 * kappa * pow<2>(sinh(b * fldS(0_c) / 2)); }
    auto dG() const { return b * (cosh(b * fldS(0_c)) - kappa * sinh(b * fldS(0_c))); }
    auto d2G() const { return pow<2>(b) * (sinh(b * fldS(0_c)) - kappa * cosh(b * fldS(0_c))); }

    auto potentialTerms(Tag<0>) const // Inflaton potential energy
    {
      return pref * pow(G(), k);
    }
    auto potentialTerms(Tag<1>) const // Interaction energy
    {
      return FloatType(0.5) * qInt * pow<2>(fldS(0_c) * fldS(1_c));
    }

    /////////
    // Derivatives of the program potential with respect to the fields
    /////////

    auto potDeriv(Tag<0>) // Derivative with respect to the inflaton
    {
      return pref * k * pow(G(), k - 1) * dG() + qInt * fldS(0_c) * pow<2>(fldS(1_c));
    }
    auto potDeriv(Tag<1>) // Derivative with respect to the daughter field
    {
      return qInt * pow<2>(fldS(0_c)) * fldS(1_c);
    }

    /////////
    // Second derivatives of the program potential with respect to the fields
    /////////

    auto potDeriv2(Tag<0>) // Second derivative with respect to the inflaton
    {
      return pref * k * ((k - 1) * pow(G(), k - 2) * pow<2>(dG()) + pow(G(), k - 1) * d2G()) +
             qInt * pow<2>(fldS(1_c));
    }
    auto potDeriv2(Tag<1>) // Second derivative with respect to the daughter field
    {
      return qInt * pow<2>(fldS(0_c));
    }
  };
} // namespace TempLat

#endif // ATTRACTORE_H
