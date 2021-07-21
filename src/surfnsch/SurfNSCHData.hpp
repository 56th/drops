#ifndef SURF_NAVIER_STOKES_CAHN_HILLIARD_DATA_HPP
#define SURF_NAVIER_STOKES_CAHN_HILLIARD_DATA_HPP

#include "surfnavierstokes/SurfNavierStokesData.hpp"
#include "surfphasesep/SurfCahnHilliardData.hpp"

namespace DROPS {

    using SurfNSCHData = std::pair<SurfNavierStokesData, SurfCahnHilliardData>;

    SurfNSCHData surfNSCHDataFactory(Surface const & surface, std::string const & name, ParamCL& params) {
        std::string funcName = __func__;
        if (!surface.isStationary()) throw std::invalid_argument(funcName + ": use stationary surface");
        SurfNavierStokesData dataNS;
        SurfCahnHilliardData dataCH;
        params.put_child("SurfNavierStokes", params.get_child("SurfNSCH.NS"));
        params.put_child("SurfCahnHilliard", params.get_child("SurfNSCH.CH"));
        if (name == "WanDerVaalsExact") {
            if (params.get<double>("SurfNSCH.NS.LineTension")) throw std::invalid_argument(funcName + ": use zero line tension for '" + name + "' test");
            if (params.get<double>("SurfNSCH.NS.rho.min") != params.get<double>("SurfNSCH.NS.rho.max") || params.get<double>("SurfNSCH.NS.rho.max") != 1.) throw std::invalid_argument(funcName + ": use rho_min = rho_max = 1 for '" + name + "' test");
            auto omega = params.get<double>("SurfNSCH.IC.Params." + name + ".AngularVelocity");
            params.put("SurfNavierStokes.IC.Params.KillingExact.AngularVelocity", omega);
            dataNS = surfNavierStokesDataFactory(surface, "KillingExact", params);
            params.put("SurfCahnHilliard.IC.Params.WanDerVaals.AngularVelocity", omega);
            params.put("SurfCahnHilliard.IC.Params.WanDerVaals.Noise", 0.);
            dataCH = surfCahnHilliardDataFactory(surface, "WanDerVaals", params);
        }
        else if (name == "KelvinHelmholtz") {
            params.put_child("SurfNavierStokes.IC.Params." + name, params.get_child("SurfNSCH.IC.Params." + name));
            params.put("SurfNavierStokes.IC.Params." + name + ".delta_0", params.get<double>("SurfNSCH.CH.Epsilon"));
            dataNS = surfNavierStokesDataFactory(surface, name, params);
            params.put("SurfCahnHilliard.IC.Params.WanDerVaals.AngularVelocity", 0.);
            params.put("SurfCahnHilliard.IC.Params.WanDerVaals.Noise", 0.);
            dataCH = surfCahnHilliardDataFactory(surface, "WanDerVaals", params);
            dataCH.exact = false;
            dataCH.f = zeroInstatScalarFunction;
        }
        else if (name == "RayleighTaylor") {
            dataNS = surfNavierStokesDataFactory(surface, "0", params);
            auto g = params.get<double>("SurfNSCH.IC.Params." + name + ".GravityScaling");
            dataNS.f_T = [=, &surface](Point3DCL const & x, double t) { return surface.P(surface.ext(x, t), t) * Point3DCL(0., 0., -g); };
            params.put("SurfCahnHilliard.IC.Params.WanDerVaals.AngularVelocity", 0.);
            params.put("SurfCahnHilliard.IC.Params.WanDerVaals.Noise", params.get<double>("SurfNSCH.IC.Params." + name + ".Noise"));
            dataCH = surfCahnHilliardDataFactory(surface, "WanDerVaals", params);
            dataCH.exact = false;
            dataCH.f = zeroInstatScalarFunction;
        }
        else if (name == "RandomBernoulli") {
            dataNS = surfNavierStokesDataFactory(surface, "0", params);
            auto raftRatio = params.get<double>("SurfNSCH.IC.Params.RandomBernoulli.RaftRatio");
            params.put("SurfCahnHilliard.IC.Params." + name + ".RaftRatio", raftRatio);
            dataCH = surfCahnHilliardDataFactory(surface, "RandomBernoulli", params);
            dataCH.exact = false;
            dataCH.f = zeroInstatScalarFunction;
        }
        else if (name == "RandomUniform") {
            dataNS = surfNavierStokesDataFactory(surface, "0", params);
            auto raftRatio = params.get<double>("SurfNSCH.IC.Params.RandomUniform.RaftRatio");
            auto raftRatioNoisePercent = params.get<double>("SurfNSCH.IC.Params.RandomUniform.RaftRatioNoiseFraction");
            params.put("SurfCahnHilliard.IC.Params." + name + ".RaftRatioNoiseFraction", raftRatioNoisePercent);
            params.put("SurfCahnHilliard.IC.Params." + name + ".RaftRatio", raftRatio);
            dataCH = surfCahnHilliardDataFactory(surface, "RandomUniform", params);
            dataCH.exact = false;
            dataCH.f = zeroInstatScalarFunction;
        }
        else if (name == "SyntheticRafts") {
            dataNS = surfNavierStokesDataFactory(surface, "0", params);
            dataCH = surfCahnHilliardDataFactory(surface, "SyntheticRafts", params);
            dataCH.exact = false;
            dataCH.f = zeroInstatScalarFunction;
        }
        else throw std::invalid_argument(funcName + ": IC '" + name + "' is not defined");
        return { dataNS, dataCH };
    }

}

#endif // SURF_NAVIER_STOKES_CAHN_HILLIARD_DATA_HPP