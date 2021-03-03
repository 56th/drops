#ifndef SURF_NAVIER_STOKES_CAHN_HILLIARD_DATA_HPP
#define SURF_NAVIER_STOKES_CAHN_HILLIARD_DATA_HPP

#include "surfnavierstokes/SurfNavierStokesData.hpp"
#include "surfphasesep/SurfCahnHilliardData.hpp"

namespace DROPS {

    using SurfNSCHData = std::pair<SurfNavierStokesData, SurfCahnHilliardData>;

    SurfNSCHData surfNSCHDataFactory(Surface const & surface, std::string const & name, ParamCL& params) {
        std::string funcName = __func__;
        SurfNavierStokesData dataNS;
        SurfCahnHilliardData dataCH;
        if (name == "WanDerVaalsExact") {
            if (params.get<double>("SurfNSCH.NS.LineTension")) throw std::invalid_argument(funcName + ": use zero line tension for '" + name + "' test");
            if (params.get<double>("SurfNSCH.NS.rho.min") != params.get<double>("SurfNSCH.NS.rho.max") || params.get<double>("SurfNSCH.NS.rho.max") != 1.) throw std::invalid_argument(funcName + ": use rho_min = rho_max = 1 for '" + name + "' test");
            if (params.get<bool>("SurfNSCH.CH.UseDegenerateMobility")) throw std::invalid_argument(funcName + ": use constant mobility for '" + name + "' test");
            auto omega = get<double>("SurfNSCH.IC.Params." + name + ".AngularVelocity");
            params.put("SurfNavierStokes.IC.Params.KillingExact.AngularVelocity", omega);
            dataNS = surfNavierStokesDataFactory(surface, "KillingExact", params);
            params.put("SurfCahnHilliard.Epsilon", params.get<double>("SurfNSCH.CH.Epsilon"));
            params.put("SurfCahnHilliard.IC.Params.WanDerVaals.AngularVelocity", omega);
            dataCH = surfCahnHilliardDataFactory(surface, "WanDerVaals", params);
        }
        else if (name == "KelvinHelmholtz") {
            params.set_child("SurfNavierStokes.IC.Params." + name, params.get_child("SurfNSCH.IC.Params." + name));
            params.put("SurfNavierStokes.IC.Params." + name + ".delta_0", .5 * params.get<double>("SurfNSCH.CH.Epsilon"));
            dataNS = surfNavierStokesDataFactory(surface, name, params);
            params.put("SurfCahnHilliard.Epsilon", params.get<double>("SurfNSCH.CH.Epsilon"));
            params.put("SurfCahnHilliard.IC.Params.WanDerVaals.AngularVelocity", 0.);
            dataCH = surfCahnHilliardDataFactory(surface, "WanDerVaals", params);
            dataCH.exact = false;
            dataCH.f = zeroInstatScalarFunction;
        }
        else if (name == "RayleighTaylor") {
            dataNS = surfNavierStokesDataFactory(surface, "0", params);
            auto g = params.get<double>("SurfNSCH.IC.Params." + name + ".GravityScaling");
            dataNS.f_T = [=, &surface](Point3DCL const & x, double t) { return surface.P(surface.ext(x, t), t) * Point3DCL(0., 0., -g); };
            params.put("SurfCahnHilliard.Epsilon", params.get<double>("SurfNSCH.CH.Epsilon"));
            params.put("SurfCahnHilliard.IC.Params.WanDerVaals.AngularVelocity", 0.);
            dataCH = surfCahnHilliardDataFactory(surface, "WanDerVaals", params);
            // perturb chi ...
            dataCH.exact = false;
            dataCH.f = zeroInstatScalarFunction;
        }
        else throw std::invalid_argument(funcName + ": IC '" + name + "' is not defined");
        return { dataNS, dataCH };
    }

}

#endif // SURF_NAVIER_STOKES_CAHN_HILLIARD_DATA_HPP