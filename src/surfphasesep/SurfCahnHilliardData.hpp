//
// Created by Alexander Zhiliakov on 2/18/20.
//

#ifndef SURF_CAHN_HILIARD_DATA_HPP
#define SURF_CAHN_HILIARD_DATA_HPP

#include <random>
#include "surfactant/Surface.hpp"
#include "num/functions.hpp"

namespace DROPS {

    struct SurfCahnHilliardData {
        bool exact;
        double raftRatio;
        InstatScalarFunction chi, f;
        std::string description;
    };

    SurfCahnHilliardData surfCahnHilliardDataFactory(Surface const & surface, std::string const & name, ParamCL const & params) {
        std::string funcName = __func__;
        InstatScalarFunction chi = zeroInstatScalarFunction, f = zeroInstatScalarFunction;
        SurfCahnHilliardData data;
        data.exact = false;
        data.raftRatio = -1;
        // cf. mathematica/rhs.wls
        if (name == "WanDerVaals") {
            data.description = "chi_0 = .5 * (1 + tanh(z / (2 sqrt(2) eps)))";
            auto eps = params.get<double>("SurfCahnHilliard.Epsilon");
            auto omega = params.get<double>("SurfCahnHilliard.IC.Params." + name + ".AngularVelocity");
            chi = [=](Point3DCL const & x, double t) { return .5 * (1. + tanh((x[2] * cos(omega * t) - x[1] * sin(omega * t))/(2. * std::sqrt(2.) * eps))); };
            auto sphere = dynamic_cast<Sphere const *>(&surface);
            if (sphere && !params.get<bool>("SurfCahnHilliard.UseDegenerateMobility")) {
                data.exact = true;
                data.description =
                    "chi = .5 * (1 + tanh((z cos(omega t) - y sin(omega t))/ (2 sqrt(2) eps)))\n"
                    "rotation of Wan der Vaals 'tanh' solution around x-axis w/ velocity u = omega (0, -z, y)";
                auto sech = [](double x) { return 1. / cosh(x); };
                auto M = params.get<double>("SurfCahnHilliard.MobilityScaling");
                f = [=](Point3DCL const & x, double t) {
                    auto arg = x[2] * cos(omega * t) - x[1] * sin(omega * t);
                    return 0.015625 * M * (pow(sech((0.35355339059327373*arg)/eps),5)*pow(arg,3)*(11.313708498984761*eps*cosh((1.0606601717798212*arg)/eps) - arg*sinh((1.0606601717798212*arg)/eps)) + 16.*pow(eps,2)*pow(sech((0.35355339059327373*arg)/eps),2)*(2.8284271247461903*eps*arg - (-3.*pow(sphere->r(t),2) + 7.*pow(arg,2))*tanh((0.35355339059327373*arg)/eps)) + pow(sech((0.35355339059327373*arg)/eps),4)*arg*(-33.941125496954285*pow(x[2],2)*eps*pow(cos(t*omega),2) + 16.970562748477143*eps*((2. - cosh((0.7071067811865475*arg)/eps))*pow(sphere->r(t),2) - 2.*x[1]*(x[1]*pow(sin(t*omega),2) - x[2]*sin(2.*t*omega))) + arg*(2.*(-5. + cosh((0.7071067811865475*arg)/eps))*pow(sphere->r(t),2) + 11.*pow(arg,2))*tanh((0.35355339059327373*arg)/eps)))/(pow(eps,2)*pow(sphere->r(t),4));
                };
            }
        }
        else if (name == "RandomUniform") {
            data.raftRatio = params.get<double>("SurfCahnHilliard.IC.Params." + name + ".RaftRatio");
            auto raftRatioNoisePercent = params.get<double>("SurfCahnHilliard.IC." + name + ".RaftRatioNoiseFraction");
            auto a = data.raftRatio - raftRatioNoisePercent * data.raftRatio;
            auto b = data.raftRatio + raftRatioNoisePercent * data.raftRatio;
            data.description = "chi_0 ~ Uniform(" + std::to_string(a) + ", " + std::to_string(b) + ")";
            chi = [=](Point3DCL const &, double) mutable {
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_real_distribution<> dis(a, b);
                return dis(gen);
            };
        }
        else if (name == "RandomBernoulli") {
            data.raftRatio = params.get<double>("SurfCahnHilliard.IC.Params." + name + ".RaftRatio");
            data.description = "chi_0 ~ Bernoulli(" + std::to_string(data.raftRatio) + ")";
            chi = [=](Point3DCL const &, double) mutable {
                std::random_device rd;
                std::mt19937 gen(rd());
                std::bernoulli_distribution dis(data.raftRatio);
                return dis(gen);
            };
        }
        else throw std::invalid_argument(funcName + ": IC '" + name + "' is not defined");
        data.chi = [=, &surface](Point3DCL const & x, double t) { return chi(surface.ext(x, t), t); };
        data.f = [=, &surface](Point3DCL const & x, double t) { return f(surface.ext(x, t), t); };
        return data;
    }

}

#endif // SURF_CAHN_HILIARD_DATA_HPP