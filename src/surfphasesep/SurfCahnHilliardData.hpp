//
// Created by Alexander Zhiliakov on 2/18/20.
//

#ifndef SURF_CAHN_HILIARD_DATA_HPP
#define SURF_CAHN_HILIARD_DATA_HPP

#include <random>
#include <string>
#include "num/functions.hpp"

namespace DROPS {

    struct SurfCahnHilliardData {
        bool exactSoln;
        double raftRatio = -1.;
        InstatScalarFunction chi, omega, f;
        Surface surface;
        std::string description;
    };

    SurfCahnHilliardData SurfCahnHilliardDataFactory(ParamCL const & param) {
        auto test = param.get<std::string>("SurfCahnHilliard.TestName");
        SurfCahnHilliardData data;
        Surface sphere;
        sphere.phi = [](Point3DCL const & p, double) {
            return std::sqrt(std::pow(p[0], 2.) + std::pow(p[1], 2.) + std::pow(p[2], 2.)) - 1.;
        };
        sphere.n = [](Point3DCL const & p, double) {
            auto den = std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
            Point3DCL v(p[0]/den, p[1]/den, p[2]/den);
            return v;
        };
        sphere.e = sphere.n;
        sphere.H = [](Point3DCL const & p, double) {
            SMatrixCL<3, 3> res;
            auto den = std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2);
            res(0, 0) = (std::pow(p[1],2) + std::pow(p[2],2)) / den;
            res(0, 1) = (-1.*p[0]*p[1]) / den;
            res(0, 2) = (-1.*p[0]*p[2]) / den;
            res(1, 0) = res(0, 1);
            res(1, 1) = (std::pow(p[0],2) + std::pow(p[2],2)) / den;
            res(1, 2) = (-1.*p[1]*p[2]) / den;
            res(2, 0) = res(0, 2);
            res(2, 1) = res(1, 2);
            res(2, 2) = (std::pow(p[0],2) + std::pow(p[1],2)) / den;
            return res;
        };
        if (test == "RandomUniform") {
            data.raftRatio = param.get<double>("SurfCahnHilliard.IC.Random.RaftRatio");
            data.exactSoln = false;
            auto raftRatioNoisePercent = param.get<double>("SurfCahnHilliard.IC.Random.RaftRatioNoiseFraction");
            auto a = data.raftRatio - raftRatioNoisePercent * data.raftRatio;
            auto b = data.raftRatio + raftRatioNoisePercent * data.raftRatio;
            data.description =
                    "phi = sqrt(x^2 + y^2 + z^2) - 1, chi ~ Uniform(" + std::to_string(a) + ", " + std::to_string(b) + ")\n"
                    "no force terms\n";
            data.chi = [=](Point3DCL const &, double) mutable {
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_real_distribution<> dis(a, b);
                return dis(gen);
            };
            data.omega = data.f = [](Point3DCL const &, double) {
                return 0.;
            };
            data.surface = sphere;
        }
        else if (test == "RandomBernoulli") {
            data.raftRatio = param.get<double>("SurfCahnHilliard.IC.Random.RaftRatio");
            data.exactSoln = false;
            data.description =
                    "phi = sqrt(x^2 + y^2 + z^2) - 1, chi ~ Bernoulli(" + std::to_string(data.raftRatio) + ")\n"
                    "no force terms\n";
            data.chi = [=](Point3DCL const &, double) mutable {
                std::random_device rd;
                std::mt19937 gen(rd());
                std::bernoulli_distribution dis(data.raftRatio);
                return dis(gen);
            };
            data.omega = data.f = [](Point3DCL const &, double) {
                return 0.;
            };
            data.surface = sphere;
        }
        else if (test == "OneZeroSphere") {
            data.exactSoln = false;
            data.description =
                    "phi = sqrt(x^2 + y^2 + z^2) - 1, chi = .5 * (1 + tanh(y / (2 sqrt(2) eps)))\n"
                    "Note: rhs is set to zero\n";
            auto eps = param.get<double>("SurfCahnHilliard.Epsilon");
            data.chi = [=](Point3DCL const & p, double t) mutable {
                auto x = sphere.e(p, t);
                return .5 * (1. + std::tanh(x[1] / (2. * std::sqrt(2.) * eps)));
            };
            data.omega = data.f = [](Point3DCL const &, double) {
                return 0.;
            };
            data.surface = sphere;
        }
        else throw std::invalid_argument("test '" + test + "' is not defined");
        data.description += data.exactSoln ? "exact solution is available and the errors will be computed\n" : "exact solution is NOT available and the errors will NOT be computed\n";
        return data;
    }

}

#endif // SURF_CAHN_HILIARD_DATA_HPP