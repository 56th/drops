//
// Created by Alexander Zhiliakov on 2/18/20.
//

#ifndef SURF_CAHN_HILIARD_DATA_HPP
#define SURF_CAHN_HILIARD_DATA_HPP

#include <random>
#include <string>
#include "../num/discretize.h"

namespace DROPS {

    struct SurfCahnHilliardData {
        bool exactSoln;
        InstatScalarFunction chi, omega, rhs3, rhs4;
        struct Surface {
            InstatScalarFunction phi;
            InstatVectorFunction n;
            InstatVectorFunction e;
            InstatMatrixFunction H;
        };
        Surface surface;
        std::string description;
    };

    SurfCahnHilliardData SurfCahnHilliardDataFactory(ParamCL const & param) {
        auto test = param.get<std::string>("SurfCahnHilliard.TestName");
        SurfCahnHilliardData data;
        SurfCahnHilliardData::Surface sphere;
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
            data.exactSoln = false;
            auto raftRatio = param.get<double>("SurfCahnHilliard.IC.Random.RaftRatio");
            auto raftRatioNoisePercent = param.get<double>("SurfCahnHilliard.IC.Random.RaftRatioNoiseFraction");
            auto a = raftRatio - raftRatioNoisePercent * raftRatio;
            auto b = raftRatio + raftRatioNoisePercent * raftRatio;
            data.description =
                    "phi = sqrt(x^2 + y^2 + z^2) - 1, chi ~ Uniform(" + std::to_string(a) + ", " + std::to_string(b) + ")\n"
                    "no force terms\n";
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(a, b);
            data.chi = [=](Point3DCL const &, double) mutable {
                return dis(gen);
            };
            data.omega = data.rhs3 = data.rhs4 = [](Point3DCL const &, double) {
                return 0.;
            };
            data.surface = sphere;
        }
        else if (test == "RandomBernoulli") {
            data.exactSoln = false;
            auto raftRatio = param.get<double>("SurfCahnHilliard.IC.Random.RaftRatio");
            data.description =
                    "phi = sqrt(x^2 + y^2 + z^2) - 1, chi ~ Bernoulli(" + std::to_string(raftRatio) + ")\n"
                    "no force terms\n";
            std::random_device rd;
            std::mt19937 gen(rd());
            std::bernoulli_distribution dis(raftRatio);
            data.chi = [=](Point3DCL const &, double) mutable {
                return dis(gen);
            };
            data.omega = data.rhs3 = data.rhs4 = [](Point3DCL const &, double) {
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
            data.omega = data.rhs3 = data.rhs4 = [](Point3DCL const &, double) {
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