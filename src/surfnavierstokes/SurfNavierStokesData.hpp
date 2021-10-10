//
// Created by Alexander Zhiliakov on 9/27/19.
//

#ifndef SURF_NAVIER_STOKES_DATA_HPP
#define SURF_NAVIER_STOKES_DATA_HPP

#include <unordered_map>
#include "surfactant/Surface.hpp"
#include "num/functions.hpp"

namespace DROPS {

    struct SurfNavierStokesData {
        bool exact;
        InstatVectorFunction u_T, f_T, w_T;
        InstatScalarFunction p, m_g; // m_g is "-g"
        std::string description;
    };

    std::unordered_map<std::string, InstatDiffFunction> pressureData = {
        { "0", { "p = 0", zeroInstatScalarFunction, zeroInstatVectorFunction } },
        { "PolynomialExact", { "p = x y^3 + z", [](Point3DCL const & x, double) { return x[0] * pow(x[1], 3.) + x[2]; }, [](Point3DCL const & x, double) { return Point3DCL(pow(x[1], 3.), 3. * x[0] * x[1] * x[1], 1.); } } }
    };

    SurfNavierStokesData surfNavierStokesDataFactory(Surface const & surface, std::string const & velName, std::string const & preName, ParamCL const & params) {
        std::string funcName = __func__;
        if (pressureData.find(preName) == pressureData.end()) throw std::invalid_argument(funcName + ": pressure data '" + preName + "' is not defined");
        auto convectionTermType = params.get<std::string>("SurfNavierStokes.ConvectionTermType");
        InstatVectorFunction u_T = zeroInstatVectorFunction;
        struct { InstatVectorFunction stokesTerm, convTerm; } f_T = { zeroInstatVectorFunction, zeroInstatVectorFunction };
        InstatScalarFunction p = zeroInstatScalarFunction, m_g = zeroInstatScalarFunction;
        SurfNavierStokesData data;
        data.description = "convection term type: " + convectionTermType + '\n';
        data.exact = false;
        // cf. mathematica/rhs.wls
        if (velName == "0") {
            if (surface.isStationary()) {
                data.exact = true;
                data.description += "u = 0";
            }
            else data.description += "IC: u = 0";
        }
        else if (velName == "KelvinHelmholtz") {
            if (!surface.rotationalInvariance()[2]) throw std::invalid_argument(funcName + ": surface must be rotational invariant around z-axis for '" + velName + "' test");
            data.description += "IC: u = Kelvin-Helmholtz, cf. https://arxiv.org/abs/1909.02990";
            auto delta_0 = params.get<double>("SurfNavierStokes.IC.Velocity.Params." + velName + ".delta_0");
            auto cn = params.get<double>("SurfNavierStokes.IC.Velocity.Params." + velName + ".cn");
            auto aa = params.get<double>("SurfNavierStokes.IC.Velocity.Params." + velName + ".aa");
            auto ab = params.get<double>("SurfNavierStokes.IC.Velocity.Params." + velName + ".ab");
            auto ma = params.get<double>("SurfNavierStokes.IC.Velocity.Params." + velName + ".ma");
            auto mb = params.get<double>("SurfNavierStokes.IC.Velocity.Params." + velName + ".mb");
            auto eta = [](Point3DCL const & x) { return -0.3183098861837907 * asin(x[2] / std::sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2))); };
            auto eXi = [](Point3DCL const & x) {
                if (pow(x[0], 2) + pow(x[1], 2) == 0.) return Point3DCL(0., 0., 0.);
                return Point3DCL(-1. * x[1] * std::sqrt(1 / (pow(x[0], 2) + pow(x[1], 2))), x[0] * std::sqrt(1 / (pow(x[0], 2) + pow(x[1], 2))), 0.);
            };
            auto Hs = [=](double eta) { return tanh(2. * eta / delta_0); };
            auto arctan = [](double x, double y) {
                if (x > 0 && y >= 0) return atan(y/x);
                if (x > 0 && y < 0)  return atan(y/x) + 2. * M_PI;
                if (x < 0)           return atan(y/x) + M_PI;
                return sign(y) * M_PI / 2.;
            };
            auto pert = [=](Point3DCL const & x) {
                return Point3DCL(
                    (0.05066059182116889*(-4. * x[1] * std::sqrt((pow(x[0], 2) + pow(x[1], 2)) / (pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2))) * std::sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2)) * asin(x[2] / std::sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2))) * (aa * cos(0.5 * ma * arctan(x[0], x[1])) + ab * cos(0.5 * mb * arctan(x[0], x[1]))) + 9.869604401089358 * x[0] * x[2] * pow(delta_0, 2) * (aa * ma * sin(0.5 * ma * arctan(x[0], x[1])) + ab * mb * sin(0.5 * mb * arctan(x[0], x[1]))))) / (pow(2.718281828459045, (0.10132118364233778 * pow(asin(x[2] / std::sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2))), 2)) / pow(delta_0, 2)) * (pow(x[0], 2) + pow(x[1], 2)) * std::sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2)) * pow(delta_0, 2)),
                    (0.05066059182116889*(4. * x[0] * std::sqrt((pow(x[0], 2) + pow(x[1], 2)) / (pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2))) * std::sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2)) * asin(x[2] / std::sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2))) * (aa * cos(0.5 * ma * arctan(x[0], x[1])) + ab * cos(0.5 * mb * arctan(x[0], x[1]))) + 9.869604401089358 * x[1] * x[2] * pow(delta_0, 2) * (aa * ma * sin(0.5 * ma * arctan(x[0], x[1])) + ab * mb * sin(0.5 * mb * arctan(x[0], x[1]))))) / (pow(2.718281828459045, (0.10132118364233778 * pow(asin(x[2] / std::sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2))), 2)) / pow(delta_0, 2)) * (pow(x[0], 2) + pow(x[1], 2)) * std::sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2)) * pow(delta_0, 2)),
                    (-0.5*(aa*ma*sin(0.5*ma*arctan(x[0], x[1])) + ab * mb * sin(0.5 * mb * arctan(x[0], x[1])))) / (pow(2.718281828459045, (0.10132118364233778 * pow(asin(x[2] / std::sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2))), 2)) / pow(delta_0, 2)) * std::sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2)))
                );
            };
            auto r = [&surface](double x, double y) { return (std::sqrt(x * x + y * y) - surface.radius()[0].min) / (surface.radius()[0].max - surface.radius()[0].min); };
            u_T = [=](Point3DCL const & x, double) { return r(x[0], x[1]) * (Hs(eta(x)) * eXi(x) + cn * pert(x)); };
        }
        else if (velName == "PolynomialExact") {
            auto sphere = dynamic_cast<Sphere const *>(&surface);
            if (!sphere) throw std::invalid_argument(funcName + ": test '" + velName + "' is defined for sphere only");
            data.exact = true;
            data.description += "u = P (-z^2, y, x)";
            u_T = [=](Point3DCL const & x, double) { return sphere->P(x) * Point3DCL(-x[2] * x[2], x[1], x[0]); };
            m_g = [=](Point3DCL const & x, double t) { return -(x[0]*x[2]*(4.*x[2] - 3.) - 3.*pow(x[1],2) + pow(sphere->r(t),2) + 2.*sphere->r(t)*sphere->r_prime(t))/pow(sphere->r(t),2); };
            auto nu = params.get<double>("SurfNavierStokes.nu");
            f_T.stokesTerm = [=](Point3DCL const & x, double t) {
                return Point3DCL(
                    (0.5*(-2.*(5.*x[0]*pow(x[1],2) + x[2]*(-5. + 11.*x[2])*(pow(x[1],2) + pow(x[2],2)))*nu + (2.*pow(x[1],2) - 5.*x[2] + 14.*pow(x[2],2))*nu*pow(sphere->r(t),2) + 2.*sphere->r(t)*(x[0]*pow(x[1],2) + (-1. + x[2])*x[2]*(pow(x[1],2) + pow(x[2],2)) + x[2]*pow(sphere->r(t),2))*sphere->r_prime(t)))/pow(sphere->r(t),4),
                    (x[1]*nu*(-5.*pow(x[1],2) + x[0]*x[2]*(-5. + 11.*x[2]) + (5. - x[0])*pow(sphere->r(t),2)) + x[1]*sphere->r(t)*(pow(x[1],2) + x[0]*x[2] - x[0]*pow(x[2],2) + pow(sphere->r(t),2))*sphere->r_prime(t))/pow(sphere->r(t),4),
                    (0.5*(2.*x[2]*(-5.*pow(x[1],2) + x[0]*x[2]*(-5. + 11.*x[2]))*nu + x[0]*(5. - 14.*x[2])*nu*pow(sphere->r(t),2) + 2.*sphere->r(t)*(pow(x[1],2)*x[2] + x[0]*pow(x[2],2) - x[0]*pow(x[2],3) + x[0]*pow(sphere->r(t),2))*sphere->r_prime(t)))/pow(sphere->r(t),4)
                );
            };
            f_T.convTerm = [=](Point3DCL const & x, double t) {
                return Point3DCL(
                    (pow(x[1],2)*x[2]*(-4. + 5.*x[2])*(pow(x[1],2) + pow(x[2],2)) + x[0]*(2.*pow(x[1],4) - 2.*pow(x[2],4) + 5.*pow(x[2],5) - 3.*pow(x[2],6) + pow(x[1],2)*pow(x[2],2)*(-2. + 5.*x[2] - 3.*pow(x[2],2))) - (x[0]*pow(x[1],2) + 2.*(-2. + x[0])*pow(x[1],2)*x[2] - 2.*(x[0] - pow(x[1],2))*pow(x[2],2) + 3.*x[0]*pow(x[2],3))*pow(sphere->r(t),2) - sphere->r(t)*(x[0]*pow(x[1],2) + x[2]*(-1. + 4.*x[2])*(pow(x[1],2) + pow(x[2],2)) + x[2]*pow(sphere->r(t),2))*sphere->r_prime(t))/pow(sphere->r(t),4),
                    (x[1]*(2.*pow(x[1],4) - 2.*pow(x[2],4) + 5.*pow(x[2],5) - 3.*pow(x[2],6) + pow(x[1],2)*x[2]*(4.*x[0] - 2.*x[2] - 5.*x[0]*x[2] + 5.*pow(x[2],2) - 3.*pow(x[2],3)) + pow(sphere->r(t),2)*(pow(x[1],2)*(-3. - 2.*x[2]) - 2.*x[0]*x[2] + pow(x[2],2)*(2. + 2.*x[0] - 6.*x[2] + 3.*pow(x[2],2)) + (1. + 2.*x[2])*pow(sphere->r(t),2)) + sphere->r(t)*(-pow(x[1],2) - x[0]*x[2] + 4.*x[0]*pow(x[2],2) + pow(sphere->r(t),2))*sphere->r_prime(t)))/pow(sphere->r(t),4),
                    (2.*pow(x[1],4)*x[2] - 2.*pow(x[2],5) + 5.*pow(x[2],6) - 3.*pow(x[2],7) + pow(x[1],2)*pow(x[2],2)*(4.*x[0] - 2.*x[2] - 5.*x[0]*x[2] + 5.*pow(x[2],2) - 3.*pow(x[2],3)) + pow(sphere->r(t),2)*(-2.*x[0]*pow(x[1],2) + 4.*pow(x[2],3) - 8.*pow(x[2],4) + 3.*pow(x[2],5) + pow(x[1],2)*(x[2] - 4.*pow(x[2],2)) + x[2]*(-2. + 3.*x[2])*pow(sphere->r(t),2)) + sphere->r(t)*(x[2]*(-pow(x[1],2) - x[0]*x[2] + 4.*x[0]*pow(x[2],2)) + x[0]*pow(sphere->r(t),2))*sphere->r_prime(t))/pow(sphere->r(t),4)
                );
            };
        }
        else if (velName == "DirectionChangeExact") {
            auto sphere = dynamic_cast<Sphere const *>(&surface);
            if (!sphere) throw std::invalid_argument(funcName + ": test '" + velName + "' is defined for sphere only");
            data.exact = true;
            data.description += "u = P (1 - 2 t, 0, 0)";
            u_T = [=](Point3DCL const & x, double t) { return sphere->P(x, t) * Point3DCL(1. - 2. * t, 0., 0.); };
            m_g = [=](Point3DCL const & x, double t) { return (-2.*((-1. + 2.*t)*x[0] + sphere->r(t)*sphere->r_prime(t)))/pow(sphere->r(t),2); };
            auto nu = params.get<double>("SurfNavierStokes.nu");
            f_T.stokesTerm = [=](Point3DCL const & x, double t) {
                return Point3DCL(
                    (-(pow(x[1],2) + pow(x[2],2))*((-1. + 2.*t)*nu + 2.*pow(sphere->r(t),2) + (-1. + 2.*t)*sphere->r(t)*sphere->r_prime(t)))/pow(sphere->r(t),4),
                    (x[0]*x[1]*((-1. + 2.*t)*nu + 2.*pow(sphere->r(t),2) + (-1. + 2.*t)*sphere->r(t)*sphere->r_prime(t)))/pow(sphere->r(t),4),
                    (x[0]*x[2]*((-1. + 2.*t)*nu + 2.*pow(sphere->r(t),2) + (-1. + 2.*t)*sphere->r(t)*sphere->r_prime(t)))/pow(sphere->r(t),4)
                );
            };
            f_T.convTerm = [=](Point3DCL const & x, double t) {
                return Point3DCL(
                    (-pow(1. - 2.*t,2)*x[0]*(pow(x[1],2) + pow(x[2],2)))/pow(sphere->r(t),4),
                    (pow(1. - 2.*t,2)*pow(x[0],2)*x[1])/pow(sphere->r(t),4),
                    (pow(1. - 2.*t,2)*pow(x[0],2)*x[2])/pow(sphere->r(t),4)
                );
            };
        }
        else if (velName == "KillingExact") {
            auto sphere = dynamic_cast<Sphere const *>(&surface);
            if (!sphere) throw std::invalid_argument(funcName + ": test '" + velName + "' is defined for sphere only");
            data.exact = true;
            auto omega = params.get<double>("SurfNavierStokes.IC.Params." + velName + ".AngularVelocity");
            data.description += "u = omega (0, -z, y), omega = " + std::to_string(omega);
            u_T = [=](Point3DCL const & x, double) { return omega * Point3DCL(0., -x[2], x[1]); };
            m_g = [=](Point3DCL const &, double t) { return -2. * sphere->r_prime(t) / sphere->r(t); };
            f_T.stokesTerm = [=](Point3DCL const & x, double t) { return Point3DCL(0., -x[2] * omega * sphere->r_prime(t) / sphere->r(t), x[1] * omega * sphere->r_prime(t) / sphere->r(t)); };
            f_T.convTerm = [=](Point3DCL const & x, double t) {
                return Point3DCL(
                    (x[0]*(pow(x[1],2) + pow(x[2],2))*pow(omega,2))/pow(sphere->r(t),2),
                    (omega*(x[1]*omega*(pow(x[1],2) + pow(x[2],2) - 1.*pow(sphere->r(t),2)) - 1.*x[2]*sphere->r(t)*sphere->r_prime(t)))/pow(sphere->r(t),2),
                    (omega*(x[2]*omega*(pow(x[1],2) + pow(x[2],2) - 1.*pow(sphere->r(t),2)) + x[1]*sphere->r(t)*sphere->r_prime(t)))/pow(sphere->r(t),2)
                );
            };
        }
        else throw std::invalid_argument(funcName + ": velocity IC '" + velName + "' is not defined");
        data.description += ", " + pressureData[preName].description;
        data.u_T = [=, &surface](Point3DCL const & x, double t) { return u_T(surface.ext(x, t), t); };
        data.p = [=, &surface](Point3DCL const & x, double t) { return pressureData[preName].f(surface.ext(x, t), t); };
        data.m_g = [=, &surface](Point3DCL const & x, double t) { return m_g(surface.ext(x, t), t); };
        if (convectionTermType == "Stokes") {
            data.w_T = zeroInstatVectorFunction;
            data.f_T = [=, &surface](Point3DCL const & y, double t) {
                auto x = surface.ext(y, t);
                return f_T.stokesTerm(x, t) + surface.P(x, t) * pressureData[preName].grad(x, t);
            };
        }
        else if (convectionTermType == "Oseen") {
            if (!data.exact) throw std::invalid_argument(funcName + ": 'Oseen' convection term type is unavailable since exact solution is not provided; use 'Stokes' or 'NavierStokes'");
            data.w_T = data.u_T;
            data.f_T = [=, &surface](Point3DCL const & y, double t) {
                auto x = surface.ext(y, t);
                return f_T.stokesTerm(x, t) + f_T.convTerm(x, t) + surface.P(x, t) * pressureData[preName].grad(x, t);
            };
        }
        else if (convectionTermType == "NavierStokes") {
            data.w_T = nullptr;
            data.f_T = [=, &surface](Point3DCL const & y, double t) {
                auto x = surface.ext(y, t);
                return f_T.stokesTerm(x, t) + f_T.convTerm(x, t) + surface.P(x, t) * pressureData[preName].grad(x, t);
            };
        }
        else throw std::invalid_argument(funcName + ": convection term type '" + convectionTermType + "' is not defined");
        return data;
    }

}

#endif // SURF_NAVIER_STOKES_DATA_HPP