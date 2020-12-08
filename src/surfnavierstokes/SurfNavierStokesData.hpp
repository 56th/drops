//
// Created by Alexander Zhiliakov on 9/27/19.
//

#ifndef SURF_NAVIER_STOKES_DATA_HPP
#define SURF_NAVIER_STOKES_DATA_HPP

#include "../num/discretize.h"
// #include "surfnavierstokes_funcs.h"
#include "hermite_cubic/hermite_cubic.hpp"

namespace DROPS {

    struct SurfNavierStokesData {
        bool exactSoln;
        InstatVectorFunction u_T, f_T, w_T = nullptr; // w_T = nullptr for the Navier-Stokes case
        InstatScalarFunction p, m_g; // m_g is "-g"
        struct Surface {
            InstatScalarFunction phi, u_N;
            InstatVectorFunction n;
            InstatVectorFunction e;
            InstatMatrixFunction H;
        };
        Surface surface;
        std::string description;
    };

    double arctan(double x, double y) {
        if (x > 0 && y >= 0) return std::atan(y/x);
        if (x > 0 && y < 0)  return std::atan(y/x) + 2. * M_PI;
        if (x < 0)           return std::atan(y/x) + M_PI;
        return sign(y) * M_PI / 2.;
    }

    SurfNavierStokesData SurfNavierStokesDataFactory(std::string const & test, double nu, ParamCL const & param) {
        SurfNavierStokesData data;
        if (test.find("Sphere") != std::string::npos) {
            data.description = "phi = x^2 + y^2 + z^2 - 1\n";
            data.surface.u_N = [](Point3DCL const &, double) {
                return 0.;
            };
            data.surface.phi = [](Point3DCL const &p, double) {
                return std::pow(p[0], 2.) + std::pow(p[1], 2.) + std::pow(p[2], 2.) - 1.;
            };
            data.surface.n = [](Point3DCL const &p, double) {
                auto den = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                Point3DCL v(p[0] / den, p[1] / den, p[2] / den);
                return v;
            };
            data.surface.e = data.surface.n;
            data.surface.H = [](Point3DCL const &p, double) {
                SMatrixCL<3, 3> res;
                auto den = std::pow(p[0], 2) + std::pow(p[1], 2) + std::pow(p[2], 2);
                res(0, 0) = (std::pow(p[1], 2) + std::pow(p[2], 2)) / den;
                res(0, 1) = (-1. * p[0] * p[1]) / den;
                res(0, 2) = (-1. * p[0] * p[2]) / den;
                res(1, 0) = res(0, 1);
                res(1, 1) = (std::pow(p[0], 2) + std::pow(p[2], 2)) / den;
                res(1, 2) = (-1. * p[1] * p[2]) / den;
                res(2, 0) = res(0, 2);
                res(2, 1) = res(1, 2);
                res(2, 2) = (std::pow(p[0], 2) + std::pow(p[1], 2)) / den;
                return res;
            };
        }
        else if (test.find("Bubble") != std::string::npos) {
            auto r0 = param.get<double>("SurfNavStokes.IC.Bubble.r0");
            auto a  = param.get<double>("SurfNavStokes.IC.Bubble.a");
            data.description = "phi = x^2 + y^2 + z^2 - r^2(t), r(t) = r0 (1 + a sin(2 Pi t))\n";
            data.surface.u_N = [=](Point3DCL const &, double t) {
                return 6.283185307179586*a*r0*cos(6.283185307179586*t);
            };
            data.surface.phi = [=](Point3DCL const & p, double t) {
                return std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2) - 1.*std::pow(r0,2)*std::pow(1. + a*sin(6.283185307179586*t),2);
            };
            data.surface.n = [](Point3DCL const &p, double) {
                auto den = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                Point3DCL v(p[0] / den, p[1] / den, p[2] / den);
                return v;
            };
        }
        else if (test.find("TorusVarTube") != std::string::npos) {
            auto r_0 = param.get<double>("SurfNavStokes.IC.TorusVarTube.r_0");
            auto r_1 = param.get<double>("SurfNavStokes.IC.TorusVarTube.r_1");
            auto R   = param.get<double>("SurfNavStokes.IC.TorusVarTube.R");
            data.description =
                    "torus w/ const distance $R$ from the center of the tube to the center of the torus and\n"
                    "      w/ variable tube radius $" + std::to_string(r_0) + " =: r_0 <= r(\\xi) <= r_1 =: " + std::to_string(r_1) + ", \\xi \\in [0, \\2 pi]\n"
                    "$\\phi = (x^2 + y^2 + z^2 + R^2 - r(\\xi)^2)^2 - 4 R^2 (x^2 + y^2)$\n";
            data.surface.u_N = [](Point3DCL const &, double) {
                return 0.;
            };
            data.surface.phi = [=](Point3DCL const & p, double) {
                // auto xi = arctan(p[0], std::fabs(p[1]));
                // auto r = r_0 * (M_PI - xi) / M_PI + r_1 * xi / M_PI;
                // auto xi = arctan(p[0],p[1]);
                // auto r = r_1 + (r_0 - r_1) / (M_PI * M_PI) * (xi - M_PI) * (xi - M_PI);
                // auto xi = arctan(p[0], p[1]);
                // auto r = (r_1 - r_0) / M_PI * std::sqrt(M_PI * M_PI - (xi - M_PI) * (xi - M_PI)) + r_0;
                // auto xi = arctan(p[0], std::fabs(p[1]));
                // double r, dr, d2r, d3r;
                // hermite_cubic_value(0., r_0, 0., M_PI, r_1, 0., 1, &xi, &r, &dr, &d2r, &d3r);
                auto r = r_0 + .5 * (r_1 - r_0) * (1. - p[0] / std::sqrt(p[0] * p[0] + p[1] * p[1]));
                return std::pow(norm_sq(p) + R * R - r * r, 2.) - 4. * R * R * (std::pow(p[0], 2.) + std::pow(p[1], 2.));
            };
            data.surface.e = [](Point3DCL const &p, double) { // TODO: correct using distance func
                return p;
            };
        }
        else if (test.find("Torus") != std::string::npos) {
            data.description = "phi = (x^2 + y^2 + z^2 + R^2 - r^2)^2 - 4 R^2 (x^2 + y^2), R = 1, r = 0.5\n";
            data.surface.u_N = [](Point3DCL const &, double) {
                return 0.;
            };
            data.surface.phi = [](Point3DCL const & p, double) {
                return -4.*(std::pow(p[0],2) + std::pow(p[1],2)) + std::pow(0.75 + std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
            };
            data.surface.e = [](Point3DCL const &p, double) { // TODO: correct using distance func
                return p;
            };
        }
        else throw std::invalid_argument("unknown surface");
        if (test == "OseenBubbleDirection") {
            data.exactSoln = true;
            data.description +=
                    "u = P (1 - 2 t, 0, 0)^e, p = 0\n"
                    "u = u^e is tangential, mean of p = p^e is zero\n"
                    "wind field is u + u_N n\n";
            data.u_T = [=](Point3DCL const & p, double t) {
                Point3DCL v;
                v[0] = ((1. - 2.*t)*(std::pow(p[1],2) + std::pow(p[2],2)))/(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2));
                v[1] = (-1.*(1. - 2.*t)*p[0]*p[1])/(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2));
                v[2] = (-1.*(1. - 2.*t)*p[0]*p[2])/(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2));
                return v;
            };
            data.w_T = data.u_T;
            data.p = [](Point3DCL const & p, double) {
                return 0.;
            };
            auto r0 = param.get<double>("SurfNavStokes.IC.Bubble.r0");
            auto a  = param.get<double>("SurfNavStokes.IC.Bubble.a");
            data.f_T = [=](Point3DCL const & p, double t) {
                Point3DCL v;
                v[0] = (-1.*(std::pow(p[1],2) + std::pow(p[2],2))*(std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*(std::pow(1. - 2.*t,2)*p[0] + 2.*std::pow(p[0],2) + 2.*std::pow(p[1],2) + 2.*std::pow(p[2],2) - 1.*nu + 2.*t*nu) + 6.283185307179586*a*r0*(-1. + 2.*t)*(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*cos(6.283185307179586*t)))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2.5);
                v[1] = (p[0]*p[1]*(std::pow(1. - 2.*t,2)*p[0] + 2.*std::pow(p[0],2) + 2.*std::pow(p[1],2) + 2.*std::pow(p[2],2) - 1.*nu + 2.*t*nu + 6.283185307179586*a*r0*(-1. + 2.*t)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*cos(6.283185307179586*t)))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
                v[2] = (p[0]*p[2]*(std::pow(1. - 2.*t,2)*p[0] + 2.*std::pow(p[0],2) + 2.*std::pow(p[1],2) + 2.*std::pow(p[2],2) - 1.*nu + 2.*t*nu + 6.283185307179586*a*r0*(-1. + 2.*t)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*cos(6.283185307179586*t)))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
                return v;
            };
            data.m_g = [=](Point3DCL const & p, double t) {
                return (-2.*(-1. + 2.*t)*p[0])/(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) - (12.566370614359172*a*r0*cos(6.283185307179586*t))/std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2));
            };
        }
        else if (test == "StokesSphereSimple" || test == "OseenSphereSimple") {
            data.exactSoln = true;
            data.description +=
                    "u = P (1, 0, 0)^e, p = 0\n"
                    "u = u^e is tangential, mean of p = p^e is zero\n";
            data.u_T = [=](Point3DCL const & p, double) {
                Point3DCL v;
                v[0] = (std::pow(p[1], 2) + std::pow(p[2], 2)) / (std::pow(p[0], 2) + std::pow(p[1], 2) + std::pow(p[2], 2));
                v[1] = (-1. * p[0] * p[1]) / (std::pow(p[0], 2) + std::pow(p[1], 2) + std::pow(p[2], 2));
                v[2] = (-1. * p[0] * p[2]) / (std::pow(p[0], 2) + std::pow(p[1], 2) + std::pow(p[2], 2));
                return v;
            };
            data.w_T = [](Point3DCL const &, double) {
                return Point3DCL(0., 0., 0.);
            };
            data.p = [](Point3DCL const & p, double) {
                return 0.;
            };
            data.f_T = [=](Point3DCL const & p, double) {
                Point3DCL v;
                v[0] = ((std::pow(p[1],2) + std::pow(p[2],2))*nu)/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
                v[1] = (-1.*p[0]*p[1]*nu)/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
                v[2] = (-1.*p[0]*p[2]*nu)/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
                return v;
            };
            data.m_g = [](Point3DCL const & p, double) {
                return (2.*p[0])/(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2));
            };
            if (test == "OseenSphereSimple") {
                data.description += "wind field is the same as soln\n";
                data.w_T = data.u_T;
                data.f_T = [=](Point3DCL const & p, double) {
                    Point3DCL v;
                    v[0] = (-1.*(std::pow(p[1],2) + std::pow(p[2],2))*(p[0] - 1.*nu))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
                    v[1] = (p[0]*p[1]*(p[0] - 1.*nu))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
                    v[2] = (p[0]*p[2]*(p[0] - 1.*nu))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
                    return v;
                };
            }
        }
        else if (test == "StokesSphere" || test == "OseenSphere") {
            data.exactSoln = true;
            data.description +=
                    "u = P (-z^2, y, x)^e, p = (x y^3 + z)^e\n"
                    "u = u^e is tangential, mean of p = p^e is zero\n";
            data.u_T = [=](Point3DCL const & p, double) {
                Point3DCL v;
                v[0] = (-1.*(std::pow(p[2],4) + std::pow(p[0],2)*p[2]*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) + std::pow(p[1],2)*(std::pow(p[2],2) + p[0]*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
                v[1] = (p[1]*(p[0]*std::pow(p[2],2) - 1.*p[0]*p[2]*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) + (std::pow(p[0],2) + std::pow(p[2],2))*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
                v[2] = (p[0]*std::pow(p[2],3) + p[0]*(std::pow(p[0],2) + std::pow(p[1],2))*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) - 1.*std::pow(p[1],2)*p[2]*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
                return v;
            };
            data.w_T = [](Point3DCL const &, double) {
                return Point3DCL(0., 0., 0.);
            };
            data.p = [](Point3DCL const & p, double) {
                return (p[0]*std::pow(p[1],3))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2) + p[2]/std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2));
            };
            data.f_T = [=](Point3DCL const & p, double) {
                Point3DCL v;
                v[0] = (0.5*(-2.*std::pow(p[0],3)*p[2]*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) + std::pow(p[0],2)*(-6.*std::pow(p[1],3) + 2.*std::pow(p[1],2)*nu + p[2]*(14.*p[2] - 5.*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))*nu) + (std::pow(p[1],2) + std::pow(p[2],2))*(2.*std::pow(p[1],3) + 2.*std::pow(p[1],2)*nu + p[2]*(-8.*p[2] + 5.*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))*nu) - 2.*p[0]*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*(std::pow(p[2],3) + std::pow(p[1],2)*(p[2] + 5.*nu))))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),3);
                v[1] = (-1.*p[1]*(p[2]*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*(std::pow(p[1],2) + p[2]*(p[2] - 5.*nu)) + std::pow(p[0],2)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*(p[2] - 5.*nu) + std::pow(p[0],3)*(-3.*p[1] + nu) + p[0]*(std::pow(p[1],3) - 3.*p[1]*std::pow(p[2],2) + std::pow(p[1],2)*nu + 5.*p[2]*(-2.*p[2] + std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))*nu)))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),3);
                v[2] = (0.5*(-8.*p[0]*std::pow(p[1],3)*p[2] + 2.*(std::pow(p[0],2) + std::pow(p[1],2))*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),1.5) - 14.*std::pow(p[0],3)*p[2]*nu - 14.*p[0]*std::pow(p[1],2)*p[2]*nu + 8.*p[0]*std::pow(p[2],3)*nu + 5.*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*(std::pow(p[0],3) - 2.*std::pow(p[1],2)*p[2] + p[0]*(p[1] - 1.*p[2])*(p[1] + p[2]))*nu))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),3);
                return v;
            };
            data.m_g = [](Point3DCL const & p, double) {
                return (-1.*(4.*p[0]*std::pow(p[2],2) - 3.*std::pow(p[1],2)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) - 3.*p[0]*p[2]*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) + std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),1.5)))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2);
            };
            if (test == "OseenSphere") {
                data.description += "wind field is the same as soln\n";
                data.w_T = data.u_T;
                data.f_T = [=](Point3DCL const & p, double) {
                    Point3DCL v;
                    v[0] = (0.5*(-2.*std::pow(p[0],5)*(std::pow(p[1],2) + p[2]*(-2.*p[2] + std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))) - 2.*std::pow(p[0],3)*(-2.*std::pow(p[2],4) + 4.*std::pow(p[1],2)*p[2]*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) + 5.*std::pow(p[2],3)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) + 5.*std::pow(p[1],2)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*nu) + std::pow(p[0],4)*(-6.*std::pow(p[1],3) + p[2]*(14.*p[2] - 5.*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))*nu + 2.*std::pow(p[1],2)*(4.*p[2] + nu)) + 2.*p[0]*(std::pow(p[1],2) + std::pow(p[2],2))*(std::pow(p[1],4) + std::pow(p[2],3)*(-3.*p[2] + std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))) - 1.*std::pow(p[1],2)*(std::pow(p[2],2) + 3.*p[2]*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) + 5.*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*nu)) + (std::pow(p[1],2) + std::pow(p[2],2))*(2.*std::pow(p[1],5) + 2.*std::pow(p[1],3)*std::pow(p[2],2) + 2.*std::pow(p[1],4)*nu + std::pow(p[2],3)*(-8.*p[2] + 5.*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))*nu + std::pow(p[1],2)*p[2]*(6.*p[2]*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) - 6.*p[2]*nu + 5.*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*nu)) + 2.*std::pow(p[0],2)*(-2.*std::pow(p[1],2)*std::pow(p[2],2)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) - 1.*(std::pow(p[1],2) + std::pow(p[2],2))*(2.*std::pow(p[1],3) - 3.*std::pow(p[2],2)*nu - 2.*std::pow(p[1],2)*(2.*p[2] + nu)))))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),4);
                    v[1] = (p[1]*(-3.*std::pow(p[2],4)*(std::pow(p[1],2) + std::pow(p[2],2)) - 5.*p[0]*std::pow(p[1],2)*std::pow(p[2],2)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) + 5.*std::pow(p[1],2)*std::pow(p[2],3)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) + 5.*std::pow(p[2],5)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) - 2.*std::pow(p[1],2)*p[2]*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),1.5) + 2.*p[0]*std::pow(p[2],2)*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),1.5) - 6.*std::pow(p[2],3)*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),1.5) + p[2]*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2.5) + std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),3) - 5.*std::pow(p[1],2)*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),1.5)*nu - 5.*p[0]*p[2]*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),1.5)*nu + 5.*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2.5)*nu - 1.*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2)*(3.*std::pow(p[1],2) - 2.*std::pow(p[2],2) + p[0]*(-3.*p[1] + 2.*p[2] + nu)) + (std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*(2.*std::pow(p[1],4) - 2.*std::pow(p[1],2)*std::pow(p[2],2) + std::pow(p[2],4) + p[0]*(-4.*std::pow(p[1],3) + 4.*std::pow(p[1],2)*p[2] + 11.*std::pow(p[2],2)*nu))))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),4);
                    v[2] = (-0.5*(6.*std::pow(p[2],5)*(std::pow(p[1],2) + std::pow(p[2],2)) + 10.*p[0]*std::pow(p[1],2)*std::pow(p[2],3)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) - 10.*std::pow(p[1],2)*std::pow(p[2],4)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) - 10.*std::pow(p[2],6)*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)) + 8.*std::pow(p[1],2)*std::pow(p[2],2)*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),1.5) + 16.*std::pow(p[2],4)*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),1.5) - 4.*std::pow(p[2],2)*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2.5) + 4.*p[2]*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),3) - 2.*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),3.5) + 10.*std::pow(p[1],2)*p[2]*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),1.5)*nu + 10.*p[0]*std::pow(p[2],2)*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),1.5)*nu - 5.*p[0]*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2.5)*nu + 2.*std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),2)*(2.*p[0]*std::pow(p[1],2) - 1.*std::pow(p[1],2)*p[2] - 4.*std::pow(p[2],3) + 7.*p[0]*p[2]*nu) - 2.*p[2]*(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*(2.*std::pow(p[1],4) - 2.*std::pow(p[1],2)*std::pow(p[2],2) + std::pow(p[2],4) + p[0]*(-4.*std::pow(p[1],3) + 4.*std::pow(p[1],2)*p[2] + 11.*std::pow(p[2],2)*nu))))/std::pow(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2),4);
                    return v;
                };
            }
        }
        else if (test == "KelvinHelmholtzSphere") {
            data.exactSoln = false;
            data.description += "u_0 = Kelvin-Helmholtz, p = 0\n";
            auto ang = [](double x, double y) {
                auto pi = 3.1415926535897932;
                if (x > 0 && y >= 0) return std::atan(y/x);
                if (x > 0 && y < 0)  return std::atan(y/x) + 2. * pi;
                if (x < 0)           return std::atan(y/x) + pi;
                return sign(y) * pi / 2.;
            };
            auto freq = param.get<double>("SurfNavStokes.IC." + test + ".Frequency");
            auto ampl = param.get<double>("SurfNavStokes.IC." + test + ".Amplitude");
            data.u_T = [=](Point3DCL const & p, double) {
                Point3DCL v(0., 0., 0.);
                auto x = p / norm(p);
                if (x[2] < ampl * std::sin(freq * ang(x[0], x[1]))) return v;
                v[0] = -x[1];
                v[1] = x[0];
                return v;
            };
            data.p = [](Point3DCL const &, double) {
                return 0.;
            };
            data.f_T = [=](Point3DCL const & p, double) {
                return Point3DCL(0., 0., 0.);
            };
            data.m_g = [](Point3DCL const & p, double) {
                return 0.;
            };
        }
        else if (test.find("KelvinHelmholtzCristoph") != std::string::npos) {
            data.description += "u_0 = Kelvin-Helmholtz, p = 0\n";
            data.exactSoln = false;
            auto eta = [=](Point3DCL const & p) {
                return -0.3183098861837907*std::asin(p[2]/std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)));
            };
            auto xi = [=](Point3DCL const & p) {
                return 0.15915494309189535*arctan(p[0],p[1]);
            };
            auto eXi = [=](Point3DCL const & p) {
                if (std::pow(p[0],2) + std::pow(p[1],2) == 0.) return Point3DCL(0., 0., 0.);
                return Point3DCL(-1.*p[1]*std::sqrt(1/(std::pow(p[0],2) + std::pow(p[1],2))), p[0]*std::sqrt(1/(std::pow(p[0],2) + std::pow(p[1],2))), 0.);
            };
            auto delta_0 = param.get<double>("SurfNavStokes.IC.KelvinHelmholtzCristoph.Delta_0");
            auto cn = param.get<double>("SurfNavStokes.IC.KelvinHelmholtzCristoph.cn");
            auto aa = param.get<double>("SurfNavStokes.IC.KelvinHelmholtzCristoph.aa");
            auto ab = param.get<double>("SurfNavStokes.IC.KelvinHelmholtzCristoph.ab");
            auto ma = param.get<double>("SurfNavStokes.IC.KelvinHelmholtzCristoph.ma");
            auto mb = param.get<double>("SurfNavStokes.IC.KelvinHelmholtzCristoph.mb");
            auto Hs = [=](double eta) {
                return std::tanh(2. * eta / delta_0);
            };
            auto pert = [=](Point3DCL const & p) {
                return Point3DCL(
                    (0.05066059182116889*(-4.*p[1]*std::sqrt((std::pow(p[0],2) + std::pow(p[1],2))/(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*std::asin(p[2]/std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))*(aa*cos(0.5*ma*arctan(p[0],p[1])) + ab*cos(0.5*mb*arctan(p[0],p[1]))) + 9.869604401089358*p[0]*p[2]*std::pow(delta_0,2)*(aa*ma*sin(0.5*ma*arctan(p[0],p[1])) + ab*mb*sin(0.5*mb*arctan(p[0],p[1])))))/(std::pow(2.718281828459045,(0.10132118364233778*std::pow(std::asin(p[2]/std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))),2))/std::pow(delta_0,2))*(std::pow(p[0],2) + std::pow(p[1],2))*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*std::pow(delta_0,2)),
                    (0.05066059182116889*(4.*p[0]*std::sqrt((std::pow(p[0],2) + std::pow(p[1],2))/(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*std::asin(p[2]/std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))*(aa*cos(0.5*ma*arctan(p[0],p[1])) + ab*cos(0.5*mb*arctan(p[0],p[1]))) + 9.869604401089358*p[1]*p[2]*std::pow(delta_0,2)*(aa*ma*sin(0.5*ma*arctan(p[0],p[1])) + ab*mb*sin(0.5*mb*arctan(p[0],p[1])))))/(std::pow(2.718281828459045,(0.10132118364233778*std::pow(std::asin(p[2]/std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))),2))/std::pow(delta_0,2))*(std::pow(p[0],2) + std::pow(p[1],2))*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))*std::pow(delta_0,2)),
                    (-0.5*(aa*ma*sin(0.5*ma*arctan(p[0],p[1])) + ab*mb*sin(0.5*mb*arctan(p[0],p[1]))))/(std::pow(2.718281828459045,(0.10132118364233778*std::pow(std::asin(p[2]/std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2))),2))/std::pow(delta_0,2))*std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)))
                );
            };
            auto zDistMin = 0.;
            if (test.find("Torus") != std::string::npos)
                zDistMin = .5;
            auto zDist = [=](Point3DCL const & p) {
                return std::sqrt(p[0] * p[0] + p[1] * p[1]) - zDistMin;
            };
            data.u_T = [=](Point3DCL const & p, double) {
                auto x = data.surface.e(p, 0.);
                return zDist(x) * (Hs(eta(x)) * eXi(x) + cn * pert(x));
            };
            data.p = [](Point3DCL const &, double) {
                return 0.;
            };
            data.f_T = [=](Point3DCL const & p, double) {
                return Point3DCL(0., 0., 0.);
            };
            data.m_g = [](Point3DCL const & p, double) {
                return 0.;
            };
        }
        else throw std::invalid_argument("test '" + test + "' is not defined");
        data.description += data.exactSoln ? "exact solution is available and the errors will be computed\n" : "exact solution is NOT available and the errors will NOT be computed\n";
        return data;
    }

}

#endif // SURF_NAVIER_STOKES_DATA_HPP