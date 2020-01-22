//
// Created by Alexander Zhiliakov on 9/27/19.
//

#ifndef SURF_NAVIER_STOKES_DATA_HPP
#define SURF_NAVIER_STOKES_DATA_HPP

#include "../num/discretize.h"
// #include "surfnavierstokes_funcs.h"

namespace DROPS {

    struct SurfNavierStokesData {
        bool exactSoln;
        InstatVectorFunction u_T, f_T, w_T = nullptr; // w_T = nullptr for the Navier-Stokes case
        InstatScalarFunction u_N, p, m_g; // m_g is "-g"
        struct Surface {
            InstatScalarFunction phi;
            InstatVectorFunction n;
            InstatVectorFunction e;
            InstatMatrixFunction H;
        };
        Surface surface;
        std::string description;
    };

    SurfNavierStokesData SurfNavierStokesDataFactory(std::string const & test, double nu, ParamCL const & param) {
        SurfNavierStokesData data;
        SurfNavierStokesData::Surface sphere;
        sphere.phi = [](Point3DCL const & p, double) {
            return std::pow(p[0], 2.) + std::pow(p[1], 2.) + std::pow(p[2], 2.) - 1.;
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
        if (test == "StokesSphereSimple" || test == "OseenSphereSimple") {
            data.exactSoln = true;
            data.description =
                    "phi = x^2 + y^2 + z^2 - 1, u = P (1, 0, 0)^e, p = 0\n"
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
            data.u_N = [](Point3DCL const &, double) {
                return 0.;
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
            data.surface = sphere;
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
            data.description =
                    "phi = x^2 + y^2 + z^2 - 1, u = P (-z^2, y, x)^e, p = (x y^3 + z)^e\n"
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
            data.u_N = [](Point3DCL const &, double) {
                return 0.;
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
            data.surface = sphere;
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
            data.description =
                    "phi = x^2 + y^2 + z^2 - 1, u_0 = Kelvin-Helmholtz, p = 0\n";
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
            data.u_N = [](Point3DCL const &, double) {
                return 0.;
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
            data.surface = sphere;
        }
        else if (test == "KelvinHelmholtzCristophSphere") {
            data.exactSoln = false;
            data.description = "phi = x^2 + y^2 + z^2 - 1, u_0 = Kelvin-Helmholtz, p = 0\n";
            auto eta = [=](Point3DCL const & p) {
                return -0.3183098861837907*std::asin(p[2]/std::sqrt(std::pow(p[0],2) + std::pow(p[1],2) + std::pow(p[2],2)));
            };
            auto arctan = [](double x, double y) {
                auto pi = 3.1415926535897932;
                if (x > 0 && y >= 0) return std::atan(y/x);
                if (x > 0 && y < 0)  return std::atan(y/x) + 2. * pi;
                if (x < 0)           return std::atan(y/x) + pi;
                return sign(y) * pi / 2.;
            };
            auto xi = [=](Point3DCL const & p) {
                return 0.15915494309189535*arctan(p[0],p[1]);
            };
            auto eXi = [=](Point3DCL const & p) {
                if (std::pow(p[0],2) + std::pow(p[1],2) == 0.) return Point3DCL(0., 0., 0.);
                return Point3DCL(-1.*p[1]*std::sqrt(1/(std::pow(p[0],2) + std::pow(p[1],2))), p[0]*std::sqrt(1/(std::pow(p[0],2) + std::pow(p[1],2))), 0.);
            };
            auto delta_0 = param.get<double>("SurfNavStokes.IC." + test + ".Delta_0");
            auto cn = param.get<double>("SurfNavStokes.IC." + test + ".cn");
            auto aa = param.get<double>("SurfNavStokes.IC." + test + ".aa");
            auto ab = param.get<double>("SurfNavStokes.IC." + test + ".ab");
            auto ma = param.get<double>("SurfNavStokes.IC." + test + ".ma");
            auto mb = param.get<double>("SurfNavStokes.IC." + test + ".mb");
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
            data.u_T = [=](Point3DCL const & p, double) {
                auto x = sphere.e(p, 0.);
                return Hs(eta(x)) * std::sqrt(1. - x[2] * x[2]) * eXi(x) + cn * pert(x);
            };
            data.u_N = [](Point3DCL const &, double) {
                return 0.;
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
            data.surface = sphere;
        }
        else throw std::invalid_argument("test '" + test + "' is not defined");
        data.description += data.exactSoln ? "exact solution is available and the errors will be computed\n" : "exact solution is NOT available and the errors will NOT be computed\n";
        return data;
    }

}

#endif // SURF_NAVIER_STOKES_DATA_HPP