//
// Created by Alexander Zhiliakov on 9/27/19.
//

#ifndef SURF_NAVIER_STOKES_DATA_HPP
#define SURF_NAVIER_STOKES_DATA_HPP

#include "../num/discretize.h"
#include "surfnavierstokes_funcs.h"

namespace DROPS {

    struct SurfNavierStokesData {
        bool exactSoln;
        instat_vector_fun_ptr u_T, f_T, w_T = nullptr; // w_T = nullptr for the Navier-Stokes case
        instat_scalar_fun_ptr u_N, p, m_g; // m_g is "-g"
        struct Surface {
            instat_scalar_fun_ptr phi;
            instat_vector_fun_ptr n;
            instat_matrix_fun_ptr H;
        };
        Surface surface;
        std::string description;
    };

    SurfNavierStokesData SurfNavierStokesDataFactory(std::string const & test, double nu) {
        SurfNavierStokesData data;
        SurfNavierStokesData::Surface sphere;
        sphere.phi = [](Point3DCL const & p, double) {
            return pow(p[0], 2.) + pow(p[1], 2.) + pow(p[2], 2.) - 1.;
        };
        sphere.n = [](Point3DCL const & p, double) {
            auto den = std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
            Point3DCL v(p[0]/den, p[1]/den, p[2]/den);
            return v;
        };
        sphere.H = [](Point3DCL const & p, double) {
            SMatrixCL<3, 3> res;
            auto den = pow(p[0],2) + pow(p[1],2) + pow(p[2],2);
            res(0, 0) = (pow(p[1],2) + pow(p[2],2)) / den;
            res(0, 1) = (-1.*p[0]*p[1]) / den;
            res(0, 2) = (-1.*p[0]*p[2]) / den;
            res(1, 0) = res(0, 1);
            res(1, 1) = (pow(p[0],2) + pow(p[2],2)) / den;
            res(1, 2) = (-1.*p[1]*p[2]) / den;
            res(2, 0) = res(0, 2);
            res(2, 1) = res(1, 2);
            res(2, 2) = (pow(p[0],2) + pow(p[1],2)) / den;
            return res;
        };
        if (test == "StokesSphereSimple" || test == "OseenSphereSimple") {
            data.exactSoln = true;
            data.description =
                    "phi = x^2 + y^2 + z^2 - 1, u = P (1, 0, 0)^e, p = 0\n"
                    "u = u^e is tangential, mean of p = p^e is zero\n";
            data.u_T = [=](Point3DCL const & p, double) {
                Point3DCL v;
                v[0] = (pow(p[1], 2) + pow(p[2], 2)) / (pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2));
                v[1] = (-1. * p[0] * p[1]) / (pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2));
                v[2] = (-1. * p[0] * p[2]) / (pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2));
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
                v[0] = ((pow(p[1],2) + pow(p[2],2))*nu)/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
                v[1] = (-1.*p[0]*p[1]*nu)/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
                v[2] = (-1.*p[0]*p[2]*nu)/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
                return v;
            };
            data.m_g = [](Point3DCL const & p, double) {
                return (2.*p[0])/(pow(p[0],2) + pow(p[1],2) + pow(p[2],2));
            };
            data.surface = sphere;
            if (test == "OseenSphereSimple") {
                data.description += "wind field is the same as soln\n";
                data.w_T = data.u_T;
                data.f_T = [=](Point3DCL const & p, double) {
                    Point3DCL v;
                    v[0] = (-1.*(pow(p[1],2) + pow(p[2],2))*(p[0] - 1.*nu))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
                    v[1] = (p[0]*p[1]*(p[0] - 1.*nu))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
                    v[2] = (p[0]*p[2]*(p[0] - 1.*nu))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
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
                v[0] = (-1.*(pow(p[2],4) + pow(p[0],2)*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + pow(p[1],2)*(pow(p[2],2) + p[0]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)))))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
                v[1] = (p[1]*(p[0]*pow(p[2],2) - 1.*p[0]*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + (pow(p[0],2) + pow(p[2],2))*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
                v[2] = (p[0]*pow(p[2],3) + p[0]*(pow(p[0],2) + pow(p[1],2))*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) - 1.*pow(p[1],2)*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
                return v;
            };
            data.w_T = [](Point3DCL const &, double) {
                return Point3DCL(0., 0., 0.);
            };
            data.u_N = [](Point3DCL const &, double) {
                return 0.;
            };
            data.p = [](Point3DCL const & p, double) {
                return (p[0]*pow(p[1],3))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2) + p[2]/std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2));
            };
            data.f_T = [=](Point3DCL const & p, double) {
                Point3DCL v;
                v[0] = (0.5*(-2.*pow(p[0],3)*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + pow(p[0],2)*(-6.*pow(p[1],3) + 2.*pow(p[1],2)*nu + p[2]*(14.*p[2] - 5.*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)))*nu) + (pow(p[1],2) + pow(p[2],2))*(2.*pow(p[1],3) + 2.*pow(p[1],2)*nu + p[2]*(-8.*p[2] + 5.*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)))*nu) - 2.*p[0]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))*(pow(p[2],3) + pow(p[1],2)*(p[2] + 5.*nu))))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),3);
                v[1] = (-1.*p[1]*(p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))*(pow(p[1],2) + p[2]*(p[2] - 5.*nu)) + pow(p[0],2)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))*(p[2] - 5.*nu) + pow(p[0],3)*(-3.*p[1] + nu) + p[0]*(pow(p[1],3) - 3.*p[1]*pow(p[2],2) + pow(p[1],2)*nu + 5.*p[2]*(-2.*p[2] + std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)))*nu)))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),3);
                v[2] = (0.5*(-8.*p[0]*pow(p[1],3)*p[2] + 2.*(pow(p[0],2) + pow(p[1],2))*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) - 14.*pow(p[0],3)*p[2]*nu - 14.*p[0]*pow(p[1],2)*p[2]*nu + 8.*p[0]*pow(p[2],3)*nu + 5.*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))*(pow(p[0],3) - 2.*pow(p[1],2)*p[2] + p[0]*(p[1] - 1.*p[2])*(p[1] + p[2]))*nu))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),3);
                return v;
            };
            data.m_g = [](Point3DCL const & p, double) {
                return (-1.*(4.*p[0]*pow(p[2],2) - 3.*pow(p[1],2)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) - 3.*p[0]*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5)))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
            };
            data.surface = sphere;
            if (test == "OseenSphere") {
                data.description += "wind field is the same as soln\n";
                data.w_T = data.u_T;
                data.f_T = [=](Point3DCL const & p, double) {
                    Point3DCL v;
                    v[0] = (0.5*(-2.*pow(p[0],5)*(pow(p[1],2) + p[2]*(-2.*p[2] + std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)))) - 2.*pow(p[0],3)*(-2.*pow(p[2],4) + 4.*pow(p[1],2)*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + 5.*pow(p[2],3)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + 5.*pow(p[1],2)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))*nu) + pow(p[0],4)*(-6.*pow(p[1],3) + p[2]*(14.*p[2] - 5.*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)))*nu + 2.*pow(p[1],2)*(4.*p[2] + nu)) + 2.*p[0]*(pow(p[1],2) + pow(p[2],2))*(pow(p[1],4) + pow(p[2],3)*(-3.*p[2] + std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))) - 1.*pow(p[1],2)*(pow(p[2],2) + 3.*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + 5.*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))*nu)) + (pow(p[1],2) + pow(p[2],2))*(2.*pow(p[1],5) + 2.*pow(p[1],3)*pow(p[2],2) + 2.*pow(p[1],4)*nu + pow(p[2],3)*(-8.*p[2] + 5.*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)))*nu + pow(p[1],2)*p[2]*(6.*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) - 6.*p[2]*nu + 5.*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))*nu)) + 2.*pow(p[0],2)*(-2.*pow(p[1],2)*pow(p[2],2)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) - 1.*(pow(p[1],2) + pow(p[2],2))*(2.*pow(p[1],3) - 3.*pow(p[2],2)*nu - 2.*pow(p[1],2)*(2.*p[2] + nu)))))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),4);
                    v[1] = (p[1]*(-3.*pow(p[2],4)*(pow(p[1],2) + pow(p[2],2)) - 5.*p[0]*pow(p[1],2)*pow(p[2],2)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + 5.*pow(p[1],2)*pow(p[2],3)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + 5.*pow(p[2],5)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) - 2.*pow(p[1],2)*p[2]*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) + 2.*p[0]*pow(p[2],2)*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) - 6.*pow(p[2],3)*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) + p[2]*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2.5) + pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),3) - 5.*pow(p[1],2)*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5)*nu - 5.*p[0]*p[2]*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5)*nu + 5.*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2.5)*nu - 1.*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2)*(3.*pow(p[1],2) - 2.*pow(p[2],2) + p[0]*(-3.*p[1] + 2.*p[2] + nu)) + (pow(p[0],2) + pow(p[1],2) + pow(p[2],2))*(2.*pow(p[1],4) - 2.*pow(p[1],2)*pow(p[2],2) + pow(p[2],4) + p[0]*(-4.*pow(p[1],3) + 4.*pow(p[1],2)*p[2] + 11.*pow(p[2],2)*nu))))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),4);
                    v[2] = (-0.5*(6.*pow(p[2],5)*(pow(p[1],2) + pow(p[2],2)) + 10.*p[0]*pow(p[1],2)*pow(p[2],3)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) - 10.*pow(p[1],2)*pow(p[2],4)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) - 10.*pow(p[2],6)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + 8.*pow(p[1],2)*pow(p[2],2)*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) + 16.*pow(p[2],4)*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) - 4.*pow(p[2],2)*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2.5) + 4.*p[2]*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),3) - 2.*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),3.5) + 10.*pow(p[1],2)*p[2]*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5)*nu + 10.*p[0]*pow(p[2],2)*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5)*nu - 5.*p[0]*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2.5)*nu + 2.*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2)*(2.*p[0]*pow(p[1],2) - 1.*pow(p[1],2)*p[2] - 4.*pow(p[2],3) + 7.*p[0]*p[2]*nu) - 2.*p[2]*(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))*(2.*pow(p[1],4) - 2.*pow(p[1],2)*pow(p[2],2) + pow(p[2],4) + p[0]*(-4.*pow(p[1],3) + 4.*pow(p[1],2)*p[2] + 11.*pow(p[2],2)*nu))))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),4);
                    return v;
                };
            }
        }
        else if (test == "KelvinHelmholtzSphere") {
            data.exactSoln = false;
            data.description =
                    "phi = x^2 + y^2 + z^2 - 1, u_0 = K-H, p = 0\n";
            data.u_T = Test_A_plus_M_vSolVectorFun18;
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