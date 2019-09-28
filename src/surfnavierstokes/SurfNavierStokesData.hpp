//
// Created by Alexander Zhiliakov on 9/27/19.
//

#ifndef SURF_NAVIER_STOKES_DATA_HPP
#define SURF_NAVIER_STOKES_DATA_HPP

#include "../num/discretize.h"

namespace DROPS {

    struct SurfNavierStokesData {
        instat_vector_fun_ptr u_T, f_T;
        instat_scalar_fun_ptr u_N, p, m_g; // m_g is "-g"
        struct {
            instat_scalar_fun_ptr phi;
            instat_vector_fun_ptr n;
            instat_matrix_fun_ptr H;
        } surface;
        std::string description;
    };

    SurfNavierStokesData SurfNavierStokesDataFactory(std::string const & test, double rho, double mu) {
        SurfNavierStokesData data;
        if (test == "85") {
            // TODO: add mu and rho!
            data.u_T = [](const DROPS::Point3DCL& p, double) {
                DROPS::Point3DCL v;
                v[0] = (-1.*(pow(p[2],4) + pow(p[0],2)*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + pow(p[1],2)*(pow(p[2],2) + p[0]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)))))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
                v[1] = (p[1]*(p[0]*pow(p[2],2) - 1.*p[0]*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + (pow(p[0],2) + pow(p[2],2))*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
                v[2] = (p[0]*pow(p[2],3) + p[0]*(pow(p[0],2) + pow(p[1],2))*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) - 1.*pow(p[1],2)*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
                return v;
            };
            data.p = [](const DROPS::Point3DCL& p, double) {
                return (p[0]*pow(p[1],3))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2) + p[2]/std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2));
            };
            data.f_T = [](const DROPS::Point3DCL& p, double) {
                DROPS::Point3DCL v;
                v[0] = (-1.*(pow(p[0],4)*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + pow(p[0],3)*(pow(p[1],2) + p[2])*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + p[0]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))*(pow(p[1],4) + pow(p[2],3) + pow(p[1],2)*(10. + p[2] + pow(p[2],2))) - 1.*(pow(p[1],2) + pow(p[2],2))*(pow(p[1],3) - 1.*pow(p[1],2)*(-2. + pow(p[2],2)) - 1.*p[2]*(8.*p[2] + pow(p[2],3) - 5.*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)))) + pow(p[0],2)*(3.*pow(p[1],3) + pow(p[1],2)*(-2. + p[2]*(p[2] + std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)))) + p[2]*(5.*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + p[2]*(-14. + p[2]*(p[2] + std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2))))))))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),3);
                v[1] = (p[1]*(-2.*pow(p[0],3) + 3.*pow(p[0],3)*p[1] - 2.*p[0]*pow(p[1],2) - 1.*p[0]*pow(p[1],3) + 20.*p[0]*pow(p[2],2) + pow(p[0],3)*pow(p[2],2) + 3.*p[0]*p[1]*pow(p[2],2) + p[0]*pow(p[1],2)*pow(p[2],2) + p[0]*pow(p[2],4) + 10.*pow(p[0],2)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) - 10.*p[0]*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + 10.*pow(p[2],2)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + pow(p[0],2)*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) - 1.*p[2]*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) - 1.*p[0]*p[2]*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) + pow(p[2],2)*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5)))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),3);
                v[2] = (-4.*p[0]*pow(p[1],3)*p[2] + 22.*p[0]*pow(p[2],3) - 10.*pow(p[1],2)*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) - 10.*p[0]*pow(p[2],2)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) - 14.*p[0]*p[2]*(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + p[0]*pow(p[2],3)*(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + 5.*p[0]*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) - 1.*pow(p[1],2)*p[2]*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) - 1.*pow(p[2],2)*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) - 1.*p[0]*pow(p[2],2)*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5) + pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2.5) + p[0]*pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2.5))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),3);
                return v;
            };
            data.m_g = [](const DROPS::Point3DCL& p, double) {
                return (-1.*(4.*p[0]*pow(p[2],2) - 3.*pow(p[1],2)*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) - 3.*p[0]*p[2]*std::sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) + pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),1.5)))/pow(pow(p[0],2) + pow(p[1],2) + pow(p[2],2),2);
            };
            data.surface.phi = [](const DROPS::Point3DCL& p, double) {
                return pow(p[0], 2.) + pow(p[1], 2.) + pow(p[2], 2.) - 1.;
            };
            data.surface.n = [](const DROPS::Point3DCL& p, double) {
                auto den = std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
                DROPS::Point3DCL v(p[0]/den, p[1]/den, p[2]/den);
                return v;
            };
            data.surface.H = [](const DROPS::Point3DCL& p, double) {
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
        }
        else throw std::invalid_argument("test '" + test + "' is not defined");
        return data;
    }

}

#endif // SURF_NAVIER_STOKES_DATA_HPP