//
// Created by Alexander Zhiliakov on 2/24/21.
//

#ifndef SURFACE_HPP
#define SURFACE_HPP

#include <cmath>
#include <array>
#include <memory>
#include "num/functions.hpp"

namespace DROPS {

    struct Surface {
    protected:
        const bool isStationary_;
        std::string description_;
        struct interval { double min, max; };
        using box = std::array<interval, 3>;
    public:
        Surface(bool isStationary) : isStationary_(isStationary), description_(isStationary_ ? "stationary surface\n" : "evolving surface\n") {}
        bool isStationary() const { return isStationary_; }
        std::string description() const { return description_; }
        virtual box bounds(double t = 0.) const {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": not implemented");
        }
        virtual std::array<bool, 3> rotationalInvariance(double t = 0.) const {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": not implemented");
        }
        virtual double u_N(Point3DCL const &, double t = 0.) const {
            std::string funcName = __func__;
            if (!isStationary_) throw std::logic_error(funcName + ": not implemented");
            return 0.;
        }
        virtual double dist(Point3DCL const &, double t = 0.) const = 0;
        virtual Point3DCL dist_grad(Point3DCL const &, double t = 0.) const {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": normal vector function is not available");
        }
        SMatrixCL<3, 3> P(Point3DCL const & x, double t = 0.) const {
            return eye<3, 3>() - outer_product(dist_grad(x, t), dist_grad(x, t));
        }
        Point3DCL ext(Point3DCL const & x, double t = 0.) const {
            try {
                return x - dist(x, t) * dist_grad(x, t);
            } catch(...) {
                return x;
            }
        }
    };

    struct Sphere : Surface {
        const double r_0, A;
        double r(double t) const { return r_0 * (1. + A * sin(2. * M_PI * t)); }
        double r_prime(double t) const { return 2. * A * M_PI * r_0 * cos(2. * M_PI * t); }
        explicit Sphere(double r_0 = 1., double A = 0.) : Surface(!A), r_0(r_0), A(A) {
            std::string funcName = __func__;
            if (r_0 <= 0.) throw std::invalid_argument(funcName + ": bounds must be positive");
            if (A < 0. || A >= 1.) throw std::invalid_argument(funcName + ": A must be in [0, 1)");
            description_ += isStationary_ ? "phi = x^2 + y^2 + z^2 - R, R = " + std::to_string(r_0) :
                "phi = x^2 + y^2 + z^2 - r^2(t), r(t) := r0 (1 + A sin(2 pi t))\n"
                "A = " + std::to_string(A) + ", r_0 = " + std::to_string(r_0);
        }
        virtual ~Sphere() {}
        box bounds(double t) const final { return {0., r(t), 0., r(t), 0., r(t) }; }
        std::array<bool, 3> rotationalInvariance(double) const final { return { true, true, true }; }
        double u_N(Point3DCL const &, double t) const final { return r_prime(t); }
        double dist(Point3DCL const & x, double t) const final { return norm(x) - r(t); }
        Point3DCL dist_grad(Point3DCL const & x, double) const final { return x / norm(x); }
    };

    struct OscillatingInextensibleSphere : Surface {
        const double r_0, eps, eps_unsym, omega, omega_unsym;
        double H_2(Point3DCL const & x) const { return .25 * std::sqrt(5. / M_PI) * (-pow(x[0],2) - pow(x[1],2) + 2.*pow(x[2],2))/norm_sq(x); };
        double H_3(Point3DCL const & x) const { return .25 * std::sqrt(7. / M_PI) * (x[2]*(-3.*pow(x[0],2) - 3.*pow(x[1],2) + 2.*pow(x[2],2)))/pow(norm_sq(x),1.5); };
        double H_31(Point3DCL const & x) const { return .25 * std::sqrt(10.5 / M_PI) * x[0] * (4. * pow(x[2], 2) - pow(x[0], 2) - pow(x[1], 2))/pow(norm_sq(x), 1.5); };
        double H_42(Point3DCL const & x) const { return .375 * std::sqrt(5. / M_PI) * (pow(x[0], 2) - pow(x[1], 2)) * (7. * pow(x[2], 2) - norm_sq(x)) / pow(norm_sq(x), 2); };
        double A_2(double t) const { return .5 * cos(omega * 2. * M_PI * t); }
        double A_2_prime(double t) const { return -omega * M_PI * sin(omega * 2. * M_PI * t); }
        double A_3(double t) const { return pow(10., -.5) * sin(omega * 2. * M_PI * t); }
        double A_3_prime(double t) const { return std::sqrt(.4) * omega * M_PI * cos(omega * 2. * M_PI * t); }
        double A_31(double t) const { return .5 * cos(omega_unsym * 2. * M_PI * t); }
        double A_31_prime(double t) const { return -omega_unsym * M_PI * sin(omega_unsym * 2. * M_PI * t); }
        double A_42(double t) const { return .5 * sin(omega_unsym * 2. * M_PI * t); }
        double A_42_prime(double t) const { return omega_unsym * M_PI * cos(omega_unsym * 2. * M_PI * t); }
        double r(Point3DCL const & x, double t) const { return r_0 + eps * (A_2(t) * H_2(x) + A_3(t) * H_3(x)) + eps_unsym * (A_42(t) * H_42(x) + A_31(t) * H_31(x)); }
        explicit OscillatingInextensibleSphere(double r_0 = 1., double eps = 0.1, double omega = 1., double eps_unsym = 0.1, double omega_unsym = 0.1) : Surface(!eps && !eps_unsym), r_0(r_0), eps(eps), eps_unsym(eps_unsym), omega(omega), omega_unsym(omega_unsym) {
            std::string funcName = __func__;
            if (eps < 0.) throw std::invalid_argument(funcName + ": eps must be >= 0");
            if (r_0 <= eps) throw std::invalid_argument(funcName + ": eps must be << r_0");
            if (eps_unsym < 0.) throw std::invalid_argument(funcName + ": eps_unsym must be >= 0");
            if (r_0 <= eps_unsym) throw std::invalid_argument(funcName + ": eps_unsym must be << r_0");
            description_ += isStationary_ ? "phi = x^2 + y^2 + z^2 - R, R = " + std::to_string(r_0) : "inextensible sphere";
        }
        virtual ~OscillatingInextensibleSphere() {}
        double u_N(Point3DCL const & x, double t) const final { return eps * (A_2_prime(t) * H_2(x) + A_3_prime(t) * H_3(x)) + eps_unsym * (A_42_prime(t) * H_42(x) + A_31_prime(t) * H_31(x)); }
        double dist(Point3DCL const & x, double t) const final { return norm(x) - r(x, t); }
    };

    struct Torus : Surface {
        const double R, r_min, r_max;
        const size_t x_0, x_1, x_2;
        double r(double x, double y) const { return r_min + .5 * (r_max - r_min) * (1. - x / std::sqrt(x * x + y * y)); }
        std::array<double, 2> r_grad(double x, double y) const { return { .5 * (r_min - r_max) * y * y / pow(x * x + y * y, 1.5), .5 * (r_max - r_min) * x * y / pow(x * x + y * y, 1.5) }; }
        explicit Torus(double R = 1., double r_min = .5, double r_max = .5, size_t axis = 0) : Surface(true), R(R), r_min(r_min), r_max(r_max), x_2(axis), x_0((axis + 1) % 3), x_1((axis + 2) % 3) {
            std::string funcName = __func__;
            if (r_min <= 0. || r_max < r_min || R <= r_max) throw std::invalid_argument(funcName + ": invalid torus params");
            if (axis > 2) throw std::invalid_argument(funcName + ": invalid torus axis");
            description_ +=
                "torus w/ const distance R from the center of the tube to the origin and\n"
                "w/ variable tube rad " + std::to_string(r_min) + " =: r_min <= r(xi) <= r_max =: " + std::to_string(r_max) + ", xi in [0, 2 pi]\n"
                "phi = (x^2 + y^2 + z^2 + R^2 - r(xi)^2)^2 - 4 R^2 (x^2 + y^2)";
        }
        virtual ~Torus() {}
        box bounds(double) const final {
            box res;
            res[x_0] = { R - r_max, R + r_max };
            res[x_1] = { R - .5 * (r_max + r_min), R + .5 * (r_max + r_min) };
            res[x_2] = { 0., r_max };
            return res;
        }
        std::array<bool, 3> rotationalInvariance(double) const final {
            std::array<bool, 3> res = { false, false, false };
            res[x_2] = r_min == r_max;
            return res;
        }
        double dist(Point3DCL const & x, double) const final { return std::sqrt(pow(std::sqrt(x[x_0] * x[x_0] + x[x_1] * x[x_1]) - R, 2.) + x[x_2] * x[x_2]) - r(x[x_0], x[x_1]); }
        Point3DCL dist_grad(Point3DCL const & x, double) const final {
            auto [r_x, r_y] = r_grad(x[x_0], x[x_1]);
            return Point3DCL(
                -r_x + (x[x_0] * (-R + std::sqrt(pow(x[x_0], 2) + pow(x[x_1], 2)))) / (std::sqrt(pow(x[x_0], 2) + pow(x[x_1], 2)) * std::sqrt(pow(R - std::sqrt(pow(x[x_0], 2) + pow(x[x_1], 2)), 2) + pow(x[x_2], 2))),
                -r_y + (x[x_1] * (-R + std::sqrt(pow(x[x_0], 2) + pow(x[x_1], 2)))) / (std::sqrt(pow(x[x_0], 2) + pow(x[x_1], 2)) * std::sqrt(pow(R - std::sqrt(pow(x[x_0], 2) + pow(x[x_1], 2)), 2) + pow(x[x_2], 2))),
                x[x_2] / std::sqrt(pow(R - std::sqrt(pow(x[x_0], 2) + pow(x[x_1], 2)), 2) + pow(x[x_2], 2))
            );
        }
    };

    std::unique_ptr<Surface> surfaceFactory(ParamCL const & params) {
        std::string funcName = __func__;
        auto name = params.get<std::string>("Surface.Name");
        if (name == "Sphere") return std::make_unique<Sphere>(params.get<double>("Surface.Params." + name + ".r_0"), params.get<double>("Surface.Params." + name + ".A"));
        if (name == "Torus") return std::make_unique<Torus>(params.get<double>("Surface.Params." + name + ".R"), params.get<double>("Surface.Params." + name + ".r_min"), params.get<double>("Surface.Params." + name + ".r_max"), params.get<size_t>("Surface.Params." + name + ".axis"));
        if (name == "OscillatingInextensibleSphere") return std::make_unique<OscillatingInextensibleSphere>(params.get<double>("Surface.Params." + name + ".r_0"), params.get<double>("Surface.Params." + name + ".eps"), params.get<double>("Surface.Params." + name + ".omega"), params.get<double>("Surface.Params." + name + ".eps_unsym"), params.get<double>("Surface.Params." + name + ".omega_unsym"));
        throw std::invalid_argument(funcName + ": unknown surface '" + name + "'");
    }

}

#endif // SURFACE_HPP