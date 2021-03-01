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
        struct radius_ { double min, max; };
    public:
        Surface(bool isStationary) : isStationary_(isStationary), description_(isStationary_ ? "stationary surface\n" : "evolving surface\n") {}
        bool isStationary() const { return isStationary_; }
        std::string description() const { return description_; }
        virtual std::array<radius_, 3> radius(double t = 0.) const {
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
        virtual double dist(Point3DCL const &, double t = 0.) const {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": distance function not available");
        }
        virtual double phi(Point3DCL const & x, double t = 0.) const { return dist(x, t); }
        virtual Point3DCL normal(Point3DCL const &, double t = 0.) const = 0;
        SMatrixCL<3, 3> P(Point3DCL const & x, double t = 0.) const { return eye<3, 3>() - outer_product(normal(x, t), normal(x, t)); }
        Point3DCL ext(Point3DCL const & x, double t = 0.) const { return x - dist(x, t) * normal(x, t); }
    };

    struct Sphere : Surface {
        const double r_0, A;
        double r(double t) const { return r_0 * (1. + A * sin(2. * M_PI * t)); }
        double r_prime(double t) const { return 2. * A * M_PI * r_0 * cos(2. * M_PI * t); }
        explicit Sphere(double r_0 = 1., double A = 0.) : Surface(!A), r_0(r_0), A(A) {
            std::string funcName = __func__;
            if (r_0 <= 0.) throw std::invalid_argument(funcName + ": radius must be positive");
            if (A < 0. || A >= 1.) throw std::invalid_argument(funcName + ": A must be in [0, 1)");
            description_ += isStationary_ ? "phi = x^2 + y^2 + z^2 - R, R = " + std::to_string(r_0) :
                "phi = x^2 + y^2 + z^2 - r^2(t), r(t) := r0 (1 + A sin(2 pi t))\n"
                "A = " + std::to_string(A) + ", r_0 = " + std::to_string(r_0);
        }
        virtual ~Sphere() {}
        std::array<radius_, 3> radius(double t) const final { return {0., r(t), 0., r(t), 0., r(t) }; }
        std::array<bool, 3> rotationalInvariance(double) const final { return { true, true, true }; }
        double u_N(Point3DCL const &, double t) const final { return r_prime(t); }
        double dist(Point3DCL const & x, double t) const final { return norm(x) - r(t); }
        double phi(Point3DCL const & x, double t) const final { return norm_sq(x) - r(t) * r(t); } // global P2 func
        Point3DCL normal(Point3DCL const & x, double) const final { return x / norm(x); }
    };

    struct Torus : Surface {
        double r(double x, double y) const { return r_min + .5 * (r_max - r_min) * (1. - x / std::sqrt(x * x + y * y)); }
        const double R, r_min, r_max;
        explicit Torus(double R = 1., double r_min = .5, double r_max = .5) : Surface(true), R(R), r_min(r_min), r_max(r_max) {
            std::string funcName = __func__;
            if (r_min <= 0. || r_max < r_min || R <= r_max) throw std::invalid_argument(funcName + ": invalid torus params");
            description_ +=
                "torus w/ const distance R from the center of the tube to the origin and\n"
                "w/ variable tube rad " + std::to_string(r_min) + " =: r_min <= r(xi) <= r_max =: " + std::to_string(r_max) + ", xi in [0, 2 pi]\n"
                "phi = (x^2 + y^2 + z^2 + R^2 - r(xi)^2)^2 - 4 R^2 (x^2 + y^2)";
        }
        virtual ~Torus() {}
        std::array<radius_, 3> radius(double) const final { return {R - r_max, R + r_max, R - .5 * (r_max + r_min), R + .5 * (r_max + r_min), 0., r_max }; }
        std::array<bool, 3> rotationalInvariance(double) const final { return {false, false, r_min == r_max }; }
        double dist(Point3DCL const & x, double) const final { return std::sqrt(pow(std::sqrt(x[0] * x[0] + x[1] * x[1]) - R, 2.) + x[2] * x[2]) - r(x[0], x[1]); }
        double phi(Point3DCL const & x, double) const final { return pow(norm_sq(x) + R * R - r(x[0], x[1]) * r(x[0], x[1]), 2.) - 4. * R * R * (pow(x[0], 2.) + pow(x[1], 2.)); }
        Point3DCL normal(Point3DCL const & x, double) const final {
            return Point3DCL(
                    .5 * (r_max - r_min) * x[1] * x[1] / pow(pow(x[0], 2) + pow(x[1], 2), 1.5) + (x[0] * (std::sqrt(pow(x[0], 2) + pow(x[1], 2)) - R)) / (std::sqrt(pow(x[0], 2) + pow(x[1], 2)) * std::sqrt(pow(R - std::sqrt(pow(x[0], 2) + pow(x[1], 2)), 2) + pow(x[2], 2))),
                    .5 * (r_min - r_max) * x[0] * x[1] / pow(pow(x[0], 2) + pow(x[1], 2), 1.5) + (x[1] * (std::sqrt(pow(x[0], 2) + pow(x[1], 2)) - R)) / (std::sqrt(pow(x[0], 2) + pow(x[1], 2)) * std::sqrt(pow(R - std::sqrt(pow(x[0], 2) + pow(x[1], 2)), 2) + pow(x[2], 2))),
                x[2] / std::sqrt(pow(R - std::sqrt(pow(x[0], 2) + pow(x[1], 2)), 2) + pow(x[2], 2))
            );
        }
    };

    std::unique_ptr<Surface> surfaceFactory(ParamCL const & params) {
        std::string funcName = __func__;
        auto name = params.get<std::string>("Surface.Name");
        if (name == "Sphere") return std::make_unique<Sphere>(params.get<double>("Surface.Params." + name + ".r_0"), params.get<double>("Surface.Params." + name + ".A"));
        else if (name == "Torus") return std::make_unique<Torus>(params.get<double>("Surface.Params." + name + ".R"), params.get<double>("Surface.Params." + name + ".r_min"), params.get<double>("Surface.Params." + name + ".r_max"));
        throw std::invalid_argument(funcName + ": unknown surface '" + name + "'");
    }

}

#endif // SURFACE_HPP