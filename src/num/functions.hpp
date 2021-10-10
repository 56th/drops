#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <functional>
#include "misc/container.h"

namespace DROPS {
    using ScalarFunction = std::function<double(Point3DCL const &)>;
    using VectorFunction = std::function<Point3DCL(Point3DCL const &)>;
    using MatrixFunction = std::function<SMatrixCL<3, 3>(Point3DCL const &)>;
    template <typename T> using InstatFunction = std::function<T(Point3DCL const &, double)>;
    using InstatScalarFunction = InstatFunction<double>;
    using InstatVectorFunction = InstatFunction<Point3DCL>;
    using InstatMatrixFunction = InstatFunction<SMatrixCL<3, 3>>;
    using MatchFunction = std::function<bool(const Point3DCL&, const Point3DCL&)>;
    inline double zeroInstatScalarFunction(Point3DCL const &, double) { return 0.; }
    inline Point3DCL zeroInstatVectorFunction(Point3DCL const &, double) { return Point3DCL(0., 0., 0.); }
    struct InstatDiffFunction { std::string description; InstatScalarFunction f; InstatVectorFunction grad; };
}

#endif // FUNCTIONS_HPP