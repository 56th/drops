#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <functional>
#include "misc/container.h"

namespace DROPS {
    using ScalarFunction = std::function<double(Point3DCL const &)>;
    using VectorFunction = std::function<Point3DCL(Point3DCL const &)>;
    using MatrixFunction = std::function<SMatrixCL<3, 3>(Point3DCL const &)>;
    using InstatScalarFunction = std::function<double(Point3DCL const &, double)>;
    using InstatVectorFunction = std::function<Point3DCL(Point3DCL const &, double)>;
    using InstatMatrixFunction = std::function<SMatrixCL<3, 3>(Point3DCL const &, double)>;
    using MatchFunction = std::function<bool(const Point3DCL&, const Point3DCL&)>;
}

#endif // FUNCTIONS_HPP