/// \file csg.h
/// \brief classes for simple constructive solid geometry (CSG) for level sets
/// \author LNM RWTH Aachen: Joerg Grande

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_CSG_H
#define DROPS_CSG_H

#include <vector>
#include <stack>
#include "misc/params.h"
#include "misc/singletonmap.h"

namespace DROPS {
namespace CSG {

/// A CSG-object is a directed acyclic graph of bodies (BodyCL-objects):
///    * Leave-nodes:
///         * Halfspace
///         * Sphere
///         * Levelset
///    * inner nodes: Operators:
///         * Complement (Body b)
///         * Intersection (Body b0, Body b1)
///         * Union (Body b0, Body b1)

/// Description of the CSG-object in json-file as stack-machine:
/// Start with empty stack S.
/// Consume a sequence of instructions:
///    * Leave-instructions create a new object and push a reference to it on S.
///    * Operator-instructions pop the required number of pointers from S, create the new operator-object and push a reference to it on S.
/// After the last instruction, the top of S is popped and returned as final body. It is an error if S is not empty after that.

/// Description of a CSG-object in json:
///    * As a json object; keys:
///             * "Type": "CSG-body (v0)", required
///             * "Instructions": array of InstructionType, required
///             * "Name": string, optional
///        * InstructionType:
///             * "Complement"
///             * "Intersection"
///             * "Union"
///             * Json-object with key: "Type": string, required; optional further key-value pairs.

///{ "Type": "CreateReference", "Name": "XXX" } /// Memoizes s.top() as XXX.
///{ "Type": "PushReference", "Name": "XXX" }   /// Push the reference to XXX on s.
/// *Do not* create cycles or you get, what you deserve.
/// Use these sparingly. They are appropriate for very complex objects to reuse already complex subobjects. E.g., name the result of a module-load.
/// { "Type": "LoadFromModule", "Path": "path/to/json" } /// calls body_builder on path/to/json; the result is pushed to s; s.top owns the so-constructed CSG-object.
/// { "Type": "LoadFromModule", "Path": "path/to/json", "Name": "XXX" } /// calls body_builder on the child with name XXX in path/to/json; the result is pushed to s; s.top owns the so-constructed CSG-object.
/// { "Type": "ApplySimilarityToDomain", "Translation": [], "RotationAxis": [], "RotationAngle": a, "Scaling" s } all parts are optional, but the Rotation-Stuff must be present or absent as a whole.

class BodyCL; ///< Base-type of all CSG-Objects: Evaluation, memory-management.

/// \brief Type to store references to BodyCL for memory menagement.
typedef std::vector<const BodyCL*> BodyVectorT;

class BodyCL
{
  private:
    mutable BodyVectorT owned_ptrs_; ///< All bodies in this container are destroyed in the destructor.

  public:
    virtual ~BodyCL () {
        for (BodyVectorT::iterator it= owned_ptrs_.begin(); it != owned_ptrs_.end(); ++it)
            delete *it;
    }
    ///\brief Transfer ownership of the bodies in refs to this object.
    void make_owner_of (const BodyVectorT& refs) const { owned_ptrs_= refs; }

    ///\brief Evaluation of the level set function.
    virtual double operator() (const Point3DCL& x, double t) const= 0;
};


/// \brief Reads the json object in p an returns the corresponding body.
/// The returned body owns all referenced subobjects.
const BodyCL* body_builder (const ParamCL& p); // :-)

} // end of namespace DROPS::CSG
} // end of namespace DROPS

#endif
