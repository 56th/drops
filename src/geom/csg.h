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
#include <string>
#include "misc/params.h"
#include "misc/singletonmap.h"

namespace DROPS {
namespace CSG {

/// A CSG-object is a directed acyclic graph (DAG) of BodyCL-objects (bodies):
///    * Leave-nodes (do not refer to other nodes):
///         * Halfspace
///         * InfiniteCone
///         * InfiniteCylinder
///         * Levelset
///         * Sphere
///         * Torus

///    * inner nodes: Operators (refer to at least one node):
///         * Complement (Body b)
///         * Intersection (Body b0, Body b1)
///         * Union (Body b0, Body b1)
///         * ApplySimilarityToDomain( Body b)
///         * Reference (Body b)
///         * Module -- refers to the Body created by loading the module.

/// All bodies support
///     virtual double operator() (const Point3DCL& x, double t) const
/// to evaluate the level set function at (x,t).

/// From am mathematical perspective, it might seem odd to use the same class for
/// the operators and the geometric primitives. However, it is convenient for
/// implementation in C++, if the data-structure is homogeneous.

/// All inner nodes take ownership of the bodies they refer to. This means, that
/// an inner node will delete the bodies it refers to, when it is destroyed itself.
///
/// The only exception is the Reference-body. It does not destroy the body it
/// refers to.
/// References allow for smaller CSG-objects, if many instances of the
/// same CSG-sub-object are used in a CSG-object. They can instead all refer
/// to the same one.

/// The sequential description of the CSG-object uses a stack-machine with a stack s.
/// All bodies are passed via s.
/// Inner nodes pop the required bodies from s. New bodies are pushed to s.
///    * Start with empty stack s.
///    * Consume the sequence of instructions.
///    * After the last instruction, the top of S is popped and returned as
///      the desired CSG-object ret.
/// It is an error, if s is not empty after popping ret.

/// Every instruction is mapped to the constructor of a body.
/// Instructions are classified, on whether they accept additional arguments, which
/// are not passed on the stack.
///    * SimpleInstructions only pass the current stack and the name, with which
///      they were invoked to the body-constructor.
///    * ComplexInstructions pass the current stack a ParamCL-reference to the body-constructor.

/// The ComplexInstruction "CreateReference" is not mapped to a body-constructor.
/// It only memoizes a name for s.top().

/// Without references, one can describe arbitrary trees with a stack machine (postorder traversal).
/// The Reference-operator allows the description of arbitrary directed graphs.
/// *Do not* create cycles -- or you will get what you deserve.

/// All bodies have the function known_keys (), although it is only used for
/// ComplexInstructions.
/// The body_builder uses the returned vector of strings to verify, that only
/// known keys were supplied to the instruction. This guards against misspellings.

/// Description of a CSG-object in json:
///    * A json-object with the following entries:
///             * key "Type", value "CSG-body (v0)", required
///             * key "Instructions" [json-array of InstructionType],
///               must be non-empty, required
///             * key "Name" [string], optional
///    * InstructionType:
///         * SimpleInstruction": a string
///         * ComplexInstruction": A json-object with the following entries:
///             * key "Type" [string], required
///             * optional further key-value pairs.

/// =The following SimpleInstructions exist=
/// ==Complement: set-theoretic complement==
/// -ls;
/// Pops one arg from s, pushes one result.
/// ComplementBodyCL

/// ==Intersection: set-theoretic intersection==
/// max(ls0, ls1);
/// Pops two args from s, pushes one result.
/// IntersectionBodyCL

/// ==Union: set-theoretic union==
/// min(ls0, ls1);
/// Pops two args from s, pushes one result.
/// UnionBodyCL

/// ==SmoothIntersection: set-theoretic intersection==
/// smooth, distance-like R-function instead of max; preserves the level-set
/// exactly; otherwise as Intersection.
/// IntersectionBodyCL

/// ==SmoothUnion: set-theoretic union==
/// smooth, distance-like R-function instead of min; preserves the level-set
/// exactly; otherwise as Union.
/// UnionBodyCL


/// =The following ComplexInstructions exist=
/// ==HalfSpace: Linear level set function for a plane==
///    * key "Normal" [Point3DCL], required;
///    * key "ValueAtOrigin" [double], optional, defaults to 0.0;
///    * key "MakeUnitNormal" [bool]; optional, defaults to false;
/// Pushes one body on the stack;
/// HalfspaceBodyCL

/// ==Sphere: Distance function for a sphere==
///    * key "Center" [Point3DCL]; optional, defaults to (0,0,0)^T;
///    * key "Radius", [double]; optional, defaults to 0.0;
/// Pushes one body on the stack;
/// SphereBodyCL

/// ==InfiniteCone: Level set function function for an infinite cone==
///    * key "Apex" [Point3DCL]; optional, defaults to (0,0,0)^T;
///    * key "Axis", [Point3DCL]; optional, defaults to (1,0,0)^T;
///    * key "SemiAperture", [double]; optional, defaults to M_PI/4.;
/// The semi-aperture is the angle between the axis and the straight lines in the
/// hull of the cone. The Axis points from the Apex into the cone.
/// Pushes one body on the stack;
/// InfiniteConeBodyCL

/// ==InfiniteCylinder: Level set function function for an infinite cylinder==
///    * key "Origin" [Point3DCL]; optional, defaults to (0,0,0)^T;
///    * key "Axis", [Point3DCL]; optional, defaults to (1,0,0)^T;
///    * key "Radius", [double]; optional, defaults to 1.;
/// The Origin defines a point on the Axis.
/// Pushes one body on the stack;
/// InfiniteCylinderBodyCL

/// ==Torus: Level set function function for a torus==
///    * key "Origin" [Point3DCL]; optional, defaults to (0,0,0)^T;
///    * key "Axis", [Point3DCL]; optional, defaults to (1,0,0)^T;
///    * key "BigRadius", [double]; optional, defaults to 0.75;
///    * key "SmallRadius", [double]; optional, defaults to 0.25;
/// The Origin defines a point on the Axis.
/// Pushes one body on the stack;
/// TorusBodyCL

/// ==Levelset: Use level set function from InScaMap==
///    * key "Function", value of type string, required;
/// Pushes one body on the stack;
/// LevelsetBodyCL

/// ==CreateReference: Create and memoize a name for current s.top()==
///    * key "Name" [string], required;
/// Does not alter the stack.
/// Implemented in body_builder.

/// ==PushReference: Push a reference to the named body on the stack.
///    * key "Name" [string], required;
/// The value of Name must have been created by CreateReference.
/// Pushes one ReferenceBodyCL on the stack.
/// ReferenceBodyCL

/// ==LoadFromModule: create a CSG-object from a file by calling body_builder==
///    * key "Path" [string], required; path to a json-file.
///    * key "Name" [string], optional
/// If Name is not present the file at Path must contain a valid CSG-object.
/// If Name is present, the json file must have a child of this name, which
/// is a valid CSG-object.
/// Pushes one body on the stack;
/// ModuleBodyCL

/// ==ApplySimilarityToDomain: Apply a similarity-transformation==
/// y= s*r*x + t is applied to the domain and  v= u/s to the image of the level set function.
///    * key "Scaling" [double], optional, defaults to s=1.0;
///    * key "Translation" [Point3DCL], optional, defaults to t=(0, 0, 0)^T;
///    * key "RotationAxis", [Point3DCL]; optional, the axis of the rotation;
///    * key "RotationAngle" [double]; optional;
/// The rotation-parameters must either both be present or none of them. In
/// the latter case, no rotation is performed.
/// Pops one arg from the stack, pushes one body on the stack.
/// SimilarityTransformBodyCL

class BodyCL; ///< Base-type of all CSG-Objects: Evaluation, memory-management.

/// \brief Type to store references to BodyCL for memory menagement.
typedef std::vector<const BodyCL*> BodyVectorT;

/// \brief Stack-type for the stack-machine.
typedef std::stack<const BodyCL*>  BodyStackT;

/// \brief Type to store known keys of complex instructions.
typedef std::vector<std::string> VecOfStringT;

class BodyCL
{
  private:
    mutable BodyVectorT owned_ptrs_; ///< All bodies in this container are destroyed in the destructor.

  protected:
    ///\brief Register s.top() for deletion, pop it from s, and return it.
    const BodyCL* pop_and_take_ownership (BodyStackT& s) const {
        if (s.empty())
            throw DROPSErrCL("pop_and_take_ownership: Cannot pop an empty stack.\n");
        const BodyCL* ref= s.top();
        s.pop();
        owned_ptrs_.push_back( ref);
        return ref;
    }

  public:
    virtual ~BodyCL () {
        for (BodyVectorT::iterator it= owned_ptrs_.begin(); it != owned_ptrs_.end(); ++it)
            delete *it;
    }

    ///\brief Return a vector of known keys in ParamCL.
    virtual const VecOfStringT& known_keys () const {
        static const VecOfStringT k;
        return k;
    }

    ///\brief Evaluation of the level set function.
    virtual double operator() (const Point3DCL& x, double t) const= 0;
};


/// \brief Reads the json object in p an returns the corresponding body.
/// The returned body recursively owns all referenced subobjects.
/// Destroy with delete.
const BodyCL* body_builder (const ParamCL& p); // :-)

} // end of namespace DROPS::CSG
} // end of namespace DROPS

#endif
