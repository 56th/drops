/// \file marking_strategy.h
/// \brief
/// \author LNM RWTH Aachen: Matthias Kirchhart

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
 * Copyright 2011 LNM/SC RWTH Aachen, Germany
*/
#ifndef MARKING_STRATEGY_H
#define MARKING_STRATEGY_H

#include "num/accumulator.h"
#include "num/discretize.h" // For instat_scalar_fun_ptr.

#include <vector>

namespace DROPS
{


enum MarkingDecisionT { CoarsenC, RefineC, KeepC, DontCareC };

/*!
 * \brief Base class for all marking strategies.
 *
 * This class represents the interface that is common to all marking strategies.
 * A marking strategy is a TetraAccumulatorCL that performs the following tasks:
 * <ol>
 * <li>Make a decision what to do about a given tetrahedron. This decision is
 *     one out of:
 * <ul>
 * <li><b>CoarsenC:</b> the tetrahedron should be "coarsened".
 * <li><b>RefineC:</b> the tetrahedron should be refined.
 * <li><b>KeepC:</b> the tetrahedron should neither be coarsened nor refined.
 * <li><b>DontCareC:</b> it doesn't matter what happens to the tetrahedron.
 * </ul>
 * </li>
 * <li>The decision is stored and can be retrieved using the "GetDecision" member.</li>
 * <li>Set the mark of the tetrahedron accordingly, if necessary.</li>
 * <li>If the mark of the tetrahedron was changed, set the status to modified.</li>
 * </ol>
 *
 * On input, it is expected that the mark of the tetrahedron is set to
 * NoRefMark.
 */
class MarkingStrategyCL: public TetraAccumulatorCL
{
public:
    virtual ~MarkingStrategyCL() {}
    virtual  void visit( const TetraCL& t ) = 0;
    virtual  TetraAccumulatorCL* clone (int clone_id) = 0;
    virtual  MarkingStrategyCL*  clone_strategy() = 0;

    virtual  bool modified() const = 0;
    virtual  void SetUnmodified() = 0;

    virtual  MarkingDecisionT GetDecision() const = 0;

    virtual  Uint GetFineLevel() const = 0;
    virtual  Uint GetCoarseLevel() const = 0;
};


/*!
 * \brief Refinement strategy for uniform refinement.
 *
 * This marking strategy marks all tetrahedra for refinement until a prescribed level is reached.
 */
class UniformMarkingStrategyCL: public MarkingStrategyCL
{
public:
    UniformMarkingStrategyCL( Uint Level) : level_(Level), decision_(DontCareC), modified_(false) {}
    ~UniformMarkingStrategyCL() {}
    void visit( const TetraCL& t );
    TetraAccumulatorCL* clone (int clone_id);
    MarkingStrategyCL*  clone_strategy();

    bool modified() const { return modified_; }
    void SetUnmodified()  { modified_= false; }

    MarkingDecisionT GetDecision() const { return decision_; }

    Uint GetFineLevel()   const { return level_; }
    Uint GetCoarseLevel() const { return 0; }

private:
    Uint level_;
    MarkingDecisionT decision_;
    bool modified_;
};




/*!
 * \brief Combines several MarkingStrategyCL into one.
 *
 * This class can be used for combining various marking strategies into a
 * single one. In order to do so, it stores *copies* of given marking strategies.
 * Based on the decisions returned by GetDecision() of these strategies, it then
 * decides what to do to a tetrahedron. It implements a "Refinement always wins"
 * approach.
 */
class StrategyCombinerCL: public MarkingStrategyCL
{
public:
    StrategyCombinerCL();
    StrategyCombinerCL( const StrategyCombinerCL &rhs );
    ~StrategyCombinerCL();

    StrategyCombinerCL& operator=( const StrategyCombinerCL& rhs );

    void visit( const TetraCL& t );
    TetraAccumulatorCL* clone( int );
    MarkingStrategyCL*  clone_strategy();

    bool modified() const;
    void SetUnmodified();

    MarkingDecisionT GetDecision() const;

    Uint GetFineLevel() const;
    Uint GetCoarseLevel() const;

    bool empty();
    void push_back( MarkingStrategyCL &s );
    void pop_back();

private:
    bool modified_;
    MarkingDecisionT decision_;
    std::vector<MarkingStrategyCL*> strategies_;
};




class LevelsetP2CL;
class EdgeCL;
class VertexCL;
class ValueGetterCL;

/*!
 * \brief Refinement strategy according to a levelset function.
 *
 * This marking strategy refines the mesh around the zero level of a given
 * levelset function.
 */
class DistMarkingStrategyCL: public MarkingStrategyCL
{
public:
    DistMarkingStrategyCL( instat_scalar_fun_ptr fct,
                           double width, Uint coarse_level, Uint fine_level,
                           double time = 0. );
    DistMarkingStrategyCL( const LevelsetP2CL &fct,
                           double width, Uint coarse_level, Uint fine_level );
    DistMarkingStrategyCL( const DistMarkingStrategyCL &rhs );
    ~DistMarkingStrategyCL();
    DistMarkingStrategyCL& operator=( const DistMarkingStrategyCL& rhs );
    TetraAccumulatorCL* clone( int );
    MarkingStrategyCL*  clone_strategy();

    void visit( const TetraCL& t );

    bool modified() const;
    void SetUnmodified();

    MarkingDecisionT GetDecision() const;

    void SetDistFct( instat_scalar_fun_ptr fct, double time = 0 );
    void SetDistFct( const LevelsetP2CL& fct );

    double GetWidth() const;
    void   SetWidth( double width );
    Uint   GetCoarseLevel() const;
    void   SetCoarseLevel( Uint level );
    Uint   GetFineLevel() const;
    void   SetFineLevel( Uint level );

private:
    ValueGetterCL *getter_;
    double width_;
    Uint c_level_;
    Uint f_level_;
    bool modified_;
    MarkingDecisionT decision_;
};


///////////////////////////////////////////
// ValueGetter for DistMarkingStrategyCL //
///////////////////////////////////////////

class ValueGetterCL
{
public:
    virtual ~ValueGetterCL() {}
    virtual double GetValue( const VertexCL& v ) const = 0;
    virtual double GetValue( const EdgeCL&   e ) const = 0;
    virtual double GetValue( const TetraCL&  t ) const = 0;
    virtual ValueGetterCL* clone() const = 0;
};

class LevelsetP2GetterCL: public ValueGetterCL
{
public:
    LevelsetP2GetterCL( const LevelsetP2CL& lset );

    double GetValue( const VertexCL& v ) const;
    double GetValue( const EdgeCL&   e ) const;
    double GetValue( const TetraCL&  t ) const;
    ValueGetterCL* clone() const;

private:
    const LevelsetP2CL& lset_;
};

class FunPtrGetterCL: public ValueGetterCL
{
public:
    FunPtrGetterCL( instat_scalar_fun_ptr fct, double time = 0. );

    double GetValue( const VertexCL& v ) const;
    double GetValue( const EdgeCL&   e ) const;
    double GetValue( const TetraCL&  t ) const;
    ValueGetterCL* clone() const;

private:
    double time_;
    instat_scalar_fun_ptr fct_;
};

}

#endif

