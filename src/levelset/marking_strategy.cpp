/// \file marking_strategy.cpp
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
#include "levelset/marking_strategy.h"

#include "levelset/levelset.h"

#include <limits>

namespace DROPS
{

TetraAccumulatorCL* UniformMarkingStrategyCL::clone( int )
{
    return new UniformMarkingStrategyCL( *this );
}

MarkingStrategyCL* UniformMarkingStrategyCL::clone_strategy()
{
    return new UniformMarkingStrategyCL( *this );
}

void UniformMarkingStrategyCL::visit( const TetraCL &t )
{
    // level should be level_.
    const Uint l= t.GetLevel();

    if ( l !=  level_ || (l == level_ && ! t.IsRegular()) )
    {
        // tetra will be marked for refinement/removement
        if ( l <= level_ )
        {
            t.SetRegRefMark();
            modified_ = true;
            decision_ = RefineC;
        }
        else
        {
            t.SetRemoveMark();
            modified_ = true;
            decision_ = CoarsenC;
        }
    }
    else if ( l == level_ )
    {
        // Tetrahedra has exactly the level that we require.
        decision_ = KeepC;
    }
}


StrategyCombinerCL::StrategyCombinerCL():
 modified_( false ), decision_( DontCareC ), strategies_( 0 )
{}


StrategyCombinerCL::StrategyCombinerCL( const StrategyCombinerCL& rhs ):
 modified_( rhs.modified_ ), decision_( rhs.decision_ ),
 strategies_( rhs.strategies_.size() )
{
    for ( Uint i = 0; i < strategies_.size(); ++i )
    {
        strategies_[ i ] = rhs.strategies_[ i ]->clone_strategy();
    }
}


StrategyCombinerCL::~StrategyCombinerCL()
{
    for ( Uint i = 0; i < strategies_.size(); ++i )
    {
        delete strategies_[ i ];
    }
}


StrategyCombinerCL& StrategyCombinerCL::operator=( const StrategyCombinerCL& rhs )
{
    if ( &rhs == this ) return *this;

    modified_ = rhs.modified_;
    decision_ = rhs.decision_;

    for ( Uint i = 0; i < strategies_.size(); ++i )
    {
        delete strategies_[ i ];
        strategies_[ i ] = 0; // Needed for exception safety.
    }

    strategies_.resize( rhs.strategies_.size() );
    for ( Uint i = 0; i < strategies_.size(); ++i )
    {
        strategies_[ i ] = rhs.strategies_[ i ]->clone_strategy();
    }

    return *this;
}


void StrategyCombinerCL::visit( const TetraCL& t )
{
    decision_ = DontCareC;

    for ( Uint i = 0; i < strategies_.size(); ++i )
    {
        // Every strategy expects Refmark to be on NoRef on input.
	t.SetNoRefMark();

        MarkingStrategyCL *s = strategies_[ i ];
        s->visit( t );

        switch ( s->GetDecision() )
        {
        case DontCareC:
            // Don't do anything.
            break;
        case KeepC:
            // Do not coarsen.
            if ( decision_ == DontCareC || decision_ == CoarsenC )
            {
                decision_ = KeepC;
            }
            break;
        case CoarsenC:
            if ( decision_ == DontCareC )
            {
                decision_ = CoarsenC;
            }
            break;
        case RefineC:
            // Refinement always wins. Return.
            decision_ = RefineC;
            modified_ = true;
            return;
        }
    }

    if ( decision_ == CoarsenC )
    {
        modified_ = true;
        t.SetRemoveMark();
    }
    else t.SetNoRefMark();
}

TetraAccumulatorCL* StrategyCombinerCL::clone( int )
{
    return new StrategyCombinerCL(*this);
}

MarkingStrategyCL* StrategyCombinerCL::clone_strategy()
{
    return new StrategyCombinerCL(*this);
}

bool StrategyCombinerCL::modified() const
{
    return modified_;
}


void StrategyCombinerCL::SetUnmodified()
{
    modified_ = false;
}


MarkingDecisionT StrategyCombinerCL::GetDecision() const
{
    return decision_;
}


Uint StrategyCombinerCL::GetFineLevel() const
{
    Uint result = 0;
    for ( Uint i = 0; i < strategies_.size(); ++i )
    {
        result = std::max( result, strategies_[ i ]->GetFineLevel() );
    }
    return result;
}


Uint StrategyCombinerCL::GetCoarseLevel() const
{
    if ( strategies_.size() == 0 ) return 0;

    Uint result = std::numeric_limits<Uint>::max();
    for ( Uint i = 0; i < strategies_.size(); ++i )
    {
        result = std::min( result, strategies_[ i ]->GetCoarseLevel() );
    }
    return result;
}


bool StrategyCombinerCL::empty()
{
    return strategies_.empty();
}


void StrategyCombinerCL::push_back( MarkingStrategyCL &s )
{
    strategies_.push_back( s.clone_strategy() );
}


void StrategyCombinerCL::pop_back()
{
    strategies_.pop_back();
}




///////////////////////////////////////////
LevelsetP2GetterCL::LevelsetP2GetterCL( const LevelsetP2CL& lset ): lset_( lset )
{}

double LevelsetP2GetterCL::GetValue( const VertexCL& v ) const
{
    return lset_.GetSolution().val(v);
}

double LevelsetP2GetterCL::GetValue( const EdgeCL& e ) const
{
    return lset_.GetSolution().val(e);
}

double LevelsetP2GetterCL::GetValue( const TetraCL& t ) const
{
    return lset_.GetSolution().val( t, 0.25, 0.25, 0.25 );
}

ValueGetterCL* LevelsetP2GetterCL::clone() const
{
    return new LevelsetP2GetterCL(*this);
}

FunPtrGetterCL::FunPtrGetterCL( instat_scalar_fun_ptr fct, double time ):
 time_( time ), fct_( fct )
{}

double FunPtrGetterCL::GetValue( const VertexCL& v ) const
{
    return fct_( v.GetCoord(), time_ );
}

double FunPtrGetterCL::GetValue( const EdgeCL& e ) const
{
    return fct_( GetBaryCenter(e), time_ );
}

double FunPtrGetterCL::GetValue( const TetraCL& t ) const
{
    return fct_( GetBaryCenter(t), time_ );
}

ValueGetterCL* FunPtrGetterCL::clone() const
{
    return new FunPtrGetterCL(*this);
}


///////////////////////////
// DistMarkingStrategyCL //
///////////////////////////

DistMarkingStrategyCL::DistMarkingStrategyCL( const LevelsetP2CL &dist,
                                              double width,
                                              Uint coarse_level, Uint fine_level ):
 getter_( new LevelsetP2GetterCL(dist) ), width_( width ), c_level_( coarse_level ),
 f_level_( fine_level ), modified_( false ), decision_( DontCareC )
{
    if ( f_level_ < c_level_ )
    {
        throw DROPSErrCL( "Refinement levels are cheesy.\n" );
    }
}

DistMarkingStrategyCL::DistMarkingStrategyCL( instat_scalar_fun_ptr fct,
                                              double width,
                                              Uint coarse_level, Uint fine_level,
                                              double time ):
 getter_( new FunPtrGetterCL(fct,time) ), width_( width ), c_level_( coarse_level ),
 f_level_( fine_level ), modified_( false ), decision_( DontCareC )
{
    if ( f_level_ < c_level_ )
    {
        throw DROPSErrCL( "Refinement levels are cheesy.\n" );
    }
}

DistMarkingStrategyCL::DistMarkingStrategyCL( const DistMarkingStrategyCL& rhs ):
 getter_( rhs.getter_->clone() ),  width_( rhs.width_ ),
 c_level_( rhs.c_level_ ), f_level_( rhs.f_level_ ),
 modified_( rhs.modified_ ), decision_( rhs.decision_ )
{}

DistMarkingStrategyCL::~DistMarkingStrategyCL()
{
    delete getter_;
}

DistMarkingStrategyCL& DistMarkingStrategyCL::operator=( const DistMarkingStrategyCL& rhs )
{
    if ( &rhs == this ) return *this;

    delete getter_;
    getter_ = 0; // Needed for exception safety.

    getter_ = rhs.getter_->clone();
    width_  = rhs.width_;
    c_level_ = rhs.c_level_;
    f_level_ = rhs.f_level_;
    modified_ = rhs.modified_;
    decision_ = rhs.decision_;

    return *this;
}

TetraAccumulatorCL* DistMarkingStrategyCL::clone( int )
{
    return new DistMarkingStrategyCL( *this );
}

MarkingStrategyCL* DistMarkingStrategyCL::clone_strategy()
{
    return new DistMarkingStrategyCL( *this );
}

void DistMarkingStrategyCL::visit( const TetraCL &t )
{
    double d = std::numeric_limits<double>::max();
    int num_pos = 0;
    for ( Uint j = 0; j < 4; ++j )
    {
        const double dist = getter_->GetValue( *t.GetVertex(j) );
        if ( dist >= 0 ) ++num_pos;
        d = std::min( d, std::abs(dist) );
    }

    for ( Uint j = 0; j < 6; ++j )
    {
        const double dist = getter_->GetValue( *t.GetEdge(j) );
        if ( dist >= 0 ) ++num_pos;
        d = std::min( d, std::abs(dist) );
    }

    d = std::min( d, std::abs(getter_->GetValue(t)) );

    const bool vzw = num_pos!=0 && num_pos!=10; // change of sign
    const Uint l = t.GetLevel();

    // In the shell:      level should be f_level_.
    // Outside the shell: level should be c_level_.
    const Uint soll_level = ( d <= width_ || vzw ) ? f_level_ : c_level_;

    if ( l !=  soll_level || (l == soll_level && ! t.IsRegular()) )
    {
        // tetra will be marked for refinement/removement
        if ( l <= soll_level )
        {
            t.SetRegRefMark();
            modified_ = true;
            decision_ = RefineC;
        }
        else
        {
            t.SetRemoveMark();
            modified_ = true;
            decision_ = CoarsenC;
        }
    }
    else if ( l == soll_level )
    {
        // Tetrahedra has exactly the level that we require.
        decision_ = KeepC;
    }
}

bool DistMarkingStrategyCL::modified() const
{
    return modified_;
}

void DistMarkingStrategyCL::SetUnmodified()
{
    modified_ = false;
}

MarkingDecisionT DistMarkingStrategyCL::GetDecision() const
{
    return decision_;
}

void DistMarkingStrategyCL::SetDistFct( instat_scalar_fun_ptr fct, double time )
{
    delete getter_;
    getter_ = 0; // Needed for exception safety.

    getter_ = new FunPtrGetterCL( fct, time );
}

void DistMarkingStrategyCL::SetDistFct( const LevelsetP2CL &fct )
{
    delete getter_;
    getter_ = 0; // Needed for exception safety.

    getter_ = new LevelsetP2GetterCL( fct );
}

double DistMarkingStrategyCL::GetWidth() const
{
    return width_;
}

void DistMarkingStrategyCL::SetWidth( double width )
{
    width_ = width;
}


Uint DistMarkingStrategyCL::GetCoarseLevel() const
{
    return c_level_;
}

void DistMarkingStrategyCL::SetCoarseLevel( Uint level )
{
    c_level_ = level;
    if ( f_level_ < c_level_ )
    {
        throw DROPSErrCL( "Refinementlevels are cheesy.\n" );
    }
}

Uint DistMarkingStrategyCL::GetFineLevel() const
{
    return f_level_;
}

void DistMarkingStrategyCL::SetFineLevel( Uint level )
{
    f_level_ = level;
    if ( f_level_ < c_level_ )
    {
        throw DROPSErrCL( "Refinementlevels are cheesy.\n" );
    }
}

}

