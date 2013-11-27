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

#include <vector>
#include <limits>

namespace DROPS
{

/*!
 * \brief The base class for all marking strategies.
 *
 * This class provides the interface for all marking strategies. On input, a
 * marking strategy's visit function expects a tetrahedron which is marked as
 * 'NoRefMark'. This is true when it is used in the AdapTriangCL. It can then
 * choose between three options: change the mark to RegRefMark, change the mark
 * to RemoveMark, or do not touch the mark.
 * 
 * If a mark was changed, the modified() function will return true afterwards.
 */
class MarkingStrategyCL: public TetraAccumulatorCL
{
public:
    virtual ~MarkingStrategyCL() {}
    virtual  void visit( const TetraCL& t ) = 0;
    virtual  TetraAccumulatorCL* clone (int clone_id) = 0;

    virtual  bool modified() const = 0;
    virtual  void SetUnmodified() = 0;

    virtual  Uint GetFineLevel() const = 0;
    virtual  Uint GetCoarseLevel() const = 0;
};


/*!
 * \brief Combine various marking strategies into one.
 *
 * AdapTriangCL always uses a single marking strategy. This class can be used
 * to combine arbitrary many marking strategies into a single one. It uses a
 * simple approach: for a given tetrahedron, mark it as RegRefMark if only one
 * strategy does so. Only mark it with RemoveMark if no stratetgy sets RegRefMark
 * and at least one strategy sets RemoveMark.
 */
class StrategyCombinerCL: public MarkingStrategyCL
{
public:
    void visit( const TetraCL& t );
    TetraAccumulatorCL* clone( int );

    bool modified() const;
    void SetUnmodified();

    Uint GetFineLevel() const;
    Uint GetCoarseLevel() const;

    void push( MarkingStrategyCL* s );
    MarkingStrategyCL* pop(); 

private:
    bool modified_;
    std::vector<MarkingStrategyCL*> strategies_;
};

inline
TetraAccumulatorCL* StrategyCombinerCL::clone( int )
{
    return new StrategyCombinerCL(*this);
}

inline
bool StrategyCombinerCL::modified() const
{
    return modified_;
}

inline
void StrategyCombinerCL::SetUnmodified()
{
    modified_ = false;
}

inline
Uint StrategyCombinerCL::GetFineLevel() const
{
    Uint result = 0;
    for ( Uint i = 0; i < strategies_.size(); ++i )
    {
        result = std::max( result, strategies_[ i ]->GetFineLevel() );
    }
    return result;
}

inline
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

inline
void StrategyCombinerCL::push( MarkingStrategyCL *s )
{
    strategies_.push_back( s );
}

inline
MarkingStrategyCL* StrategyCombinerCL::pop()
{
    if ( strategies_.size() == 0 ) return 0;

    MarkingStrategyCL *result = strategies_.back();
    strategies_.pop_back();
    return result;
}

}

#endif

