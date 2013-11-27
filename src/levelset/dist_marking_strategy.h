/// \file dist_marking_strategy.h
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
#ifndef DIST_MARKING_STRATEGY_H
#define DIST_MARKING_STRATEGY_H

#include "levelset/marking_strategy.h"

#include "levelset/levelset.h"
#include <limits>

namespace DROPS
{

/*!
 * \brief Marking strategy that refines around the interface.
 *
 * This Marking strategy marks the functions around the zero-level of a given
 * level-set function for refinement and the remaining ones for coarsening.
 */
template <class DistFctT>
class DistMarkingStrategyCL: public MarkingStrategyCL
{
public:
    DistMarkingStrategyCL( DistFctT *dist, double width,
                           Uint coarse_level, Uint fine_level ):
    dist_( dist ), width_( width ), c_level_( coarse_level ),
    f_level_( fine_level ) {}

    void visit( const TetraCL& t );
    TetraAccumulatorCL* clone( int );

    bool modified() const;
    void SetUnmodified();

    DistFctT* GetDistFct() const;
    void      SetDistFct( DistFctT* fct );
    double GetWidth() const;
    void   SetWidth( double width );
    Uint   GetCoarseLevel() const;
    void   SetCoarseLevel( Uint level );
    Uint   GetFineLevel() const;
    void   SetFineLevel( Uint level ); 

private:
    double GetValue( const VertexCL& v ) const;
    double GetValue( const EdgeCL&   e ) const;
    double GetValue( const TetraCL&  t ) const;

private:
    DistFctT* dist_;
    double width_;
    Uint   c_level_;
    Uint   f_level_;

    bool modified_;

};

template <class DistFctT> inline
TetraAccumulatorCL* DistMarkingStrategyCL<DistFctT>::clone( int )
{
    return new DistMarkingStrategyCL<DistFctT>( *this );
}

template <class DistFctT> inline
bool DistMarkingStrategyCL<DistFctT>::modified() const
{
    return modified_;
}

template <class DistFctT> inline
void DistMarkingStrategyCL<DistFctT>::SetUnmodified()
{
    modified_ = false;
}

template <class DistFctT> inline
DistFctT* DistMarkingStrategyCL<DistFctT>::GetDistFct() const
{
    return dist_;
}

template <class DistFctT> inline
void DistMarkingStrategyCL<DistFctT>::SetDistFct( DistFctT* fct )
{
    dist_ = fct;
}

template <class DistFctT> inline
double DistMarkingStrategyCL<DistFctT>::GetWidth() const
{
    return width_;
}

template <class DistFctT> inline
void DistMarkingStrategyCL<DistFctT>::SetWidth( double width )
{
    width_ = width;
}

template <class DistFctT> inline
Uint DistMarkingStrategyCL<DistFctT>::GetCoarseLevel() const
{
    return c_level_;
}

template <class DistFctT> inline
void DistMarkingStrategyCL<DistFctT>::SetCoarseLevel( Uint level )
{
    c_level_ = level;
}

template <class DistFctT> inline
Uint DistMarkingStrategyCL<DistFctT>::GetFineLevel() const
{
    return f_level_;
}

template <class DistFctT> inline
void DistMarkingStrategyCL<DistFctT>::SetFineLevel( Uint level )
{
    f_level_ = level;
}

template <class DistFctT> inline
double DistMarkingStrategyCL<DistFctT>::GetValue( const VertexCL& v ) const
{
    return dist_->GetSolution().val(v);
}

template <class DistFctT> inline
double DistMarkingStrategyCL<DistFctT>::GetValue( const EdgeCL& e ) const
{
    return dist_->GetSolution().val(e);
}

template <class DistFctT> inline
double DistMarkingStrategyCL<DistFctT>::GetValue( const TetraCL& t ) const
{
    return dist_->GetSolution().val(t,0.25,0.25,0.25);
}



/*!
 * \brief A specialisation for function pointers.
 */
template <>
class DistMarkingStrategyCL<instat_scalar_fun_ptr>: public MarkingStrategyCL
{
public:
    DistMarkingStrategyCL( instat_scalar_fun_ptr dist, double width,
                           Uint coarse_level, Uint fine_level,
                           double time = 0. ):
    dist_( dist ), width_( width ), c_level_( coarse_level ),
    f_level_( fine_level ), time_( time ) {}

    void visit( const TetraCL& t );
    TetraAccumulatorCL* clone( int );

    bool modified() const;
    void SetUnmodified();

    instat_scalar_fun_ptr GetDistFct() const;
    void   SetDistFct( instat_scalar_fun_ptr fct );
    double GetWidth() const;
    void   SetWidth( double width );
    Uint   GetCoarseLevel() const;
    void   SetCoarseLevel( Uint level );
    Uint   GetFineLevel() const;
    void   SetFineLevel( Uint level ); 
    double GetTime() const;
    void   SetTime( double time );

private:
    double GetValue( const VertexCL& v ) const;
    double GetValue( const EdgeCL&   e ) const;
    double GetValue( const TetraCL&  t ) const;

private:
    instat_scalar_fun_ptr dist_;
    double width_;
    Uint   c_level_;
    Uint   f_level_;
    double time_;

    bool modified_;
};


inline
TetraAccumulatorCL* DistMarkingStrategyCL<instat_scalar_fun_ptr>::clone( int )
{
    return new DistMarkingStrategyCL<instat_scalar_fun_ptr>( *this );
}

inline
bool DistMarkingStrategyCL<instat_scalar_fun_ptr>::modified() const
{
    return modified_;
}

inline
void DistMarkingStrategyCL<instat_scalar_fun_ptr>::SetUnmodified()
{
    modified_ = false;
}

inline
instat_scalar_fun_ptr DistMarkingStrategyCL<instat_scalar_fun_ptr>::GetDistFct() const
{
    return dist_;
}

inline
void DistMarkingStrategyCL<instat_scalar_fun_ptr>::SetDistFct( instat_scalar_fun_ptr fct )
{
    dist_ = fct;
}

inline
double DistMarkingStrategyCL<instat_scalar_fun_ptr>::GetWidth() const
{
    return width_;
}

inline
void DistMarkingStrategyCL<instat_scalar_fun_ptr>::SetWidth( double width )
{
    width_ = width;
}

inline
Uint DistMarkingStrategyCL<instat_scalar_fun_ptr>::GetCoarseLevel() const
{
    return c_level_;
}

inline
void DistMarkingStrategyCL<instat_scalar_fun_ptr>::SetCoarseLevel( Uint level )
{
    c_level_ = level;
}

inline
Uint DistMarkingStrategyCL<instat_scalar_fun_ptr>::GetFineLevel() const
{
    return f_level_;
}

inline
void DistMarkingStrategyCL<instat_scalar_fun_ptr>::SetFineLevel( Uint level )
{
    f_level_ = level;
}

inline
double DistMarkingStrategyCL<instat_scalar_fun_ptr>::GetTime() const
{
    return time_;
}

inline
void DistMarkingStrategyCL<instat_scalar_fun_ptr>::SetTime( double time )
{
    time_ = time;
}

inline
double DistMarkingStrategyCL<instat_scalar_fun_ptr>::GetValue( const VertexCL& v ) const
{
    return dist_( v.GetCoord(), time_ );
}

inline
double DistMarkingStrategyCL<instat_scalar_fun_ptr>::GetValue( const EdgeCL&   e ) const
{
    return dist_( GetBaryCenter(e), time_ );
}

inline
double DistMarkingStrategyCL<instat_scalar_fun_ptr>::GetValue( const TetraCL&  t ) const
{
    return dist_( GetBaryCenter(t), time_ );
}

}

#include "levelset/dist_marking_strategy.tpp"
#endif

