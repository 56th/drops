/// \file indicator.h
/// \brief classes that estimates/indicates error for mass transport
/// \author LNM RWTH Aachen: Christoph Lehrenfeld

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
 * Copyright 2013 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_INDICATOR_H
#define DROPS_INDICATOR_H

#include "levelset/marking_strategy.h"
#include "levelset/levelset.h"

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "misc/problem.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "misc/params.h"
#include "levelset/levelset.h"
#include "spacetimetransp/stxfem.h"


namespace DROPS
{

typedef DROPS::BndDataCL<> cBndDataCL;
typedef cBndDataCL::bnd_val_fun  c_bnd_val_fun;

/*!
 * \brief Refinement strategy according to a concentration field.
 *
 * This marking strategy does only heuristic stuff right now.. (TODO!)
 */
class ConcentrationMarkingStrategyCL: public MarkingStrategyCL
{
public:
    ConcentrationMarkingStrategyCL( const LevelsetP2CL &fct,
                                    P1EvalCL<double,cBndDataCL,VecDescCL > & solneg, 
                                    P1EvalCL<double,cBndDataCL,VecDescCL > & solpos, 
                                    double * thresholdlist, Uint coarse_level, Uint fine_level,
                                    bool hacked,
                                    double hacked_width,
                                    std::ofstream * & distributionout,
                                    double time = 0. );
    ConcentrationMarkingStrategyCL( const ConcentrationMarkingStrategyCL &rhs );
    ~ConcentrationMarkingStrategyCL();
    ConcentrationMarkingStrategyCL& operator=( const ConcentrationMarkingStrategyCL& rhs );
    TetraAccumulatorCL* clone( int );
    MarkingStrategyCL*  clone_strategy();

    Point3DCL GetGradientOfTetra( const TetraCL &t, bool is_pos ); 
    double GetMeanConcentrationOfTetra( const TetraCL &t, bool is_pos ); 

    void visit( const TetraCL& t );

    bool modified() const;
    void SetUnmodified();

    MarkingDecisionT GetDecision() const;

    /* void SetDistFct( instat_scalar_fun_ptr fct, double time = 0 ); */
    void SetDistFct( const LevelsetP2CL& fct );

    /* double GetWidth() const; */
    /* void   SetWidth( double width ); */
    Uint   GetCoarseLevel() const;
    void   SetCoarseLevel( Uint level );
    Uint   GetFineLevel() const;
    void   SetFineLevel( Uint level );

    void   ResetOutput(int step){ 
        delete distributionout_; 
        std::ostringstream name;
        name << "distrib" << step << ".out";
        distributionout_ = new std::ofstream(name.str().c_str()); 
    }

private:
    ValueGetterCL *lsetgetter_;
    P1EvalCL<double,cBndDataCL,VecDescCL > solneg;
    P1EvalCL<double,cBndDataCL,VecDescCL > solpos;
    double * thresholdlist_;
    Uint c_level_;
    Uint f_level_;
    bool hacked_;
    double hacked_width_;
    bool modified_;
    MarkingDecisionT decision_;
    std::ofstream *& distributionout_;
    bool owndistout_;

};

}

#endif 
