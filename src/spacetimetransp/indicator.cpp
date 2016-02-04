/// \file indicator.cpp
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

#include "spacetimetransp/indicator.h"

namespace DROPS
{


///////////////////////////
// ConcentrationMarkingStrategyCL //
///////////////////////////

ConcentrationMarkingStrategyCL::ConcentrationMarkingStrategyCL
                                  ( const LevelsetP2CL &lset,
                                    P1EvalCL<double,cBndDataCL,VecDescCL > & a_solneg, 
                                    P1EvalCL<double,cBndDataCL,VecDescCL > & a_solpos, 
                                    double * thresholdlist, Uint coarse_level, Uint fine_level,
                                    bool a_hacked,
                                    double a_hacked_width,
                                    std::ofstream * & distributionout,
                                    double )
: lsetgetter_( new LevelsetP2GetterCL(lset) ), solneg(a_solneg), solpos(a_solpos), 
  thresholdlist_( thresholdlist ), c_level_( coarse_level ),
  f_level_( fine_level ), hacked_(a_hacked), hacked_width_(a_hacked_width), modified_( false ), decision_( DontCareC ), distributionout_( distributionout), owndistout_(true)
{
    if (distributionout_ != NULL) delete distributionout_;
    distributionout_ = new std::ofstream("distrib.out");

    if ( f_level_ < c_level_ )
    {
        throw DROPSErrCL( "Refinement levels are cheesy.\n" );
    }
}

ConcentrationMarkingStrategyCL::ConcentrationMarkingStrategyCL( const ConcentrationMarkingStrategyCL& rhs ):
 lsetgetter_( rhs.lsetgetter_->clone() ), solneg(rhs.solneg), solpos(rhs.solpos), thresholdlist_( rhs.thresholdlist_ ),
 c_level_( rhs.c_level_ ), f_level_( rhs.f_level_ ), hacked_( rhs.hacked_), hacked_width_( rhs.hacked_width_),
 modified_( rhs.modified_ ), decision_( rhs.decision_ ), distributionout_( rhs.distributionout_), owndistout_(false)
{
}

ConcentrationMarkingStrategyCL::~ConcentrationMarkingStrategyCL()
{
    delete lsetgetter_;
    if (owndistout_)
        delete distributionout_;
}

ConcentrationMarkingStrategyCL& ConcentrationMarkingStrategyCL::operator=( const ConcentrationMarkingStrategyCL& rhs )
{
    if ( &rhs == this ) return *this;

    delete lsetgetter_;
    lsetgetter_ = 0; // Needed for exception safety.

    lsetgetter_ = rhs.lsetgetter_->clone();
    solneg = rhs.solneg;
    solpos = rhs.solpos;
    thresholdlist_  = rhs.thresholdlist_;
    c_level_ = rhs.c_level_;
    f_level_ = rhs.f_level_;
    modified_ = rhs.modified_;
    decision_ = rhs.decision_;

    return *this;
}

TetraAccumulatorCL* ConcentrationMarkingStrategyCL::clone( int )
{
    return new ConcentrationMarkingStrategyCL( *this );
}

MarkingStrategyCL* ConcentrationMarkingStrategyCL::clone_strategy()
{
    return new ConcentrationMarkingStrategyCL( *this );
}

Point3DCL ConcentrationMarkingStrategyCL::GetGradientOfTetra( const TetraCL &t, bool is_pos) 
{
    LocalP1CL<> c(t, is_pos ? solpos : solneg);
    double det = 0;
    Point3DCL gradsi[4];
    P1DiscCL::GetGradients ( gradsi, det, t);
    Point3DCL grad = MakePoint3D(0,0,0);
    for (Uint i = 0; i < 4; ++i)
        grad += c[i] * gradsi[i];

    return grad;
}

double ConcentrationMarkingStrategyCL::GetMeanConcentrationOfTetra( const TetraCL &t, bool is_pos) 
{
    LocalP1CL<> c(t, is_pos ? solpos : solneg);
    double res = 0.0;
    for (Uint i = 0; i < 4; ++i)
        res += 0.25 * c[i];
    return res;
}

void ConcentrationMarkingStrategyCL::visit( const TetraCL &t )
{

    double lsetmean = 0.0;

    Uint soll_level = c_level_;
    const Uint l = t.GetLevel();

    double d = std::numeric_limits<double>::max();
    int num_pos = 0;
    for ( Uint j = 0; j < 4; ++j )
    {
        const double dist = lsetgetter_->GetValue( *t.GetVertex(j) );
        lsetmean += dist/10.0;
        if ( dist >= 0 ) ++num_pos;
        d = std::min( d, std::abs(dist) );
    }

    for ( Uint j = 0; j < 6; ++j )
    {
        const double dist = lsetgetter_->GetValue( *t.GetEdge(j) );
        lsetmean += dist/10.0;
        if ( dist >= 0 ) ++num_pos;
        d = std::min( d, std::abs(dist) );
    }

    d = std::min( d, std::abs(lsetgetter_->GetValue(t)) );
    const bool vzw = num_pos!=0 && num_pos!=10; // change of sign


    bool centerline = false;
    const bool hacked_centerline = hacked_;
    if (hacked_centerline)
    {
        const double time = solneg.GetSolution()->t;
        for ( Uint j = 0; j < 4; ++j )
        {
            const double x = (t.GetVertex(j)->GetCoord())[0];
            const double y = (t.GetVertex(j)->GetCoord())[1];
            const double z = (t.GetVertex(j)->GetCoord())[2];
            const double r = std::sqrt(x*x+y*y);
            if ( r < hacked_width_ + 1e-16)
                if ( (z >= -0.002 ) && (z <= 0.1 * time - 0.002 ))
                    centerline = true;
        }
    }

    if ((hacked_centerline && centerline) || vzw)
    {
        soll_level = f_level_;
    }
    else
    {
        double errormeas = 0;
// #ifdef _PAR
        // double res = GetMeanConcentrationOfTetra(t,num_pos>0);
        // errormeas = res; // 1e-4 * grad.norm() * grad.norm();
// #else
        Point3DCL grad = GetGradientOfTetra(t,num_pos>0);
        for (Uint f = 0; f < NumFacesC; ++f)
        {
            const TetraCL * neighb = t.GetNeighbor(f);
            if (neighb != NULL)
            {
                const FaceCL * face = t.GetFace(f);
                
                Point3DCL points[3];
                for (Uint v = 0; v < 3; ++v)
                {
                    const VertexCL * verti = face->GetVertex(v);
                    points[v] = verti->GetCoord();
                }

                Point3DCL diffa = points[1] - points[0];
                Point3DCL diffb = points[2] - points[0];
                Point3DCL diffcross;
                cross_product(diffcross,diffa,diffb);
                const double area = diffcross.norm();

                Point3DCL normal;
                t.GetOuterNormal(f, normal);
                double gradn = inner_prod(normal, grad);
                
                Point3DCL neighbgrad = GetGradientOfTetra(*neighb,num_pos>0);
                double neighbgradn = inner_prod(normal, neighbgrad);
                
                errormeas += area * (neighbgradn - gradn) * (neighbgradn - gradn);
            }
        }
// #endif

        // double threshold * = {1e-3, 1e-4, 1e-5};
        
#pragma omp critical(distout)
        *distributionout_ << lsetmean << "\t" << errormeas << "\n";

        const Uint nl = f_level_ - c_level_;
        for (Uint i = 0; i < nl; ++i)
        {
            // if( i==nl-1 && num_pos > 0) continue;
            Uint tl = nl - 1 - i;
            // if (thresholdlist_[tl] > 0)
            {
                if (errormeas > 0.25*thresholdlist_[tl] && l == f_level_ - tl) soll_level = f_level_-tl; //stay
                if (errormeas > thresholdlist_[tl] && l < f_level_ - tl) soll_level = f_level_-tl; //refine
                if (errormeas < 0.25*thresholdlist_[tl] && l >= f_level_ - tl) soll_level = f_level_- (tl+1);//crsn
            }
            // else // concentration below 1 - thresholdlist ? (hack for parallel)
            // {
            //     if (errormeas < 1+0.25*thresholdlist_[tl] && l == f_level_ - tl) soll_level = f_level_-tl; //stay
            //     if (errormeas < 1+thresholdlist_[tl] && l < f_level_ - tl) soll_level = f_level_-tl; //refine
            //     if (errormeas < 1+0.25*thresholdlist_[tl] && l >= f_level_ - tl) soll_level = f_level_- (tl+1);//crsn

            // }
        }

    }

    //avoid levels out of bounds
    soll_level = std::max(soll_level, c_level_);
    soll_level = std::min(soll_level, f_level_);

    // In the shell:      level should be f_level_.
    // Outside the shell: level should be c_level_.

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

bool ConcentrationMarkingStrategyCL::modified() const
{
    return modified_;
}

void ConcentrationMarkingStrategyCL::SetUnmodified()
{
    modified_ = false;
}

MarkingDecisionT ConcentrationMarkingStrategyCL::GetDecision() const
{
    return decision_;
}

// void ConcentrationMarkingStrategyCL::SetDistFct( instat_scalar_fun_ptr fct, double time )
// {
//     delete lsetgetter_;
//     lsetgetter_ = 0; // Needed for exception safety.

//     lsetgetter_ = new FunPtrGetterCL( fct, time );
// }

void ConcentrationMarkingStrategyCL::SetDistFct( const LevelsetP2CL &fct )
{
    delete lsetgetter_;
    lsetgetter_ = 0; // Needed for exception safety.

    lsetgetter_ = new LevelsetP2GetterCL( fct );
}

// double ConcentrationMarkingStrategyCL::GetWidth() const
// {
//     return width_;
// }

// void ConcentrationMarkingStrategyCL::SetWidth( double width )
// {
//     width_ = width;
// }


Uint ConcentrationMarkingStrategyCL::GetCoarseLevel() const
{
    return c_level_;
}

void ConcentrationMarkingStrategyCL::SetCoarseLevel( Uint level )
{
    c_level_ = level;
    if ( f_level_ < c_level_ )
    {
        throw DROPSErrCL( "Refinementlevels are cheesy.\n" );
    }
}

Uint ConcentrationMarkingStrategyCL::GetFineLevel() const
{
    return f_level_;
}

void ConcentrationMarkingStrategyCL::SetFineLevel( Uint level )
{
    f_level_ = level;
    if ( f_level_ < c_level_ )
    {
        throw DROPSErrCL( "Refinementlevels are cheesy.\n" );
    }
}


}
