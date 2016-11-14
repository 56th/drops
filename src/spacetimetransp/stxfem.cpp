/// \file stxfem.cpp
/// \brief space time xfem classes: IdxDescCL
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
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

#include "spacetimetransp/stxfem.h"
#include "num/interfacePatch.h"

#include "misc/params.h" // just for development time (for ParamCL to access development options..)
#include "misc/problem.h"

#ifdef _PAR
#  include "parallel/parallel.h"
#  include "parallel/exchange.h"
#endif

extern DROPS::ParamCL P;

namespace DROPS
{

namespace STXFEM
{

/// not yet sure if the old criteria is a good approach also in space-time
void UpdateSTXNumbering( IdxDescCL* Idx, const MultiGridCL& mg, 
                         const VecDescCL& lset_old, const VecDescCL& lset_new, 
                         instat_scalar_fun_ptr lset_fpt, const TimeInterval ti,
                         const BndDataCL<>& lsetbnd, bool NumberingChanged, double vmax)
{
    std::cout << " starting STXNumbering ... ";
            
    const int ints_per_space_edge = P.get<double>("Transp.Quadrature.SubIntervalsPerEdge");
    const int subtimeintervals = P.get<double>("Transp.Quadrature.SubTimeIntervals");
    ExtIdxDescCL& extidxdesc = Idx->GetXidx();
    const Uint sysnum= Idx->GetIdx(),
        level= Idx->TriangLevel(),
        stride= Idx->IsScalar() ? 1: 3; // scalar or vector-valued FE
                                        // (!) 2 unknowns per (spatial) vertex due to space time
    IdxT extIdx = NumberingChanged ? Idx->NumUnknowns() : extidxdesc.GetNumUnknownsStdFE();
    extidxdesc.SwapExtIdxOldNew(); // stores the least "current" as "old" before changing new "current"
    extidxdesc.ResetExtIdx(extIdx,NoIdx);

    const double omit_bound = extidxdesc.GetBound();

    LocalP2CL<> ps[4]; // values of phi_i*phi_i
    LocalP2CL<> p0; // values of phi_i*phi_i

    for (int i=0; i<4; ++i) //initialize hat_ii
        ps[i][i]=1.;

    LocalP2CL<> locPhi_old, locPhi_new;

    Point4DContainer pcontref;
    GeneralizedPrism4CL refprism4(pcontref);

    const double dt = lset_new.t - lset_old.t;
    Ulint no_skipped_ints = 0;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, level, it)
    {
        locPhi_old.assign( *it, lset_old, lsetbnd);
        locPhi_new.assign( *it, lset_new, lsetbnd);

        Point3DCL cc; // unused; circumcenter of T dummy
        double hT; // radius of circumcircle of T
        circumcircle(*it, cc, hT);

        bool inside_outer_band = false;
        bool outside_inner_band = false;
    
        for (int i = 0; i < 10; ++i)
        {
            if (locPhi_old[i] < vmax*dt || locPhi_new[i] < vmax*dt) 
                inside_outer_band = true;
            if (locPhi_old[i] > -vmax*dt || locPhi_new[i] >- vmax*dt) 
                outside_inner_band = true;
        }
        if (!inside_outer_band || !outside_inner_band) //not close to interface
            continue;

        CompositeSTQuadCL<QuadRule> * p_cstquad;
        if (lset_fpt ==NULL)
            p_cstquad = new CompositeSTQuadCL<QuadRule>(refprism4,locPhi_old, locPhi_new, 
                                                                      ints_per_space_edge, subtimeintervals);
        else
            p_cstquad = new CompositeSTQuadCL<QuadRule>(*it, ti, lset_fpt,
                                                                      ints_per_space_edge, subtimeintervals);
        CompositeSTQuadCL<QuadRule> & cstquad(*p_cstquad);

        if (!cstquad.HasInterface()){
            delete p_cstquad;
            continue;
        }

        const double hs3= it->GetVolume()*6.0,
            hs= cbrt( hs3), //, hs5= hs*hs*hs3, // hs^5 (hs: spatial h)
            limit= hs*hs*omit_bound;

        for (Uint b = 0; b < 2; b++)
        {
            const bool newtime = (b==1);
            for (Uint i = 0; i < 4; ++i)
            {
                const bool sign = cstquad.GetVertexSign(i,/* newtime? */ newtime);
                GridFunctionCL<double> pi = newtime ? cstquad.EvalLinearOnPart( p0, ps[i], !sign) 
                    : cstquad.EvalLinearOnPart( ps[i], p0, !sign);
                
                GridFunctionCL<double> pi2 = pi;
                pi2 *= pi;
                // GridFunctionCL<double> pi2 = pi*pi;
                //GridFunctionCL<double> pi2 (1.0);
                const double locint = cstquad.QuadOnPart( pi2, !sign);

                // - ? - ASK SVEN - ? - => if (cut.IntersectsInterior())
                // std::cout << "locint = " << locint[i] << " | limit = " << limit << std::endl;
                if (locint < limit) {
                    no_skipped_ints++;
                    continue; // omit DoFs of minor importance, which lead to unstable solutions
                }
                if (!it->GetVertex(i)->Unknowns.Exist( sysnum)) continue;
                const IdxT nr= it->GetVertex(i)->Unknowns(sysnum);
                if (extidxdesc[nr+b*stride]==NoIdx) // not extended yet
                    for (Uint k=0; k<stride; ++k)
                        extidxdesc[nr+k+b*stride]= extIdx++;
            }
        }
        delete p_cstquad;


    }
#ifdef _PAR
    // communicate extended dofs on vertices
    ExtIdxDescCL::CommunicateXFEMNumbCL comm( Idx);
    comm.Call();

    // number all extended dofs from other procs (where extended dof is flagged by NoIdx-1)
    for (size_t i=0; i<extidxdesc.GetNumUnknownsStdFE(); ++i)
        if (extidxdesc[i] == NoIdx-1) // remote is extended
            extidxdesc[i]= extIdx++;
#endif

    Idx->ChangeNumUnknowns(extIdx);
#ifdef _PAR
    Idx->GetEx().CreateList(mg, Idx, true, true);
#endif
    std::cout << " finished. " << std::endl;

    // return extIdx;
    // for (size_t i=0; i<extidxdesc.GetNumUnknownsStdFE(); ++i)
    //          std::cout << " i = " << i << ": "  << extidxdesc[i] << std::endl;
    std::cout << " UpdateSTXNumbering: \n"; 
    std::cout << " number of std-functions: " << extidxdesc.GetNumUnknownsStdFE() << std::endl;
    std::cout << " number of x-functions: " << extIdx - extidxdesc.GetNumUnknownsStdFE() << std::endl;
#ifdef _PAR
    std::cout << " Idx->GetNumOwnedUnknowns() = " << Idx->GetNumOwnedUnknowns() << std::endl;
    std::cout << " Idx->GetGlobalNumUnknowns() = " << Idx->GetGlobalNumUnknowns() << std::endl;
#endif
    std::cout << " skipped " << no_skipped_ints << " parts (not dof!) " << std::endl;

}

template <TIMEDIRECTION timedir>
void GetTrace( const VecDescCL& stsol, VecDescCL& sol, const MultiGridCL& mg)
{
    const Uint lvl = stsol.GetLevel();
    const Uint stidx = stsol.RowIdx->GetIdx();
    const Uint idx = sol.RowIdx->GetIdx();
    const ExtIdxDescCL& STXidx= stsol.RowIdx->GetXidx();
    const ExtIdxDescCL& Xidx= sol.RowIdx->GetXidx();
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {    
        if (it->Unknowns.Exist( stidx)){
            const IdxT traceunkn = it->Unknowns( stidx) + timedir; //first dof is past - second future
            const IdxT tracexunkn = STXidx[traceunkn];
            const IdxT spaceunkn = it->Unknowns( idx);
            const IdxT spacexunkn = Xidx[spaceunkn];
            sol.Data[spaceunkn] = stsol.Data[traceunkn];
            if (spacexunkn != NoIdx && tracexunkn != NoIdx)
                sol.Data[spacexunkn] = stsol.Data[tracexunkn];
        }
    }
    
}


void GetFutureTrace( const VecDescCL& stsol, VecDescCL& sol, const MultiGridCL& mg) 
{
    GetTrace<FUTURE>(stsol,sol,mg);
}

void GetPastTrace( const VecDescCL& stsol, VecDescCL& sol, const MultiGridCL& mg)
{
    GetTrace<PAST>(stsol,sol,mg);
}



}//end of namespace DROPS::STXFEM
}//end of namespace DROPS
