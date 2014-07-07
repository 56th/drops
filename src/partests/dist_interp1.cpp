/// \file dist_interp1.cpp
/// \brief tests implementation of RepairP1CL
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "num/fe.h"
#include "num/fe_repair.h"
#include "misc/problem.h"
#include "num/discretize.h"
#include <tr1/unordered_set>


using namespace DROPS;

typedef double (*fun_ptr)(const SVectorCL<3>&);

enum  OutputModeT { SILENT, NOISY };


double f(const SVectorCL<3>& p)
//{ return 42; }
{ return p[0] +10.*p[1] +100.*p[2]+1; }

void MarkDrop(DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}

void UnMarkDrop(DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte( 0.5);

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It( mg.GetTriangTetraBegin(maxLevel)),
             ItEnd( mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It) {
        if ( (GetBaryCenter( *It)-Mitte).norm() <= std::max( 0.2, 1.5*std::pow( It->GetVolume(), 1.0/3.0)) )
            It->SetRemoveMark();
    }
}

typedef NoBndDataCL<double> BndCL;
BndCL Bnd;

void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg, fun_ptr f)
{
    vd.Data.resize( vd.RowIdx->NumUnknowns());
    P2EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun( &vd, &Bnd, &mg);
    const Uint lvl= vd.GetLevel();
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(lvl),
         theend= mg.GetTriangVertexEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( sit->GetCoord()));
    }
}


int CheckResult(DROPS::P2EvalCL<double, BndCL,
                const DROPS::VecDescCL>& fun, fun_ptr f, OutputModeT om, double eps= 0.)
{
    std::string filename("simplex.chk");
    DROPS::ProcCL::AppendProcNum(filename);
    std::ofstream sanityfile( filename.c_str());
    int ret= 0;
    const VertexCL* v= 0;
    const DROPS::MultiGridCL& mg= fun.GetMG();
    const DROPS::Uint trilevel= fun.GetLevel();
    if (om!=SILENT) std::cout << "Verts:" << std::endl;
    double diff, vmaxdiff= 0., vval= 1.;
    for (MultiGridCL::const_TriangVertexIteratorCL sit=mg.GetTriangVertexBegin( trilevel),
         theend= mg.GetTriangVertexEnd( trilevel); sit!=theend; ++sit) {
        diff= fun.val( *sit) - f( sit->GetCoord());
        if (std::abs( diff) > eps && std::abs(diff) > vmaxdiff) {
            ++ret;
            vmaxdiff= std::abs(diff);
            v= &*sit;
            vval= fun.val( *sit);
        }
    }
    if (om!=SILENT) {
        if (vval == 0.) vval= 1.;
        sanityfile << "maximale Differenz Vertices: " << vmaxdiff << " (exakter Wert: " << vval << ") auf\n";
        if (v) v->DebugInfo( sanityfile);
        sanityfile << std::endl;
    }
    return ret;
}


int TestRepairUniform()
{
    BndDataCL<> bnd( 6);
    int ret= 0;
    BrickBuilderCL builder( DROPS::std_basis<3>( 0), DROPS::std_basis<3>( 1),
                            DROPS::std_basis<3>( 2), DROPS::std_basis<3>( 3),
                            1, 1, 1);
    DROPS::MultiGridCL mg(builder);
    DROPS::LoadBalCL lb( mg);  // loadbalancing
    lb.DoMigration();          // distribute initial grid

    DROPS::IdxDescCL i0( P1_FE, Bnd), i1( P1_FE, Bnd);
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair for uniform refinement with quadratic function:\n";
    for (DROPS::Uint i=0; i<5; ++i) {
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        DROPS::VecDescCL v0, v1;
        v0.SetIdx( &i0);
        SetFun( v0, mg, f);
        RepairP1CL<double>::type repairp1( mg, v0, bnd);
        MarkAll( mg);
        mg.Refine();
        i1.CreateNumbering( i0.TriangLevel(), mg);
        v1.SetIdx( &i1);
        repairp1.repair( v1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1, f, NOISY, 1e-10);
        i0.DeleteNumbering( mg);
        i1.DeleteNumbering( mg);
        lb.DoMigration();
    }

    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair for uniform coarsening with quadratic function:\n";
    for (DROPS::Uint i=0; i<5; ++i) {
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        DROPS::VecDescCL v0, v1;
        v0.SetIdx(&i0);
        SetFun(v0, mg, f);
        RepairP1CL<double>::type repairp1( mg, v0, bnd);
        UnMarkAll( mg);
        mg.Refine();

        Uint i1_Level= i0.TriangLevel();
        if (mg.GetLastLevel() < i1_Level) {
            std::cout << "Letztes Level entfernt!" << std::endl;
            --i1_Level;
        }
        else {
            std::cout << "Letztes Level behalten!" << std::endl;
        }
        i1.CreateNumbering( i1_Level, mg);
        v1.SetIdx( &i1);
        repairp1.repair( v1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1, f, NOISY, 1e-10);
        if (mg.GetLastLevel() < i0.TriangLevel()) {
            Uint level= mg.GetLastLevel();
            DROPS::DeleteNumbOnSimplex( i0.GetIdx(), mg.GetAllVertexBegin( level),
                                        mg.GetAllVertexEnd( level));
            DROPS::DeleteNumbOnSimplex( i0.GetIdx(), mg.GetAllEdgeBegin( level),
                                        mg.GetAllEdgeEnd( level));
        }
        DROPS::DeleteNumbOnSimplex( i1.GetIdx(), mg.GetAllVertexBegin( i1.TriangLevel()),
                                    mg.GetAllVertexEnd( i1.TriangLevel()));
        DROPS::DeleteNumbOnSimplex( i1.GetIdx(), mg.GetAllEdgeBegin( i1.TriangLevel()),
                                    mg.GetAllEdgeEnd( i1.TriangLevel()));
        lb.DoMigration();

    }
    return ret;
}

int TestRepair()
{
    BndDataCL<> bnd( 6);
    int ret= 0;
    BrickBuilderCL builder( DROPS::std_basis<3>( 0), DROPS::std_basis<3>( 1),
                            DROPS::std_basis<3>( 2), DROPS::std_basis<3>( 3),
                            1, 1, 1);
    DROPS::MultiGridCL mg(builder);
    DROPS::LoadBalCL lb( mg);  // loadbalancing
    lb.DoMigration();          // distribute initial grid
    DROPS::IdCL<DROPS::VertexCL>::ResetCounter();
    DROPS::IdxDescCL i0( P1_FE, Bnd), i1( P1_FE, Bnd);
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair for drop refinement with quadratic function:\n";
    for (DROPS::Uint i=0; i<8; ++i) {
        std::cout << "i: " << i;
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        std::cout << " i0.TriangLevel(): " << i0.TriangLevel();
        DROPS::VecDescCL v0, v1;
        v0.SetIdx( &i0);
        SetFun( v0, mg, f);
        RepairP1CL<double>::type repairp1( mg, v0, bnd);
        MarkDrop( mg, mg.GetLastLevel());
        mg.Refine();
        std::cout << " mg.GetLastLevel() after refine: " << mg.GetLastLevel();
        i1.CreateNumbering( i0.TriangLevel(), mg);
        std::cout << " i1.TriangLevel(): " << i1.TriangLevel();
        v1.SetIdx( &i1);
        repairp1.repair( v1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1, f, NOISY, 1e-10);
        i0.DeleteNumbering( mg);
        i1.DeleteNumbering( mg);
        std::cout << '\n';
        lb.DoMigration();
    }

    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair for drop coarsening with quadratic function:\n";
    for (DROPS::Uint i=0; i<8; ++i) {
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        DROPS::VecDescCL v0, v1;
        v0.SetIdx(&i0);
        SetFun(v0, mg, f);
        RepairP1CL<double>::type repairp1( mg, v0, bnd);
        UnMarkDrop( mg, mg.GetLastLevel());
        mg.Refine();
        Uint i1_Level= i0.TriangLevel();
        if (mg.GetLastLevel() < i1_Level) {
            std::cout << "Letztes Level entfernt!" << std::endl;
            --i1_Level;
        }
        else {
            std::cout << "Letztes Level behalten!" << std::endl;
        }
        i1.CreateNumbering( i1_Level, mg);
        v1.SetIdx( &i1);
        repairp1.repair( v1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1, f, NOISY, 1e-10);
        if (mg.GetLastLevel() < i0.TriangLevel())
            i0.DeleteNumbering( mg);
        i1.DeleteNumbering( mg);
        lb.DoMigration();

    }
    return ret;
}


int main ( int argc, char **argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
#endif
  int ret = 0;
  try {
    ret +=TestRepairUniform();
    int sum = ProcCL::GlobalSum(ret);
    std::cout << "ret: " << sum << std::endl;
    ret += TestRepair();
    sum = ProcCL::GlobalSum(ret);
    std::cout << "ret: " << sum << std::endl;
  }
  catch (DROPS::DROPSErrCL& err) { err.handle(); }

  return ret;
}
