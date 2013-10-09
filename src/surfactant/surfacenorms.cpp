/// \file
/// \brief Evaluation-tool for solutions produced by surfactant/surfactant.cpp
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
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "surfactant/ifacetransp.h"
#include "misc/params.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "levelset/surfacetension.h"
#include "out/ensightOut.h"
#include "misc/dynamicload.h"

#include <fstream>

DROPS::ParamCL P;

DROPS::Point3DCL u_func (const DROPS::Point3DCL&, double)
{
    return P.get<DROPS::Point3DCL>("Exp.Velocity");
}

double sphere_2 (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL x( p - P.get<DROPS::Point3DCL>("Exp.PosDrop"));

    return x.norm() - P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
}

double sphere_2move (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x( p - (P.get<DROPS::Point3DCL>("Exp.PosDrop") + t*u_func(p, t)));

    return x.norm() - P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
}

typedef double (*dist_funT) (const DROPS::Point3DCL&, double);

void LSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, dist_funT d, double t)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( 0.5*(it->GetVertex( 0)->GetCoord() + it->GetVertex( 1)->GetCoord()), t);
}

const double a( -13./8.*std::sqrt( 35./M_PI));

double sol0t (const DROPS::Point3DCL& p, double t)
{
    const DROPS::Point3DCL q( p - (P.get<DROPS::Point3DCL>("Exp.PosDrop") + t*u_func(p, t)));
    const double val( a*(3.*q[0]*q[0]*q[1] - q[1]*q[1]*q[1]));

    return q.norm_sq()/(12. + q.norm_sq())*val;
}

/// \brief Short-hand for simple loops over the interface.
/// \param t  - Reference to a tetra
/// \param ls - Levelset-reference: Something that can be handed to InterfacePatchCL::Init as 2nd argument.
/// \param bnd  - ???
/// \param p  - The InterfacePatchCL that should be used.
/// \param n  - Name of the integer to reference the interface-triangles
#define DROPS_FOR_TETRA_INTERFACE_BEGIN( t, ls, bnd, p, n) \
    (p).Init( (t), (ls), (bnd)); \
    if ((p).Intersects()) { /*We are at the phase boundary.*/ \
        for (int ch__= 0; ch__ < 8; ++ch__) { \
            (p).ComputeForChild( ch__); \
            for (int n= 0; n < (p).GetNumTriangles(); ++n) \

#define DROPS_FOR_TETRA_INTERFACE_END }}


template <typename DiscP1FunT>
double Integral_Gamma (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const DROPS::BndDataCL<>& lsbnd,
    const DiscP1FunT& discsol)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<> qdiscsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        DROPS_FOR_TETRA_INTERFACE_BEGIN( *it, ls, lsbnd, triangle, tri) {
            qdiscsol.assign(  *it, &triangle.GetBary( tri), discsol);
            d+= qdiscsol.quad( triangle.GetAbsDet( tri));
        }
        DROPS_FOR_TETRA_INTERFACE_END
    }
    return d;
}


template <typename DiscP1FunT>
double L2_norm (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const DROPS::BndDataCL<>& lsbnd,
    const DiscP1FunT& discsol)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<> qdiscsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        DROPS_FOR_TETRA_INTERFACE_BEGIN( *it, ls, lsbnd, triangle, tri) {
            qdiscsol.assign(  *it, &triangle.GetBary( tri), discsol);
            d+= DROPS::Quad5_2DCL<>( qdiscsol*qdiscsol).quad( triangle.GetAbsDet( tri));
        }
        DROPS_FOR_TETRA_INTERFACE_END
    }
    return std::sqrt( d);
}

template <typename DiscP1FunT>
double L2_norm_Omega (const DROPS::MultiGridCL& mg, const DiscP1FunT& discsol)
{
    double d( 0.);
    const DROPS::Uint lvl = discsol.GetLevel();
    DROPS::Quad5CL<> qdiscsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        qdiscsol.assign(  *it, discsol);
        d+= DROPS::Quad5CL<>( qdiscsol*qdiscsol).quad( it->GetVolume()*6.);
    }
    return std::sqrt( d);
}

template <typename DiscP1FunT>
double L2_err (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const DROPS::BndDataCL<>& lsbnd,
    const DiscP1FunT& discsol, dist_funT sol, double t)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<> qdiscsol, qsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        DROPS_FOR_TETRA_INTERFACE_BEGIN( *it, ls, lsbnd, triangle, tri) {
            qdiscsol.assign(  *it, &triangle.GetBary( tri), discsol);
            qsol.assign( *it, &triangle.GetBary( tri), sol, t);
            qdiscsol-= qsol;
            d+= DROPS::Quad5_2DCL<>( qdiscsol*qdiscsol).quad( triangle.GetAbsDet( tri));
        }
        DROPS_FOR_TETRA_INTERFACE_END
    }
    return std::sqrt( d);
}

namespace DROPS // for Strategy
{

const int N= 11;
const int L= 4;

void Strategy (DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
{
    using namespace DROPS;

    MultiGridCL& mg= adap.GetMG();
    IdxDescCL fullidx( P1_FE);
    fullidx.CreateNumbering( mg.GetLastLevel(), mg);
    NoBndDataCL<> bnd;

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);

    VecDescCL DV[L][N];
    InterfaceP1RepairCL* ifrep[L][N];
    const char* lvlstr[]= { "-lvl2/", "-lvl3/", "-lvl4/", "-lvl5/" };

  for (int l= 0; l < L; ++l) {
    for (int i= 0; i < N; ++i) {
        DV[l][i].SetIdx( &fullidx);
        DV[l][i].t = 1.0;
    }

    std::string ensdir = P.get<std::string>("EnsightDir");
    std::string enscase = P.get<std::string>("EnsightCase");
    ReadEnsightP2SolCL reader( mg);
    reader.ReadScalar( ensdir + lvlstr[l] + enscase + ".sur1",   DV[l][0], bnd);
    reader.ReadScalar( ensdir + lvlstr[l] + enscase + ".sur2",   DV[l][1], bnd);
    reader.ReadScalar( ensdir + lvlstr[l] + enscase + ".sur4",   DV[l][2], bnd);
    reader.ReadScalar( ensdir + lvlstr[l] + enscase + ".sur8",   DV[l][3], bnd);
    reader.ReadScalar( ensdir + lvlstr[l] + enscase + ".sur16",  DV[l][4], bnd);
    reader.ReadScalar( ensdir + lvlstr[l] + enscase + ".sur32",  DV[l][5], bnd);
    reader.ReadScalar( ensdir + lvlstr[l] + enscase + ".sur64",  DV[l][6], bnd);
    reader.ReadScalar( ensdir + lvlstr[l] + enscase + ".sur128", DV[l][7], bnd);
    reader.ReadScalar( ensdir + lvlstr[l] + enscase + ".sur256", DV[l][8], bnd);
    reader.ReadScalar( ensdir + lvlstr[l] + enscase + ".sur512", DV[l][9], bnd);
    reader.ReadScalar( ensdir + lvlstr[l] + enscase + ".sur1024", DV[l][10], bnd);

    if (l == L - 1) continue;
    for (int i= 0; i < N; ++i) {
        ifrep[l][i]= new InterfaceP1RepairCL( mg, lset.Phi, lset.GetBndData(), DV[l][i]);
        adap.push_back( ifrep[l][i]);
    }
    adap.SetFineLevel( mg.GetLastLevel() + 1);
    adap.UpdateTriang( lset);
    LSInit( mg, lset.Phi, &sphere_2move, 1.);
  }

    // LSInit( mg, lset.Phi, &sphere_2move, 0.);
    // VecDescCL DVinitial( &fullidx);
    // reader.ReadScalar( "ensight/surfactant.sur0",   DVinitial, bnd);
    // std::cout << "\ninitial discretization error: "
    //           << L2_err( mg, lset.Phi, make_P1Eval( mg, bnd, DVinitial, 0.), &sol0t, 0.) << std::endl;

    LSInit( mg, lset.Phi, &sphere_2move, 1.);
    // std::cout << "\nfinal discretization error:" << std::endl;
    // for (int i= 0; i < N - 1; ++i)
    //     std::cout << L2_err( mg, lset.Phi, make_P1Eval( mg, bnd, DV[0][i], 1.), &sol0t, C.tm_StepSize*C.tm_NumSteps) << ", ";

    std::cout << std::endl;
  for (int l= 0; l < L; ++l) {
    std::cout << " Level: " << l + 2 << ":\n";
    std::cout << "Integral on \\Gamma: ";
    for (int i= 0; i < N - 1; ++i) {
        std::cout << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval( mg, bnd, DV[l][i])) << ", ";
    }
    std::cout << "\nL_2 error: ";
    for (int i= 0; i < N - 1; ++i) {
        DV[l][i].Data-=DV[L-1][N - 1].Data;
        // std::cout << norm( DV[i].Data) << ", ";
        std::cout << L2_norm( mg, lset.Phi, lset.GetBndData(), make_P1Eval( mg, bnd, DV[l][i])) << ", ";
    }
    // std::cout << std::endl;
    // for (int i= 0; i < N - 1; ++i) {
    //     std::cout << L2_norm_Omega( mg, make_P1Eval( mg, bnd, DV[i], 1.)) << ", ";
    // }
    std::cout << std::endl;
  }
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "surfactant.json");
    std::cout << P << std::endl;

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

    std::cout << "Setting up interface-PDE:\n";
    DROPS::BrickBuilderCL brick( DROPS::MakePoint3D( -2., -2., -2.),
                                 4.*DROPS::std_basis<3>( 1),
                                 4.*DROPS::std_basis<3>( 2),
                                 4.*DROPS::std_basis<3>( 3),
                                 P.get<int>("InitialDivisions"), P.get<int>("InitialDivisions"), P.get<int>("InitialDivisions"));
    DROPS::MultiGridCL mg( brick);
    DROPS::AdapTriangCL adap( mg, P.get<double>("AdaptRef.Width"), 0, P.get<int>("AdaptRef.FinestLevel"));
    adap.MakeInitialTriang( sphere_2);

    const DROPS::BndCondT bcls[6]= { DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC };
    const DROPS::LsetBndDataCL::bnd_val_fun bfunls[6]= { 0,0,0,0,0,0};
    DROPS::LsetBndDataCL lsbnd( 6, bcls, bfunls);

    DROPS::instat_scalar_fun_ptr sigma (0);
    DROPS::SurfaceTensionCL sf( sigma, 0);
    DROPS::LevelsetP2CL & lset( * DROPS::LevelsetP2CL::Create( mg, lsbnd, sf) );

    lset.idx.CreateNumbering( mg.GetLastLevel(), mg);
    lset.Phi.SetIdx( &lset.idx);
    Strategy( adap, lset);
    delete &lset;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
