/// \file f_Gamma.cpp
/// \brief test surface force term
/// \author LNM RWTH Aachen: Sven Gross, Joerg Peters, Volker Reichelt; SC RWTH Aachen:

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

#include "geom/multigrid.h"
#include "out/output.h"
#include "out/vtkOut.h"
#include "geom/builder.h"
#include "stokes/instatstokes2phase.h"
#include "num/krylovsolver.h"
#include "num/precond.h"
#include "levelset/coupling.h"
#include "levelset/adaptriang.h"
#include "misc/params.h"
#include "levelset/surfacetension.h"
#include "misc/dynamicload.h"
#include <fstream>
#include <sstream>


DROPS::ParamCL P;

std::auto_ptr<DROPS::VTKOutCL> vtkwriter;

int FctCode=9;

double SurfTension;
double sigmaf( const DROPS::Point3DCL& p, double)
{
//     return SurfTension;
    return 1 + std::cos( 2.*M_PI*p[0]);
}

/*
double DistanceFct( const DROPS::Point3DCL& p)
{ // cube
    const DROPS::Point3DCL d= C.Mitte-p;
    double max=-1;
    for (int i=0; i<3; ++i)
        if (fabs(d[i])>max) max= fabs(d[i]);
    return max-C.Radius;
}

double DistanceFct( const DROPS::Point3DCL& p)
{ // plane perpendicular to n=PosDrop with distance Radius from origin.
    return inner_prod( C.Mitte/norm(C.Mitte), p) - C.Radius;
}

*/

DROPS::Point3DCL c;
double r;
double DistanceFct( const DROPS::Point3DCL& p, double)
{ // ball
//     return (p-c).norm()-r;
    return inner_prod( p - c, p - c) - r*r;
}

// double DistanceFct( const DROPS::Point3DCL& p, double)
// { // plane... XXX The line integral on \pa\iface must be implemented.
//     return p[2] - 1./3.;
// }

DROPS::Point3DCL Null( const DROPS::Point3DCL&, double)
{ return DROPS::Point3DCL(0.); }

static DROPS::RegisterVectorFunction regveczero("VecZero",  Null);


void WriteFct( std::ostream& os)
{
    const int idx= FctCode%10,
        fct= FctCode/10;
    os << "(";
    std::string term;
    switch (fct%10)
    {
        case 1:
            term= "x"; break;
        case 2:
            term= "y"; break;
        case 3:
            term= "z"; break;
        case 0: // 1
            term= "1"; break;
        default:
            std::cout << "Fehler in WriteFct, fct = " << fct << std::endl;
    }
    if (fct>=9)
        switch (fct/10)
        {
            case 1:
                term+= "x"; break;
            case 2:
                term+= "y"; break;
            case 3:
                term+= "z"; break;
            case 0: // 1
                term+= "1"; break;
            default:
                std::cout << "Fehler in WriteFct (quadratisch), fct = " << fct << std::endl;
        }
    for (int i=0; i<3; ++i)
        os << ' ' << (i==idx-1 ? term : "0");
    os << " )";
}


DROPS::Point3DCL TestFct( const DROPS::Point3DCL& p, double)
{
    const int idx= FctCode%10,
        fct= FctCode/10,
        fct1= fct%10, fct2= fct/10;
    double val= 0;
    switch (fct1)
    {
        case 1: case 2: case 3: // x/y/z
            val= p[fct1-1]; break;
        case 0: // 1
            val= 1; break;
        default:
            std::cout << "Fehler in TestFct, fct = " << fct << std::endl;
    }
    if (fct>=9)
        switch (fct2)
        {
            case 1: case 2: case 3: // x/y/z
                val*= p[fct2-1]; break;
            case 0: // 1
                val*= 1; break;
            default:
                std::cout << "Fehler in TestFct (quadratisch), fct = " << fct << std::endl;
        }

    DROPS::Point3DCL ret;
    if (idx<1 || idx>3) std::cout << "Fehler in TestFct, idx = " << idx << std::endl;
    ret[idx-1]= val;
    return ret;
}

namespace DROPS // for Strategy
{

void ApplyToTestFct( InstatStokes2PhaseP2P1CL& Stokes, const LsetBndDataCL& lsbnd)
// program for testing the approximation order of
// the discretization of the surface force term f_Gamma
// using simple trial functions:
//     | f_Gamma(v) - f_{Gamma_h}(v) |  =   O(h^p)
{
    MultiGridCL& MG= Stokes.GetMG();
    const double curv= 2/P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
    SurfaceTensionCL sf( sigmaf, 0);

    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsbnd, sf, P.get_child("Levelset")) );

//    lset.SetSurfaceForce( SF_Const);

    MLIdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    lset.CreateNumbering( MG.GetLastLevel(), lidx);

    lset.Phi.SetIdx( lidx);
    lset.Init( DistanceFct);

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, NULL, &lset);

    VecDescCL f_Gamma( vidx), v( vidx);
    MG.SizeInfo( std::cout);
    Stokes.b.SetIdx( vidx);
    Stokes.c.SetIdx( pidx);
    Stokes.p.SetIdx( pidx);
    Stokes.v.SetIdx( vidx);

    std::cout << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cout << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cout << lset.Phi.Data.size() << " levelset unknowns.\n";

    const double Vol= 4./3*M_PI*std::pow( P.get<DROPS::Point3DCL>("Exp.RadDrop")[0], 3);
    std::cout << "Volumen = " << Vol << "\tKruemmung = " << curv << "\n\n";
    typedef std::map<int,double> LsgMapT;
    LsgMapT reflsg;
    reflsg[1]= 0;
//    reflsg[2]= 0;
//    reflsg[3]= 0;
    reflsg[11]= curv*Vol;
    reflsg[12]= 0;
//    reflsg[13]= 0;
//    reflsg[21]= 0;
    reflsg[22]= curv*Vol;
//    reflsg[23]= 0;
//    reflsg[31]= 0;
//    reflsg[32]= 0;
    reflsg[33]= curv*Vol;
    reflsg[111]= curv*2*Vol*P.get<DROPS::Point3DCL>("Exp.PosDrop")[0];
    reflsg[121]= curv*Vol*P.get<DROPS::Point3DCL>("Exp.PosDrop")[1];
    reflsg[131]= curv*Vol*P.get<DROPS::Point3DCL>("Exp.PosDrop")[2];
//    reflsg[231]= 0;
//    reflsg[221]= 0;
//    reflsg[331]= 0;

    f_Gamma.Clear( Stokes.v.t);
    lset.AccumulateBndIntegral( f_Gamma);

    std::vector<double> errVec, refVec;
    for (LsgMapT::iterator it= reflsg.begin(), end= reflsg.end(); it!=end; ++it)
    {
        FctCode= it->first;
        const double Lsg= it->second;
        Stokes.InitVel( &v, TestFct);
        const double fv= dot( f_Gamma.Data, v.Data);
        const double err= std::abs(fv - Lsg);
        WriteFct( std::cout);
        std::cout << "\nfct " << FctCode << ":\tval = " << fv << "\tref = " << Lsg << "\nerror = " << err << "\n\n";
        errVec.push_back( err);
        refVec.push_back( Lsg);
    }

    std::cout << "\n\n\"f\",";
    for (LsgMapT::iterator it= reflsg.begin(), end= reflsg.end(); it!=end; ++it)
    { std::cout << "\""; FctCode= it->first; WriteFct( std::cout); std::cout << "\","; }

    std::cout << "\n\n\"Referenz\",";
    for (size_t i=0; i<refVec.size(); ++i)
        std::cout << refVec[i] << ",\t";

    std::cout << "\n\n" << P.get<int>("AdaptRef.FinestLevel") << ",\t";
    for (size_t i=0; i<errVec.size(); ++i)
        std::cout << errVec[i] << ",\t";
    std::cout << "\n\n";
    delete &lset;
}

void Compare_LaplBeltramiSF_ConstSF( InstatStokes2PhaseP2P1CL& Stokes, const LsetBndDataCL& lsbnd)
// computation of order of LB discretization, cf. paper
// S. Gross, A. Reusken: Finite element discretization error analysis of a surface tension force in two-phase incompressible flows,
// SIAM J. Numer. Anal. 45, 1679--1700 (2007)

{
    MultiGridCL& MG= Stokes.GetMG();
    // Levelset-Disc.: Crank-Nicholson
    const double curv= 2/P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
    SurfaceTensionCL sf( sigmaf, 0);
    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsbnd, sf, P.get_child("Levelset")) );


    MLIdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    lset.CreateNumbering( MG.GetLastLevel(), lidx);

    lset.Phi.SetIdx( lidx);
    lset.Init( DistanceFct);

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, NULL, &lset);

    VecDescCL f_Const( vidx), f_LaplBeltrami( vidx), v( vidx);
    MG.SizeInfo( std::cout);
    Stokes.b.SetIdx( vidx);
    Stokes.c.SetIdx( pidx);
    Stokes.p.SetIdx( pidx);
    Stokes.v.SetIdx( vidx);
    Stokes.A.SetIdx( vidx, vidx);
    Stokes.M.SetIdx( vidx, vidx);

    std::cout << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cout << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cout << lset.Phi.Data.size() << " levelset unknowns.\n";

    const double Vol= 4./3*M_PI*std::pow( P.get<DROPS::Point3DCL>("Exp.RadDrop")[0], 3);
    std::cout << "Volumen = " << Vol << "\tKruemmung = " << curv << "\n\n";

    Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &Stokes.b, lset, 0.);

    f_LaplBeltrami.Clear( Stokes.v.t);
    lset.SetSurfaceForce( SF_ImprovedLB);
//     lset.SetSurfaceForce( SF_LB);
    lset.AccumulateBndIntegral( f_LaplBeltrami);

    f_Const.Clear( Stokes.v.t);
    lset.SetSurfaceForce( SF_Const);
    lset.AccumulateBndIntegral( f_Const);

    VectorCL d( curv*f_Const.Data - f_LaplBeltrami.Data);
    std::cout << "|d| = \t\t" << norm(d) << std::endl;
    VectorCL A_inv_d( d.size());
    SSORPcCL pc;
    typedef PCGSolverCL<SSORPcCL>     PCG_SsorCL;
    PCG_SsorCL cg( pc, 1000, 1e-18);
    std::cout << "Solving system with stiffness matrix:\t";
    cg.Solve( Stokes.A.Data, A_inv_d, d, vidx->GetEx());
    std::cout << cg.GetIter() << " iter,\tresid = " << cg.GetResid();
    const double sup= std::sqrt(dot( A_inv_d, d));

    std::cout << "\n\nsup |f1(v)-f2(v)|/|v|_1 = \t\t" << sup
              << "\n|A^-1 d| = \t\t" << norm( A_inv_d) << std::endl;
    MLMatrixCL MA;
    MA.LinComb( 1, Stokes.M.Data, 1, Stokes.A.Data);
    VectorCL MA_inv_d( A_inv_d);
    std::cout << "Solving system with MA matrix:\t";
    cg.Solve( MA, MA_inv_d, d, vidx->GetEx());
    std::cout << cg.GetIter() << " iter,\tresid = " << cg.GetResid();
    const double sup2= std::sqrt(dot( MA_inv_d, d));
    std::cout << "\n\nsup |f1(v)-f2(v)|/||v||_1 = \t\t" << sup2
              << "\n|(MA)^-1 d| = " << norm( MA_inv_d) << std::endl;
    delete &lset;
}



/// \brief Accumulator for the interfacial tension term with variable interfacial tension coefficient.
class SphereObliqueLaplaceBeltramiAccuCL : public TetraAccumulatorCL
{
  private:
    VecDescCL& f;
    const LevelsetP2CL& ls_;
    LocalNumbP2CL n_;

    const PrincipalLatticeCL& lat_;
    LocalP2CL<> loc_ls_;
    std::valarray<double> ls_val_;
    SurfacePatchCL p_;
    QuadDomain2DCL q_;

    LocalP1CL<Point3DCL> gradref_[10],
                         grad_[10];
    GridFunctionCL<Point3DCL> w_[10],
                              qnt_,
                              qnh_;
    SMatrixCL<3, 3> T_;
    GridFunctionCL<double> qsigma;

    void
    resize_and_scatter_piecewise_spatial_normal (const SPatchCL<3>& surf, const QuadDomainCodim1CL<3>& qdom, std::valarray<Point3DCL>& spatial_normal);
   template <class ResultIteratorT>
      void
      evaluate_sigma_on_vertexes (const TetraCL& tet, const QuadDomain2DCL& q, double t, ResultIteratorT res);


  public:
    SphereObliqueLaplaceBeltramiAccuCL (const LevelsetP2CL& ls, VecDescCL& f_Gamma)
     : f( f_Gamma), ls_( ls), lat_( PrincipalLatticeCL::instance( 2)),
       ls_val_( 10)
    {
        P2DiscCL::GetGradientsOnRef( gradref_);
    }

    void begin_accumulation () {
        f.Data= 0.;
        std::cerr << "SphereObliqueLaplaceBeltramiAccuCL::begin_accumulation: " << f.Data.size() << "dof.\n";
    }
    void finalize_accumulation() {
        std::cerr << "SphereObliqueLaplaceBeltramiAccuCL::finalize_accumulation.\n";
    }

    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new SphereObliqueLaplaceBeltramiAccuCL ( *this); };
};

void
SphereObliqueLaplaceBeltramiAccuCL::resize_and_scatter_piecewise_spatial_normal (const SPatchCL<3>& surf, const QuadDomainCodim1CL<3>& qdom, std::valarray<Point3DCL>& spatial_normal)
{
    spatial_normal.resize( qdom.vertex_size());
    if (spatial_normal.size() == 0)
        return;
    if (surf.normal_empty()) // As qdom has vertexes, the must be facets, i.e. normals.
        throw DROPSErrCL( "resize_and_scatter_piecewise_spatial_normal: normals were not precomputed.\n");

    const Uint NodesPerFacet= qdom.vertex_size()/surf.facet_size();
    if (qdom.vertex_size()%surf.facet_size() != 0)
        throw DROPSErrCL( "resize_and_scatter_piecewise_spatial_normal: qdom.vertex_size is not a multiple of surf.facet_size.\n");

    const SPatchCL<3>::const_normal_iterator n= surf.normal_begin();
    for (Uint i= 0; i < surf.facet_size(); ++i) {
        std::fill_n( &spatial_normal[i*NodesPerFacet], NodesPerFacet, n[i]);
    }
}

template <class ResultIteratorT>
void SphereObliqueLaplaceBeltramiAccuCL::evaluate_sigma_on_vertexes (const TetraCL& tet, const QuadDomain2DCL& q, double t, ResultIteratorT res)
{
    Bary2WorldCoordCL mapper( tet);
    for (typename QuadDomain2DCL::const_vertex_iterator v= q.vertex_begin(), e= q.vertex_end(); v != e; ++v) {
        const Point3DCL& p= mapper( *v);
        *res++= sigmaf( (r/p.norm())*p, t);
    }
}

void SphereObliqueLaplaceBeltramiAccuCL::visit (const TetraCL& t)
{
    evaluate_on_vertexes( ls_.GetSolution(), t, lat_, Addr( ls_val_));
    if (equal_signs( ls_val_))
        return;

    loc_ls_.assign( t, ls_.GetSolution());
    double det; // dummy
    SMatrixCL<3,3> T;
    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( grad_, gradref_, T);
    LocalP1CL<Point3DCL> nt; // gradient of quadratic level set function
    for (Uint i= 0; i < 10 ; ++i)
        nt+= grad_[i]*loc_ls_[i];

    p_.make_patch<MergeCutPolicyCL>( lat_, ls_val_);
    p_.compute_normals( t);
    make_CompositeQuad5Domain2D( q_, p_, t);
    resize_and_scatter_piecewise_spatial_normal( p_, q_, qnh_); // unit-length normal to linear interface
    resize_and_evaluate_on_vertexes( nt, q_, qnt_); // normal to quadratic interface
    for (Uint i= 0; i < q_.vertex_size(); ++i)
        qnt_[i]/= norm( qnt_[i]);
    GridFunctionCL<> qalpha( dot( qnh_, qnt_));
    // If a triangle has zero area, its normal is returned as 0; we avoid 
    // division by zero... the value will not matter later as the func-det in quad_2D is also 0.
    for (size_t i=0; i < qalpha.size(); ++i)
        if (std::fabs( qalpha[i]) == 0.)
            qalpha[i]= 1.;

    // \alpha (1+d/R)
    GridFunctionCL<double> qf;
    resize_and_evaluate_on_vertexes( DistanceFct, t, q_, 0.,  qf);
    qf/= r;
    qf+= 1.;
    qf*=std::abs( qalpha); // For the projectors, the orientation of the normals cancels, but not here.

    GridFunctionCL<Point3DCL> qgradi( q_.vertex_size());
    for (Uint i= 0; i < 10; ++i) {
        evaluate_on_vertexes( grad_[i], q_, Addr( qgradi));
        w_[i].resize( q_.vertex_size());
        w_[i]= qgradi - GridFunctionCL<double>( dot( qnt_, qgradi)/qalpha)*qnh_; // \bQt D b^i, where b^i is the ith scalar P2-basis-function.
    }
    qsigma.resize( q_.vertex_size());
    this->evaluate_sigma_on_vertexes( t, q_, /*time*/ 0., Addr( qsigma)); // interfacial tension

    n_.assign_indices_only( t, *f.RowIdx);
    for (Uint i= 0; i < 10; ++i) {
        const GridFunctionCL<double>& tmp= GridFunctionCL<>( qsigma/qf);
        add_to_global_vector( f.Data, -quad_2D(tmp*w_[i], q_), n_.num[i]);
    }
}

void LSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, double (*d)(const Point3DCL&, double), double t= 0.)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( 0.5*(it->GetVertex( 0)->GetCoord() + it->GetVertex( 1)->GetCoord()), t);
    ls.t= t;
}

void Compare_Oblique_Improved (DROPS::AdapTriangCL& adap, InstatStokes2PhaseP2P1CL& Stokes, LevelsetP2CL& lset)
{
    MultiGridCL& MG= Stokes.GetMG();

    MLIdxDescCL* vidx= &Stokes.vel_idx;
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);

    MG.SizeInfo( std::cout);
    Stokes.b.SetIdx( vidx);
    Stokes.v.SetIdx( vidx);
    Stokes.A.SetIdx( vidx, vidx);
    Stokes.M.SetIdx( vidx, vidx);

    std::cout << vidx->NumUnknowns() << " velocity unknowns,\n";
    std::cout << lset.Phi.Data.size() << " levelset unknowns.\n";

    const double Vol= 4./3*M_PI*std::pow( P.get<DROPS::Point3DCL>("Exp.RadDrop")[0], 3);
    const double curv= 2./P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
    std::cout << "Volumen = " << Vol << "\tKruemmung = " << curv << "\n\n";

    VecDescCL f_improved( vidx);
    lset.SetSurfaceForce( SF_ImprovedLBVar);
    lset.AccumulateBndIntegral( f_improved);

    VecDescCL f_oblique( vidx);
    lset.SetSurfaceForce( SF_ObliqueLBVar);
    lset.AccumulateBndIntegral( f_oblique);

    vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
    vtkwriter->Register( make_VTKVector( make_P2Eval( MG, Stokes.GetBndData().Vel, f_oblique),  "f_oblique"));
    vtkwriter->Register( make_VTKVector( make_P2Eval( MG, Stokes.GetBndData().Vel, f_improved), "f_improved"));
    BndDataCL<> nobnd( 0);
    VecDescCL vd_sigmaf( &lset.idx);
    LSInit( MG, vd_sigmaf, sigmaf, 0.);
    vtkwriter->Register( make_VTKScalar( make_P2Eval( MG, nobnd, vd_sigmaf),                    "sigmaf"));
    vtkwriter->Write( 0.);

    VectorCL d( f_oblique.Data - f_improved.Data);
    std::cout << "|d| = \t\t" << norm(d) << std::endl;
// std::cerr << "oblique: " << f_oblique.Data.size() << ", improved: " << f_improved.Data.size() << std::endl;
// for (size_t i= 0; i < f_oblique.Data.size(); ++i)
//     if (isnan( f_oblique.Data[i]))
//         std::cerr << "i: " << i << std::endl;
// for (size_t i= 0; i < f_improved.Data.size(); ++i)
//     if (isnan( f_improved.Data[i]))
//         std::cerr << "Hallo i: " << i << std::endl;
// for (size_t i= 0; i < d.size(); ++i)
//     if (isnan( d[i]))
//         std::cerr << "Hallo d i: " << i << std::endl;

    Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &Stokes.b, lset, 0.);
    MLMatrixCL MA;
    MA.LinComb( 1, Stokes.M.Data, 1, Stokes.A.Data);
    VectorCL MA_inv_d( d);
    std::cout << "Solving system with MA matrix:\t";
    SSORPcCL pc;
    typedef PCGSolverCL<SSORPcCL> PCG_SsorCL;
    PCG_SsorCL cg( pc, 1000, 1e-18);
    cg.Solve( MA, MA_inv_d, d, vidx->GetEx());
    std::cout << cg.GetIter() << " iter,\tresid = " << cg.GetResid();
    const double sup2= std::sqrt(dot( MA_inv_d, d));
    std::cout << "\n\nsup |f1(v)-f2(v)|/||v||_1 = \t\t" << sup2
              << "\n|(MA)^-1 d| = " << norm( MA_inv_d) << std::endl;
}

void Compare_Oblique_Helper (DROPS::AdapTriangCL& adap, InstatStokes2PhaseP2P1CL& Stokes, LevelsetP2CL& lset)
{
    MultiGridCL& MG= Stokes.GetMG();

    MLIdxDescCL* vidx= &Stokes.vel_idx;
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);

    MG.SizeInfo( std::cout);
    Stokes.b.SetIdx( vidx);
    Stokes.v.SetIdx( vidx);
    Stokes.A.SetIdx( vidx, vidx);
    Stokes.M.SetIdx( vidx, vidx);

    std::cout << vidx->NumUnknowns() << " velocity unknowns,\n";
    std::cout << lset.Phi.Data.size() << " levelset unknowns.\n";

    const double Vol= 4./3*M_PI*std::pow( P.get<DROPS::Point3DCL>("Exp.RadDrop")[0], 3);
    const double curv= 2./P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
    std::cout << "Volumen = " << Vol << "\tKruemmung = " << curv << "\n\n";

    VecDescCL f_improved( vidx);
    TetraAccumulatorTupleCL accus;
    accus.push_back_acquire( new SphereObliqueLaplaceBeltramiAccuCL( lset, f_improved));
    accumulate( accus, MG, vidx->TriangLevel(), vidx->GetMatchingFunction(), vidx->GetBndInfo());

    VecDescCL f_oblique( vidx);
    lset.SetSurfaceForce( SF_ObliqueLBVar);
    lset.AccumulateBndIntegral( f_oblique);

    vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
    vtkwriter->Register( make_VTKVector( make_P2Eval( MG, Stokes.GetBndData().Vel, f_oblique),  "f_oblique"));
    vtkwriter->Register( make_VTKVector( make_P2Eval( MG, Stokes.GetBndData().Vel, f_improved), "f_improved"));
    BndDataCL<> nobnd( 0);
    VecDescCL vd_sigmaf( &lset.idx);
    LSInit( MG, vd_sigmaf, sigmaf, 0.);
    vtkwriter->Register( make_VTKScalar( make_P2Eval( MG, nobnd, vd_sigmaf),                    "sigmaf"));
    vtkwriter->Write( 0.);

    VectorCL d( f_oblique.Data - f_improved.Data);
    std::cout << "|d| = \t\t" << norm(d) << std::endl;
// std::cerr << "oblique: " << f_oblique.Data.size() << ", improved: " << f_improved.Data.size() << std::endl;
// for (size_t i= 0; i < f_oblique.Data.size(); ++i)
//     if (isnan( f_oblique.Data[i]))
//         std::cerr << "i: " << i << std::endl;
// for (size_t i= 0; i < f_improved.Data.size(); ++i)
//     if (isnan( f_improved.Data[i]))
//         std::cerr << "Hallo i: " << i << std::endl;
// for (size_t i= 0; i < d.size(); ++i)
//     if (isnan( d[i]))
//         std::cerr << "Hallo d i: " << i << std::endl;

    Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &Stokes.b, lset, 0.);
    MLMatrixCL MA;
    MA.LinComb( 1, Stokes.M.Data, 1, Stokes.A.Data);
    VectorCL MA_inv_d( d);
    std::cout << "Solving system with MA matrix:\t";
    SSORPcCL pc;
    typedef PCGSolverCL<SSORPcCL> PCG_SsorCL;
    PCG_SsorCL cg( pc, 1000, 1e-18);
    cg.Solve( MA, MA_inv_d, d, vidx->GetEx());
    std::cout << cg.GetIter() << " iter,\tresid = " << cg.GetResid();
    const double sup2= std::sqrt(dot( MA_inv_d, d));
    std::cout << "\n\nsup |f1(v)-f2(v)|/||v||_1 = \t\t" << sup2
              << "\n|(MA)^-1 d| = " << norm( MA_inv_d) << std::endl;
}

} // end of namespace DROPS


/// \brief Set Default parameters here s.t. they are initialized.
/// The result can be checked when Param-list is written to the output.
void SetMissingParameters(DROPS::ParamCL& P){
    P.put_if_unset<std::string>("Exp.VolForce", "ZeroVel");
    P.put_if_unset<double>("SurfTens.ShearVisco", 0.0);
    P.put_if_unset<double>("SurfTens.DilatationalVisco", 0.0);
    P.put_if_unset<double>("Mat.DensDrop", 1);
    P.put_if_unset<double>("Mat.ViscDrop", 1);
    P.put_if_unset<double>("Mat.DensFluid", 1);
    P.put_if_unset<double>("Mat.ViscFluid", 1);
    P.put_if_unset<double>("Mat.SmoothZone", 1e-05);
    P.put_if_unset<DROPS::Point3DCL>("Exp.Gravity", DROPS::Point3DCL());
}

int main (int argc, char** argv)
{
  try
  {
    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "f_Gamma.json");
    SetMissingParameters(P);
    std::cout << P << std::endl;

    SurfTension= P.get<double>( "SurfTens.SurfTension");
    c= P.get<DROPS::Point3DCL>("Exp.PosDrop");
    r= P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

    std::auto_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P.get_child( "Domain")));
    DROPS::MultiGridCL mg( *builder);
    DROPS::AdapTriangCL adap(
        mg,
        P.get<double>("AdaptRef.Width"),
        P.get<int>("AdaptRef.CoarsestLevel"),
        P.get<int>("AdaptRef.FinestLevel")
    );
    adap.MakeInitialTriang( DistanceFct);

    DROPS::SurfaceTensionCL sf( sigmaf);
    DROPS::LsetBndDataCL lsbnd( 0);
    read_BndData( lsbnd, mg, P.get_child( "Levelset.BndData"));
    DROPS::LevelsetP2CL& lset( *DROPS::LevelsetP2CL::Create( mg, lsbnd, sf, P.get_child("Levelset")));
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    lset.Init( DistanceFct);

    DROPS::StokesBndDataCL::VelBndDataCL velbnd( 0);
    read_BndData( velbnd, mg, P.get_child( "Stokes.VelocityBndData"));
    DROPS::StokesBndDataCL::PrBndDataCL  prbnd( 0);
    read_BndData( prbnd,  mg, P.get_child( "Stokes.PressureBndData"));

    typedef DROPS::InstatStokes2PhaseP2P1CL MyStokesCL;
    MyStokesCL prob( mg, DROPS::TwoPhaseFlowCoeffCL(P), DROPS::StokesBndDataCL( velbnd, prbnd));

    vtkwriter= std::auto_ptr<DROPS::VTKOutCL>( new DROPS::VTKOutCL(
        mg,
        "DROPS data",
        1,
        P.get<std::string>("VTK.VTKDir"),
        P.get<std::string>("VTK.VTKName"),
        P.get<std::string>("VTK.TimeFileName"),
        P.get<int>("VTK.Binary"), 
        P.get<bool>("VTK.UseOnlyP1"),
        false, /* <- P2DG */
        -1,    /* <- level */
        P.get<bool>("VTK.ReUseTimeFile")
    ));
//     Compare_LaplBeltramiSF_ConstSF( prob, lsbnd);
//     Compare_Oblique_Improved( adap, prob, lset);
    Compare_Oblique_Helper( adap, prob, lset);

    delete &lset;
    return 0;
  }
  catch (DROPS::DROPSErrCL& err) { err.handle(); }
}

