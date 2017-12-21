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
#include "levelset/marking_strategy.h"
#include "levelset/levelsetmapper.h"
#include "misc/params.h"
#include "levelset/surfacetension.h"
#include "levelset/mzelle_hdr.h"
#include "misc/dynamicload.h"
#include "num/gradient_recovery.h"
#include "surfactant/ifacetransp.h"
#include <fstream>
#include <sstream>
#include <sys/resource.h>


DROPS::ParamCL P;

std::unique_ptr<DROPS::VTKOutCL> vtkwriter;

int FctCode=9;

double SurfTension;

double const_sigmaf( const DROPS::Point3DCL&, double)
{
    return -1;
}

double sigmaf( const DROPS::Point3DCL& p, double)
{
//     return SurfTension;
    return 1 + std::cos( 2.*M_PI*p[0]);
}

DROPS::SMatrixCL<3,3> sigmaf_matrix( const DROPS::Point3DCL& p, double)
{
    using namespace DROPS;
    SMatrixCL<3,3> ret;
    const double normp= p.norm();
    if (normp == 0.)
        return ret;

    const double tau= sigmaf( p, 0.);
    const SMatrixCL<3,3> P= eye<3,3>() - outer_product( p/normp, p/normp);
    SMatrixCL<3,3> M1;
    M1( 0, 0)= std::sin( 2.*M_PI*p[1]);
    M1( 0, 1)= 2.*std::sin( M_PI*p[1]);
    M1( 1, 0)= std::cos( 3.*M_PI*p[2]);
    M1( 1, 1)= std::cos( M_PI*p[2]);
    M1( 2, 2)= -1.;
    return tau*P*M1*P;
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
    return (p - c).norm() - r;
//     return inner_prod( p - c, p - c) - r*r;
}

DROPS::Point3DCL GradDistanceFct( const DROPS::Point3DCL& p, double)
{ // ball
    if ((p - c).norm() < 1e-12)
        return DROPS::Point3DCL();
    return (p - c)/(p - c).norm();
}

DROPS::SMatrixCL<3,3> dp_sphere (const DROPS::Point3DCL& x, double)
{
    const double normx= x.norm();
    return normx == 0. ? DROPS::SMatrixCL<3,3>() : r/normx*(DROPS::eye<3,3>() - outer_product( x/normx, x/normx));
}

double h_iface;

// f(x)= exp(x0)*( cos(pi*x0), cos(2*pi*x1), cos(3*pi*x2))
// Df(x)= e0*f(x)^T - exp(x0)*pi*diag(sin(pi*x0), 2*sin(2*pi*x1), 3*sin(3*pi*x2))
DROPS::SMatrixCL<3,3> exp_test_function (const DROPS::Point3DCL& x, double)
{
    using std::cos;
    using std::sin;
    const DROPS::SVectorCL<3> f= std::exp( x[0])*DROPS::MakePoint3D( cos( M_PI*x[0]), cos( 2.*M_PI*x[1]), cos( 3.*M_PI*x[2]));
    DROPS::SMatrixCL<3,3> ret= DROPS::outer_product( DROPS::std_basis<3>( 1), f);
    for (int i= 0; i < 3; ++i)
        ret( i, i)-= std::exp(x[0])*(i + 1)*M_PI*sin( (i + 1.)*M_PI*x[i]);
    return ret;
}

// double DistanceFct( const DROPS::Point3DCL& p, double)
// { // plane... XXX The line integral on \pa\iface must be implemented.
//     return p[2] - 1./3.;
// }

DROPS::Point3DCL Null( const DROPS::Point3DCL&, double)
{ return DROPS::Point3DCL(0.); }

static DROPS::RegisterVectorFunction regveczero("VecZero",  Null);

DROPS::Point3DCL TestFunLowerBound (const DROPS::Point3DCL& p, double)
{
//     return std::abs( p[2])*p;
//     return DROPS::MakePoint3D( std::abs( p[1])*p[1], std::abs( p[0])*p[0], std::abs(p[2])*p[2])/(r*r);
//  f(x)= exp(x0)*( cos(pi*x0), cos(2*pi*x1), cos(3*pi*x2)):
    using std::cos;
//     return std::exp( p[0])*DROPS::MakePoint3D( cos(M_PI*p[0]), cos( 2.*M_PI*p[1]), cos( 3.*M_PI*p[2]));

    DROPS::Point3DCL ret= std::exp( p[0])*DROPS::MakePoint3D( cos(M_PI*p[0]), cos( 2.*M_PI*p[1]), cos( 3.*M_PI*p[2]));
    double drel= DistanceFct( p, 0.)/h_iface;
    return (1. - drel*drel)*ret;
}

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
    const double curv= 2/P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0];
    SurfaceTensionCL sf( sigmaf, 0);

    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsbnd, sf, P.get_child("Levelset")) );

//    lset.SetSurfaceForce( SF_Const);

    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    lset.CreateNumbering( MG.GetLastLevel());

    lset.Init( DistanceFct);

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, &lset);

    VecDescCL f_Gamma( vidx), v( vidx);
    MG.SizeInfo( std::cout);
    Stokes.b.SetIdx( vidx);
    Stokes.c.SetIdx( pidx);
    Stokes.p.SetIdx( pidx);
    Stokes.v.SetIdx( vidx);

    std::cout << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cout << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cout << lset.Phi.Data.size() << " levelset unknowns.\n";

    const double Vol= 4./3*M_PI*std::pow( P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0], 3);
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
    reflsg[111]= curv*2*Vol*P.get<DROPS::Point3DCL>("Levelset.PosDrop")[0];
    reflsg[121]= curv*Vol*P.get<DROPS::Point3DCL>("Levelset.PosDrop")[1];
    reflsg[131]= curv*Vol*P.get<DROPS::Point3DCL>("Levelset.PosDrop")[2];
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

    std::cout << "\n\n" << P.get<int>("Mesh.AdaptRef.FinestLevel") << ",\t";
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
    const double curv= 2/P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0];
    SurfaceTensionCL sf( const_sigmaf, 0);
    LevelsetP2CL & lset( * LevelsetP2CL::Create( MG, lsbnd, sf, P.get_child("Levelset")) );


    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    lset.CreateNumbering( MG.GetLastLevel());

    lset.Init( DistanceFct);

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, &lset);

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

    const double Vol= 4./3*M_PI*std::pow( P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0], 3);
    std::cout << "Volumen = " << Vol << "\tKruemmung = " << curv << "\n\n";

    Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &Stokes.b, lset, 0.);

    f_LaplBeltrami.Clear( Stokes.v.t);
    lset.SetSurfaceForce( SF_ImprovedLBVar);
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

    const PrincipalLatticeCL* latp_;
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
     : f( f_Gamma), ls_( ls), latp_( &PrincipalLatticeCL::instance( 2)),
       ls_val_( latp_->vertex_size())
    {
        P2DiscCL::GetGradientsOnRef( gradref_);
    }

    /// Set the sublevel of the principal-lattice for integration.
    void set_sublevel (Uint lvl) {
        std::cout << "set_sublevel: sublevel: " << lvl << ", subdivisions: " << (1u << lvl) << ".\n";
        latp_= &PrincipalLatticeCL::instance( 1u << lvl);
        ls_val_.resize( latp_->vertex_size());
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
    evaluate_on_vertexes( ls_.GetSolution(), t, *latp_, Addr( ls_val_));
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

    p_.make_patch<MergeCutPolicyCL>( *latp_, ls_val_);
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

void MyInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, Point3DCL (*d)(const Point3DCL&, double), double t= 0.)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        if (!it->Unknowns.Exist( idx))
            continue;
        const Point3DCL& v= d( it->GetCoord(), t);
        ls.Data[it->Unknowns( idx)  ]= v[0];
        ls.Data[it->Unknowns( idx)+1]= v[1];
        ls.Data[it->Unknowns( idx)+2]= v[2];
    }

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
        if (!it->Unknowns.Exist( idx))
            continue;
        const Point3DCL& v= d( 0.5*(it->GetVertex( 0)->GetCoord() + it->GetVertex( 1)->GetCoord()), t);
        ls.Data[it->Unknowns( idx)  ]= v[0];
        ls.Data[it->Unknowns( idx)+1]= v[1];
        ls.Data[it->Unknowns( idx)+2]= v[2];
    }
    ls.t= t;
}

/// \brief Accumulator for the interfacial tension term with variable interfacial tension coefficient.
/// Uses a higher order approximation of the normal field to obtain order 2.5.
class VarObliqueLaplaceBeltrami2AccuCL : public SurfTensAccumulatorCL
{
  private:
    const SurfaceTensionCL& sf_;
    LocalNumbP2CL n_;

    const InterfaceCommonDataP2CL& cdata_;
    const DROPS_STD_UNORDERED_MAP<const TetraCL*, QuadDomain2DCL>& qmap_;

    NoBndDataCL<Point3DCL> nobndvec_;
    const VecDescCL& nt_; // Point3D-valued, h^3-approximation of the exact normal to \Gamma.

    LocalP1CL<Point3DCL> gradref_[10],
                         grad_[10];
    GridFunctionCL<Point3DCL> w_[10],
                              qnq_, // normal to quadratic interface
                              qnt_, // higher order normal to quadratic interface
                              qnl;  // normal to the linear interface
    SMatrixCL<3, 3> T_;
    GridFunctionCL<double> qsigma;

    instat_matrix_fun_ptr sigmaf_matrix;
    GridFunctionCL< SMatrixCL<3,3> > qsigma_matrix;

    double area;

    bool use_mapped_fe_space_;
    bool use_linear_subsampling_;

    double test_result;
    instat_matrix_fun_ptr test_fun;
    GridFunctionCL< SMatrixCL<3,3> > qtest_fun;

    void visit_mapped_P2 (const TetraCL&);
    void visit_P2 (const TetraCL&);
    void visit_P2_fun (const TetraCL&);

  public:
    VarObliqueLaplaceBeltrami2AccuCL (const LevelsetP2CL& ls, VecDescCL& f_Gamma, const SurfaceTensionCL& sf, const InterfaceCommonDataP2CL& cdata, const VecDescCL& nt, const DROPS_STD_UNORDERED_MAP<const TetraCL*, QuadDomain2DCL>& qmap, bool use_mapped_fe_space)
     : SurfTensAccumulatorCL( ls, f_Gamma), sf_( sf), cdata_( cdata), qmap_( qmap), nt_( nt), sigmaf_matrix( 0), use_mapped_fe_space_( use_mapped_fe_space), use_linear_subsampling_( false), test_fun( 0)
    {
        P2DiscCL::GetGradientsOnRef( gradref_);
    }

    void set_test_function (instat_matrix_fun_ptr f) { test_fun= f; }
    void use_linear_subsampling (bool use) { use_linear_subsampling_= use; }
    void set_matrix_tension (instat_matrix_fun_ptr p) { sigmaf_matrix= p; }

    void begin_accumulation () {
        f.Data= 0.;
        area= 0.;
        test_result= 0.;
        std::cerr << "VarObliqueLaplaceBeltrami2AccuCL::begin_accumulation: " << f.Data.size() << "dof.\n";
    }
    void finalize_accumulation() {
        std::cerr << "VarObliqueLaplaceBeltrami2AccuCL::finalize_accumulation.\n";
        std::cerr << "    area: " << area << ".\n";
        if (test_fun) {
            const double exact= -1.796234565533891; // Maple: 16 digits.
            std::cout.precision( 15);
            std::cout << "    test_result: " << test_result << " error: " << test_result - exact << ".\n";
        }
    }

    void visit (const TetraCL& t) {
        if (test_fun)
            visit_P2_fun( t);
        else if (use_mapped_fe_space_ || use_linear_subsampling_)
            visit_mapped_P2( t);
        else
            visit_P2( t);
    }

    TetraAccumulatorCL* clone (int /*tid*/) { return new VarObliqueLaplaceBeltrami2AccuCL ( *this); };
};

void VarObliqueLaplaceBeltrami2AccuCL::visit_mapped_P2 (const TetraCL& t)
{
    const InterfaceCommonDataP2CL& cdata= cdata_.get_clone();
    if (cdata.empty())
        return;

    double det; // dummy
    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( grad_, cdata_.gradrefp2, T);
    LocalP1CL<Point3DCL> nq; // gradient of quadratic level set function
    for (Uint i= 0; i < 10 ; ++i)
        nq+= grad_[i]*cdata.locp2_ls[i];

    const QuadDomain2DCL&       qdom=           cdata.qdom;
    const TetraBaryPairVectorT& qdom_projected= cdata.qdom_projected.vertexes();

    // resize_and_evaluate_on_vertexes( nq, qdom_projected, qnq_);
    qnq_.resize( qdom.vertex_size());
    for (Uint i= 0; i < qdom.vertex_size(); ++i) {
        if (qdom_projected[i].first == &t)
            qnq_[i]= nq( qdom_projected[i].second);
        else {
            GetTrafoTr( T, det, *qdom_projected[i].first);
            LocalP1CL<Point3DCL> mygrad[10];
            P2DiscCL::GetGradients( mygrad, cdata_.gradrefp2, T);
            LocalP2CL<> mylocp2_ls( cdata.get_local_p2_ls( *qdom_projected[i].first));
            LocalP1CL<Point3DCL> mynq; // gradient of quadratic level set function
            for (Uint j= 0; j < 10 ; ++j)
                mynq+= mygrad[j]*mylocp2_ls[j];
            qnq_[i]= mynq( qdom_projected[i].second);
        }
    }

    for (Uint i= 0; i < qdom.vertex_size(); ++i)
        qnq_[i]/= norm( qnq_[i]);

    P2EvalCL<Point3DCL, const NoBndDataCL<Point3DCL>, const VecDescCL> nteval( &nt_, &nobndvec_, 0);
    resize_and_evaluate_on_vertexes( nteval, qdom_projected, qnt_);
    for (Uint i= 0; i < qdom.vertex_size(); ++i)
        qnt_[i]/= norm( qnt_[i]);

    GridFunctionCL<> qalpha( dot( qnq_, qnt_));
    // ??? If a triangle has zero area, its normal is returned as 0; we avoid
    // division by zero... the value will not matter later as the func-det in quad_2D is also 0.
    for (size_t i= 0; i < qalpha.size(); ++i)
        if (std::fabs( qalpha[i]) == 0.)
            qalpha[i]= 1.;

    GridFunctionCL<Point3DCL> qgradi( qdom.vertex_size());
    if (use_linear_subsampling_) {
        for (int i= 0; i < 10; ++i) {
            evaluate_on_vertexes ( grad_[i], qdom, Addr( qgradi));
            w_[i].resize( qdom.vertex_size());
            w_[i]= qgradi - GridFunctionCL<double>( dot( qnt_, qgradi)/qalpha)*qnq_; // \bQt D b^i, where b^i is the ith scalar P2-basis-function.
        }
    }
    else {
        if (cdata.surf.normal_empty())
            cdata.surf.compute_normals( t);
        GridFunctionCL<SMatrixCL<3,3> > Winv( cdata.qdom.vertex_size());
        QRDecompCL<3,3> qr;
        SVectorCL<3> tmp;
        for (Uint i= 0; i < cdata.qdom.vertex_size(); ++i) {
            gradient_trafo( t, cdata.qdom.vertex_begin()[i], cdata.quaqua, cdata.surf, qr.GetMatrix());
            qr.prepare_solve();
            for (Uint j= 0; j < 3; ++j) {
                tmp= std_basis<3>( j + 1);
                qr.Solve( tmp);
                Winv[i].col( j, tmp);
            }
        }
        resize_and_scatter_piecewise_normal( cdata.surf, cdata.qdom, qnl);
        for (int i= 0; i < 10; ++i) {
            evaluate_on_vertexes ( grad_[i], cdata.qdom, Addr( qgradi));
            for (Uint j= 0; j < qgradi.size(); ++j) {
                tmp=  qgradi[j] - inner_prod( qnl[j], qgradi[j])*qnl[j];
                qgradi[j]= Winv[j]*tmp;
            }
            w_[i].resize( qdom.vertex_size());
            w_[i]= qgradi - GridFunctionCL<double>( dot( qnt_, qgradi)/qalpha)*qnq_; // \bQt D b^i, where b^i is the ith scalar P2-basis-function.
        }
    }

    if (!sigmaf_matrix) {
        qsigma.resize( qdom.vertex_size());
        evaluate_on_vertexes( sf_.GetSigma(), qdom_projected, /*time*/ 0., Addr( qsigma)); // interfacial tension
    }
    else {
        qsigma_matrix.resize( qdom.vertex_size());
        evaluate_on_vertexes( sigmaf_matrix, qdom_projected, /*time*/ 0., Addr( qsigma_matrix)); // interfacial tension
    }
    if (!use_linear_subsampling_) {
        if (!sigmaf_matrix)
            qsigma*= cdata.qdom_projected.absdets();
        else
            for (Uint i= 0; i < qsigma_matrix.size(); ++i)
                qsigma_matrix[i]*= cdata.qdom_projected.absdets()[i];
    }

    n_.assign_indices_only( t, *f.RowIdx);

    if (!sigmaf_matrix)
        for (Uint i= 0; i < 10; ++i)
            add_to_global_vector( f.Data, -quad_2D( qsigma*w_[i], qdom), n_.num[i]);
    else
        for (Uint i= 0; i < 10; ++i)
            add_to_global_vector( f.Data, -quad_2D( qsigma_matrix*w_[i], qdom), n_.num[i]);
}

void VarObliqueLaplaceBeltrami2AccuCL::visit_P2 (const TetraCL& t)
{
    if (qmap_.count( &t) == 0) // Consider only tetras which meet the piecewise quadratic interface.
        return;

    const InterfaceCommonDataP2CL& cdata= cdata_.get_clone();
    const QuadDomain2DCL& qdom= qmap_.find( &t)->second;

std::valarray<double> ones( 1., qdom.vertex_size());
area+= quad_2D( ones, qdom);

    double det; // dummy
    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( grad_, cdata_.gradrefp2, T);
    LocalP1CL<Point3DCL> nq; // gradient of quadratic level set function
    for (Uint i= 0; i < 10 ; ++i)
        nq+= grad_[i]*cdata.locp2_ls[i];

    resize_and_evaluate_on_vertexes( nq, qdom, qnq_);
    for (Uint i= 0; i < qdom.vertex_size(); ++i)
        qnq_[i]/= norm( qnq_[i]);

    P2EvalCL<Point3DCL, const NoBndDataCL<Point3DCL>, const VecDescCL> nteval( &nt_, &nobndvec_, 0);
    resize_and_evaluate_on_vertexes( nteval, t, qdom, qnt_);
    for (Uint i= 0; i < qdom.vertex_size(); ++i)
        qnt_[i]/= norm( qnt_[i]);

    GridFunctionCL<> qalpha( dot( qnq_, qnt_));
    // ??? If a triangle has zero area, its normal is returned as 0; we avoid
    // division by zero... the value will not matter later as the func-det in quad_2D is also 0.
    for (size_t i= 0; i < qalpha.size(); ++i)
        if (std::fabs( qalpha[i]) == 0.)
            qalpha[i]= 1.;
    GridFunctionCL<Point3DCL> qgradi( qdom.vertex_size());
    for (Uint i= 0; i < 10; ++i) {
        evaluate_on_vertexes( grad_[i], qdom, Addr( qgradi));
        w_[i].resize( qdom.vertex_size());
        w_[i]= qgradi - GridFunctionCL<double>( dot( qnt_, qgradi)/qalpha)*qnq_; // \bQt D b^i, where b^i is the ith scalar P2-basis-function.
    }
    qsigma.resize( qdom.vertex_size());
    sf_.evaluate_on_vertexes( t, qdom, /*time*/ 0., Addr( qsigma)); // interfacial tension

    n_.assign_indices_only( t, *f.RowIdx);
    for (Uint i= 0; i < 10; ++i)
        add_to_global_vector( f.Data, -quad_2D( qsigma*w_[i], qdom), n_.num[i]);
}


void VarObliqueLaplaceBeltrami2AccuCL::visit_P2_fun (const TetraCL& t)
{
    if (qmap_.count( &t) == 0) // Consider only tetras which meet the piecewise quadratic interface.
        return;

    const InterfaceCommonDataP2CL& cdata= cdata_.get_clone();
    const QuadDomain2DCL& qdom= qmap_.find( &t)->second;

std::valarray<double> ones( 1., qdom.vertex_size());
area+= quad_2D( ones, qdom);

    double det; // dummy
    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( grad_, cdata_.gradrefp2, T);
    LocalP1CL<Point3DCL> nq; // gradient of quadratic level set function
    for (Uint i= 0; i < 10 ; ++i)
        nq+= grad_[i]*cdata.locp2_ls[i];

    resize_and_evaluate_on_vertexes( nq, qdom, qnq_);
    for (Uint i= 0; i < qdom.vertex_size(); ++i)
        qnq_[i]/= norm( qnq_[i]);

    P2EvalCL<Point3DCL, const NoBndDataCL<Point3DCL>, const VecDescCL> nteval( &nt_, &nobndvec_, 0);
    resize_and_evaluate_on_vertexes( nteval, t, qdom, qnt_);
    for (Uint i= 0; i < qdom.vertex_size(); ++i)
        qnt_[i]/= norm( qnt_[i]);

    GridFunctionCL<> qalpha( dot( qnq_, qnt_));
    // ??? If a triangle has zero area, its normal is returned as 0; we avoid
    // division by zero... the value will not matter later as the func-det in quad_2D is also 0.
    for (size_t i= 0; i < qalpha.size(); ++i)
        if (std::fabs( qalpha[i]) == 0.)
            qalpha[i]= 1.;
    resize_and_evaluate_on_vertexes( test_fun, t, qdom, 0., qtest_fun);
    // \bQt Dv, where v is the (vector-valued) test-function.
    qtest_fun-=  outer_product( GridFunctionCL<>(1./qalpha)*qnq_, transp_mul( qtest_fun, qnt_));
    qsigma.resize( qdom.vertex_size());
    sf_.evaluate_on_vertexes( t, qdom, /*time*/ 0., Addr( qsigma)); // interfacial tension

    test_result+= -quad_2D( qsigma*trace( qtest_fun), qdom);
}

void AccumulateBndIntegral (LevelsetP2CL& lset, const PrincipalLatticeCL& lat, VecDescCL& f_Gamma)
{
    // SF_ObliqueLBVar2:
    ScopeTimerCL scope("AccumulateBndIntegral2");

    MultiGridCL& mg= lset.GetMG();
    const Uint lvl= lset.Phi.GetLevel();

    // Recover the gradient of the level set function
    IdxDescCL vecp2idx( vecP2_FE);
    vecp2idx.CreateNumbering( lvl, mg);
    VecDescCL nt( &vecp2idx);
    MyInit( mg, nt, &GradDistanceFct, 0.);

    // Compute neighborhoods of the tetras at the interface
    TetraToTetrasT tetra_neighborhoods;
    compute_tetra_neighborhoods( mg, lset.Phi, lset.GetBndData(), lat, tetra_neighborhoods);

    VecDescCL lsgradrec( &vecp2idx);
    averaging_P2_gradient_recovery( mg, lset.Phi, lset.GetBndData(), lsgradrec);

    QuaQuaMapperCL quaqua( mg, lset.Phi, lsgradrec, &tetra_neighborhoods,
                           P.get<int>( "LevelsetMapper.Iter"),
                           P.get<double>( "LevelsetMapper.Tol"),
                           P.get<std::string>( "LevelsetMapper.Method") == "FixedPointWithLineSearch",
                           P.get<double>( "LevelsetMapper.ArmijoConstant"));

    InterfaceCommonDataP2CL cdatap2( lset.Phi, lset.GetBndData(), quaqua, lat);
    QuaQuaQuadDomainMapperAccuCL hoqdom_accu( cdatap2);
    if (P.get<std::string>( "TestCase.ComparisonSource") == "ObliqueLBVar2") {
        TetraAccumulatorTupleCL hoaccus;
        hoaccus.push_back( &cdatap2);
        hoaccus.push_back( &hoqdom_accu);
        accumulate( hoaccus, lset.GetMG(), lvl, lset.Phi.RowIdx->GetBndInfo());
    }

    TetraAccumulatorTupleCL accus;
    accus.push_back( &cdatap2);
    VarObliqueLaplaceBeltrami2AccuCL accu( lset, f_Gamma, lset.GetSF(), cdatap2, nt, hoqdom_accu.qmap, P.get<bool>( "NavStokes.Coeff.SurfTens.UseMappedFESpace"));
    if (P.get<std::string>( "NavStokes.Coeff.SurfTens.TestFunction") == "exp_test_function")
        accu.set_test_function( &exp_test_function);
    if (P.get<bool>( "NavStokes.Coeff.SurfTens.UseMatrixTension") == true)
        accu.set_matrix_tension( &sigmaf_matrix);
    if (P.get<std::string>( "TestCase.ComparisonSource") == "ObliqueLBVar3") {
        const Uint subsampling= std::ceil( P.get<double>( "TestCase.SubsamplingFactor")*std::pow( 2., P.get<double>( "TestCase.SubsamplingExponent")*lvl));
//         const Uint subsampling= 1u << (lvl == 0 ? 0 : lvl - 1);
        std::cout << "ObliqueLBVar3: subsampling: " << subsampling << ".\n";
        accu.use_linear_subsampling( true);
        cdatap2.set_lattice( PrincipalLatticeCL::instance( subsampling));
        cdatap2.compute_absdet( false);
    }
    accus.push_back( &accu);
    accumulate( accus, lset.GetMG(), lvl, lset.Phi.RowIdx->GetBndInfo());
    if (P.get<std::string>( "TestCase.ComparisonSource") == "ObliqueLBVar3") {
        cdatap2.set_lattice( PrincipalLatticeCL::instance( 1)); // XXX: Enhancement: Computation of absdets now also works for general lattices.
        cdatap2.compute_absdet( true);
    }

//     InterfaceDebugP2CL debug_accu( cdatap2);
//     debug_accu.set_true_area( 4.*M_PI*r*r);
//     debug_accu.set_ref_dp( &dp_sphere);
// //         p2debugaccu.set_ref_abs_det( &abs_det_sphere);
//     TetraAccumulatorTupleCL debug_accus;
//     debug_accus.push_back( &cdatap2);
//     debug_accus.push_back( &debug_accu);
//     accumulate( debug_accus, lset.GetMG(), lvl, lset.Phi.RowIdx->GetMatchingFunction(), lset.Phi.RowIdx->GetBndInfo());
}

void Compare_Oblique_Coarse (DROPS::AdapTriangCL&, InstatStokes2PhaseP2P1CL& Stokes, LevelsetP2CL& lset, std::string comparison_source)
{
    MultiGridCL& mg= Stokes.GetMG();
    if (mg.GetLastLevel() == 0)
        throw DROPSErrCL( "Compare_Oblique: Need level > 0.\n");

    const Uint flvl= mg.GetLastLevel();

    const double Vol= 4./3*M_PI*std::pow( r, 3);
    const double curv= 2./r;
    std::cout << "Volumen = " << Vol << "\tKruemmung = " << curv << "\n\n";

    MLIdxDescCL mlidx( vecP2_FE, flvl+1, Stokes.GetBndData().Vel);
    mlidx.CreateNumbering( flvl, mg);
    MLVecDescCL ff_oblique( &mlidx);
    std::cout << "ff_oblique: " << ff_oblique.size() << " levels.\n";
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 1);

    if (comparison_source == "ObliqueLBVar") {
        lset.SetSurfaceForce( SF_ObliqueLBVar);
        lset.AccumulateBndIntegral( ff_oblique.GetFinest());
    }
    else if (comparison_source == "ObliqueLBVar2" || comparison_source == "ObliqueLBVar3")
        AccumulateBndIntegral( lset, lat, ff_oblique.GetFinest());
    else
        throw DROPSErrCL( "Compare_Oblique_Coarse: Unknown source: " + comparison_source + ".\n");
    WriteFEToFile( ff_oblique.GetFinest(), mg, "ff_oblique.txt");

    {
        MLDataCL<ProlongationCL<Point3DCL> > mlprolongation;
        SetupProlongationMatrix( mg, mlprolongation, &mlidx, &mlidx);
        std::cout << "mlprolongation: " << mlprolongation.size() << " levels.\n";
        MLVecDescCL::iterator vdit= ff_oblique.GetFinestIter();
        MLDataCL<ProlongationCL<Point3DCL> >::const_iterator pit= mlprolongation.GetFinestIter();
        for (Uint l= flvl; l > 0; --vdit, --pit, --l) {
            --vdit;
            VectorCL& vc= vdit->Data;
            ++vdit;
            vc= transp_mul( *pit, vdit->Data);
        }
    }

    DROPS::SurfaceTensionCL csf( sigmaf);
    std::unique_ptr<DROPS::LevelsetP2CL> clsetp( DROPS::LevelsetP2CL::Create( mg, lset.GetBndData(), csf, P.get_child("Levelset")));
    DROPS::LevelsetP2CL& clset= *clsetp;
    clset.SetSurfaceForce( SF_ObliqueLBVar);
    MLVecDescCL::const_iterator ffit= ff_oblique.GetCoarsestIter();
    for (Uint l= 0; l < flvl; ++l, ++ffit) {
        if (l > 0) {
            clset.ClearMat();
            clset.DeleteNumbering( &clset.idx);
            Stokes.ClearMat();
            Stokes.DeleteNumbering( &Stokes.vel_idx);
        }
        clset.CreateNumbering( l, &clset.idx);
        std::cout << clset.idx.NumUnknowns() << " levelset unknowns.\n";
        clset.Phi.SetIdx( &clset.idx);
        clset.Init( DistanceFct);

        MLIdxDescCL* vidx= &Stokes.vel_idx;
        Stokes.CreateNumberingVel( l, vidx);
        std::cout << vidx->NumUnknowns() << " velocity unknowns,\n";
        Stokes.b.SetIdx( vidx);
        Stokes.v.SetIdx( vidx);
        Stokes.A.SetIdx( vidx, vidx);
        Stokes.M.SetIdx( vidx, vidx);
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &Stokes.b, clset, 0.);
        MLMatrixCL MA;
        MA.LinComb( 1, Stokes.M.Data, 1, Stokes.A.Data);

        VecDescCL oblique( vidx);
        if (comparison_source == "ObliqueLBVar")
            clset.AccumulateBndIntegral( oblique);
        else if (comparison_source == "ObliqueLBVar2" || comparison_source == "ObliqueLBVar3")
            AccumulateBndIntegral( clset, lat, oblique);
        else
            throw DROPSErrCL( "Compare_Oblique_Coarse: Unknown source: " + comparison_source + " (2).\n");
        VectorCL d( oblique.Data - ffit->Data);
        std::cout << "|d| = \t\t" << norm(d) << std::endl;
        VectorCL MA_inv_d( d);
        std::cout << "Solving system with MA matrix:\t";
        SSORPcCL pc;
        typedef PCGSolverCL<SSORPcCL> PCG_SsorCL;
//         PCG_SsorCL cg( pc, 1000, 1e-14);
        PCG_SsorCL cg( pc, 1000, 1e-7);
        cg.Solve( MA, MA_inv_d, d, vidx->GetEx());
        std::cout << cg.GetIter() << " iter,\tresid = " << cg.GetResid();
        const double sup2= std::sqrt(dot( MA_inv_d, d));
        std::cout.precision( 12);
        std::cout << "\n\nsup |f1(v)-f2(v)|/||v||_1 = \t\t" << sup2
                  << "\n|(MA)^-1 d| = " << norm( MA_inv_d) << std::endl;
    }
}

void Compare_Oblique (DROPS::AdapTriangCL&, InstatStokes2PhaseP2P1CL& Stokes, LevelsetP2CL& lset, std::string comparison_source, std::string comparison_target)
{
    MultiGridCL& mg= Stokes.GetMG();

    MLIdxDescCL* vidx= &Stokes.vel_idx;
    Stokes.CreateNumberingVel( mg.GetLastLevel(), vidx);

    mg.SizeInfo( std::cout);
    Stokes.b.SetIdx( vidx);
    Stokes.v.SetIdx( vidx);
    Stokes.A.SetIdx( vidx, vidx);
    Stokes.M.SetIdx( vidx, vidx);

    std::cout << vidx->NumUnknowns() << " velocity unknowns,\n";
    std::cout << lset.Phi.Data.size() << " levelset unknowns.\n";

    const double Vol= 4./3*M_PI*std::pow( r, 3);
    const double curv= 2./r;
    std::cout << "Volumen = " << Vol << "\tKruemmung = " << curv << "\n\n";

    VecDescCL f_oblique( vidx);
    if (comparison_source == "ObliqueLBVar") {
        lset.SetSurfaceForce( SF_ObliqueLBVar);
        lset.AccumulateBndIntegral( f_oblique);
    }
    else if (comparison_source == "ObliqueLBVar2" || comparison_source == "ObliqueLBVar3") {
        const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 1);
        AccumulateBndIntegral( lset, lat, f_oblique);
    }
    else
        throw DROPSErrCL( "Compare_Oblique: Unknown source: " + comparison_source + ".\n");

    VecDescCL f_improved( vidx);
    VecDescCL funvd( vidx);
    if (comparison_target == "Helper" || comparison_target == "VariableHelper") {
        TetraAccumulatorTupleCL accus;
        SphereObliqueLaplaceBeltramiAccuCL accu( lset, f_improved);
        accus.push_back( &accu);
        if (comparison_target == "VariableHelper")
            accu.set_sublevel( mg.GetLastLevel());
        accumulate( accus, mg, vidx->TriangLevel(), vidx->GetBndInfo());
    }
    else if (comparison_target == "Improved") {
        lset.SetSurfaceForce( SF_ImprovedLBVar);
        lset.AccumulateBndIntegral( f_improved);
    }
    else if (comparison_target == "LowerBound") {
        MyInit( mg, funvd, TestFunLowerBound);
//         const double exact= -4.*M_PI*std::pow( r, 3);
//         const double exact= -2.*M_PI*r;
        const double exact= -1.796234565533891; // Maple: 16 digits.
        std::cout << "LowerBound dot( f_oblique, funvd): " << dot( f_oblique.Data, funvd.Data) << "\n"
                     "fh(v) - f(v): " << dot( f_oblique.Data, funvd.Data) - exact << ".\n";
        return;
    }
    else if (comparison_target == "Norm") {
        f_improved.Data= 0.;
    }
    else
        throw DROPSErrCL( "Compare_Oblique: Unknown target: " + comparison_target + ".\n");

    vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "Levelset") );
    vtkwriter->Register( make_VTKVector( make_P2Eval( mg, Stokes.GetBndData().Vel, f_oblique),  "f_oblique"));
    vtkwriter->Register( make_VTKVector( make_P2Eval( mg, Stokes.GetBndData().Vel, f_improved), "f_improved"));
    BndDataCL<> nobnd( 0);
    VecDescCL vd_sigmaf( &lset.idx);
    LSInit( mg, vd_sigmaf, sigmaf, 0.);
    vtkwriter->Register( make_VTKScalar( make_P2Eval( mg, nobnd, vd_sigmaf),                    "sigmaf"));
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
    if (comparison_target == "LowerBound") {
        const double h1norm= std::sqrt( dot( funvd.Data, MA*funvd.Data));
        std::cout << "|f(v)|/||v||_1 = \t\t" << dot( f_oblique.Data, funvd.Data)/h1norm << ".\n";
    }
    VectorCL MA_inv_d( d);
    std::cout << "Solving system with MA matrix:\t";
    SSORPcCL pc;
    typedef PCGSolverCL<SSORPcCL> PCG_SsorCL;
//     PCG_SsorCL cg( pc, 1000, 1e-18);
    PCG_SsorCL cg( pc, 1000, 1e-14);
    cg.Solve( MA, MA_inv_d, d, vidx->GetEx());
    std::cout << cg.GetIter() << " iter,\tresid = " << cg.GetResid();
    const double sup2= std::sqrt(dot( MA_inv_d, d));
    std::cout.precision( 12);
    std::cout << "\n\nsup |f1(v)-f2(v)|/||v||_1 = \t\t" << sup2
              << "\n|(MA)^-1 d| = " << norm( MA_inv_d) << std::endl;
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    ScopeTimerCL scope("main");
    std::cout.precision( 15);


    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/tests/f_Gamma/f_Gamma.json");
    std::cout << P << std::endl;

    SurfTension= P.get<double>( "NavStokes.Coeff.SurfTens.SurfTension");
    c= P.get<DROPS::Point3DCL>("Levelset.PosDrop");
    r= P.get<DROPS::Point3DCL>("Levelset.RadDrop")[0];
    h_iface= P.get<DROPS::Point3DCL>( "Mesh.E1")[0]/P.get<double>( "Mesh.N1")/std::pow(2., P.get<DROPS::Uint>( "Mesh.AdaptRef.FinestLevel"));

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"), P.get<std::vector<std::string> >("General.DynamicLibs") );

    std::unique_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P));
    DROPS::MultiGridCL mg( *builder);
    typedef DROPS::DistMarkingStrategyCL MarkerT;
    MarkerT InitialMarker( DistanceFct,
                           P.get<double>("Mesh.AdaptRef.Width"),
                           P.get<double>("Mesh.AdaptRef.CoarsestLevel"), P.get<double>("Mesh.AdaptRef.FinestLevel") );

    DROPS::AdapTriangCL adap( mg, &InitialMarker);
    adap.MakeInitialTriang();

    DROPS::SurfaceTensionCL sf( sigmaf);
    DROPS::LsetBndDataCL lsbnd( 0);
    read_BndData( lsbnd, mg, P.get_child( "Levelset.BndData"));
    std::unique_ptr<DROPS::LevelsetP2CL> lsetp( DROPS::LevelsetP2CL::Create( mg, lsbnd, sf, P.get_child("Levelset")));
    DROPS::LevelsetP2CL& lset= *lsetp;
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    lset.Init( DistanceFct);

    DROPS::StokesBndDataCL::VelBndDataCL velbnd( 0);
    DROPS::read_BndData( velbnd, mg, P.get_child( "NavStokes.BoundaryData.Velocity"));
    DROPS::StokesBndDataCL::PrBndDataCL  prbnd( 0);
    DROPS::read_BndData( prbnd,  mg, P.get_child( "NavStokes.BoundaryData.Pressure"));

    typedef DROPS::InstatStokes2PhaseP2P1CL MyStokesCL;
    MyStokesCL prob( mg, DROPS::TwoPhaseFlowCoeffCL(P), DROPS::StokesBndDataCL( velbnd, prbnd));

    vtkwriter= std::unique_ptr<DROPS::VTKOutCL>( new DROPS::VTKOutCL(
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
    Compare_LaplBeltramiSF_ConstSF( prob, lsbnd);
    if (P.get<std::string>( "TestCase.ComparisonTarget") == "CoarseLevel")
        Compare_Oblique_Coarse( adap, prob, lset, P.get<std::string>( "TestCase.ComparisonSource"));
    else
        Compare_Oblique( adap, prob, lset, P.get<std::string>( "TestCase.ComparisonSource"), P.get<std::string>( "TestCase.ComparisonTarget"));

    rusage usage;
    getrusage( RUSAGE_SELF, &usage);
    printf( "ru_maxrss: %li kB.\n", usage.ru_maxrss);
    return 0;
  }
  catch (DROPS::DROPSErrCL& err) { err.handle(); }
}

