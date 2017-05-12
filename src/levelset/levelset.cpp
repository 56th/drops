/// \file levelset.cpp
/// \brief levelset equation for two phase flow problems
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

#include "levelset/levelset.h"
#include "levelset/volume_adjustment.h"
#include "num/krylovsolver.h"
#include "num/precond.h"
#include "levelset/fastmarch.h"
#include "num/lattice-eval.h"
#include "num/quadrature.h"
#include <fstream>

namespace DROPS
{

inline double SmoothedSign( double x, double alpha)
{
    return x/std::sqrt(x*x+alpha);
}

/// \brief Base class for all surface tension accumulators
class SurfTensAccumulatorCL : public TetraAccumulatorCL
{
  protected:
    VecDescCL  SmPhi_;
    const BndDataCL<>& lsetbnd_;
    VecDescCL& f;
    SMatrixCL<3,3> T;
    InterfaceTriangleCL triangle;

  public:
    SurfTensAccumulatorCL( const LevelsetP2CL& ls, VecDescCL& f_Gamma)
     : SmPhi_(*ls.PhiC), lsetbnd_(ls.GetBndData()), f(f_Gamma)
    { ls.MaybeSmooth( SmPhi_.Data); }

    void begin_accumulation ()
    {
        // uncomment for Geomview output
        //std::ofstream fil("surf.off");
        //fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";
    }
    void finalize_accumulation()
    {
        // uncomment for Geomview output
        //fil << "}\n";
    }
    ///\brief Do setup of f_Gamma on given tetra
    virtual void visit (const TetraCL&)= 0;

};

/// \brief Accumulator for the (artificial) constant surface force, mainly used for numerical test cases.
///
/// Computes the integral
///         \f[ \sigma \int_\Gamma v \textbf n ds. \f]
class ConstSurfTensAccumulatorCL : public SurfTensAccumulatorCL
{
  private:
    const double sigma_;
    const RefRuleCL RegRef_;

    IdxT Numb[14];
    LocalP2CL<> p1abs_p[4][8], p1abs_n[4][8]; // extended basis functions on pos./neg. part, resp., for each of the 8 regular children
    LocalP2CL<> loc_phi;

  public:
    ConstSurfTensAccumulatorCL( const LevelsetP2CL& ls, VecDescCL& f_Gamma, double sigma)
     : SurfTensAccumulatorCL( ls, f_Gamma), sigma_(sigma), RegRef_(GetRefRule( RegRefRuleC)) {}

    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new ConstSurfTensAccumulatorCL ( *this); };
};

void ConstSurfTensAccumulatorCL::visit (const TetraCL& t)
{
    const Uint idx_f= f.RowIdx->GetIdx();
    const bool velXfem= f.RowIdx->IsExtended();

    loc_phi.assign( t, SmPhi_, lsetbnd_);
    triangle.Init( t, loc_phi);
    if (!triangle.Intersects())
        return;

    for (int v=0; v<10; ++v)
    { // collect data on all DoF
        const UnknownHandleCL& unk= v<4 ? t.GetVertex(v)->Unknowns : t.GetEdge(v-4)->Unknowns;
        Numb[v]= unk.Exist(idx_f) ? unk(idx_f) : NoIdx;
    }
    for (int xv=0; xv<4; ++xv)
    { // collect data on all extended DoF
        Numb[xv+10]= velXfem && Numb[xv]!=NoIdx ? f.RowIdx->GetXidx()[Numb[xv]] : NoIdx;
    }

    if (velXfem)
        P2RidgeDiscCL::GetExtBasisOnChildren( p1abs_p, p1abs_n, loc_phi);

    for (int ch=0; ch<8; ++ch)
    {
        if (!triangle.ComputeForChild(ch)) // no patch for this child
            continue;

//patch.WriteGeom( fil);

        BaryCoordCL BaryPQR, BarySQR;
        for (int i=0; i<3; ++i)
        {
            // addiere baryzentrische Koordinaten von P,Q,R bzw. S,Q,R
            BaryPQR+= triangle.GetBary(i);
            BarySQR+= triangle.GetBary(i+1);
        }
        BaryPQR/= 3; BarySQR/= 3;

        Point3DCL n; // senkrecht auf PQ, PR, nach aussen gerichtet...
        cross_product( n, triangle.GetPoint(1)-triangle.GetPoint(0), triangle.GetPoint(2)-triangle.GetPoint(0));
        n/= n.norm();

        const ChildDataCL data= GetChildData( RegRef_.Children[ch]);
        const int find_sign= triangle.GetNumSign( 1) ? 1 : -1;
        Point3DCL pos_dir;
        for (Uint i=0; i<4; ++i) // compute vector with positive direction
        {
            const Uint vert= data.Vertices[i];
            if (triangle.GetSign(vert)==find_sign)
            {
                const Point3DCL signedPoint= vert<4 ? t.GetVertex(vert)->GetCoord() : GetBaryCenter( *t.GetEdge(vert-4));
                pos_dir= signedPoint - triangle.GetPoint(0);
                if (find_sign == -1) pos_dir= -pos_dir;
                break;
            }
        }
        if (inner_prod( n, pos_dir) < 0) n= -n;

        double val_hat[4];
        for (int v=0; v<14; ++v)
        {
            if (Numb[v]==NoIdx) continue;

            for (Uint k=0; k<triangle.GetNumPoints(); ++k)
                // values of basis function in P,Q,R,S. Note: p1abs_p==p1abs_n on \f$Gamma_h\f$
                val_hat[k]= v<10 ? FE_P2CL::H(v,triangle.GetBary(k)) : p1abs_p[v-10][ch](triangle.GetBary(k));

            double v_Bary= v<10 ? FE_P2CL::H(v,BaryPQR) : p1abs_p[v-10][ch](BaryPQR),
                sum_v= 0;
            for (int k=0; k<3; ++k)
                 sum_v+= val_hat[k];

            if (triangle.IsQuadrilateral())
            {
                double sum_vSQR= 0;
                for (int k=1; k<4; ++k)
                    sum_vSQR+= val_hat[k];
                sum_v+= triangle.GetAreaFrac() * sum_vSQR;
                v_Bary+= triangle.GetAreaFrac() * (v<10 ? FE_P2CL::H(v,BarySQR) : p1abs_p[v-10][ch](BarySQR));
            }

            // Quadraturformel auf Dreieck, exakt bis zum Grad 2
            const double int_v= (1./12)*sum_v + 0.75 * v_Bary;
            const double C= sigma_*triangle.GetAbsDet()/2;
            for (int i=0; i<3; ++i)
                f.Data[Numb[v]+i]-= C * int_v*n[i];
        }
    } // Ende der for-Schleife ueber die Kinder
}


/// \brief Accumulator for the naive Laplace-Beltrami discretization of the CSF term.
///
/// Computes the integral
///         \f[ \sigma \int_\Gamma \kappa v \textbf n_h ds = \sigma \int_\Gamma P_h\nabla id P_h\nabla v ds \f]
/// with \f$P_h = I - \textbf n_h\cdot \textbf n_h^T\f$.
/// Discretization error in \f$ H^1(\Omega)\f$ dual norm has only order 1/2 w.r.t. the grid size at the interface.
/// Better use the improved Laplace-Beltrami discretization with better discretization order.
class NaiveLaplaceBeltramiAccuCL : public SurfTensAccumulatorCL
{
  private:
    const double sigma_;

    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    IdxT Numb[10];

  public:
    NaiveLaplaceBeltramiAccuCL( const LevelsetP2CL& ls, VecDescCL& f_Gamma, double sigma)
     : SurfTensAccumulatorCL( ls, f_Gamma), sigma_(sigma)
    { P2DiscCL::GetGradientsOnRef( GradRef); }

    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new NaiveLaplaceBeltramiAccuCL ( *this); };
};

void NaiveLaplaceBeltramiAccuCL::visit( const TetraCL& t)
// computes the integral
//         sigma \int_\Gamma \kappa v n ds = sigma \int_\Gamma grad id grad v ds
{
    const Uint  idx_f= f.RowIdx->GetIdx();
    double det;

    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder

    for (int v=0; v<10; ++v)
    { // collect data on all DoF
        const UnknownHandleCL& unk= v<4 ? t.GetVertex(v)->Unknowns : t.GetEdge(v-4)->Unknowns;
        Numb[v]= unk.Exist(idx_f) ? unk(idx_f) : NoIdx;
    }

    triangle.Init( t, SmPhi_, lsetbnd_);

    for (int ch=0; ch<8; ++ch)
    {
        if (!triangle.ComputeForChild(ch)) // no patch for this child
            continue;

//patch.WriteGeom( fil);

        BaryCoordCL BaryPQR, BarySQR;
        for (int i=0; i<3; ++i)
        {
            // addiere baryzentrische Koordinaten von P,Q,R bzw. S,Q,R
            BaryPQR+= triangle.GetBary(i);
            BarySQR+= triangle.GetBary(i+1);
        }

        const double C= triangle.GetAbsDet()/6*sigma_;  // 1/6 for quad. rule

        for (int v=0; v<10; ++v)
        {
            if (Numb[v]==NoIdx) continue;

            LocalP1CL<Point3DCL> gradv; // gradv = Werte von grad Hutfunktion fuer DoF v in den vier vertices
            for (int node=0; node<4; ++node)
                gradv[node]= Grad[v][node];

            Point3DCL gr= gradv( BaryPQR); // gr= grad v(P) + grad v(Q) + grad v(R)

            if (triangle.IsQuadrilateral())
                gr+= triangle.GetAreaFrac() * gradv( BarySQR);
            // nun gilt:
            // gr = [grad v(P)+...] + (a+b-1)[grad v(S)+...]

            for (int i=0; i<3; ++i)
            {
                const double val= inner_prod( gr, triangle.GetGradId(i));
                f.Data[Numb[v]+i]-= C *val;
            }
        }
    } // Ende der for-Schleife ueber die Kinder
}

void SF_ImprovedLaplBeltramiOnTriangle( const TetraCL& t, const BaryCoordCL * const p,
                                        const InterfaceTriangleCL&  triangle, const LocalP1CL<Point3DCL> Grad_f[10], const IdxT Numb[10],
                                        instat_scalar_fun_ptr sigma, const Quad5_2DCL<Point3DCL> e[3],
                                        double det, VectorCL& f)
{
    Quad5_2DCL<Point3DCL> Grad[10]; // Gradients of the P2-basis-functions
    Quad5_2DCL<Point3DCL> n;

    for (int v=0; v<10; ++v)
    {
        Grad[v].assign( Grad_f[v], p);
        n+= triangle.GetPhi(v)*Grad[v];
    }
    for (int i =0; i<Quad5_2DDataCL::NumNodesC; i++) if (n[i].norm()>1e-8) n[i]/= n[i].norm();

    Quad5_2DCL<> qsigma( t, p, sigma),  // surface tension
                 q1;                    // Term 1

    Quad5_2DCL<Point3DCL> qPhPhte,      // Common term in Term 1 and Term 2
                          qsigmaPhPhte; // for Term 1

    for (int i= 0; i < 3; ++i)
    {
        qPhPhte= (e[i] - dot(e[i],n)*n);
        qPhPhte.apply( triangle, &InterfaceTriangleCL::ApplyProj);
        qsigmaPhPhte= qsigma*qPhPhte;
        for (int v= 0; v < 10; ++v)
        {
            if (Numb[v]==NoIdx) continue;
            q1= dot (qsigmaPhPhte, Grad[v]);
            f[Numb[v]+i]-= q1.quad( det);
        }
    }
}

void SF_ImprovedLaplBeltramiOnTriangle( const TetraCL& t, const BaryCoordCL * const p,
    const InterfaceTriangleCL&  triangle, const LocalP1CL<Point3DCL> Grad_f[10], const IdxT Numb[10],
    const Quad5_2DCL<Point3DCL> e[3], double det, VectorCL& f, const SurfaceTensionCL& sf)
{
    Quad5_2DCL<>          p2[10];   // P2-Hat-Functions...
    Quad5_2DCL<Point3DCL> Grad[10]; // and their gradients
    Quad5_2DCL<Point3DCL> n;
    P2DiscCL::GetP2Basis( p2, p);
    for (int v=0; v<10; ++v)
    {
        Grad[v].assign( Grad_f[v], p);
        n+= triangle.GetPhi(v)*Grad[v];
    }
    for (int i =0; i<Quad5_2DDataCL::NumNodesC; i++) if (n[i].norm()>1e-8) n[i]/= n[i].norm();

    Quad5_2DCL<> qsigma, q1;
    Quad5_2DCL<Point3DCL> qPhPhte,                         // Common term in Term 1 and Term 2
                          qsigmaPhPhte;                    // for Term 1
    sf.ComputeSF(t, p, qsigma);

    for (int i= 0; i < 3; ++i)
    {
        qPhPhte= (e[i] - dot(e[i],n)*n);
        qPhPhte.apply( triangle, &InterfaceTriangleCL::ApplyProj);
        qsigmaPhPhte= qsigma*qPhPhte;
        for (int v= 0; v < 10; ++v)
        {
            if (Numb[v]==NoIdx) continue;
            q1= dot (qsigmaPhPhte, Grad[v]);
            f[Numb[v]+i]-= q1.quad( det);
        }
    }
}

/// \brief Accumulator for the improved Laplace-Beltrami discretization of the CSF term with variable surface tension coefficient.
class VarImprovedLaplaceBeltramiAccuCL : public SurfTensAccumulatorCL
{
  private:
    const SurfaceTensionCL& sf_;

    LocalP1CL<Point3DCL> GradRef[10], Grad[10];
    IdxT Numb[10];
    Quad5_2DCL<Point3DCL> e[3];

  public:
    VarImprovedLaplaceBeltramiAccuCL( const LevelsetP2CL& ls, VecDescCL& f_Gamma, const SurfaceTensionCL& sigma)
     : SurfTensAccumulatorCL( ls, f_Gamma), sf_(sigma)
    {
        P2DiscCL::GetGradientsOnRef( GradRef);
        for (int i= 0; i<3; ++i)
            e[i]= std_basis<3>( i + 1);
    }

    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new VarImprovedLaplaceBeltramiAccuCL ( *this); };
};

void VarImprovedLaplaceBeltramiAccuCL::visit( const TetraCL& t)
// computes the integral sigma \int_\Gamma \kappa v n ds = sigma
// \int_\Gamma grad id grad v ds
{
    const Uint idx_f= f.RowIdx->GetIdx();
    double det;

    triangle.Init( t, SmPhi_, lsetbnd_);

    for (int v= 0; v < 10; ++v)
    { // collect data on all DoF
        const UnknownHandleCL& unk= v<4 ? t.GetVertex(v)->Unknowns : t.GetEdge(v-4)->Unknowns;
        Numb[v]= unk.Exist( idx_f) ? unk( idx_f) : NoIdx;
    }
    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( Grad, GradRef, T);

    for (int ch= 0; ch < 8; ++ch)
    {
        triangle.ComputeForChild( ch);
        for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri)
            SF_ImprovedLaplBeltramiOnTriangle( t, &triangle.GetBary( tri),
                    triangle, Grad,  Numb, e, triangle.GetAbsDet( tri), f.Data, sf_);
    } // Ende der for-Schleife ueber die Kinder
}

bool MarkInterface ( instat_scalar_fun_ptr DistFct, double width, MultiGridCL& mg, Uint f_level, Uint c_level, double t)
{
    bool marked= false;
    DROPS_FOR_TRIANG_TETRA( mg, /*default-level*/-1, it)
    {
        double d= 1e99;
        int num_pos= 0;
        for (Uint j=0; j<10; ++j)
        {
            const double dist= j<4 ? DistFct( it->GetVertex( j)->GetCoord(), t)
                                   : DistFct( GetBaryCenter( *it->GetEdge(j-4)), t);
            if (dist>=0) ++num_pos;
            d= std::min( d, std::abs( dist));
        }

        const bool vzw= num_pos!=0 && num_pos!=10; // change of sign
        if ( f_level==(Uint)(-1) || c_level==(Uint)(-1)){
            if (d<=width || vzw){
                marked= true;
                it->SetRegRefMark();
            }
        }
        else{
            const Uint l= it->GetLevel();
            // In the shell:      level should be f_level_.
            // Outside the shell: level should be c_level_.
            const Uint soll_level= (d<=width || vzw) ? f_level : c_level;
            if (l !=  soll_level || (l == soll_level && !it->IsRegular()) )
            { // tetra will be marked for refinement/removement
                if (l <= soll_level)
                    it->SetRegRefMark();
                else // l > soll_level
                    it->SetRemoveMark();
                marked= true;
            }
        }
    }
#ifdef _PAR
    return ProcCL::GlobalOr(marked);
#else
    return marked;
#endif
}

void MarkInterface ( const LevelsetP2CL::const_DiscSolCL& lset, double width, MultiGridCL& mg)
{
    DROPS_FOR_TRIANG_TETRA( mg, /*default-level*/-1, it)
    {
        double d= 1e99;
        int num_pos= 0;
        for (Uint j=0; j<10; ++j)
        {
            const double dist= j<4 ? lset.val( *it->GetVertex( j))
                                   : lset.val( *it->GetEdge(j-4));
            if (dist>=0) ++num_pos;
            d= std::min( d, std::abs( dist));
        }

        const bool vzw= num_pos!=0 && num_pos!=10; // change of sign
        if (d<=width || vzw)
            it->SetRegRefMark();
    }
}


//*****************************************************************************
//                               LevelsetP2CL
//*****************************************************************************

void LevelsetP2ContCL::UpdateContinuous( ){ ; }
void LevelsetP2ContCL::UpdateDiscontinuous( ) { ; }

// setting inital values for level set function in case of continuous P2 FE
// via interpolation
void LevelsetP2ContCL::Init( instat_scalar_fun_ptr phi0, double t)
{
    const Uint lvl= Phi.GetLevel(),
               idx= Phi.RowIdx->GetIdx();

    for (MultiGridCL::TriangVertexIteratorCL it= MG_.GetTriangVertexBegin(lvl),
        end= MG_.GetTriangVertexEnd(lvl); it!=end; ++it)
    {
        if ( it->Unknowns.Exist(idx))
            Phi.Data[it->Unknowns(idx)]= phi0( it->GetCoord(), t);
    }
    for (MultiGridCL::TriangEdgeIteratorCL it= MG_.GetTriangEdgeBegin(lvl),
        end= MG_.GetTriangEdgeEnd(lvl); it!=end; ++it)
    {
        if ( it->Unknowns.Exist(idx))
            Phi.Data[it->Unknowns(idx)]= phi0( GetBaryCenter( *it), t);
    }
}


void LevelsetP2DiscontCL::UpdateContinuous( )
{
    ApplyClementInterpolation();
}


void LevelsetP2DiscontCL::UpdateDiscontinuous( )
{
    ProjectContinuousToDiscontinuous();
}


void LevelsetP2DiscontCL::ProjectContinuousToDiscontinuous()
{

    Phi.t = PhiC->t;

    const Uint lvl= Phi.GetLevel(), idx= PhiC->RowIdx->GetIdx();
    const Uint didx = Phi.RowIdx->GetIdx();

    // Phi.Data =0.;

    DROPS_FOR_TRIANG_TETRA(MG_,lvl,tet){
        for (int i=0; i <4; ++i){
            if (tet->GetVertex(i)->Unknowns.Exist(idx))
            {
                Phi.Data[tet->Unknowns(didx)+i] = PhiC->Data[(*tet).GetVertex(i)->Unknowns(idx)];
            }
            else
                throw DROPSErrCL("Projections not implemented for levelset non-trivial-bnds");
        }
        for (int i =0; i<6; ++i)
        {
            if (tet->GetEdge(i)->Unknowns.Exist(idx))
            {
                Phi.Data[tet->Unknowns(didx)+4+i] = PhiC->Data[(*tet).GetEdge(i)->Unknowns(idx)];
            }
            else
                throw DROPSErrCL("Projections not implemented for levelset non-trivial-bnds");
        }

    }
}

// setting inital values for level set function in case of discontinuous P2 FE
// via elementwise L2-projection
void LevelsetP2DiscontCL::InitProjection( instat_scalar_fun_ptr phi0, double t)
{
// local mass matrix and its QR decomposition are computet on reference element only
    Uint idx_num = idx.GetIdx();
    Uint lvl = idx.TriangLevel();
    QRDecompCL<10> localmass;
    SMatrixCL<10,10> & localmassm = localmass.GetMatrix();
    SVectorCL<10> sol;
    LocalP2CL<> lp2[10];
    Quad5CL<> lp2q[10];

    for(int i=0; i<10; ++i)
    {
        lp2[i][i] = 1.0;
        lp2q[i].assign(lp2[i]);
    }
    for(int i=0; i<10; ++i)
        for(int j=0; j<=i; ++j)
        {
            localmassm(i,j) = Quad5CL<>(lp2q[i]*lp2q[j]).quad(1.0);
            localmassm(j,i) = localmassm(i,j);
        }
    localmass.prepare_solve();

    DROPS_FOR_TRIANG_TETRA(MG_,lvl,sit)
    {
        Uint first = sit->Unknowns(idx_num);
        Quad5CL<> f_q5(*sit, phi0, t);
        for(int i=0; i<10; ++i)
            sol[i] = Quad5CL<>(f_q5*lp2q[i]).quad(1.0);
        localmass.Solve(sol);
        for(int i=0; i<10; ++i){
            Phi.Data[first++] = sol[i];
        }
    }

    UpdateContinuous();

}

void LevelsetP2DiscontCL::Init( instat_scalar_fun_ptr phi0, double t)
{
    Uint idx_num = idx.GetIdx();
    Uint lvl = idx.TriangLevel();
    SVectorCL<10> sol;
    DROPS_FOR_TRIANG_TETRA(MG_,lvl,sit) // sit is tetra
    {
        for (int i=0; i <4; ++i)
        {
          sol[i] = phi0((*sit).GetVertex(i)-> GetCoord(), t);
        }
        for (int i =0; i<6; ++i)
        {
            sol[i+4] = phi0(GetBaryCenter(*(*sit).GetEdge(i)), t);
        }
        Uint first = sit->Unknowns(idx_num);
        for(int i=0; i<10; ++i){
            Phi.Data[first++] =sol[i];
        }
    }

    ApplyClementInterpolation();
}


void LevelsetP2DiscontCL::ApplyZeroOrderClementInterpolation()
{
    // PhiC->Data.resize(
    PhiC->t = Phi.t;
    const Uint lvl= Phi.GetLevel(), idx= PhiC->RowIdx->GetIdx();
    double tetvol;
    const Uint didx = Phi.RowIdx->GetIdx();
    PhiC->Data =0.;
    VectorCL vols(PhiC->Data);
    vols =0.;

    std::cout << PhiC->Data.size() << std::endl;
    DROPS_FOR_TRIANG_TETRA(MG_,lvl,tet){
        tetvol = (*tet).GetVolume();
        for (int i=0; i <4; ++i){
            if (tet->GetVertex(i)->Unknowns.Exist(idx))
            {
                PhiC->Data[(*tet).GetVertex(i)->Unknowns(idx)] += Phi.Data[tet->Unknowns(didx)+i]*tetvol;
                vols[(*tet).GetVertex(i)->Unknowns(idx)] += tetvol;
            }
        }
        for (int i =0; i<6; ++i)
        {
            if (tet->GetEdge(i)->Unknowns.Exist(idx))
            {
                PhiC->Data[(*tet).GetEdge(i)->Unknowns(idx)] += Phi.Data[tet->Unknowns(didx)+i+4]*tetvol;
                vols[(*tet).GetEdge(i)->Unknowns(idx)] += tetvol;
            }
        }

    }
    for (Uint i=0; i<vols.size(); ++i){
        PhiC->Data[i] /= vols[i];
    }
}

void evaluate_polys(std::valarray<double>& q, const DROPS::QuadDomainCL& qdom, const Point3DCL& v, int j, const TetraCL& tet)
{
    Uint E[10][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{2,0,0},{1,1,0},{0,2,0},{1,0,1},{0,1,1},{0,0,2}};
    q.resize(qdom.vertex_size());
    QuadDomainCL::const_vertex_iterator qit = qdom.vertex_begin();
    Point3DCL quadp;

    for (Uint k=0; k< qdom.vertex_size(); ++k, ++qit){
        quadp = GetWorldCoord(tet, *qit)-v;
        q[k] = pow(quadp[0],E[j][0])*pow(quadp[1],E[j][1])*pow(quadp[2],E[j][2]);
    }
}

void LevelsetP2DiscontCL::ApplyClementInterpolation() //LevelsetP2DiscontCL& discontls
{
    PhiC->SetIdx( idxC );

    PhiC->t = Phi.t;
    const Uint lvl= PhiC->GetLevel(), idx= PhiC->RowIdx->GetIdx();
    PhiC->Data =0.;
    const IdxT num_unks= PhiC->RowIdx->NumUnknowns();
    LocalP2CL<> phiD;
    std::vector< SMatrixCL<10,10> > mat(num_unks);
    std::vector< SVectorCL<10> > rhs(num_unks);
    std::valarray<double> q[10], qphiD;
    QuadDomainCL qdom;
    double absdet;
    Point3DCL v;
    const TetraSignEnum s= AllTetraC;

    DROPS_FOR_TRIANG_TETRA(MG_,lvl,tet){// set up local matrices and right hand sides for each dof of the continuous P2-fct
        absdet = 6.*(*tet).GetVolume();
        phiD.assign(*tet,Phi,GetBndData());
        make_SimpleQuadDomain<Quad5DataCL>(qdom, s);
        resize_and_evaluate_on_vertexes (phiD,qdom,qphiD); //<LocalP2CL, QuadDomainCL, std::valarray<double> >

        for (int i=0; i<10; ++i){ // Dofs
            v = i<4? (tet->GetVertex(i))->GetCoord() : GetBaryCenter(*(*tet).GetEdge(i-4));
            const Uint dofi = i<4 ?  (*tet).GetVertex(i)->Unknowns(idx) : (*tet).GetEdge(i-4)->Unknowns(idx);
            for (int j=0; j<10; ++j) // polys
                evaluate_polys(q[j], qdom, v, j, *tet);
            for (int k=0; k<10; ++k){
                rhs[dofi][k] += quad(qphiD*q[k],absdet,qdom);
                for (int l=0; l<10; ++l)
                    mat[dofi](l,k) += quad(q[l]*q[k],absdet,qdom);

            }
        }
    }
    // solve local 10x10-systems
    QRDecompCL<10> localmat;
    SMatrixCL<10,10> & localm = localmat.GetMatrix();
    SVectorCL<10> sol;
    for (Uint i=0; i<num_unks; ++i){
        localm = mat[i];
        sol = rhs[i];
        localmat.prepare_solve();
        localmat.Solve(sol);
        PhiC->Data[i] = sol[0];
    }
}

void LevelsetP2CL::CreateNumbering( Uint level)
{
    idx.CreateNumbering( level, MG_, BndData_);
    Phi.SetIdx(&idx);
}



void LevelsetP2CL::CreateNumbering(Uint level, MLIdxDescCL* idx)
{
    idx->CreateNumbering( level, MG_, BndData_);
}


void LevelsetP2CL::Reparam( int method, bool Periodic)
/** \param method How to perform the reparametrization (see description of ReparamFactoryCL for details)
    \param Periodic: If true, a special variant of the algorithm for periodic boundaries is used.
*/
{
    std::unique_ptr<ReparamCL> reparam= ReparamFactoryCL::GetReparam( MG_, *PhiC, method, Periodic, &BndData_, perDirections);
    reparam->Perform();
    UpdateDiscontinuous();
}

void LevelsetP2CL::AccumulateBndIntegral( VecDescCL& f) const
{
    ScopeTimerCL scope("AccumulateBndIntegral");
    SurfTensAccumulatorCL* accu;

    switch (SF_)
    {
      case SF_LB:
          accu= new NaiveLaplaceBeltramiAccuCL( *this, f, sf_.GetSigma()(std_basis<3>(0), 0.)); break;
      case SF_Const:
          accu= new ConstSurfTensAccumulatorCL( *this, f, sf_.GetSigma()(std_basis<3>(0), 0.)); break;
      case SF_ImprovedLBVar:
          accu= new VarImprovedLaplaceBeltramiAccuCL( *this, f, sf_); break;
      default:
        throw DROPSErrCL("LevelsetP2CL::AccumulateBndIntegral not implemented for this SurfaceForceT");
    }
    TetraAccumulatorTupleCL accus;
    ProgressBarTetraAccumulatorCL accup(MG_, "SurfTension Setup", Phi.RowIdx->TriangLevel());
    accus.push_back( &accup);
    accus.push_back( accu);
    accumulate( accus, MG_, Phi.RowIdx->TriangLevel(), Phi.RowIdx->GetBndInfo());

    delete accu;
}

double LevelsetP2CL::GetInterfaceArea() const
{
    InterfaceTriangleCL triangle;
    const DROPS::Uint lvl = idx.TriangLevel();
    BndDataCL lsetbnd = GetBndData();
    double area = 0;
    DROPS_FOR_TRIANG_TETRA( MG_, lvl, it){
        triangle.Init( *it, Phi, lsetbnd);
        for(int ch=0;ch<8;++ch)
        {
            if (!triangle.ComputeForChild(ch)) // no patch for this child
                continue;
            for(int v=0;v<triangle.GetNumTriangles();v++)
                area += triangle.GetAbsDet(v);
        }
    }
#ifdef _PAR
    area = ProcCL::GlobalSum(area);
#endif
    return area*0.5;
}

double LevelsetP2CL::GetWetArea() const
{
    InterfaceTriangleCL triangle;
    const DROPS::Uint lvl = idx.TriangLevel();
    BndDataCL lsetbnd=GetBndData();
    BndTriangPartitionCL      bndpartition_;
    QuadDomainCL              bndq5dom_;
    PrincipalLatticeCL lat= PrincipalLatticeCL::instance( 2);
    std::valarray<double>     ls_loc(lat.vertex_size());
    double area = 0;
    GridFunctionCL<> qpr;
    LocalP2CL<> ls_loc0;
    DROPS_FOR_TRIANG_TETRA( MG_, lvl, it) {
        for(Uint v=0; v<4; v++)
        {
            if(lsetbnd.IsOnSlipBnd(*it->GetFace(v)))  // Do not use lsetbnd
            {
                ls_loc0.assign( *it, Phi, BndData_);
                const bool noCut= equal_signs(ls_loc0);
                if(noCut)
                {
                    if(ls_loc0[0]>0) continue;
                    const FaceCL& face = *it->GetFace(v);
                    double absdet = FuncDet2D(face.GetVertex(1)->GetCoord()-face.GetVertex(0)->GetCoord(),
                                              face.GetVertex(2)->GetCoord()-face.GetVertex(0)->GetCoord());
                    area += absdet/2;
                }
                else
                {
                    evaluate_on_vertexes( GetSolution(), *it, lat, Addr( ls_loc));
                    //Does this partition work for no-cut situations??
                    bndpartition_.make_partition2D<SortedVertexPolicyCL, MergeCutPolicyCL>( lat, v, ls_loc);
                    make_CompositeQuad5BndDomain2D( bndq5dom_, bndpartition_,*it);

                    LocalP1CL<double> fun;
                    for (Uint i= 0; i<4; ++i) fun[i]=1.0;
                    resize_and_evaluate_on_vertexes(fun, bndq5dom_, qpr);
                    area += quad( qpr, 1., bndq5dom_, NegTetraC);
                }
            }
        }
    }
#ifdef _PAR
    area = ProcCL::GlobalSum(area);
#endif
    return area;
}

double LevelsetP2CL::GetVolume( double translation, int l) const
{
    if (l==0)
        ++l;
    double Volume= l > 0 ? GetVolume_Composite( translation, l)
                         : GetVolume_Extrapolation( translation, -l);

#ifdef _PAR
    Volume = ProcCL::GlobalSum(Volume);
#endif
    return Volume;
}

double LevelsetP2CL::GetVolume_Extrapolation( double translation, int l) const
{
    double vol = 0.;
    QuadDomainCL qdom;
    DROPS::ExtrapolationToZeroCL extra( l, DROPS::RombergSubdivisionCL());
    // DROPS::ExtrapolationToZeroCL extra( l, DROPS::HarmonicSubdivisionCL());
    LocalP2CL<> loc_phi;
    DROPS_FOR_TRIANG_TETRA( MG_, idx.TriangLevel(), it) {
        loc_phi.assign(*it,Phi,GetBndData());
        loc_phi+= translation;
        make_ExtrapolatedQuad5Domain( qdom, loc_phi, extra);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        vol+=quad( integrand, it->GetVolume()*6., qdom, NegTetraC);
    }
    return vol;
}

double LevelsetP2CL::GetVolume_Composite( double translation, int l) const
{
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance ( l);
    double vol = 0.;
    std::valarray<double> ls_values (lat.vertex_size());
    QuadDomainCL qdom;
    LocalP2CL<> loc_phi;
    TetraPartitionCL partition;
    DROPS_FOR_TRIANG_TETRA( MG_, idx.TriangLevel(), it) {
        loc_phi.assign(*it,Phi,GetBndData());
        loc_phi+= translation;
        evaluate_on_vertexes (loc_phi, lat, Addr(ls_values));
        partition.make_partition< SortedVertexPolicyCL,MergeCutPolicyCL>(lat, ls_values);
        make_CompositeQuad5Domain( qdom, partition);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        vol+=quad( integrand, it->GetVolume()*6., qdom, NegTetraC);
    }
    return vol;
}

void LevelsetP2CL::AdjustVolume () const
{
    volume_adjuster_->AdjustVolume();
}

void LevelsetP2CL::InitVolume (double vol)
{
    volume_adjuster_->InitVolume (vol);
}

void LevelsetP2CL::SmoothPhi( VectorCL& SmPhi, double diff) const
{
    Comment("Smoothing for curvature calculation\n", DebugDiscretizeC);
    MatrixCL M, A, C;
    SetupSmoothSystem( M, A);
    C.LinComb( 1, M, diff, A);
#ifndef _PAR
    typedef SSORPcCL PcT;
#else
    typedef JACPcCL PcT;
#endif
    PcT pc;
    PCGSolverCL<PcT> pcg( pc, 500, 1e-10);
    pcg.Solve( C, SmPhi, M*Phi.Data, idx.GetEx());
#ifndef _PAR
    __UNUSED__ double inf_norm= supnorm( SmPhi-Phi.Data);
#else
    __UNUSED__ const double inf_norm= ProcCL::GlobalMax(supnorm( SmPhi-Phi.Data));
#endif
    Comment("||SmPhi - Phi||_oo = " <<inf_norm<< std::endl, DebugDiscretizeC);
}

void LevelsetP2CL::SetupSmoothSystem( MatrixCL& M, MatrixCL& A) const
// used for smoothing of Phi before computing curvature term
//
// M = mass matrix for P2 elements
// A = stiffness matrix for P2 elements
{
    const IdxT num_unks= Phi.RowIdx->NumUnknowns();
    const Uint lvl= Phi.GetLevel();

#ifndef _PAR
    __UNUSED__ const IdxT allnum_unks= num_unks;
#else
    __UNUSED__ const IdxT allnum_unks= ProcCL::GlobalSum(num_unks);
#endif
    Comment("entering Levelset::SetupSmoothSystem: " << allnum_unks << " levelset unknowns.\n", DebugDiscretizeC);

    SparseMatBuilderCL<double> Mb(&M, num_unks, num_unks);
    SparseMatBuilderCL<double> Ab(&A, num_unks, num_unks);

    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3>     T;
    P2DiscCL::GetGradientsOnRef( GradRef);

    IdxT         Numb[10];
    double det, absdet;

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);

        GetLocalNumbP2NoBnd( Numb, *sit, *Phi.RowIdx);

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
        {
            for(int j=0; j<10; ++j)
            {
                // M_ij = ( v_j, v_i),    A_ij = ( grad v_j, grad v_i)
                Mb( Numb[i], Numb[j])+= P2DiscCL::GetMass(i,j)*absdet;
                Ab( Numb[i], Numb[j])+= Quad2CL<>(dot( Grad[j], Grad[i])).quad( absdet);
            }
        }
    }
    Mb.Build();
    Ab.Build();
}

void LevelsetP2CL::GetMaxMinGradPhi(double& maxGradPhi, double& minGradPhi) const
{
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det;
    InterfacePatchCL patch;

    P2DiscCL::GetGradientsOnRef( GradRef);
    maxGradPhi= -1.;
    minGradPhi= 1e99;

    double maxNorm;
    double minNorm;

    DROPS_FOR_TRIANG_TETRA( MG_, idx.TriangLevel(), it)
    {
        GetTrafoTr( T, det, *it);
        P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder
        patch.Init( *it, *PhiC, BndData_);

        // compute maximal norm of grad *PhiC
        Quad2CL<Point3DCL> gradPhi;
        for (int v=0; v<10; ++v) // init gradPhi, Coord
            gradPhi+= patch.GetPhi(v)*Grad[v];

        VectorCL normGrad( 5);
        for (int v=0; v<5; ++v) // init normGrad
            normGrad[v]= norm( gradPhi[v]);
        maxNorm= normGrad.max();
        minNorm= normGrad.min();
        if (maxNorm > maxGradPhi) maxGradPhi= maxNorm;
        if (minNorm < minGradPhi && patch.Intersects()) minGradPhi= minNorm;
    }
#ifdef _PAR
    maxGradPhi= ProcCL::GlobalMax( maxNorm= maxGradPhi);
    minGradPhi= ProcCL::GlobalMin( minNorm= minGradPhi);
#endif
}

// Creates new Levelset-Object, has to be cleaned manually
LevelsetP2CL * LevelsetP2CL::Create(  MultiGridCL& MG, const LsetBndDataCL& lsetbnddata, SurfaceTensionCL& sf, const ParamCL & P)
{
    LevelsetP2CL * plset;
    if (P.get<int>("Discontinuous") <= 0)
        plset = new LevelsetP2ContCL ( MG, lsetbnddata, sf, P.get<double>("SD"), P.get<double>("CurvDiff"));
    else
        plset = new LevelsetP2DiscontCL ( MG, lsetbnddata, sf, P.get<double>("SD"), P.get<double>("CurvDiff"));

    plset->volume_adjuster_= VolumeAdjustmentCL::Create (plset, P);

    return plset;
}

// Creates new Levelset-Object, has to be cleaned manually
LevelsetP2CL * LevelsetP2CL::Create(  MultiGridCL& MG, const LsetBndDataCL& lsetbnddata, SurfaceTensionCL& sf, bool discontinuous, double SD, double curvdiff)
{
    LevelsetP2CL * plset;
    if (!discontinuous)
        plset = new LevelsetP2ContCL ( MG, lsetbnddata, sf, SD, curvdiff);
    else
        plset = new LevelsetP2DiscontCL ( MG, lsetbnddata, sf, SD, curvdiff);
    plset->volume_adjuster_ = std::unique_ptr<VolumeAdjustmentCL> ( new GlobalVolumeAdjustmentCL ( plset));
    return plset;
}

void LevelsetP2CL::SetNumLvl( size_t n)
{
    idx.resize( n, P2_FE, BndData_);
    MLPhi.resize(n);
}

//*****************************************************************************
//                               LevelsetRepairCL
//*****************************************************************************

void LevelsetRepairCL::pre_refine()
{
    p2repair_= std::unique_ptr<RepairP2CL<double>::type >(
        new RepairP2CL<double>::type( ls_.GetMG(), ls_.Phi, ls_.GetBndData()));
}

void
LevelsetRepairCL::post_refine ()
/// Do all things to complete the repairing of the FE level-set function
{
    VecDescCL loc_phi;
    MLIdxDescCL loc_lidx( P2_FE, ls_.idxC->size());
    VecDescCL& phiC= *ls_.PhiC;

    loc_lidx.CreateNumbering( ls_.GetMG().GetLastLevel(), ls_.GetMG(), ls_.GetBndData());
    loc_phi.SetIdx( &loc_lidx);

    if (ls_.IsDiscontinuous())
    {
        // LevelsetP2DiscontCL & lsD = dynamic_cast<LevelsetP2DiscontCL&>(ls_);
        VecDescCL& phiD= ls_.Phi;
        VecDescCL loc_phiD;
        MLIdxDescCL loc_ldidx( P2D_FE, ls_.idx.size());
        ls_.CreateNumbering( ls_.GetMG().GetLastLevel(), &loc_ldidx);
        loc_phiD.SetIdx( &loc_ldidx);
        p2repair_->repair( loc_phiD);

        phiD.Clear( phiD.t);
        ls_.DeleteNumbering( &ls_.idx);
        ls_.idx.swap( loc_ldidx);

        phiD.SetIdx( &ls_.idx);
        phiD.Data= loc_phiD.Data;

    }
    else
    {
        p2repair_->repair( loc_phi);
    }

    phiC.Clear( phiC.t);
    ls_.DeleteNumbering( ls_.idxC );
    ls_.idxC->swap( loc_lidx);

    phiC.SetIdx( ls_.idxC);

    if (!ls_.IsDiscontinuous())
        phiC.Data= loc_phi.Data;

    ls_.UpdateContinuous();
}

void
LevelsetRepairCL::pre_refine_sequence ()
{}

void
LevelsetRepairCL::post_refine_sequence ()
{
    ls_.GetVolumeAdjuster()->Repair();
}

//*****************************************************************************
//                               LevelsetModifyCL
//*****************************************************************************
void LevelsetModifyCL::maybeDoReparam( LevelsetP2CL& lset)
{
    bool doReparam= rpm_Freq_ && step_%rpm_Freq_ == 0;
    double lsetmaxGradPhi, lsetminGradPhi;

    if (doReparam) {
        lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
        doReparam = (lsetmaxGradPhi > rpm_MaxGrad_ || lsetminGradPhi < rpm_MinGrad_);
    }
    // reparam levelset function
    if (doReparam) {
        std::cout << "before reparametrization: minGradPhi " << lsetminGradPhi << "\tmaxGradPhi " << lsetmaxGradPhi << '\n';
        lset.Reparam( rpm_Method_, per_);
        lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
        std::cout << "after  reparametrization: minGradPhi " << lsetminGradPhi << "\tmaxGradPhi " << lsetmaxGradPhi << '\n';
        // volume correction after reparametrization
        lset.AdjustVolume();
        lset.GetVolumeAdjuster()->DebugOutput (std::cout);
    }
}

void LevelsetModifyCL::maybeDoVolCorr( LevelsetP2CL& lset) {
    lset.AdjustVolume();
    lset.GetVolumeAdjuster()->DebugOutput (std::cout);
}

} // end of namespace DROPS

