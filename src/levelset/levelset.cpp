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

/// \brief Accumulator for the improved Laplace-Beltrami discretization of the CSF term.
///
/// Computes the integral
///         \f[ \sigma \int_\Gamma \kappa v \textbf n_h ds = \sigma \int_\Gamma \hat P_h \nabla id \hat P_h\nabla v ds \f]
/// with \f$\hat P_h = \tilde P_h P_h, P_h = I - \textbf n_h\cdot \textbf n_h^T, \tilde P_h = I - \tilde\textbf n_h\cdot \tilde\textbf n_h^T \f$.
/// Discretization error in \f$ H^1(\Omega)\f$ dual norm has order 1 w.r.t. the grid size at the interface, numerical experiments even indicate order 1.5.
class ImprovedLaplaceBeltramiAccuCL : public SurfTensAccumulatorCL
{
  private:
    const double sigma_;

    LocalP1CL<Point3DCL> Grad[10], GradRef[10];
    IdxT Numb[10];
    LocalP2CL<> velR_p[4][8], velR_n[4][8]; // for P2R basis on children
    LocalP2CL<> loc_phi;

  public:
    ImprovedLaplaceBeltramiAccuCL( const LevelsetP2CL& ls, VecDescCL& f_Gamma, double sigma)
     : SurfTensAccumulatorCL( ls, f_Gamma), sigma_(sigma)
    { P2DiscCL::GetGradientsOnRef( GradRef); }

    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new ImprovedLaplaceBeltramiAccuCL ( *this); };
};

void ImprovedLaplaceBeltramiAccuCL::visit ( const TetraCL& t)
{
    const Uint idx_f=   f.RowIdx->GetIdx();
    const bool velXfem= f.RowIdx->IsExtended();
    double det;

    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder
    LocalP1CL<Point3DCL> n;

    loc_phi.assign( t, SmPhi_, lsetbnd_);
    triangle.Init( t, loc_phi);
    for (int v=0; v<10; ++v)
    { // collect data on all DoF
        const UnknownHandleCL& unk= v<4 ? t.GetVertex(v)->Unknowns : t.GetEdge(v-4)->Unknowns;
        Numb[v]= unk.Exist(idx_f) ? unk(idx_f) : NoIdx;
        for (int k=0; k<4; ++k)
            n[k]+= triangle.GetPhi(v)*Grad[v][k];
    }

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
        BaryPQR/= 3.;    BarySQR/= 3.;

        typedef SArrayCL<Point3DCL,3> ProjT;
        GridFunctionCL<ProjT> GradId( ProjT(), 6);  // values in P, Q, R, S, BaryPQR, BarySQR
        for (int p=0; p<6; ++p)
        {
            Point3DCL np= n( p<4 ? triangle.GetBary(p) : p==4 ? BaryPQR : BarySQR);
            if (np.norm()>1e-8) np/= np.norm();
            for (int i=0; i<3; ++i)
                GradId[p][i]= triangle.ApplyProj( std_basis<3>(i+1) - np[i]*np);
            //                     GradId[p][i]= std_basis<3>(i+1) - np[i]*np;
        }
        const double C= triangle.GetAbsDet()*sigma_/2.;
        if (velXfem)
            P2RidgeDiscCL::GetExtBasisOnChildren( velR_p, velR_n, loc_phi);
        for (int v=0; v<(velXfem ? 14 : 10); ++v)
        {
            const IdxT Numbv= v<10 ? Numb[v] : (velXfem && Numb[v-10]!=NoIdx ? f.RowIdx->GetXidx()[Numb[v-10]] : NoIdx);
            if (Numbv==NoIdx) continue;

            LocalP1CL<Point3DCL> gradv; // gradv = gradient of hat function for dof v
            if (v<10) // std basis function
                for (int node=0; node<4; ++node)
                    gradv[node]= Grad[v][node];
            else // extended basis function: tangential derivative is the same for pos./neg. part, ie., P_h grad(vx_p) == P_h grad(vx_n). W.l.o.g. take pos. part for computation.
                P2DiscCL::GetFuncGradient( gradv, velR_p[v-10][ch], Grad);

            for (int i=0; i<3; ++i)
            {
                double intSum= 0; // sum of the integrand in PQR, SQR
                for (int k=0; k<3; ++k)
                {
                    intSum+= inner_prod( GradId[k][i], gradv(triangle.GetBary(k)));
                    if (triangle.IsQuadrilateral())
                        intSum+= triangle.GetAreaFrac() * inner_prod( GradId[k+1][i], gradv(triangle.GetBary(k+1)));
                }
                double intBary= inner_prod( GradId[4][i], gradv(BaryPQR));
                if (triangle.IsQuadrilateral())
                    intBary+= triangle.GetAreaFrac() * inner_prod( GradId[5][i], gradv(BarySQR));
                f.Data[Numbv+i]-= C *(intSum/12. + 0.75*intBary);
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
    if (sf.GetInputMethod() == Sigma_X)
    {
        SF_ImprovedLaplBeltramiOnTriangle( t, p, triangle, Grad_f, Numb, sf.GetSigma(), e, det, f);
        return;
    }
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




/// \brief Accumulator for the Young's force on the three-phase contact line.
///
/// Computes the integral
///         \f[ \sigma \int_{MCL} \cos(\theta)v\cdot \tau ds \f]
/// with \f$\tau \f$ being the normal direction of the moving contact line in tangential plain of domain boundary.
class YoungForceAccumulatorCL : public  TetraAccumulatorCL
{
 private:
    VecDescCL  SmPhi_;
    const BndDataCL<>& lsetbnd_;
    VecDescCL& f;
   // SMatrixCL<3,3> T;
    InterfaceTriangleCL triangle;

    const double sigma_;
    instat_scalar_fun_ptr angle_;	//Young's contact angle
    instat_vector_fun_ptr outnormal_;//outnormal of the domain boundary

    bool SpeBnd; //special boundary condition

   // LocalP1CL<Point3DCL> Grad[10], GradRef[10];
    IdxT Numb[10];
    LocalP2CL<> velR_p[4][8], velR_n[4][8]; // for P2R basis on children
 //   LocalP2CL<> loc_phi;

  public:
    YoungForceAccumulatorCL( const LevelsetP2CL& ls, VecDescCL& f_Gamma, double sigma,instat_scalar_fun_ptr cangle,instat_vector_fun_ptr outnormal)
     :  SmPhi_(ls.Phi),lsetbnd_(ls.GetBndData()),f(f_Gamma), sigma_(sigma),angle_(cangle),outnormal_(outnormal)
    { ls.MaybeSmooth( SmPhi_.Data);
    //P2DiscCL::GetGradientsOnRef( GradRef);
    }

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
    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new YoungForceAccumulatorCL ( *this); };
};

void YoungForceAccumulatorCL::visit ( const TetraCL& t)
{
	SpeBnd=false;
	for(Uint v=0; v<4; v++)
	   	if(lsetbnd_.GetBC(*t.GetFace(v))==Slip0BC||lsetbnd_.GetBC(*t.GetFace(v))==SlipBC ||lsetbnd_.GetBC(*t.GetFace(v))==SymmBC )
	   	{
	   		SpeBnd=true; break;
	   	}
	if(!SpeBnd)
	{
		for(Uint v=0; v<6; v++)
			if(lsetbnd_.GetBC(*t.GetEdge(v))==Slip0BC||lsetbnd_.GetBC(*t.GetEdge(v))==SlipBC ||lsetbnd_.GetBC(*t.GetEdge(v))==SymmBC)
			{
				SpeBnd=true; break;
			}
		if(!SpeBnd)
		return;
	}
    const Uint idx_f=   f.RowIdx->GetIdx();
    const bool velXfem= f.RowIdx->IsExtended();
    if (velXfem)
    	throw DROPSErrCL("WARNING: YoungForceAccumulatorCL : not implemented for velocity XFEM method yet!");
  //  double det;

   // GetTrafoTr( T, det, t);

  //  loc_phi.assign( t, SmPhi_, lsetbnd_);
    triangle.BInit( t, SmPhi_,lsetbnd_); //we have to use this init function!!!!!!!!!
    triangle.SetBndOutNormal(outnormal_);
    for (int v=0; v<10; ++v)
    {   const UnknownHandleCL& unk= v<4 ? t.GetVertex(v)->Unknowns : t.GetEdge(v-4)->Unknowns;
        Numb[v]= unk.Exist(idx_f) ? unk(idx_f) : NoIdx;
    }
	LocalP2CL<double> phi[10];
	for(Uint i=0; i<10; ++i)
		phi[i][i] = 1;
	Point3DCL normal_mcl;
	double costheta[5];
	BaryCoordCL quadBarys[5];
	Point3DCL outnormalOnMcl[5]; //outnormal of the domain boundary
	double sintheta_D;//sin\theta_d
	double weight[5]={0.568888889, 0.47862867,0.47862867,0.236926885,0.236926885};
	//integral in [-1,1]
	double qupt[5]={0,-0.53846931,0.53846931,-0.906179846,0.906179846};
	for (int ch=0; ch<8; ++ch)
//	for(int ch=8;ch<9;++ch)
    {

        if (!triangle.ComputeMCLForChild(ch)) // no patch for this child
            continue;

        BaryCoordCL Barys[2];
        Point3DCL pt0,pt1;
        Point3DCL midpt;
        double length;
        Uint ncl=triangle.GetNumMCL();
        for(Uint i=0;i<ncl;i++)
        {
        	length = triangle.GetInfoMCL(i,Barys[0],Barys[1],pt0,pt1);
        	normal_mcl = triangle.GetMCLNormal(i);
        	sintheta_D= triangle.IsSymmType(i) ? 0 :sin(triangle.GetActualContactAngle(i));
        	for(Uint j=0;j<5;j++)
        	{
        		quadBarys[j]=(Barys[0]+Barys[1])/2+qupt[j]*(Barys[1]-Barys[0])/2;
        		midpt=(pt0+pt1)/2 + qupt[j]*(pt1-pt0)/2;
        		costheta[j]= triangle.IsSymmType(i) ? 1 : cos(angle_(midpt,0));
        		outnormalOnMcl[j]=outnormal_(midpt,0);
        	}


        	for (int v=0; v<10; ++v)
        	{
        		const IdxT Numbv= v<10 ? Numb[v] : (velXfem && Numb[v-10]!=NoIdx ? f.RowIdx->GetXidx()[Numb[v-10]] : NoIdx);
        		if (Numbv==NoIdx) continue;

        		Point3DCL value = ((phi[v](quadBarys[0])*costheta[0]*weight[0]+ phi[v](quadBarys[1])*costheta[1]*weight[1]
        		              + phi[v](quadBarys[2])*costheta[2]*weight[2]+ phi[v](quadBarys[3])*costheta[3]*weight[3]
        		              + phi[v](quadBarys[4])*costheta[4]*weight[4])*length/2)*normal_mcl;
        		value+=(phi[v](quadBarys[0])*weight[0]*outnormalOnMcl[0]+ phi[v](quadBarys[1])*weight[1]*outnormalOnMcl[1]
        		              + phi[v](quadBarys[2])*weight[2]*outnormalOnMcl[2]+ phi[v](quadBarys[3])*weight[3]*outnormalOnMcl[3]
        		        	  + phi[v](quadBarys[4])*weight[4]*outnormalOnMcl[4])*length/2*sintheta_D;
        		//higher order quadrature is used!!
        		for (int j=0; j<3; ++j)
        		{
        			f.Data[Numbv+j] += sigma_*value[j];
        		}
        	}
        }
    } // Ende der for-Schleife ueber die Kinder
}

/// \brief Impoved Accumulator for the Young's force on the three-phase contact line.
///
/// Computes the integral
///         \f[ \sigma \int_{MCL} \cos(\theta)v\cdot \tau ds \f]
/// with \f$\tau \f$ being the normal direction of the moving contact line in tangential plain of domain boundary.
/// we compute the force using the outnormal of the levelset instead of the straight contact line
class ImprovedYoungForceAccumulatorCL : public  TetraAccumulatorCL
{
 private:
    VecDescCL  SmPhi_;
    const BndDataCL<>& lsetbnd_;
    VecDescCL& f;
   // SMatrixCL<3,3> T;
    InterfaceTriangleCL triangle;

    const double sigma_;
    instat_scalar_fun_ptr angle_;	//Young's contact angle
    instat_vector_fun_ptr outnormal_;//outnormal of the domain boundary

    bool SpeBnd; //special boundary condition

   // LocalP1CL<Point3DCL> Grad[10], GradRef[10];
    IdxT Numb[10];
    LocalP2CL<> velR_p[4][8], velR_n[4][8]; // for P2R basis on children
 //   LocalP2CL<> loc_phi;

  public:
    ImprovedYoungForceAccumulatorCL( const LevelsetP2CL& ls, VecDescCL& f_Gamma, double sigma,instat_scalar_fun_ptr cangle,instat_vector_fun_ptr outnormal)
     :  SmPhi_(ls.Phi),lsetbnd_(ls.GetBndData()),f(f_Gamma), sigma_(sigma),angle_(cangle),outnormal_(outnormal)
    { ls.MaybeSmooth( SmPhi_.Data);
    //P2DiscCL::GetGradientsOnRef( GradRef);
    }

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
    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new ImprovedYoungForceAccumulatorCL ( *this); };
};

void ImprovedYoungForceAccumulatorCL::visit ( const TetraCL& t)
{
	SpeBnd=false;
	//check if the tetra contains one face or one edge on slip or symmetric boundary.
	for(Uint v=0; v<4; v++)
	   	if(lsetbnd_.GetBC(*t.GetFace(v))==Slip0BC||lsetbnd_.GetBC(*t.GetFace(v))==SlipBC||lsetbnd_.GetBC(*t.GetFace(v))==SymmBC )
	   	{
	   		SpeBnd=true;
	   		break;
	   	}
	if(!SpeBnd)
	{
		for(Uint v=0; v<6; v++)
			if(lsetbnd_.GetBC(*t.GetEdge(v))==Slip0BC||lsetbnd_.GetBC(*t.GetEdge(v))==SlipBC||lsetbnd_.GetBC(*t.GetEdge(v))==SymmBC )
			{
				SpeBnd=true;
				break;
			}
		if(!SpeBnd)
		return;
	}
    const Uint idx_f=   f.RowIdx->GetIdx();
    const bool velXfem= f.RowIdx->IsExtended();
    if (velXfem)
    	throw DROPSErrCL("WARNING: ImprovedYoungForceAccumulatorCL : not implemented for velocity XFEM method yet!");
    //Initialize one interface patch
    triangle.BInit( t, SmPhi_,lsetbnd_); //we have to use this init function!!!!!!!!!
    triangle.SetBndOutNormal(outnormal_);
    for (int v=0; v<10; ++v)
    {   const UnknownHandleCL& unk= v<4 ? t.GetVertex(v)->Unknowns : t.GetEdge(v-4)->Unknowns;
        Numb[v]= unk.Exist(idx_f) ? unk(idx_f) : NoIdx;
    }
	LocalP2CL<double> phi[10];
	for(Uint i=0; i<10; ++i)
		phi[i][i] = 1;
	Point3DCL normal_mcl[5];     //normal of moving contact lines in tangential surface
	double costheta[5];          //cos\theta_s
	Point3DCL outnormalOnMcl[5]; //outnormal of the domain boundary
	double sintheta_D[5];        //sin\theta_d
	BaryCoordCL quadBarys[5];
	double weight[5]={0.568888889, 0.47862867,0.47862867,0.236926885,0.236926885};
	//integral in [-1,1]
	double qupt[5]={0,-0.53846931,0.53846931,-0.906179846,0.906179846};
	Point3DCL value;
	for (int ch=0; ch<8; ++ch)
//	for(int ch=8;ch<9;++ch)
    {
        if (!triangle.ComputeMCLForChild(ch)) // no patch for this child
            continue;
        BaryCoordCL Barys[2];
        Point3DCL pt0,pt1;
        Point3DCL midpt;
        double length;
        Uint ncl=triangle.GetNumMCL();
        for(Uint i=0;i<ncl;i++)
        {
        	length = triangle.GetInfoMCL(i,Barys[0],Barys[1],pt0,pt1);
        	for(Uint j=0;j<5;j++)
        	{
        		normal_mcl[j] = triangle.GetImprovedMCLNormal(i,(qupt[j]+1)/2);
        		quadBarys[j]=(Barys[0]+Barys[1])/2+qupt[j]*(Barys[1]-Barys[0])/2;
        		midpt=(pt0+pt1)/2 + qupt[j]*(pt1-pt0)/2;
        		costheta[j]=triangle.IsSymmType(i) ? 0 : cos(angle_(midpt,0));
        		sintheta_D[j]=triangle.IsSymmType(i) ? 1 : sin(triangle.GetImprovedActualContactAngle(i,(qupt[j]+1)/2));
        		outnormalOnMcl[j]=outnormal_(midpt,0);
        	}

        	for (int v=0; v<10; ++v)
        	{
        		const IdxT Numbv= v<10 ? Numb[v] : (velXfem && Numb[v-10]!=NoIdx ? f.RowIdx->GetXidx()[Numb[v-10]] : NoIdx);
        		if (Numbv==NoIdx) continue;
				//could be written in a loop
        		value = (phi[v](quadBarys[0])*costheta[0]*weight[0]*normal_mcl[0]+ phi[v](quadBarys[1])*costheta[1]*weight[1]*normal_mcl[1]
        		                 + phi[v](quadBarys[2])*costheta[2]*weight[2]*normal_mcl[2]+ phi[v](quadBarys[3])*costheta[3]*weight[3]*normal_mcl[3]
        		                 + phi[v](quadBarys[4])*costheta[4]*weight[4]*normal_mcl[4])*length/2;
        		value+=(phi[v](quadBarys[0])*sintheta_D[0]*weight[0]*outnormalOnMcl[0]+ phi[v](quadBarys[1])*sintheta_D[1]*weight[1]*outnormalOnMcl[1]
        		                 + phi[v](quadBarys[2])*sintheta_D[2]*weight[2]*outnormalOnMcl[2]+ phi[v](quadBarys[3])*sintheta_D[3]*weight[3]*outnormalOnMcl[3]
        		                 + phi[v](quadBarys[4])*sintheta_D[4]*weight[4]*outnormalOnMcl[4])*length/2;
        		//5 points Gauss–Legendre quadrature is used.
        		for (int j=0; j<3; ++j)
        		{
        			f.Data[Numbv+j] += sigma_*value[j];
        		}
        	}
        }
    } // Ende der for-Schleife ueber die Kinder
}

//*****************************************************************************
//                               LevelsetP2CL
//*****************************************************************************

void LevelsetP2ContCL::UpdateContinuous( ){ ; }
void LevelsetP2ContCL::UpdateDiscontinuous( ) { ; }

// setting inital values for level set function in case of continuous P2 FE
// via interpolation
void LevelsetP2ContCL::Init( instat_scalar_fun_ptr phi0)
{   
    double t= 0.;
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
void LevelsetP2DiscontCL::InitProjection( instat_scalar_fun_ptr phi0)
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
        Quad5CL<> f_q5(*sit, phi0);
        for(int i=0; i<10; ++i)
            sol[i] = Quad5CL<>(f_q5*lp2q[i]).quad(1.0);
        localmass.Solve(sol);
        for(int i=0; i<10; ++i){
            Phi.Data[first++] = sol[i];
        }
    }

    UpdateContinuous();

}

void LevelsetP2DiscontCL::Init( instat_scalar_fun_ptr phi0)
{
    double t= 0.;
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
 
void LevelsetP2CL::CreateNumbering( Uint level, MLIdxDescCL* idx, match_fun match)
{
    idx->CreateNumbering( level, MG_, BndData_, match);
}


void LevelsetP2CL::Reparam( int method, bool Periodic)
/** \param method How to perform the reparametrization (see description of ReparamFactoryCL for details)
    \param Periodic: If true, a special variant of the algorithm for periodic boundaries is used.
*/
{
    std::auto_ptr<ReparamCL> reparam= ReparamFactoryCL::GetReparam( MG_, *PhiC, method, Periodic, &BndData_, perDirections);
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
      case SF_ImprovedLB:
          accu= new ImprovedLaplaceBeltramiAccuCL( *this, f, sf_.GetSigma()(std_basis<3>(0), 0.)); break;
      case SF_ImprovedLBVar:
          accu= new VarImprovedLaplaceBeltramiAccuCL( *this, f, sf_); break;
      default:
        throw DROPSErrCL("LevelsetP2CL::AccumulateBndIntegral not implemented for this SurfaceForceT");
    }
    TetraAccumulatorTupleCL accus;
    ProgressBarTetraAccumulatorCL accup(MG_, "SurfTension Setup", Phi.RowIdx->TriangLevel());
    accus.push_back( &accup);
    accus.push_back( accu);
    accumulate( accus, MG_, Phi.RowIdx->TriangLevel(), Phi.RowIdx->GetMatchingFunction(), Phi.RowIdx->GetBndInfo());

    delete accu;
}


void LevelsetP2CL::AccumulateYoungForce( VecDescCL& f) const
{
	 ScopeTimerCL scope("AccumulateYoungForce");
	TetraAccumulatorCL *accu;
    switch (SF_)
    {
    	//accu= new YoungForceAccumulatorCL( *this, f, sf_.GetSigma()(std_basis<3>(0), 0.),CA_,Bndoutnormal_);
    case SF_ImprovedLB:
       	accu= new ImprovedYoungForceAccumulatorCL( *this, f, sf_.GetSigma()(std_basis<3>(0), 0.),CA_,Bndoutnormal_); break;
    default:
    	throw DROPSErrCL("LevelsetP2CL::AccumulateYoungForce not implemented for non-constant surface tension");
    }
    TetraAccumulatorTupleCL accus;
    accus.push_back( accu);
    accumulate( accus, MG_, Phi.RowIdx->TriangLevel(), Phi.RowIdx->GetMatchingFunction(), Phi.RowIdx->GetBndInfo());
}

//>to do for parallel programe, we need add all values in different process
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
	return area*0.5;
}

//>to do for parallel programe, we need add all values in different process
double LevelsetP2CL::GetWetArea() const
{
	InterfaceTriangleCL triangle;
	const DROPS::Uint lvl = idx.TriangLevel();
	BndDataCL lsetbnd=GetBndData();
	BndTriangPartitionCL 	  bndpartition_;
    QuadDomainCL              bndq5dom_;
    PrincipalLatticeCL lat= PrincipalLatticeCL::instance( 2);
    std::valarray<double>     ls_loc(lat.vertex_size());
	double area = 0;
	GridFunctionCL<> qpr;
	LocalP2CL<> ls_loc0;
	DROPS_FOR_TRIANG_TETRA( MG_, lvl, it){

		for(Uint v=0; v<4; v++)
		{
			if(lsetbnd.GetBC(*it->GetFace(v))==Slip0BC||lsetbnd.GetBC(*it->GetFace(v))==SlipBC)
		  	{
				ls_loc0.assign( *it, Phi, BndData_);
				const bool noCut= equal_signs(ls_loc0);
				if(noCut)
				{
					if(ls_loc0[0]>0) continue;
					const FaceCL& face = *it->GetFace(v);
				    double absdet = FuncDet2D(	face.GetVertex(1)->GetCoord()-face.GetVertex(0)->GetCoord(),
				                                face.GetVertex(2)->GetCoord()-face.GetVertex(0)->GetCoord());
				    area += absdet/2;
				}
				else
				{
					evaluate_on_vertexes( GetSolution(), *it, lat, Addr( ls_loc));
					//Does this partition work for no cut situations??
					bndpartition_.make_partition2D<SortedVertexPolicyCL, MergeCutPolicyCL>( lat, v, ls_loc);
					make_CompositeQuad5BndDomain2D( bndq5dom_, bndpartition_,*it);

					LocalP1CL<double> fun;
					for (Uint i= 0; i<3; ++i)	fun[i]=1.0;
					resize_and_evaluate_on_vertexes(fun, bndq5dom_, qpr);
					area += quad( qpr, bndq5dom_, NegTetraC);
				}
			}
		}
	}
	return area;
}
//>to do for parallel programe, we need add all values in different process
double LevelsetP2CL::GetSurfaceEnergy() const
{
	InterfaceTriangleCL triangle;
	const DROPS::Uint lvl = idx.TriangLevel();
	BndDataCL lsetbnd = GetBndData();
	instat_scalar_fun_ptr surface_tension = sf_.GetSigma();
	DROPS::Quad5_2DCL<> sfdensity0,sfdensity1;

	BndTriangPartitionCL 	  bndpartition_;
    QuadDomainCL              bndq5dom;
    PrincipalLatticeCL lat= PrincipalLatticeCL::instance( 2);
    std::valarray<double>     ls_loc(lat.vertex_size());

	GridFunctionCL<> qpr;
	BaryCoordCL bary[3];
	double sftn = surface_tension(std_basis<3>(0),0); //>to do: compute standard surface for general case
	double total_energy1=0,total_energy2=0;
	DROPS_FOR_TRIANG_TETRA( MG_, lvl, it){
		triangle.Init( *it, Phi, lsetbnd);
		for(int ch=0;ch<8;++ch)
		{
			if (!triangle.ComputeForChild(ch)) // no patch for this child
		            continue;
			for(int tri=0;tri<triangle.GetNumTriangles();tri++)
			{
				sfdensity0.assign(  *it, &triangle.GetBary( tri), surface_tension,0);
				total_energy1 += sfdensity0.quad( triangle.GetAbsDet( tri));
			}
		}

		for(Uint v=0; v<4; v++)
		{
			if(lsetbnd.GetBC(*it->GetFace(v))==Slip0BC||lsetbnd.GetBC(*it->GetFace(v))==SlipBC)
			{
						//ls_loc0.assign( *it, Phi, BndData_);
						evaluate_on_vertexes( GetSolution(), *it, lat, Addr( ls_loc));
						const bool noCut= equal_signs(ls_loc);
						if(noCut)
						{
							if(ls_loc[0]>0) continue;
							const FaceCL& face = *it->GetFace(v);
							for (Uint i= 0; i<3; ++i)
							{
								bary[i][VertOfFace(v, i)]=1;
							}
							sfdensity1.assign(*it, bary, CA_,0 );
							for(Uint s=0;s<sfdensity1.size();s++)
								sfdensity1[s]=std::cos(sfdensity1[s]);//??

						    double absdet = FuncDet2D(	face.GetVertex(1)->GetCoord()-face.GetVertex(0)->GetCoord(),
						                                face.GetVertex(2)->GetCoord()-face.GetVertex(0)->GetCoord());
						    total_energy2 += -sftn*sfdensity1.quad(absdet);
						}
						else
						{

							//Does this partition work for no cut situations??
							bndpartition_.make_partition2D<SortedVertexPolicyCL, MergeCutPolicyCL>( lat, v, ls_loc);
							make_CompositeQuad5BndDomain2D( bndq5dom, bndpartition_,*it);

							resize_and_evaluate_on_vertexes(CA_,*it, bndq5dom,0, qpr);//???
							for(Uint s=0;s<qpr.size();s++)
								qpr[s]=std::cos(qpr[s]);//??
							total_energy2 += -sftn*quad( qpr, bndq5dom, NegTetraC);
						}
			}
		}

	}
	//std::cout<<total_energy1<<"  "<<total_energy2<<std::endl;
	return total_energy1+total_energy2;
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

double LevelsetP2CL::AdjustVolume (double vol, double tol, double surface, int l) const
{
    tol*=vol;

    double v0=GetVolume(0., l)-vol;
    if (std::abs(v0)<=tol) return 0;

    double d0=0, d1=v0*(surface != 0. ? 1.1/surface : 0.23/std::pow(vol,2./3.));
    // Hinweis: surf(Kugel) = [3/4/pi*vol(Kugel)]^(2/3) * 4pi
    double v1=GetVolume(d1, l)-vol;
    if (std::abs(v1)<=tol) return d1;

    // Sekantenverfahren fuer Startwert
    while (v1*v0 > 0) // gleiches Vorzeichen
    {
        const double d2=d1-1.2*v1*(d1-d0)/(v1-v0);
        d0=d1; d1=d2; v0=v1; v1=GetVolume(d1, l)-vol;
        if (std::abs(v1)<=tol) return d1;
    }

    // Anderson-Bjoerk fuer genauen Wert
    while (true)
    {
        const double d2=(v1*d0-v0*d1)/(v1-v0),
                     v2=GetVolume(d2,l)-vol;
        if (std::abs(v2)<=tol) return d2;

        if (v2*v1 < 0) // ungleiches Vorzeichen
          { d0=d1; d1=d2; v0=v1; v1=v2; }
        else
          { const double c=1.0-v2/v1; d1=d2; v1=v2; v0*= c>0 ? c : 0.5; }
    }
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
    return plset;
}




void LevelsetP2CL::SetNumLvl( size_t n)
{
    match_fun match= MG_.GetBnd().GetMatchFun();
    idx.resize( n, P2_FE, BndData_, match);
    MLPhi.resize(n);
}

//*****************************************************************************
//                               LevelsetRepairCL
//*****************************************************************************

void LevelsetRepairCL::pre_refine()
{
    p2repair_= std::auto_ptr<RepairP2CL<double>::type >(
        new RepairP2CL<double>::type( ls_.GetMG(), ls_.Phi, ls_.GetBndData()));
}

void
LevelsetRepairCL::post_refine ()
/// Do all things to complete the repairing of the FE level-set function
{
    VecDescCL loc_phi;
    MLIdxDescCL loc_lidx( P2_FE, ls_.idxC->size());
    VecDescCL& phiC= *ls_.PhiC;
    match_fun match= ls_.GetMG().GetBnd().GetMatchFun();

    loc_lidx.CreateNumbering( ls_.GetMG().GetLastLevel(), ls_.GetMG(), ls_.GetBndData(), match);
    loc_phi.SetIdx( &loc_lidx);

    if (ls_.IsDiscontinuous())
    {
        // LevelsetP2DiscontCL & lsD = dynamic_cast<LevelsetP2DiscontCL&>(ls_);
        VecDescCL& phiD= ls_.Phi;
        VecDescCL loc_phiD;
        MLIdxDescCL loc_ldidx( P2D_FE, ls_.idx.size());
        ls_.CreateNumbering( ls_.GetMG().GetLastLevel(), &loc_ldidx, match);
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

} // end of namespace DROPS

