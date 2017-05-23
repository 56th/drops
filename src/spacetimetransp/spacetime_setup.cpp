/// \file spacetime_setup.cpp
/// \brief space time setup
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
#include "spacetimetransp/spacetime_setup.h"
// #include "num/interfacePatch.h"


namespace DROPS
{

using namespace STXFEM;

void debug_kappa(LocalP2CL<> & lsold, LocalP2CL<> & lsnew, const GridFunctionCL<Point4DCL> & points, 
                   GridFunctionCL<> & kappa_neg, GridFunctionCL<> & kappa_pos)
{
    ScopeTimer scopetiming("debug_kappa");
    const Uint ints_per_space_edge = 1;
    QuadDomainCL q2dom_;
    PrincipalLatticeCL lat( PrincipalLatticeCL::instance( ints_per_space_edge ));
    TetraPartitionCL partition_;
    LocalP1CL<> one(1.0);
    GridFunctionCL<> ones;

    kappa_neg.resize(points.size());
    kappa_pos.resize(points.size());
    for (Uint k = 0; k < points.size(); ++k)
    {
        // kappa_neg[k] = 0.5;
        // kappa_pos[k] = 0.5;
        
        const double time = points[k][3];
        LocalP2CL<> lsetcur = lsnew;
        lsetcur *= time;
        lsetcur += (1.0-time) * lsold;
        
        std::valarray<double> ls_loc_( ((ints_per_space_edge+1) * (ints_per_space_edge+2) * (ints_per_space_edge+3)) / 6);
        evaluate_on_vertexes(lsetcur, lat, Addr( ls_loc_));
        partition_.make_partition<SortedVertexPolicyCL, MergeCutPolicyCL>( lat, ls_loc_);
        make_CompositeQuad2Domain( q2dom_, partition_);
        resize_and_evaluate_on_vertexes( one, q2dom_, ones);
        kappa_neg[k] = quad( ones , 6.0 , q2dom_ , NegTetraC);
        kappa_pos[k] = 1.0 - kappa_neg[k];
    }
    
}


void STVolumeAccumulator_P1SP1TXCL::transform_comp_elmats_to_xelmats()
{
    if (b_ != NULL)
        composite_elcontrib_to_xelcontrib<8>(coup_s_s_neg, coup_s_s_pos,
                                             elvec_neg, elvec_pos, sign,
                                             intf_coup_s_s, intf_coup_s_x, intf_coup_x_s, intf_coup_x_x,
                                             coup_s_s, coup_s_x, coup_x_s, coup_x_x,
                                             elvec, elvec_x);
    else
        composite_elcontrib_to_xelcontrib<8>(coup_s_s_neg, coup_s_s_pos,
                                             NULL, NULL, sign,
                                             intf_coup_s_s, intf_coup_s_x, intf_coup_x_s, intf_coup_x_x,
                                             coup_s_s, coup_s_x, coup_x_s, coup_x_x,
                                             NULL, NULL);
}

void STVolumeAccumulator_P1SP1TXCL::update_global_matrix_without_x()
{
    for(Uint i=0; i<4; ++i)          // assemble row i
        if (UnknownIdx[i]!= NoIdx)  // vertex i is not on a Dirichlet boundary
        {
            for(Uint ti = PAST; ti<= FUTURE; ti++)
            {
                const IdxT idx_i= UnknownIdx[i]+ti;
                
                if (b_ != NULL)
                    b_->Data[idx_i] += elvec[i+4*ti];

                for(Uint j=0; j<4;++j)
                    if (UnknownIdx[j]!= NoIdx) // vertex j is not on a Dirichlet boundary
                        for(Uint tj = PAST; tj<=FUTURE; tj++)
                        {
                            {
                                const IdxT idx_j= UnknownIdx[j]+tj;
                                (*A_)( idx_i,  idx_j)+=coup_s_s[i+4*ti][j+4*tj];    
                            }
                        }
            }
        } 
}

void STVolumeAccumulator_P1SP1TXCL::update_global_matrix_with_x()
{
    for(Uint ti = PAST; ti<= FUTURE; ti++)
        for(Uint i=0; i<4; ++i)          // assemble row i
            if (UnknownIdx[i]!= NoIdx)  // vertex i is not on a Dirichlet boundary
            {
                const IdxT idx_i= UnknownIdx[i]+ti;
                const IdxT xidx_i= Xidx[idx_i];

                if (b_ != NULL)
                {
                    b_->Data[idx_i] += elvec[i+4*ti];
                    if (xidx_i != NoIdx)
                        b_->Data[xidx_i] += elvec_x[i+4*ti];
                }

                for(Uint tj = PAST; tj<=FUTURE; tj++)
                    for(Uint j=0; j<4;++j)
                    {
                        if (UnknownIdx[j]!= NoIdx) // vertex j is not on a Dirichlet boundary
                        {
                            const IdxT idx_j= UnknownIdx[j]+tj;
                            const IdxT xidx_j= Xidx[idx_j];
                            (*A_)( idx_i,  idx_j)+=coup_s_s[i+4*ti][j+4*tj];    
                            if (xidx_i != NoIdx)
                                (*A_)(xidx_i,  idx_j)+=coup_x_s[i+4*ti][j+4*tj];    
                            if (xidx_j != NoIdx)
                                (*A_)( idx_i, xidx_j)+=coup_s_x[i+4*ti][j+4*tj];    
                            if (xidx_i != NoIdx && xidx_j != NoIdx)
                                (*A_)(xidx_i, xidx_j)+=coup_x_x[i+4*ti][j+4*tj];    
                        }
                    }
            } 
}

void STVolumeAccumulator_P1SP1TXCL::update_coupling(const TetraCL& sit, bool with_x)
{
    for(Uint ti = PAST; ti<= FUTURE; ti++)
        for(Uint tj = PAST; tj<=FUTURE; tj++)
            for(Uint i=0; i<4; ++i)          // assemble row i
                if (UnknownIdx[i]!= NoIdx)  // vertex i is not on a Dirichlet boundary
                {
                    const IdxT idx_i= UnknownIdx[i]+ti;
                    const IdxT xidx_i= Xidx[idx_i];
                    for(Uint j=0; j<4;++j)
                    {
                        if (UnknownIdx[j]== NoIdx) // vertex j is on a Dirichlet boundary
                        {
                            const BndDataCL<> * BndData_cur = sign[j+4*tj] ? BndData_pos : BndData_neg;
                            const BndDataCL<> * BndData_nei = sign[j+4*tj] ? BndData_neg : BndData_pos;
                            const double beta_cur = sign[j+4*tj] ? beta_pos : beta_neg;
                            const double beta_nei = sign[j+4*tj] ? beta_neg : beta_pos;
                            const double timej = tj == PAST ? told : tnew;

                            b_->Data[ idx_i]-= coup_s_s[ti*4+i][tj*4+j] * BndData_cur->GetDirBndValue(*sit.GetVertex(j),timej) * beta_cur;
                            if (with_x)
                            {
                                b_->Data[ idx_i]-= coup_s_x[ti*4+i][tj*4+j] * BndData_nei->GetDirBndValue(*sit.GetVertex(j),timej) * beta_nei;
                            }


                            if (with_x && xidx_i != NoIdx)
                            {
                                b_->Data[xidx_i]-= coup_x_x[ti*4+i][tj*4+j] * BndData_nei->GetDirBndValue(*sit.GetVertex(j),timej) * beta_nei;
                                b_->Data[xidx_i]-= coup_x_s[ti*4+i][tj*4+j] * BndData_cur->GetDirBndValue(*sit.GetVertex(j),timej) * beta_cur;
                            }
                        }
                    }
                } 
}


STVolumeAccumulator_P1SP1TXCL::STVolumeAccumulator_P1SP1TXCL(const MultiGridCL& MG, const BndDataCL<> * BndData_neg_in, 
                                                             const BndDataCL<> * BndData_pos_in, 
                                                             const LevelsetP2CL * lsetp2old_in,
                                                             const LevelsetP2CL * lsetp2new_in,
                                                             MatrixCL* Amat, VecDescCL* b, 
                                                             const IdxDescCL& RowIdx, const IdxDescCL& ColIdx, 
                                                             const double t1, const double t2,
                                                             const ParamCL::ptree_type & P):
MG_(MG), BndData_neg(BndData_neg_in), BndData_pos(BndData_pos_in), 
    Amat_(Amat), b_(b), RowIdx_(RowIdx), ColIdx_(ColIdx), 
    lsetp2old(lsetp2old_in),lsetp2new(lsetp2new_in),
    A_(0),
    Xidx(RowIdx_.GetXidx()),
    lvl(RowIdx.TriangLevel()),
    idx(RowIdx.GetIdx()),
    ints_per_space_edge(P.get<int>("Quadrature.SubIntervalsPerEdge")),
    subtimeintervals(P.get<int>("Quadrature.SubTimeIntervals")),
    told(t1), tnew(t2),
    beta_neg(P.get<double>("HNeg",1.0)),
    beta_pos(P.get<double>("HPos",1.0))
{
    
    for(int i=0; i<4; i++)
    {
        LocalP1CL<> test = 0.0;
        test[i] = 1.0;
        phi[i] = LocalP2CL<>(test);
        q_phi[i] = Quad5CL<>(phi[i]);
    }

    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            phi_i_phi_j[i][j] = phi[i]*phi[j];
   
}

void STVolumeAccumulator_P1SP1TXCL::begin_accumulation ()
{
    if (b_ != 0) b_->Clear( tnew);
    if (Amat_)    
        A_ = new MatrixBuilderCL( Amat_, RowIdx_.NumUnknowns(), ColIdx_.NumUnknowns());
}

void STVolumeAccumulator_P1SP1TXCL::finalize_accumulation ()
{
    if (A_ != 0){
        A_->Build();
        delete A_;   
    }     
}

void STVolumeAccumulator_P1SP1TXCL::visit (const TetraCL& sit)
{
    local_setup(sit);
    ScopeTimer scopetiming("STVolumeAccumulator_P1SP1TXCL::visit without local_setup");
    if (A_ != 0)
    {
        if (iscut)  
        {
            transform_comp_elmats_to_xelmats();
            update_global_matrix_with_x();
        }
        else
            update_global_matrix_without_x();          
        // output_elmats();
        // getchar();
    }
    if (b_ != 0 && BndData_neg != 0 && BndData_pos != 0) 
    {
        if (iscut)  
            update_coupling(sit, /* with_x: */ true);
        else
            update_coupling(sit, /* with_x: */ false);
    }


}

void STVolumeAccumulator_P1SP1TXCL::output_elmats ()
{
    if (iscut){
   
        std::cout << " coup_s_s_neg = \n";
        for (Uint i = 0; i < 8; ++i)
        {
            double sumi = 0.0;
            for (Uint j = 0; j < 8; ++j)
            {
                std::cout << std::setw(8) << coup_s_s_neg[i][j] << "  ";
                sumi += coup_s_s_neg[i][j];
            }
            std::cout << " | sum = " << sumi << "\n";
        }
        std::cout << std::endl;
        std::cout << " coup_s_s_pos = \n";
        for (Uint i = 0; i < 8; ++i)
        {
            double sumi = 0.0;
            for (Uint j = 0; j < 8; ++j)
            {
                std::cout << std::setw(8) << coup_s_s_pos[i][j] << "  ";
                sumi += coup_s_s_pos[i][j];
            }
            std::cout << " | sum = " << sumi << "\n";
        }
        std::cout << std::endl;
    }
    // for (Uint i = 0; i < 8; ++i)
    //     for (Uint j = 0; j < 8; ++j)
    //         coup_s_s[i][j] =  coup_s_s_neg[i][j] + coup_s_s_pos[i][j];
    std::cout << " coup_s_s = \n";
    for (Uint i = 0; i < 8; ++i)
    {
        double sumi = 0.0;
        for (Uint j = 0; j < 8; ++j)
        {
            std::cout << std::setw(8) << coup_s_s[i][j] << "  ";
            sumi += coup_s_s[i][j];
        }
        std::cout << " | sum = " << sumi << "\n";
    }
    std::cout << std::endl;
    if (iscut)
    {
        std::cout << " coup_x_s = \n";
        for (Uint i = 0; i < 8; ++i)
        {
            double sumi = 0.0;
            for (Uint j = 0; j < 8; ++j)
            {
                std::cout << std::setw(8) << coup_x_s[i][j] << "  ";
                sumi += coup_x_s[i][j];
            }
            std::cout << " | sum = " << sumi << "\n";
        }
        std::cout << std::endl;
        std::cout << " coup_s_x = \n";
        for (Uint i = 0; i < 8; ++i)
        {
            double sumi = 0.0;
            for (Uint j = 0; j < 8; ++j)
            {
                std::cout << std::setw(8) << coup_s_x[i][j] << "  ";
                sumi += coup_s_x[i][j];
            }
            std::cout << " | sum = " << sumi << "\n";
        }
        std::cout << std::endl;
        std::cout << " coup_x_x = \n";
        for (Uint i = 0; i < 8; ++i)
        {
            double sumi = 0.0;
            for (Uint j = 0; j < 8; ++j)
            {
                std::cout << std::setw(8) << coup_x_x[i][j] << "  ";
                sumi += coup_x_x[i][j];
            }
            std::cout << " | sum = " << sumi << "\n";
        }
        std::cout << std::endl;
    }

    if (b_ != NULL)
    {
        std::cout << " elvec = \n";
        for (Uint i = 0; i < 8; ++i)
            std::cout << std::setw(8) << elvec[i] << "  ";
        std::cout << std::endl;
        if (iscut)
        {
            std::cout << " elvec_x = \n";
            for (Uint i = 0; i < 8; ++i)
                std::cout << std::setw(8) << elvec_x[i] << "  ";
            std::cout << std::endl;
        }
    }
}

/* --------------- Mass test -------------- */
void MassTestAccumulator_P1SP1TXCL::local_setup (const TetraCL& sit)
{    
    ScopeTimer scopetiming("MassTestAccumulator_P1SP1TXCL::local_setup - total");



    for(Uint i=0; i<8; ++i)
        for(Uint j=0; j<8; ++j)
        {
            intf_coup_s_s[i][j] = 0.0;
            intf_coup_s_x[i][j] = 0.0;
            intf_coup_x_s[i][j] = 0.0;
            intf_coup_x_x[i][j] = 0.0;
            coup_s_s_pos [i][j] = 0.0;
            coup_s_s_neg [i][j] = 0.0;
            
        }


    P1DiscCL::GetGradients(G,det,sit);
    absdet= std::fabs(det);
    const double st_absdet= absdet * (tnew-told);
    //quad_a.assign( sit, &Coeff::DiffusionCoeff, 0.0);                  //for variable diffusion coefficient
    //const double int_a= quad_a.quad( absdet);
    
    LocalP2CL<> locPhi_old, locPhi_new;

    locPhi_old.assign( sit, lsetp2old->Phi, lsetp2old->GetBndData());
    locPhi_new.assign( sit, lsetp2new->Phi, lsetp2new->GetBndData());

    Point4DContainer pcontref;
    GeneralizedPrism4CL refprism4(pcontref);

    // get signs of vertices
    for (Uint i = 0; i < 4; ++i)
    {
        sign[i] = locPhi_old[i] > 0;
        sign[i+4] = locPhi_new[i] > 0;
    }
    
    LocalP2CL<> p0; // values for other time level are zero

    CompositeSTQuadCL<QuadRule> 
        cstquad(refprism4,locPhi_old, locPhi_new, ints_per_space_edge, subtimeintervals);
    
    iscut = cstquad.HasInterface();

    TPSpaceTimeMapping stmap(sit, TimeInterval(told,tnew));
    GridFunctionCL<double> rhspos = cstquad.EvalOnPart( rhs_pos, true, &stmap); 
    GridFunctionCL<double> rhsneg = cstquad.EvalOnPart( rhs_neg, false, &stmap); 

    if (!iscut)
    {
        ScopeTimer scopetiming("MassTestAccumulator_P1SP1TXCL::local_setup - uncut");
        for(Uint ti = PAST; ti<= FUTURE; ti++)
            for(Uint i=0; i<4; ++i)
            {
                GridFunctionCL<double> pi = ti == FUTURE ? cstquad.EvalLinearOnPart( p0, phi[i], sign[0]) 
                    : cstquad.EvalLinearOnPart( phi[i], p0, sign[0]);
                if (b_ != NULL)
                {
                    GridFunctionCL<double> rhspi (sign[0]?rhspos*pi:rhsneg*pi);
                    elvec[ti*4+i] = cstquad.QuadOnPart( rhspi, sign[0]) * st_absdet;
                }
                    
                for(Uint tj = PAST; tj<=FUTURE; tj++)
                    for(Uint j=0; j<4; ++j){
                        GridFunctionCL<double> pj = tj == FUTURE ? cstquad.EvalLinearOnPart( p0, phi[j], sign[0]) 
                            : cstquad.EvalLinearOnPart( phi[j], p0, sign[0]);
                        GridFunctionCL<double> pij (pi * pj);
                        coup_s_s[ti*4+i][tj*4+j] = cstquad.QuadOnPart( pij, sign[0]) * st_absdet;
                    }
            }
    }
    else
    {
        ScopeTimer scopetiming("MassTestAccumulator_P1SP1TXCL::local_setup - cut");
        for(Uint s = 0; s < 2; s++) // both cases: sign = false and sign = true 
        {
            const bool csign = s==1;
            double (* coup_s_s_cur)[8][8] = csign ? &coup_s_s_pos : &coup_s_s_neg;
            double (* elvec_cur)[8] = csign ? &elvec_pos : &elvec_neg;
            for(Uint ti = PAST; ti<= FUTURE; ti++)
                for(Uint i=0; i<4; ++i)
                {
                    GridFunctionCL<double> pi = ti == FUTURE ? cstquad.EvalLinearOnPart( p0, phi[i], csign) 
                        : cstquad.EvalLinearOnPart( phi[i], p0, csign);

                    if (b_ != NULL)
                    {
                        GridFunctionCL<double> rhspi (csign?rhspos*pi:rhsneg*pi);
                        (*elvec_cur)[ti*4+i] = cstquad.QuadOnPart( rhspi, csign) * st_absdet;
                    }

                    for(Uint tj = PAST; tj<=FUTURE; tj++)
                        for(Uint j=0; j<4; ++j){
                            GridFunctionCL<double> pj = tj == FUTURE ? cstquad.EvalLinearOnPart( p0, phi[j], csign) 
                                : cstquad.EvalLinearOnPart( phi[j], p0, csign);
                            GridFunctionCL<double> pij (pi * pj);
                            (*coup_s_s_cur)[ti*4+i][tj*4+j] = cstquad.QuadOnPart( pij, csign) * st_absdet;
                        }
                }
        }
    }
    for(Uint i=0; i<4; ++i)
        UnknownIdx[i]= sit.GetVertex(i)->Unknowns.Exist(idx) ? sit.GetVertex(i)->Unknowns(idx) : NoIdx;

}

/* --------------- Spatial Laplace -------------- */
void SpatialLaplaceAccumulator_P1SP1TXCL::local_setup (const TetraCL& sit)
{    
    ScopeTimer scopetiming("SpatialLaplaceAccumulator_P1SP1TXCL::local_setup - total");

    P1DiscCL::GetGradients(G,det,sit);
    absdet= std::fabs(det);
    const double st_absdet= absdet * (tnew-told);
    //quad_a.assign( sit, &Coeff::DiffusionCoeff, 0.0);                  //for variable diffusion coefficient
    //const double int_a= quad_a.quad( absdet);
    
    LocalP2CL<> locPhi_old, locPhi_new;

    locPhi_old.assign( sit, lsetp2old->Phi, lsetp2old->GetBndData());
    locPhi_new.assign( sit, lsetp2new->Phi, lsetp2new->GetBndData());

    Point4DContainer pcontref;
    GeneralizedPrism4CL refprism4(pcontref);

    // get signs of vertices
    for (Uint i = 0; i < 4; ++i)
    {
        sign[i] = locPhi_old[i] > 0;
        sign[i+4] = locPhi_new[i] > 0;
    }
    
    CompositeSTQuadCL<QuadRule> 
        cstquad(refprism4,locPhi_old, locPhi_new, ints_per_space_edge, subtimeintervals);
    
    iscut = cstquad.HasInterface();

    TPSpaceTimeMapping stmap(sit, TimeInterval(told,tnew));
    GridFunctionCL<double> rhspos = cstquad.EvalOnPart( rhs_pos, true, &stmap); 
    GridFunctionCL<double> rhsneg = cstquad.EvalOnPart( rhs_neg, false, &stmap); 

    const double alpha = 1.0;
    
    if (!iscut)
    {
        ScopeTimer scopetiming("SpatialLaplaceAccumulator_P1SP1TXCL::local_setup - uncut");

        for(Uint i=0; i<4; ++i)
        {
            for(Uint j=0; j<i; ++j){
                coup_s_s_space[i][j] = alpha*inner_prod( G[i], G[j])/6.0; //diffusion
                coup_s_s_space[j][i] = coup_s_s_space[i][j];
            }
            coup_s_s_space[i][i] = alpha*inner_prod( G[i], G[i])/6.0; //diffusion
        }

        for(Uint ti = PAST; ti<= FUTURE; ti++)
            for(Uint i=0; i<4; ++i)
                for(Uint tj = PAST; tj<=FUTURE; tj++)
                    for(Uint j=0; j<4; ++j)
                        coup_s_s[ti*4+i][tj*4+j] = coup_s_s_space[i][j] * mass_1D_p1p1[ti][tj] * st_absdet; 
    }
    else
    {
        ScopeTimer scopetiming("SpatialLaplaceAccumulator_P1SP1TXCL::local_setup - cut");

        // LocalP2CL<Point3DCL> gradi[4];
        // for (Uint i = 0; i < 4; ++i)
        // {
        //     gradi[i] = G[i];
        // }
    
        LocalP2CL<> p0; // values for other time level are zero
        LocalP2CL<Point3DCL> vecp0; // values for other time level are zero
        LocalP2CL<> constone(1.0);

        for(Uint s = 0; s < 2; s++) // both cases: sign = false and sign = true 
        {
            const bool csign = s==1;
            const double alpha_cur = 1.0; //csign ? alpha_pos : alpha_neg;
            double (* coup_s_s_cur)[8][8] = csign ? &coup_s_s_pos : &coup_s_s_neg;
            // double (* elvec_cur)[8] = csign ? &elvec_pos : &elvec_neg;
            for(Uint ti = PAST; ti<= FUTURE; ti++)
            {
                GridFunctionCL<double> pti = ti == FUTURE ? cstquad.EvalLinearOnPart( p0, constone, csign) 
                    : cstquad.EvalLinearOnPart( constone, p0, csign);

                // if (b_ != NULL)
                // {
                //     GridFunctionCL<double> rhspi (csign?rhspos*pi:rhsneg*pi);
                //     (*elvec_cur)[ti*4+i] = cstquad.QuadOnPart( rhspi, csign) * st_absdet;
                // }

                for(Uint tj = PAST; tj<=FUTURE; tj++)
                {
                    GridFunctionCL<double> ptj = tj == FUTURE ? cstquad.EvalLinearOnPart( p0, constone, csign) 
                        : cstquad.EvalLinearOnPart( constone, p0, csign);
                    GridFunctionCL<double> ptitj (pti * ptj);
                    const double stint_t = cstquad.QuadOnPart( ptitj, csign) * st_absdet;
                    for(Uint i=0; i<4; ++i)
                        for(Uint j=0; j<4; ++j)
                            (*coup_s_s_cur)[ti*4+i][tj*4+j] =  stint_t * alpha_cur * inner_prod( G[i], G[j]) ;
                }
            }
        }
    }
    for(Uint i=0; i<4; ++i)
        UnknownIdx[i]= sit.GetVertex(i)->Unknowns.Exist(idx) ? sit.GetVertex(i)->Unknowns(idx) : NoIdx;

}

STTransportVolumeAccumulator_P1SP1TXCL::
STTransportVolumeAccumulator_P1SP1TXCL (const MultiGridCL& MG, 
                                        const BndDataCL<> * BndData_neg_in, 
                                        const BndDataCL<> * BndData_pos_in, 
                                        const LevelsetP2CL * lsetp2old_in,
                                        const LevelsetP2CL * lsetp2new_in,
                                        instat_scalar_fun_ptr lset_fpt_in,
                                        instat_scalar_fun_ptr rhs_neg_in,
                                        instat_scalar_fun_ptr rhs_pos_in,
                                        STVelocityContainer& convection_in,
                                        MatrixCL* Amat, VecDescCL* b, 
                                        const IdxDescCL& RowIdx, const IdxDescCL& ColIdx, 
                                        const double t1, const double t2,
                                        const VecDescCL & oldsol_neg, const VecDescCL & oldsol_pos,
                                        const ParamCL::ptree_type & P)
: base_(MG,BndData_neg_in,BndData_pos_in,lsetp2old_in,lsetp2new_in,Amat,b,RowIdx,ColIdx,t1,t2,P),
    lset_fpt(lset_fpt_in), rhs_neg(rhs_neg_in),rhs_pos(rhs_pos_in),convection(convection_in),
    alpha_neg(P.get<double>("DiffNeg")),alpha_pos(P.get<double>("DiffPos")),
    // beta_neg(P.get<double>("HNeg")),beta_pos(P.get<double>("HPos")),
    delta_beta(beta_pos-beta_neg),
    lambda_stab(P.get<double>("NitschePenalty")),
    vmax(P.get<double>("MaxVelocity")),
    oldsolution_neg(oldsol_neg),oldsolution_pos(oldsol_pos),
    use_stkappa_rule(P.get<int>("KappaRule")==1),
    use_space_kappa_rule(P.get<int>("KappaRule")==2),
    antisym_conv(P.get<int>("AntiSymmetricConvection")==1),
    SD_stab(P.get<double>("SD",0.0)),
    quad1d_rhs_order(P.get<int>("RhsApproximationOrderInTime")),
    quad1d_conv_order(P.get<int>("ConvectionApproximationOrderInTime")),
    adjoint(P.get<int>("Adjoint",0)==1)
{
    if (antisym_conv)
        throw DROPSErrCL("no antisym_conv version right now...");
}


/* --------------- STVolume Transport term (time deriv., convection, diffusion)  -------------- */
void STTransportVolumeAccumulator_P1SP1TXCL::local_setup (const TetraCL& sit)
{    
//    ScopeTimer scopetiming("STTranspVolAccu_P1SP1TXCL::local_setup - total");
    
    GetTrafoTr( T, det, sit);
    P1DiscCL::GetGradients(G,det,sit);

    absdet= std::fabs(det);

    const double h = std::pow(absdet,1.0/3.0);
    const double dt = tnew-told;
    const double st_absdet= absdet * dt;
    
    LocalP2CL<> locPhi_old, locPhi_new;

    locPhi_old.assign( sit, lsetp2old->Phi, lsetp2old->GetBndData());
    locPhi_new.assign( sit, lsetp2new->Phi, lsetp2new->GetBndData());


    // get signs of vertices
    for (Uint i = 0; i < 4; ++i)
    {
        sign[i] = locPhi_old[i] > 0;
        sign[i+4] = locPhi_new[i] > 0;
    }

    CompositeSTQuadCL<QuadRule> * p_cstquad = NULL;

    bool inside_outer_band = false;
    bool outside_inner_band = false;
    
    for (int i = 0; i < 10; ++i)
    {
        if (locPhi_old[i] < dt*vmax || locPhi_new[i] < dt*vmax) 
            inside_outer_band = true;
        if (locPhi_old[i] > -dt*vmax || locPhi_new[i] >- dt*vmax) 
            outside_inner_band = true;
    }

    if (inside_outer_band && outside_inner_band)
    {

        Point4DContainer pcontref;
        GeneralizedPrism4CL refprism4(pcontref);

        if (lset_fpt == NULL)
            p_cstquad = new CompositeSTQuadCL<QuadRule>(refprism4,locPhi_old, locPhi_new,
                                                                      ints_per_space_edge, subtimeintervals);
        else
            p_cstquad = new CompositeSTQuadCL<QuadRule>(sit, TimeInterval(told,tnew), lset_fpt, 
                                                                      ints_per_space_edge, subtimeintervals);
        iscut = p_cstquad->HasInterface();
    }
    else
    {
        iscut = false;
    }

    if (!iscut)
    {
        ScopeTimer scopetiming("STTranspVolAccu_P1SP1TXCL::local_setup - uncut");

        const double & alpha (sign[0] ? alpha_pos : alpha_neg);
        const double & beta (sign[0] ? beta_pos : beta_neg);
        
        // only spatial integrals 
        double coup_s_s_space_lap[4][4];
        for(Uint i=0; i<4; ++i)
        {
            for(Uint j=0; j<i; ++j){
                coup_s_s_space[i][j] = 1.0/beta * P1DiscCL::GetMass(i,j)*absdet;
                coup_s_s_space[j][i] = coup_s_s_space[i][j];
                coup_s_s_space_lap[i][j] = alpha/beta * inner_prod( G[i], G[j])/6.0; //mass
                coup_s_s_space_lap[j][i] = coup_s_s_space_lap[i][j];
            }
            coup_s_s_space[i][i] = 1.0/beta * P1DiscCL::GetMass(i,i)*absdet; //mass
            coup_s_s_space_lap[i][i] = alpha/beta*inner_prod( G[i], G[i])/6.0; //diffusion
        }

        //values for  \int_T f^(t_tl) phi_i dx where t_tl is told, tmid or tnew
        double f_tl_phi_i [quad1d_max_n_points][4]; 
        //values for  \int_T b^(t_tl) \nabla phi_j phi_i dx where t_tl is told, tmid or tnew
        double b_tl_phi_i_grad_phi_j [quad1d_max_n_points][4][4]; 

        //values for  bstar^(t_tl) where t_tl is told, tmid or tnew 
        Quad5CL<double> dirdev_star [quad1d_max_n_points][8]; 
        
        for (int k = 0; k < quad1d_rhs_order+1; k++)
        {
            const double time = told + MomentsQuad1D::GetPoint(quad1d_rhs_order,k)*(tnew-told);
            Quad5CL<> q_f_now (sit,  sign[0] ? rhs_pos : rhs_neg, time);
            for (Uint i = 0; i < 4; ++i)
                f_tl_phi_i[k][i] = Quad5CL<>(q_f_now*q_phi[i]).quad(absdet);
        }


        Quad5CL<Point3DCL> q_b_now;
        Quad5CL<Point3DCL> q_b_old, q_b_new;
        if (convection.VelocityAsP2())
        {
            q_b_old.assign(sit, convection.GetVelocityAsP2_Old());
            q_b_new.assign(sit, convection.GetVelocityAsP2_New());
        }

        double SD_stab_param = 0.0;
        double vmax = 0.0;
        for (int k = 0; k < quad1d_conv_order+1 ; k++)
        {
            const double w = MomentsQuad1D::GetPoint(quad1d_conv_order,k);
            const double time = told + w*(tnew-told);
            if (convection.VelocityAsFunctionPointer())
                q_b_now.assign(sit, convection.GetVelocityAsFunctionPointer(), time);
            else
                q_b_now = (1.0-w)* q_b_old + w * q_b_new;
            
            for (Uint ip = 0; ip < Quad5DataCL::NumNodesC; ++ip)
                vmax = std::max(vmax, q_b_now[ip].norm());

            for (Uint j = 0; j < 4; ++j)
            {
                Quad5CL<double> q_b_now_grad_phi_j = Quad5CL<double>(dot(G[j],q_b_now));
                for (Uint i = 0; i < 4; ++i)
                    if (adjoint)
                        b_tl_phi_i_grad_phi_j[k][j][i] = -Quad5CL<>(q_b_now_grad_phi_j*q_phi[i]).quad(absdet);
                    else
                        b_tl_phi_i_grad_phi_j[k][i][j] = Quad5CL<>(q_b_now_grad_phi_j*q_phi[i]).quad(absdet);

                if (SD_stab>0.0)
                    for(Uint tj = PAST; tj<= FUTURE; tj++)
                    {
                        double psi_tj =  tj == PAST ? 1.0 - w : w;
                        double dpsidt_tj =  tj == PAST ? -1.0/dt : 1.0/dt;
                        dirdev_star[k][tj*4+j] = q_b_now_grad_phi_j * psi_tj + q_phi[j] * dpsidt_tj;
                    }
            }

            
        }

        //SD param calculation
        if (SD_stab > 0)
        {
            const double h = std::pow(absdet,1.0/3.0);
            // const double g1 = (2.0/dt)*(2.0/dt) + (2*vmax/h)*(2*vmax/h)+9*(4*alpha/(h*h))*(4*alpha/(h*h));
            // const double g1 = (2*vmax/h)*(2*vmax/h); // unused
            const double gT = 2* h/vmax;

            if (2 * vmax  * h > alpha)
            {
                // std::cout << " gT = " << gT << std::endl;
                SD_stab_param = SD_stab * gT;
            }
            else 
                SD_stab_param = 0.0;
        }


        // only space-time integrals 
        for(Uint ti = PAST; ti<= FUTURE; ti++)
        {
            for(Uint i=0; i<4; ++i)
            {
                for(Uint tj = PAST; tj<=FUTURE; tj++)
                {
                    for(Uint j=0; j<4; ++j)
                    {
                        if (adjoint)
                            coup_s_s[ti*4+i][tj*4+j] = -coup_s_s_space[j][i] * dudtv_1D_p1p1[tj][ti]; 
                        else
                            coup_s_s[ti*4+i][tj*4+j] = coup_s_s_space[i][j] * dudtv_1D_p1p1[ti][tj]; 
                        coup_s_s[ti*4+i][tj*4+j] += coup_s_s_space_lap[i][j] * mass_1D_p1p1[ti][tj]  * st_absdet;
                        
                        for (int k = 0; k < quad1d_conv_order+1; ++k)
                            coup_s_s[ti*4+i][tj*4+j] += dt * 1.0/beta * b_tl_phi_i_grad_phi_j[k][i][j] * 
                                MomentsQuad1D::GetIntegralPkP1P1(quad1d_conv_order,k,ti,tj);
                        
                        // STREAMLIN DIFFUSION STABILIZATION
                        if (SD_stab > 0 && SD_stab_param > 0)
                            for (int k = 0; k < quad1d_conv_order+1; ++k)
                            {
                                const double spaceint = Quad5CL<double>(dirdev_star[k][ti*4+i]*dirdev_star[k][tj*4+j]).quad(absdet);
                                coup_s_s[ti*4+i][tj*4+j] += SD_stab_param * dt * 1.0/beta * MomentsQuad1D::GetIntegrationWeight(quad1d_conv_order,k) * spaceint;
                            }
                        
                    }
                }

                // source term r.h.s.
                if (b_ != NULL)
                {
                    elvec[ti*4+i] = 0;
                    for (int k = 0; k < quad1d_rhs_order+1; ++k)
                        elvec[ti*4+i] += dt * f_tl_phi_i[k][i] * 
                            MomentsQuad1D::GetIntegralPkP1(quad1d_rhs_order,k,ti);
                }

            }
        }

        
        // upwind term l.h.s.
        for(Uint i=0; i<4; ++i)
            for(Uint j=0; j<4; ++j)
                if (adjoint)
                    coup_s_s[4+i][4+j] += coup_s_s_space[i][j]; 
                else
                    coup_s_s[i][j] += coup_s_s_space[i][j]; 
        
        // upwind term r.h.s.
        if (b_ != NULL)
        {
            double elvec_space_old[4];
            LocalP1CL<> oldtr(sit,sign[0]?oldsolution_pos:oldsolution_neg,sign[0]?*BndData_pos:*BndData_neg);
            Quad5CL<> fold( oldtr);
            for(Uint i=0; i<4; ++i)
            {
                LocalP1CL<> pip1; pip1[i] = 1.0;
                Quad5CL<> pi(pip1);
                elvec_space_old[i] = Quad5CL<>(pi*fold).quad(absdet); 
            }

            for(Uint i=0; i<4; ++i)
                elvec[0*4+i] += elvec_space_old[i];
        }
    }
    else
    {
        ScopeTimer scopetiming("STTranspVolAccu_P1SP1TXCL::local_setup - cut");

        TPSpaceTimeMapping stmap(sit, TimeInterval(told,tnew));

        LocalP2CL<Point3DCL> vecp0(MakePoint3D(0.0,0.0,0.0)); // values for other time level are zero
        LocalP2CL<> constone(1.0);
        LocalP2CL<> p0(0.0); // values for other time level are zero

        CompositeSTQuadCL<QuadRule> & cstquad = *p_cstquad;

        GridFunctionCL<double> rhspos = cstquad.EvalOnPart( rhs_pos, true, &stmap); 
        GridFunctionCL<double> rhsneg = cstquad.EvalOnPart( rhs_neg, false, &stmap);
 
        GridFunctionCL<double> gfone_neg = cstquad.EvalLinearOnPart( constone, constone, false); 
        GridFunctionCL<double> gfone_pos = cstquad.EvalLinearOnPart( constone, constone, true); 
        
        GridFunctionCL<double> kappa_neg;
        GridFunctionCL<double> kappa_pos;

        double const_kappa_neg = 0.5;
        double const_kappa_pos = 0.5;
        if(use_stkappa_rule)
        {
            const_kappa_neg = cstquad.QuadOnPart(gfone_neg,false);
            const_kappa_pos = cstquad.QuadOnPart(gfone_pos,true);
            const double sum = const_kappa_pos + const_kappa_neg;
            const_kappa_neg /= sum;
            const_kappa_pos /= sum;
        }


        if (use_space_kappa_rule)
            debug_kappa(locPhi_old,locPhi_new,cstquad.GetInterfaceIntegrationPoints(),kappa_neg,kappa_pos);

        const GridFunctionCL<double> av_alpha (alpha_neg/beta_neg*kappa_neg + alpha_pos/beta_pos*kappa_pos);
        // this is used for { a \nabla u \cdot \bn } = {a} \nabla u \cdot n for continuous functions!
        
        const double alphamean = 0.5 * alpha_neg/beta_neg + 0.5 * alpha_pos/beta_pos;

        // integration of space-time 4D-measure
        {

            GridFunctionCL<Point3DCL> velneg (rhsneg.size());
            GridFunctionCL<Point3DCL> velpos (rhspos.size());

            if (convection.VelocityAsFunctionPointer())
            {
                velneg = cstquad.EvalOnPart( convection.GetVelocityAsFunctionPointer(), false, &stmap);
                velpos = cstquad.EvalOnPart( convection.GetVelocityAsFunctionPointer(), true, &stmap); 
            }
            else
            {
                LocalP2CL<Point3DCL> veloldp2 (sit, convection.GetVelocityAsP2_Old());
                LocalP2CL<Point3DCL> velnewp2 (sit, convection.GetVelocityAsP2_New());
                velneg = cstquad.EvalLinearOnPart( veloldp2, velnewp2, false); 
                velpos = cstquad.EvalLinearOnPart( veloldp2, velnewp2, true); 
            }


            ScopeTimer scopetiming("STTranspVolAccu_P1SP1TXCL::local_setup - cut - 4Dvols");
            for(Uint s = 0; s < 2; s++) // both cases: sign = false and sign = true 
            {
                const bool csign = s==1;

                double (* coup_s_s_cur)[8][8] = csign ? &coup_s_s_pos : &coup_s_s_neg;
                double (* elvec_cur)[8] = csign ? &elvec_pos : &elvec_neg;
                const GridFunctionCL<Point3DCL>& vel_cur (csign ? velpos : velneg);
                const double & alpha (csign ? alpha_pos : alpha_neg);
                const double & beta (csign ? beta_pos : beta_neg);

                GridFunctionCL<double> psi_PAST = cstquad.EvalLinearOnPart( constone, p0, csign);
                GridFunctionCL<double> psi_FUTURE = cstquad.EvalLinearOnPart( p0, constone, csign); 


                GridFunctionCL<double> stphi[8] =
                    {cstquad.EvalLinearOnPart( phi[0], p0, csign),
                     cstquad.EvalLinearOnPart( phi[1], p0, csign),
                     cstquad.EvalLinearOnPart( phi[2], p0, csign),
                     cstquad.EvalLinearOnPart( phi[3], p0, csign),
                     cstquad.EvalLinearOnPart( p0, phi[0], csign),
                     cstquad.EvalLinearOnPart( p0, phi[1], csign),
                     cstquad.EvalLinearOnPart( p0, phi[2], csign),
                     cstquad.EvalLinearOnPart( p0, phi[3], csign)};

                GridFunctionCL<double> b_gradphi[4] =
                    {GridFunctionCL<double>(dot(G[0],vel_cur)),
                     GridFunctionCL<double>(dot(G[1],vel_cur)),
                     GridFunctionCL<double>(dot(G[2],vel_cur)),
                     GridFunctionCL<double>(dot(G[3],vel_cur))};

                GridFunctionCL<double> const_in_time_phi[4] =
                    {cstquad.EvalLinearOnPart( phi[0], phi[0], csign),
                     cstquad.EvalLinearOnPart( phi[1], phi[1], csign),
                     cstquad.EvalLinearOnPart( phi[2], phi[2], csign),
                     cstquad.EvalLinearOnPart( phi[3], phi[3], csign)}; 


                double vmax = 0.0;
                double SD_stab_param = 0.0;
                if (SD_stab > 0.0)
                {
                    for (Uint ip = 0; ip < vel_cur.size(); ++ip)
                        vmax = std::max(vmax, vel_cur[ip].norm());
                    const double h = std::pow(absdet,1.0/3.0);
                    // const double g1 = (2.0/dt)*(2.0/dt) + (2*vmax/h)*(2*vmax/h)+9*(4*alpha/(h*h))*(4*alpha/(h*h));
                    const double gT = 2*h / vmax ;

                    if (2 * vmax  * h > alpha)
                    {
                        // std::cout << " gT = " << gT << std::endl;
                        SD_stab_param = 100.0 * gT;
                    }
                    else 
                        SD_stab_param = 0.0;
                }



                for(Uint ti = PAST; ti<= FUTURE; ti++)
                {
                    const double timeder_factor_i = ((ti == FUTURE) ? 1.0/dt : -1.0/dt );
                    for(Uint i=0; i<4; ++i)
                    {
                        GridFunctionCL<double> & phi_i = stphi[4*ti+i];
                        GridFunctionCL<double> & psi_ti = ti == FUTURE ? psi_FUTURE : psi_PAST;

                        if (b_ != NULL)
                        {
                            GridFunctionCL<double> rhspi (csign?rhspos*phi_i:rhsneg*phi_i);
                            (*elvec_cur)[ti*4+i] = cstquad.QuadOnPart( rhspi, csign) * st_absdet;
                        }

                        for(Uint tj = PAST; tj<=FUTURE; tj++)
                        {

                            const double timeder_factor_j = ((tj == FUTURE) ? 1.0/dt : -1.0/dt );
                            GridFunctionCL<double> & psi_tj = tj == FUTURE ? psi_FUTURE : psi_PAST;
                            for(Uint j=0; j<4; ++j){
                                GridFunctionCL<double> & phi_j = stphi[4*tj+j];
                                GridFunctionCL<double> gradugradv ( (alpha * inner_prod( G[i], G[j])) * psi_tj * psi_ti);
                                GridFunctionCL<double> volparts ( psi_ti.size());

                                if (adjoint)
                                {
                                    GridFunctionCL<double> dudtv ( -timeder_factor_i * phi_j * const_in_time_phi[i]);
                                    GridFunctionCL<double> b_gradu_v ( -psi_ti*b_gradphi[i] * phi_j);
                                    volparts = GridFunctionCL<double>(dudtv+gradugradv+b_gradu_v);
                                }
                                else
                                {
                                    GridFunctionCL<double> dudtv ( timeder_factor_j * phi_i * const_in_time_phi[j]);
                                    GridFunctionCL<double> b_gradu_v ( psi_tj*b_gradphi[j]*phi_i);
                                    volparts = GridFunctionCL<double>(dudtv+gradugradv+b_gradu_v);
                                }

                                if (SD_stab > 0 && SD_stab_param > 0)
                                {
                                    GridFunctionCL<double> dirdev_star_j ( timeder_factor_j * const_in_time_phi[j] +  psi_tj*b_gradphi[j]);
                                    GridFunctionCL<double> dirdev_star_i ( timeder_factor_i * const_in_time_phi[i] +  psi_ti*b_gradphi[i]);
                                    volparts += SD_stab_param * dirdev_star_j * dirdev_star_i;
                                }
                                (*coup_s_s_cur)[ti*4+i][tj*4+j] = 1.0/beta * cstquad.QuadOnPart( volparts, csign) * st_absdet;
                            }
                        }
                    } // end for i
                } // end for ti
            } // end for s (neg, pos)
        }//end ( integration of space-time 4D-measure)


        {
            ScopeTimer scopetiming("STTranspVolAccu_P1SP1TXCL::local_setup - cut - pasttrace");

            GridFunctionCL<> pp;
            GridFunctionCL<> pi, fneg, fpos;
            QuadDomainCL q2dom_;
            PrincipalLatticeCL lat( PrincipalLatticeCL::instance( ints_per_space_edge ));
            TetraPartitionCL partition_;
            std::valarray<double> ls_loc_( ((ints_per_space_edge+1) * (ints_per_space_edge+2) * (ints_per_space_edge+3)) / 6);

            double pipj_space_neg[4][4];
            double pipj_space_pos[4][4];
            double fi_space_neg[4];
            double fi_space_pos[4];
            for (Uint timedir = PAST; timedir <= FUTURE; timedir++) // => timedir = PAST
            {
                if (!adjoint && timedir == FUTURE) continue;

                // std::cout << " timedir = " << timedir << std::endl;

                if (timedir == FUTURE)
                    evaluate_on_vertexes( lsetp2new->GetSolution(), sit, lat, Addr( ls_loc_));
                else
                    evaluate_on_vertexes( lsetp2old->GetSolution(), sit, lat, Addr( ls_loc_));

                partition_.make_partition<SortedVertexPolicyCL, MergeCutPolicyCL>( lat, ls_loc_);
                make_CompositeQuad5Domain( q2dom_, partition_);
                LocalP1CL<> oldtr_neg(sit,oldsolution_neg,*BndData_neg);
                LocalP1CL<> oldtr_pos(sit,oldsolution_pos,*BndData_pos);
                for(Uint i=0; i<4; ++i){

                    for(Uint j=0; j<i; ++j){
                        resize_and_evaluate_on_vertexes( phi_i_phi_j[i][j], q2dom_, pp);
                        pipj_space_neg[i][j] = quad( pp , absdet , q2dom_ , NegTetraC);
                        pipj_space_neg[j][i] = pipj_space_neg[i][j];
                        pipj_space_pos[i][j] = quad( pp , absdet , q2dom_ , PosTetraC);
                        pipj_space_pos[j][i] = pipj_space_pos[i][j];
                    }
                    resize_and_evaluate_on_vertexes( phi_i_phi_j[i][i], q2dom_, pp);
                    pipj_space_neg[i][i] = quad( pp , absdet , q2dom_ , NegTetraC);
                    pipj_space_pos[i][i] = quad( pp , absdet , q2dom_ , PosTetraC);

                    if (timedir == PAST)
                    {
                        resize_and_evaluate_on_vertexes( oldtr_pos, q2dom_, fpos);
                        resize_and_evaluate_on_vertexes( oldtr_neg, q2dom_, fneg);
                        resize_and_evaluate_on_vertexes( phi[i], q2dom_, pi);
                        fi_space_neg[i] = quad(fneg*pi, absdet, q2dom_, NegTetraC);
                        fi_space_pos[i] = quad(fpos*pi, absdet, q2dom_, PosTetraC);
                    }
                }

                for(Uint i=0; i<4; ++i){
                    if (!adjoint || timedir == FUTURE)
                        for(Uint j=0; j<4; ++j){
                            coup_s_s_neg[4*timedir+i][4*timedir+j] += 1.0/beta_neg * pipj_space_neg[i][j];
                            coup_s_s_pos[4*timedir+i][4*timedir+j] += 1.0/beta_pos * pipj_space_pos[i][j];
                        }
                    if (timedir == PAST)
                    {
                        elvec_neg[4*timedir+i] += fi_space_neg[i]; //note that uold is in "untransf." state
                        elvec_pos[4*timedir+i] += fi_space_pos[i]; //note that uold is in "untransf." state
                    }
                }
            }
        }
        // Space-time interface terms (Nitsche)
        {
            ScopeTimer scopetiming("STTranspVolAccu_P1SP1TXCL::local_setup - cut - stintface");

            //space time normals on reference prism pointing fr. neg to pos :
            const GridFunctionCL<Point4DCL> & n_st_ref (cstquad.GetNormalsOnInterface()); 
            //not correctly scaled space time normals on physical prism pointing fr. neg to pos :
            GridFunctionCL<Point4DCL> n_st = eval(transform_normals_in_spacetime,T,absdet,dt,n_st_ref);
            //after next call space time normal is now correctly scaled (||n||=1) 
            //and the ratio between the measure on the physical prism to reference prism 
            //calculated. Additionally the space-measure-correction-factor nu is mulitplied.
            GridFunctionCL<double> meas_nu = apply_and_eval(calc_nu_measure_and_normalize,n_st);
            //The corresponding space normal is extracted from the (phys.) space time normals
            const GridFunctionCL<Point3DCL> spacenormals = eval(calc_spacenormal_from_stnormal,n_st);

            for(Uint i=0; i<8; ++i)
                for(Uint j=0; j<8; ++j)
                {
                    intf_coup_s_s[i][j] = 0.0;
                    intf_coup_s_x[i][j] = 0.0;
                    intf_coup_x_s[i][j] = 0.0;
                    intf_coup_x_x[i][j] = 0.0;
                }
            
            const GridFunctionCL<double> phi_tr[8] =
                { cstquad.EvalLinearOnInterface( phi[0], p0),
                  cstquad.EvalLinearOnInterface( phi[1], p0),
                  cstquad.EvalLinearOnInterface( phi[2], p0),
                  cstquad.EvalLinearOnInterface( phi[3], p0),
                  cstquad.EvalLinearOnInterface( p0, phi[0]),
                  cstquad.EvalLinearOnInterface( p0, phi[1]),
                  cstquad.EvalLinearOnInterface( p0, phi[2]),
                  cstquad.EvalLinearOnInterface( p0, phi[3])};

            const GridFunctionCL<double> gradn[4] =
                { GridFunctionCL<double>(dot(G[0],spacenormals)),
                  GridFunctionCL<double>(dot(G[1],spacenormals)),
                  GridFunctionCL<double>(dot(G[2],spacenormals)),
                  GridFunctionCL<double>(dot(G[3],spacenormals))};


            const GridFunctionCL<double> psi_tr_PAST =  cstquad.EvalLinearOnInterface( constone, p0);
            const GridFunctionCL<double> psi_tr_FUTURE = cstquad.EvalLinearOnInterface( p0, constone);

            const double const_av_alpha = const_kappa_neg*alpha_neg/beta_neg + const_kappa_pos*alpha_pos/beta_pos;

            for(Uint ti = PAST; ti<= FUTURE; ti++)
                for(Uint i=0; i<4; ++i)
                {
                    const GridFunctionCL<double> & gradni = gradn[i];
                    const GridFunctionCL<double> & phi_i_tr = phi_tr[ti*4+i];
                    const GridFunctionCL<double> & psi_ti_tr = ti == FUTURE ? psi_tr_FUTURE : psi_tr_PAST;
                
                    const double beta_not_i = sign[i+4*ti] ? beta_neg : beta_pos;
                    const double alpha_not_i = sign[i+4*ti] ? -alpha_neg : alpha_pos;
                    const double const_kappa_not_i = sign[i+4*ti] ? const_kappa_neg : const_kappa_pos;
                    
                    for(Uint tj = PAST; tj<=FUTURE; tj++)
                        for(Uint j=0; j<4; ++j)
                        {

                            //exploiting symmetry of nitsche term
                            if (tj*4+j > ti*4+i)
                                continue;

                            const GridFunctionCL<double> & gradnj = gradn[j];
                            const GridFunctionCL<double> & phi_j_tr = phi_tr[tj*4+j];
                            const GridFunctionCL<double> & psi_tj_tr = tj == FUTURE ? psi_tr_FUTURE : psi_tr_PAST;

                            const double beta_not_j = sign[j+4*tj] ? beta_neg : beta_pos;
                            const double alpha_not_j = sign[j+4*tj] ? -alpha_neg : alpha_pos;

                            const double const_kappa_not_j = sign[j+4*tj] ? const_kappa_neg : const_kappa_pos;

                            double A;

                            const GridFunctionCL<double> intfparts_A (meas_nu*phi_i_tr*phi_j_tr);
                            A = (lambda_stab*alphamean/h) * cstquad.QuadOnInterface(intfparts_A);

                            //const GridFunctionCL<double> intfparts_B (gradnj*meas_nu*phi_i_tr*psi_tj_tr);
                        
                            //\int nu * alpha/beta * lambda/h [u] * [v] ds 
                            
                            // the following four lines are only interesting for the "untransformed" case
                            // intf_coup_s_s[4*ti+i][4*tj+j] += delta_beta*delta_beta*A;
                            // intf_coup_s_x[4*ti+i][4*tj+j] += delta_beta*beta_not_j*A;
                            // intf_coup_x_s[4*ti+i][4*tj+j] += beta_not_i*delta_beta*A;
                            // intf_coup_x_x[4*ti+i][4*tj+j] += beta_not_i*beta_not_j*A;

                            intf_coup_x_x[4*ti+i][4*tj+j] += A;

                            if (!use_space_kappa_rule){

                                const GridFunctionCL<double> intfparts_BT (gradni*meas_nu*psi_ti_tr*phi_j_tr);
                                const GridFunctionCL<double> intfparts_B (gradnj*meas_nu*phi_i_tr*psi_tj_tr);

                                const double B = cstquad.QuadOnInterface(intfparts_B);
                                const double BT = cstquad.QuadOnInterface(intfparts_BT);

                                //\int nu * {a/b dudn} * [v] ds 
                                intf_coup_x_s[4*ti+i][4*tj+j] += const_av_alpha * B;
                                intf_coup_x_x[4*ti+i][4*tj+j] += const_kappa_not_j*alpha_not_j/beta_not_j * B;

                                // //\int nu * {a/b dvdn} * [u] ds
                                intf_coup_s_x[4*ti+i][4*tj+j] += const_av_alpha * BT ;
                                intf_coup_x_x[4*ti+i][4*tj+j] += const_kappa_not_i*alpha_not_i/beta_not_i * BT;

                                
                            }
                            else // should be standard (sadly)
                            {

                                const GridFunctionCL<double> & kappa_not_j = sign[j+4*tj] ? kappa_neg : kappa_pos;
                                const GridFunctionCL<double> & kappa_not_i = sign[i+4*ti] ? kappa_neg : kappa_pos;

                                const GridFunctionCL<double> intfparts_B_av (av_alpha*gradnj*meas_nu*phi_i_tr*psi_tj_tr);
                                const GridFunctionCL<double> intfparts_B_knj (kappa_not_j*gradnj*meas_nu*phi_i_tr*psi_tj_tr);
                                const GridFunctionCL<double> intfparts_BT_av (av_alpha*gradni*meas_nu*phi_j_tr*psi_ti_tr);
                                const GridFunctionCL<double> intfparts_BT_kni (kappa_not_i*gradni*meas_nu*phi_j_tr*psi_ti_tr);

                                const double B_av = cstquad.QuadOnInterface(intfparts_B_av);
                                const double B_knj = cstquad.QuadOnInterface(intfparts_B_knj);
                                const double BT_av = cstquad.QuadOnInterface(intfparts_BT_av);
                                const double BT_kni = cstquad.QuadOnInterface(intfparts_BT_kni);
                                //\int nu * {a/b dudn} * [v] ds 
                                intf_coup_x_s[4*ti+i][4*tj+j] += B_av;
                                intf_coup_x_x[4*ti+i][4*tj+j] += alpha_not_j/beta_not_j * B_knj;

                                //\int nu * {a/b dvdn} * [u] ds
                                intf_coup_s_x[4*ti+i][4*tj+j] += BT_av;
                                intf_coup_x_x[4*ti+i][4*tj+j] += alpha_not_i/beta_not_i * BT_kni;
                            }
                            //exploiting symmetry of nitsche term
                            if (tj*4+j < ti*4+i)
                            {
                                intf_coup_s_s[4*tj+j][4*ti+i] = intf_coup_s_s[4*ti+i][4*tj+j] ;
                                intf_coup_x_s[4*tj+j][4*ti+i] = intf_coup_s_x[4*ti+i][4*tj+j] ;
                                intf_coup_s_x[4*tj+j][4*ti+i] = intf_coup_x_s[4*ti+i][4*tj+j] ;
                                intf_coup_x_x[4*tj+j][4*ti+i] = intf_coup_x_x[4*ti+i][4*tj+j] ;
                            }
                        }
                }
        }
    }

    for(Uint i=0; i<4; ++i)
        UnknownIdx[i]= sit.GetVertex(i)->Unknowns.Exist(idx) ? sit.GetVertex(i)->Unknowns(idx) : NoIdx;


    if (inside_outer_band && outside_inner_band)
        delete p_cstquad;
}


// -----------------------------------------------------------------------------

STGeomApproxTestAccumulatorCL::STGeomApproxTestAccumulatorCL(const MultiGridCL& MG, 
                                                             const LevelsetP2CL * lsetp2old_in,
                                                             const LevelsetP2CL * lsetp2new_in, 
                                                             instat_scalar_fun_ptr lset_fpt_in, 
                                                             const double t1, const double t2,
                                                             const ParamCL::ptree_type & P):
    MG_(MG), lsetp2old(lsetp2old_in),lsetp2new(lsetp2new_in),
    lset_fpt(lset_fpt_in),
    ints_per_space_edge(P.get<int>("Quadrature.SubIntervalsPerEdge")),
    subtimeintervals(P.get<int>("Quadrature.SubTimeIntervals")),
    told(t1), tnew(t2)
{
    
}

void STGeomApproxTestAccumulatorCL::begin_accumulation ()
{
    surfarea = new double(0.0);
    stsurfarea = new double(0.0);
    posvol = new double(0.0);
    negvol = new double(0.0);
    intnorm = new Point3DCL(0.0);
    intnorm4d = new Point4DCL(0.0);
}


void STGeomApproxTestAccumulatorCL::finalize_accumulation ()
{
    std::cout << " current stsurfarea = " << *stsurfarea << std::endl;
    static double summed_stsurf = 0.0;
    summed_stsurf += *stsurfarea;
    std::cout << " summed stsurfarea = " << summed_stsurf << std::endl;
    std::cout << " surfarea = " << *surfarea/(tnew-told) << std::endl;
    std::cout << " negvol = " << *negvol/(tnew-told) << std::endl;
    std::cout << " posvol = " << *posvol/(tnew-told) << std::endl;
    std::cout << " intnorm = " << *intnorm << std::endl;
    std::cout << " intnorm4d = " << *intnorm4d << std::endl;
    delete surfarea;
    delete stsurfarea;
    delete negvol;
    delete posvol;
    delete intnorm;
    delete intnorm4d;
}


void STGeomApproxTestAccumulatorCL::visit (const TetraCL& sit)
{
    ScopeTimer scopetiming("STGeomApproxTestAccumulatorCL::visit");

    GetTrafoTr( T, det, sit);
    absdet= std::fabs(det);

    const double h = std::pow(absdet,1.0/3.0);
    const double dt = tnew-told;
    const double st_absdet= absdet * dt;
    
    LocalP2CL<> locPhi_old, locPhi_new;
    locPhi_old.assign( sit, lsetp2old->Phi, lsetp2old->GetBndData());
    locPhi_new.assign( sit, lsetp2new->Phi, lsetp2new->GetBndData());

    Point4DContainer pcontref;
    GeneralizedPrism4CL refprism4(pcontref);

    // get signs of vertices
    for (Uint i = 0; i < 4; ++i)
    {
        sign[i] = locPhi_old[i] > 0;
        sign[i+4] = locPhi_new[i] > 0;
    }
    
    bool inside_outer_band = false;
    bool outside_inner_band = false;
    
    for (int i = 0; i < 10; ++i)
    {
        if (locPhi_old[i] < h || locPhi_new[i] < h) 
            inside_outer_band = true;
        if (locPhi_old[i] > -h || locPhi_new[i] >- h) 
            outside_inner_band = true;
    }

    CompositeSTQuadCL<QuadRule> * p_cstquad = NULL;


    if (inside_outer_band && outside_inner_band)
    {
        if (lset_fpt == NULL)
            p_cstquad = new CompositeSTQuadCL<QuadRule>(refprism4,locPhi_old, locPhi_new,
                                                                      ints_per_space_edge, subtimeintervals);
        else
            p_cstquad = new CompositeSTQuadCL<QuadRule>(sit, TimeInterval(told,tnew), lset_fpt, 
                                                                      ints_per_space_edge, subtimeintervals);
        iscut = p_cstquad->HasInterface();
    }
    else
        iscut = false;
    
    if (!iscut)
    {
        ScopeTimer scopetiming("GeomApproxTestAccumulatorCL::visit - uncut");

        if (!sign[0])
#pragma omp critical(negvol)            
            *negvol += 1.0/6.0*st_absdet;
        else
#pragma omp critical(posvol)            
            *posvol += 1.0/6.0*st_absdet;
    }
    else
    {
        CompositeSTQuadCL<QuadRule> & cstquad = *p_cstquad;

        static LocalP2CL<Point3DCL> vecp0; // values for other time level are zero
        static LocalP2CL<> constone(1.0);

        LocalP2CL<> p0; // values for other time level are zero
    

        {
            ScopeTimer scopetiming("GeomApproxTestAccumulatorCL::visit_setup - cut");
            GridFunctionCL<double> gfone_neg = cstquad.EvalLinearOnPart( constone, constone, false); 
            GridFunctionCL<double> gfone_pos = cstquad.EvalLinearOnPart( constone, constone, true); 
#pragma omp critical(negvol)            
            *negvol += cstquad.QuadOnPart(gfone_neg,false)*st_absdet;
#pragma omp critical(posvol)            
            *posvol += cstquad.QuadOnPart(gfone_pos,true)*st_absdet;

        }
        // Space-time interface terms (Nitsche)
        {
            ScopeTimer scopetiming("GeomApproxTestAccumulatorCL::visit_setup - cut - stintface");
            //space time normals on reference prism pointing fr. neg to pos :
            const GridFunctionCL<Point4DCL> & n_st_ref (cstquad.GetNormalsOnInterface()); 
            //not correctly space time normals on physical prism pointing fr. neg to pos :
            GridFunctionCL<Point4DCL> n_st = eval(transform_normals_in_spacetime,T,absdet,dt,n_st_ref);
            //after next call space time normal is now correctly scaled (||n||=1) 
            //and the ratio between the measure on the physical prism to reference prism 
            //calculated. Additionally the space-measure-correction-factor nu is mulitplied.
            const GridFunctionCL<double> meas = eval(calc_norm,n_st);
            const GridFunctionCL<double> meas_nu = apply_and_eval(calc_nu_measure_and_normalize,n_st);
            //The corresponding space normal is extracted from the (phys.) space time normals
            const GridFunctionCL<Point3DCL> spacenormals = eval(calc_spacenormal_from_stnormal,n_st);
            const GridFunctionCL<Point3DCL> normnu = spacenormals * meas_nu;
            const GridFunctionCL<Point4DCL> meas_norm4d = n_st * meas;
#pragma omp critical(intnorm4d)            
            *intnorm4d += cstquad.QuadOnInterface(meas_norm4d);
#pragma omp critical(intnorm)            
            *intnorm += cstquad.QuadOnInterface(normnu);
#pragma omp critical(surfarea)            
            *surfarea += cstquad.QuadOnInterface(meas_nu);
#pragma omp critical(stsurfarea)            
            *stsurfarea += cstquad.QuadOnInterface(meas);
        }
    }
    if (inside_outer_band && outside_inner_band)
        delete p_cstquad;
}

//-----------------------------------------------------------------------------

InterfaceJumpAccumulatorCL::InterfaceJumpAccumulatorCL(const MultiGridCL& MG, 
                                                       const LevelsetP2CL * lsetp2old_in,
                                                       const LevelsetP2CL * lsetp2new_in, 
                                                       instat_scalar_fun_ptr lset_fpt_in, 
                                                       const double t1, const double t2,
                                                       const VecDescCL & old_sol_neg_in,
                                                       const VecDescCL & old_sol_pos_in,
                                                       const VecDescCL & new_sol_neg_in,
                                                       const VecDescCL & new_sol_pos_in,
                                                       const BndDataCL<> & Bnd_neg_in,
                                                       const BndDataCL<> & Bnd_pos_in,
                                                       const ParamCL::ptree_type & P):
    MG_(MG), lsetp2old(lsetp2old_in),lsetp2new(lsetp2new_in),
    lset_fpt(lset_fpt_in),
    ints_per_space_edge(P.get<int>("Quadrature.SubIntervalsPerEdge")),
    subtimeintervals(P.get<int>("Quadrature.SubTimeIntervals")),
    told(t1), tnew(t2),
    old_sol_neg(old_sol_neg_in),
    old_sol_pos(old_sol_pos_in),
    new_sol_neg(new_sol_neg_in),
    new_sol_pos(new_sol_pos_in),
    Bnd_neg(Bnd_neg_in),
    Bnd_pos(Bnd_pos_in),
    beta_neg(P.get<double>("HNeg")),beta_pos(P.get<double>("HPos")),
    lambda_stab(P.get<double>("NitschePenalty"))
{
    
}

void InterfaceJumpAccumulatorCL::begin_accumulation ()
{
    ifjump = new double(0.0);
    ifjump_past = new double(0.0);
    st_ifjump = new double(0.0);
    st_ifjump_nu = new double(0.0);
}


void InterfaceJumpAccumulatorCL::finalize_accumulation ()
{
    // std::cout << "                ifjump = " << *ifjump << std::endl;
    std::cout << "          sqrt(ifjump_past) = " << std::sqrt(*ifjump_past) << std::endl;
    std::cout << "               sqrt(ifjump) = " << std::sqrt(*ifjump) << std::endl;
    // std::cout << "             st_ifjump = " << *st_ifjump << std::endl;
    // std::cout << "       sqrt(st_ifjump) = " << std::sqrt(*st_ifjump) << std::endl;
    std::cout << "sqrt(st_ifjump_nu/(tnew-told)) = " << std::sqrt(*st_ifjump_nu/(tnew-told)) << std::endl;
    static double summed_st_ifjump = 0.0;
    summed_st_ifjump += *st_ifjump;
    std::cout << "           summed_st_ifjump = " << std::sqrt(summed_st_ifjump) << std::endl;
    static double summed_st_ifjump_nu = 0.0;
    summed_st_ifjump_nu += *st_ifjump_nu;
    std::cout << "           summed_st_ifjump_nu = " << std::sqrt(summed_st_ifjump_nu) << std::endl;
    delete st_ifjump;
    delete st_ifjump_nu;
    delete ifjump;
}


void InterfaceJumpAccumulatorCL::visit (const TetraCL& sit)
{
    ScopeTimer scopetiming("InterfaceJumpAccumulatorCL::visit");

    GetTrafoTr( T, det, sit);
    absdet= std::fabs(det);

    const double h = std::pow(absdet,1.0/3.0);
    const double dt = tnew-told;
    //const double st_absdet= absdet * dt;
    
    LocalP2CL<> locPhi_old, locPhi_new;
    locPhi_old.assign( sit, lsetp2old->Phi, lsetp2old->GetBndData());
    locPhi_new.assign( sit, lsetp2new->Phi, lsetp2new->GetBndData());

    Point4DContainer pcontref;
    GeneralizedPrism4CL refprism4(pcontref);

    // get signs of vertices
    for (Uint i = 0; i < 4; ++i)
    {
        sign[i] = locPhi_old[i] > 0;
        sign[i+4] = locPhi_new[i] > 0;
    }

    bool inside_outer_band = false;
    bool outside_inner_band = false;
    
    for (int i = 0; i < 10; ++i)
    {
        if (locPhi_old[i] < h || locPhi_new[i] < h) 
            inside_outer_band = true;
        if (locPhi_old[i] > -h || locPhi_new[i] >- h) 
            outside_inner_band = true;
    }

    if (!inside_outer_band || !outside_inner_band)
        return;

    CompositeSTQuadCL<QuadRule> * p_cstquad;
    if (lset_fpt == NULL)
        p_cstquad = new CompositeSTQuadCL<QuadRule>(refprism4,locPhi_old, locPhi_new,
                                                                  ints_per_space_edge, subtimeintervals);
    else
        p_cstquad = new CompositeSTQuadCL<QuadRule>(sit, TimeInterval(told,tnew), lset_fpt, 
                                                                  ints_per_space_edge, subtimeintervals);
    
    CompositeSTQuadCL<QuadRule> & cstquad = *p_cstquad;
    
    iscut = cstquad.HasInterface();

    static LocalP2CL<Point3DCL> vecp0; // values for other time level are zero
    static LocalP2CL<> constone(1.0);

    LocalP2CL<> p0; // values for other time level are zero
    

    LocalP1CL<> oldtrnegp1(sit,old_sol_neg,Bnd_neg);
    LocalP2CL<> oldtrneg(oldtrnegp1);
    oldtrneg *= beta_neg;
    LocalP1CL<> oldtrposp1(sit,old_sol_pos,Bnd_pos);
    LocalP2CL<> oldtrpos(oldtrposp1);
    oldtrpos *= beta_pos;

    LocalP2CL<> oldtrdiff(oldtrpos);
    oldtrdiff -= oldtrneg;

    LocalP1CL<> newtrnegp1(sit,new_sol_neg,Bnd_neg);
    LocalP2CL<> newtrneg(newtrnegp1);
    newtrneg *= beta_neg;
    LocalP1CL<> newtrposp1(sit,new_sol_pos,Bnd_pos);
    LocalP2CL<> newtrpos(newtrposp1);
    newtrpos *= beta_pos;

    LocalP2CL<> newtrdiff(newtrpos);
    newtrdiff -= newtrneg;


    if (iscut)
    {
        for(Uint ti = PAST; ti<= FUTURE; ti++)
        {
            ScopeTimer scopetiming("InterfaceJumpAccumulatorCL::visit_setup - cut - fut. intface");
            QuadDomainCL q2dom_;
            PrincipalLatticeCL lat( PrincipalLatticeCL::instance( ints_per_space_edge ));
            TetraPartitionCL partition_;
            std::valarray<double> ls_loc_( ((ints_per_space_edge+1) * (ints_per_space_edge+2) * (ints_per_space_edge+3)) / 6);
            if (ti == PAST)
                evaluate_on_vertexes( lsetp2old->GetSolution(), sit, lat, Addr( ls_loc_));
            else
                evaluate_on_vertexes( lsetp2new->GetSolution(), sit, lat, Addr( ls_loc_));

            SurfacePatchCL spatch;
            spatch.make_patch<MergeCutPolicyCL>( lat, ls_loc_);

            QuadDomain2DCL  q2Ddomain;
            make_CompositeQuad5Domain2D ( q2Ddomain, spatch, sit);

            GridFunctionCL<> diff;
            if (ti == PAST)
                resize_and_evaluate_on_vertexes( oldtrdiff, q2Ddomain, diff);
            else
                resize_and_evaluate_on_vertexes( newtrdiff, q2Ddomain, diff);

            for (Uint k = 0; k < diff.size(); k++)
                diff[k] = diff[k]*diff[k];

            const double add = quad_2D( diff ,q2Ddomain);

            if (ti == FUTURE)
#pragma omp critical(ifjump)            
                *ifjump += add;
            else
#pragma omp critical(ifjump_past)            
                *ifjump_past += add;

        }
        // Space-time interface terms
        {
            ScopeTimer scopetiming("InterfaceJumpAccumulatorCL::visit_setup - cut - stintface");
            //space time normals on reference prism pointing fr. neg to pos :
            const GridFunctionCL<Point4DCL> & n_st_ref (cstquad.GetNormalsOnInterface()); 
            //not correctly space time normals on physical prism pointing fr. neg to pos :
            GridFunctionCL<Point4DCL> n_st = eval(transform_normals_in_spacetime,T,absdet,dt,n_st_ref);
            //after next call space time normal is now correctly scaled (||n||=1) 
            //and the ratio between the measure on the physical prism to reference prism 
            //calculated. Additionally the space-measure-correction-factor nu is mulitplied.

            const GridFunctionCL<double> meas = eval(calc_norm,n_st);
            const GridFunctionCL<double> meas_nu = apply_and_eval(calc_nu_measure_and_normalize,n_st);

            GridFunctionCL<double> gfdiff = cstquad.EvalLinearOnInterface( oldtrdiff, newtrdiff);
            GridFunctionCL<double> gfdiff2(gfdiff.size());
            for (Uint k = 0; k < gfdiff2.size(); k++)
                gfdiff2[k] = gfdiff[k]*gfdiff[k] * meas[k];

            GridFunctionCL<double> gfdiff2_nu(gfdiff.size());
            for (Uint k = 0; k < gfdiff2.size(); k++)
                gfdiff2_nu[k] = gfdiff[k]*gfdiff[k] * meas_nu[k];
            
            const double add = cstquad.QuadOnInterface(gfdiff2);
#pragma omp critical(st_ifjump)            
            *st_ifjump += add;

            const double add2 = cstquad.QuadOnInterface(gfdiff2_nu);
#pragma omp critical(st_ifjump_nu)            
            *st_ifjump_nu += add2;
            
        }
    }

    delete p_cstquad;
}



}//end of namespace DROPS
