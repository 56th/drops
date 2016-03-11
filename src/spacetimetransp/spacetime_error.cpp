/// \file spacetime_error.cpp
/// \brief setup routines for error calculation
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

#include "spacetimetransp/stxfem.h"
#include "spacetimetransp/spacetime_error.h"
#include "misc/funcmap.h"

// #include "num/interfacePatch.h"


namespace DROPS
{

using namespace STXFEM;

EnergyNormErrorAccumulatorCL::EnergyNormErrorAccumulatorCL(const MultiGridCL& MG, 
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
    sol_neg(InScaMap::getInstance()[P.get<std::string>("SolNeg","Zero")]),
    sol_pos(InScaMap::getInstance()[P.get<std::string>("SolPos","Zero")]),
    sol_dt_neg(InScaMap::getInstance()[P.get<std::string>("SolDtNeg","Zero")]),
    sol_dt_pos(InScaMap::getInstance()[P.get<std::string>("SolDtPos","Zero")]),
    sol_grad_neg(InVecMap::getInstance()[P.get<std::string>("SolGradNeg","ZeroVel")]),
    sol_grad_pos(InVecMap::getInstance()[P.get<std::string>("SolGradPos","ZeroVel")]),
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
    alpha_neg(P.get<double>("DiffNeg")),alpha_pos(P.get<double>("DiffPos")),
    lambda_stab(P.get<double>("NitschePenalty")),
    vmax(P.get<double>("MaxVelocity"))
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

void EnergyNormErrorAccumulatorCL::begin_accumulation ()
{
    h10norm2 = new double(0.0);
    h01norm2 = new double(0.0);
    l2norm2 = new double(0.0);
}


void EnergyNormErrorAccumulatorCL::finalize_accumulation ()
{
    std::cout << "          sqrt(h10norm2) = " << std::sqrt(*h10norm2) << std::endl;
    std::cout << "          sqrt(h01norm2) = " << std::sqrt(*h01norm2) << std::endl;
    std::cout << "          sqrt(l2norm2)  = " << std::sqrt(*l2norm2) << std::endl;

    static double summed_l2norm2 = 0.0;
    summed_l2norm2 += *l2norm2;
    std::cout << "           summed_l2norm2  = " << std::sqrt(summed_l2norm2) << std::endl;
    static double summed_h10norm2 = 0.0;
    summed_h10norm2 += *h10norm2;
    std::cout << "           summed_h10norm2 = " << std::sqrt(summed_h10norm2) << std::endl;
    static double summed_h01norm2 = 0.0;
    summed_h01norm2 += *h01norm2;
    std::cout << "           summed_h01norm2 = " << std::sqrt(summed_h01norm2) << std::endl;

    delete l2norm2;
    delete h10norm2;
    delete h01norm2;

}


void EnergyNormErrorAccumulatorCL::visit (const TetraCL& sit)
{
    ScopeTimer scopetiming("EnergyNormErrorAccumulatorCL::visit");

    GetTrafoTr( T, det, sit);
    P1DiscCL::GetGradients(G,det,sit);

    absdet= std::fabs(det);

    // const double h = std::pow(absdet,1.0/3.0); //unused
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

    // double add=0.0; //unused

    if (iscut)
    {
        ScopeTimer scopetiming("EnergyNormErrorAccumulator::local_setup - cut");
        TPSpaceTimeMapping stmap(sit, TimeInterval (told,tnew));
        LocalP2CL<Point3DCL> vecp0(MakePoint3D(0.0,0.0,0.0)); // values for other time level are zero
        LocalP2CL<> constone(1.0);
        LocalP2CL<> p0(0.0); // values for other time level are zero

        CompositeSTQuadCL<QuadRule> & cstquad = *p_cstquad;

        GridFunctionCL<double> solpos = cstquad.EvalOnPart( sol_pos, true, &stmap); 
        GridFunctionCL<double> solneg = cstquad.EvalOnPart( sol_neg, false, &stmap);
        GridFunctionCL<double> soldtpos = cstquad.EvalOnPart( sol_dt_pos, true, &stmap); 
        GridFunctionCL<double> soldtneg = cstquad.EvalOnPart( sol_dt_neg, false, &stmap);
        GridFunctionCL<Point3DCL> solgradpos = cstquad.EvalOnPart( sol_grad_pos, true, &stmap); 
        GridFunctionCL<Point3DCL> solgradneg = cstquad.EvalOnPart( sol_grad_neg, false, &stmap);
 
        GridFunctionCL<double> gfone_neg = cstquad.EvalLinearOnPart( constone, constone, false); 
        GridFunctionCL<double> gfone_pos = cstquad.EvalLinearOnPart( constone, constone, true); 
        
        ScopeTimer scopetiming2("STTranspVolAccu_P1SP1TXCL::local_setup - cut - 4Dvols");
        for(Uint s = 0; s < 2; s++) // both cases: sign = false and sign = true 
        {
            const bool csign = s==1;

            LocalP1CL<> oldtrp1(sit,csign?old_sol_pos:old_sol_neg,csign?Bnd_pos:Bnd_neg);
            LocalP1CL<> newtrp1(sit,csign?new_sol_pos:new_sol_neg,csign?Bnd_pos:Bnd_neg);

            LocalP1CL<double> oldsolp1 (sit, csign?sol_pos:sol_neg, told);
            LocalP1CL<double> newsolp1 (sit, csign?sol_pos:sol_neg, tnew);

            LocalP2CL<> oldtr(oldtrp1);
            LocalP2CL<> newtr(newtrp1);

            LocalP2CL<> oldsoltr(oldsolp1);
            LocalP2CL<> newsoltr(newsolp1);
            
            // double (* coup_s_s_cur)[8][8] = csign ? &coup_s_s_pos : &coup_s_s_neg;
            // double (* elvec_cur)[8] = csign ? &elvec_pos : &elvec_neg;
            // const GridFunctionCL<Point3DCL>& vel_cur (csign ? velpos : velneg);

            /* unused
            const double & alpha (csign ? alpha_pos : alpha_neg);
            */
            const double & beta (csign ? beta_pos : beta_neg);

            /* unused
            GridFunctionCL<double> psi_PAST = cstquad.EvalLinearOnPart( constone, p0, csign);
            GridFunctionCL<double> psi_FUTURE = cstquad.EvalLinearOnPart( p0, constone, csign); 
            GridFunctionCL<double> & gfone (csign? gfone_pos: gfone_neg);
            */

            GridFunctionCL<double> uh = cstquad.EvalLinearOnPart( oldtr, newtr, csign);
            LocalP2CL<> dtr = newtr;
            newtr -= oldtr;
            dtr *= 1.0/dt;

            LocalP2CL<> soldtr = newsoltr;
            newsoltr -= oldsoltr;
            soldtr *= 1.0/dt;

            GridFunctionCL<double> duhdt = cstquad.EvalLinearOnPart( dtr, dtr, csign);
            GridFunctionCL<double> dudt = cstquad.EvalLinearOnPart( soldtr, soldtr, csign);

            Point3DCL gradold = MakePoint3D(0,0,0);
            Point3DCL gradnew = MakePoint3D(0,0,0);
            for (Uint i = 0; i < 4; ++i)
                gradold += G[i]*oldtrp1[i];
            for (Uint i = 0; i < 4; ++i)
                gradnew += G[i]*newtrp1[i];

            LocalP2CL<Point3DCL> gradoldp2(gradold);
            LocalP2CL<Point3DCL> gradnewp2(gradnew);
            GridFunctionCL<Point3DCL> graduh = cstquad.EvalLinearOnPart( gradoldp2, gradnewp2, csign);

            GridFunctionCL<double> &u (csign?solpos:solneg);

            // GridFunctionCL<double> &dudt (csign?soldtpos:soldtneg);
            GridFunctionCL<Point3DCL> &gradu (csign?solgradpos:solgradneg);

            GridFunctionCL<double> err = uh;
            err -= u;
            GridFunctionCL<double> err2 = eval(square, err);


            GridFunctionCL<double> errdt = duhdt;
            errdt -= dudt;
            GridFunctionCL<double> err2dt = eval(square, errdt);

            GridFunctionCL<Point3DCL> errgrad = graduh;
            errgrad -= gradu;
            GridFunctionCL<double> err2grad = eval(calc_norm_sqr, errgrad);

            // std::cout << " u  = \n" << u  << std::endl;
            // std::cout << " u2  = \n" << u2  << std::endl;
            // std::cout << " uh = \n" << uh << std::endl;
            // std::cout << " gradu  = \n" << gradu  << std::endl;
            // std::cout << " graduh = \n" << graduh << std::endl;
            // std::cout << " errgrad  = \n" << errgrad  << std::endl;
            // std::cout << " err2grad  = \n" << err2grad  << std::endl;


            double tmp = beta * cstquad.QuadOnPart( err2, csign) * st_absdet;
#pragma omp critical(l2norm2)
            *l2norm2 += tmp;

            tmp = beta * cstquad.QuadOnPart( err2grad, csign) * st_absdet;
#pragma omp critical(l2norm2)
            *h10norm2 += tmp;

            tmp = beta * cstquad.QuadOnPart( err2dt, csign) * st_absdet;
#pragma omp critical(l2norm2)
            *h01norm2 += tmp;

        }
    }
    else
    {
        // std::cout << " is not cut " << std::endl; 

        ScopeTimer scopetiming("EnergyNormErrorAccumulator::local_setup - uncut");

        const double & alpha (sign[0] ? alpha_pos : alpha_neg);
        const double & beta (sign[0] ? beta_pos : beta_neg);

        // only spatial integrals 
        double coup_s_s_space_lap[4][4];
        double coup_s_s_space_mass[4][4];
        for(Uint i=0; i<4; ++i)
        {
            for(Uint j=0; j<i; ++j){
                coup_s_s_space_mass[i][j] = beta * P1DiscCL::GetMass(i,j);
                coup_s_s_space_mass[j][i] = coup_s_s_space_mass[i][j];
                coup_s_s_space_lap[i][j] = beta * alpha*inner_prod( G[i], G[j])/6.0; //mass
                coup_s_s_space_lap[j][i] = coup_s_s_space_lap[i][j];
            }
            coup_s_s_space_mass[i][i] = beta * P1DiscCL::GetMass(i,i); //mass
            coup_s_s_space_lap[i][i] = beta*alpha*inner_prod( G[i], G[i])/6.0; //diffusion
        }

        // iterated integrals
        SMatrixCL<8,8> coup_s_s_lap;
        SMatrixCL<8,8> coup_s_s_mass;
        SMatrixCL<8,8> coup_s_s_dtdt;
        for(Uint ti = PAST; ti<= FUTURE; ti++)
            for(Uint i=0; i<4; ++i)
                for(Uint tj = PAST; tj<=FUTURE; tj++)
                    for(Uint j=0; j<4; ++j)
                    {
                        coup_s_s_lap(ti*4+i,tj*4+j) += coup_s_s_space_lap[i][j] * mass_1D_p1p1[ti][tj]  * st_absdet;
                        coup_s_s_mass(ti*4+i,tj*4+j) += coup_s_s_space_mass[i][j] * mass_1D_p1p1[ti][tj]  * st_absdet;
                        coup_s_s_dtdt(ti*4+i,tj*4+j) += coup_s_s_space_mass[i][j] * dudtdvdt_1D_p1p1[ti][tj]  * absdet / dt;
                    }

        SVectorCL<8> err_elvec_st;
        LocalP1CL<> oldtr(sit,sign[0]?old_sol_pos:old_sol_neg,sign[0]?Bnd_pos:Bnd_neg);
        LocalP1CL<> newtr(sit,sign[0]?new_sol_pos:new_sol_neg,sign[0]?Bnd_pos:Bnd_neg);

        LocalP1CL<double> oldsolp1 (sit, sign[0]?sol_pos:sol_neg, told);
        LocalP1CL<double> newsolp1 (sit, sign[0]?sol_pos:sol_neg, tnew);


        // std::cout << "oldtr =\n " << oldtr << std::endl;
        // std::cout << "newtr =\n " << newtr << std::endl;
        // std::cout << "oldsolp1 =\n " << oldsolp1 << std::endl;
        // std::cout << "newsolp1 =\n " << newsolp1 << std::endl;
        // getchar();

        for(Uint i=0; i<4; ++i)
            err_elvec_st[PAST*4+i] = oldtr[i]-oldsolp1[i];
        for(Uint i=0; i<4; ++i)
            err_elvec_st[FUTURE*4+i] = newtr[i]-newsolp1[i];

        SVectorCL<8> tmp = coup_s_s_lap * err_elvec_st;
#pragma omp critical(h10norm2)
        *h10norm2 += inner_prod(tmp, err_elvec_st);

        tmp = coup_s_s_mass * err_elvec_st;
#pragma omp critical(l2norm2)
        *l2norm2 += inner_prod(tmp, err_elvec_st);

        tmp = coup_s_s_dtdt * err_elvec_st;
#pragma omp critical(h01norm2)
        *h01norm2 += inner_prod(tmp, err_elvec_st);

    }

                
    delete p_cstquad;
}



}//end of namespace DROPS
