/// \file spacetime_surf.cpp
/// \brief space time setup test st interface
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

#include "misc/problem.h"
#include "misc/params.h"
#include "misc/progressaccu.h"
#include "num/spacetime_geom.h"
#include "num/spacetime_quad.h"
// #include "spacetimetransp/stxfem.h"
// #include "num/interfacePatch.h"
#include "num/discretize.h"
#include "num/accumulator.h"
// #include "num/quadrature.h"
#include "levelset/levelset.h"

namespace DROPS
{

// Accumulator for space time integrals
class STSurfaceTestAccumulatorCL : public TetraAccumulatorCL
{
    protected:
    const MultiGridCL& MG_;

    const LevelsetP2CL * lsetp2old;
    const LevelsetP2CL * lsetp2new;
    instat_scalar_fun_ptr lset_fpt;

    // - sharable (for future changes)
    Point3DCL G[4];
    double det;
    double absdet;
    SMatrixCL<3,3> T;
    IdxT UnknownIdx[4];  // the std. unknowns on one time level (the other follows by +1 (stride.. ))

    bool sign[8]; //sign of levelset at vertices (4(space) x 2(time))

//    const Uint lvl;


    const double told;
    const double tnew;

    double iscut;

    const Uint ints_per_space_edge;
    const Uint subtimeintervals;

    double * surfarea;

    double * minlst;
    double * maxlst;

public:
    STSurfaceTestAccumulatorCL (const MultiGridCL& MG,
                                instat_scalar_fun_ptr lset_fpt_in,
                                const LevelsetP2CL * lsetp2old_in,
                                const LevelsetP2CL * lsetp2new_in,
                                const double t1, const double t2,
                                const int subdivspace_in,
                                const int subdivtim_in,
                                double* surfarea_in,
                                double* minlst_in,
                                double* maxlst_in);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    virtual void local_setup (const TetraCL& sit);

    virtual void visit (const TetraCL&);

    virtual TetraAccumulatorCL* clone (int /*tid*/) { return new STSurfaceTestAccumulatorCL ( *this); }

};


STSurfaceTestAccumulatorCL::STSurfaceTestAccumulatorCL (const MultiGridCL& MG,
                                                        instat_scalar_fun_ptr lset_fpt_in,
                                                        const LevelsetP2CL * lsetp2old_in,
                                                        const LevelsetP2CL * lsetp2new_in,
                                                        const double t1, const double t2,
                                                        const int subdivspace_in,
                                                        const int subdivtim_in,
                                                        double* surfarea_in,
                                                        double* minlst_in,
                                                        double* maxlst_in):
    MG_(MG), lsetp2old(lsetp2old_in),lsetp2new(lsetp2new_in), lset_fpt(lset_fpt_in),
    told(t1), tnew(t2),
    ints_per_space_edge(subdivspace_in), subtimeintervals(subdivtim_in), surfarea(surfarea_in),
    minlst(minlst_in), maxlst(maxlst_in)
{
}

void STSurfaceTestAccumulatorCL::begin_accumulation ()
{
    *surfarea = 0.0;
}

void STSurfaceTestAccumulatorCL::finalize_accumulation ()
{
    // std::cout << " surfarea = " << *surfarea << std::endl;
    // std::cout << " | surfarea - PI | = " << std::fabs(*surfarea - M_PI) << std::endl;
}

void STSurfaceTestAccumulatorCL::visit (const TetraCL& sit)
{
  local_setup(sit);
}

/* --------------- STVolume Transport term (time deriv., convection, diffusion)  -------------- */
void STSurfaceTestAccumulatorCL::local_setup (const TetraCL& sit)
{
    ScopeTimer scopetiming("STSurfaceTestAccumulatorCL::local_setup - total");

    GetTrafoTr( T, det, sit);
    P1DiscCL::GetGradients(G,det,sit);

    absdet= std::fabs(det);

    // const double h = std::pow(absdet,1.0/3.0);
    const double dt = tnew-told;
    // const double st_absdet= absdet * dt;


    LocalP2CL<> locPhi_old, locPhi_new;
    if (lset_fpt == NULL)
    {
        locPhi_old.assign( sit, lsetp2old->Phi, lsetp2old->GetBndData());
        locPhi_new.assign( sit, lsetp2new->Phi, lsetp2new->GetBndData());
    }
    Point4DContainer pcont;
    GeneralizedPrism4CL refprism4(pcont);

    // get signs of vertices
    for (Uint i = 0; i < 4; ++i)
    {
        sign[i] = locPhi_old[i] > 0;
        sign[i+4] = locPhi_new[i] > 0;
    }

    TPSpaceTimeMapping stmap(sit, TimeInterval(told,tnew));

    CompositeSTQuadCL<QuadRule> * p_cstquad;
    if (lset_fpt == NULL)
        p_cstquad = new CompositeSTQuadCL<QuadRule>(refprism4,locPhi_old, locPhi_new,
                                                                  ints_per_space_edge, subtimeintervals);
    else
        p_cstquad = new CompositeSTQuadCL<QuadRule>(sit, TimeInterval(told,tnew), lset_fpt,
                                                                  ints_per_space_edge, subtimeintervals);

    CompositeSTQuadCL<QuadRule> & cstquad = *p_cstquad;

    iscut = cstquad.HasInterface();

    LocalP2CL<> p0; // values for other time level are zero

    if (iscut)
    {
        ScopeTimer scopetiming("STSurfaceTestAccumulatorCL::local_setup - cut - stintface");

        const GridFunctionCL<Point4DCL> & n_st_ref (cstquad.GetNormalsOnInterface());
        GridFunctionCL<Point4DCL> n_st = eval(transform_normals_in_spacetime,T,absdet,dt,n_st_ref);
        const GridFunctionCL<double> meas_nu = apply_and_eval(calc_nu_measure_and_normalize,n_st);
        // const GridFunctionCL<double> meas_nu = apply_and_eval(calc_norm_and_normalize,n_st);
        const GridFunctionCL<Point3DCL> spacenormals = eval(calc_spacenormal_from_stnormal,n_st);

#pragma omp critical(surfarea)
        *surfarea += cstquad.QuadOnInterface(meas_nu);


        const GridFunctionCL<Point4DCL> if_intpts = cstquad.GetInterfaceIntegrationPoints();
        GridFunctionCL<Point3DCL> if_tr_pts(if_intpts.size());
        GridFunctionCL<Point4DCL> if_tr_pts2(if_intpts.size());

        Point3DCL v[4];
        for (Uint i = 0; i < 4; ++i)
            v[i] = sit.GetVertex(i)->GetCoord();
        for (Uint i = 0; i < if_intpts.size(); ++i)
        {
            if_tr_pts[i] = MakePoint3D(0.0,0.0,0.0);
            BaryCoordCL b(MakeBaryCoord(RestrictToPoint3D(if_intpts[i])));
            for (Uint j = 0; j < 4; ++j)
                if_tr_pts[i] += b[j] * v[j];

            if_tr_pts2[i] = Point4DCL(if_tr_pts[i], if_intpts[i][3]);
        }

        // static int first = true;
        // if (first)
        // std::cout << " if_tr_pts2 = " << if_tr_pts2 << std::endl;
        // first = false;

        const GridFunctionCL<double> lsetvals(lset_fpt ? cstquad.EvalOnInterface(lset_fpt, &stmap) : cstquad.EvalLinearOnInterface(locPhi_old, locPhi_new));
        // std::cout << " lsetvals = " << lsetvals << std::endl;
        for (Uint i = 0; i < lsetvals.size(); ++i)
        {
            *maxlst = std::max(*maxlst,lsetvals[i]);
            *minlst = std::min(*minlst,lsetvals[i]);
        }
        // std::cout << " lsetvals = " << minlst << " \t " << maxlst << std::endl;


    }
    delete p_cstquad;
}



double lset_test (const Point3DCL& p, double t)
{

    // Point3DCL d = p - MakePoint3D(-0.25+t*0.5,0,0);
    // return d.norm() - 0.5;

    Point3DCL d = p - MakePoint3D(-0.25+t*0.5,0,0);
    return std::sqrt(d[0]*d[0]+d[1]*d[1]) - 0.5;

    // return p[0] - 0.5*t + 0.25;
}

double Run (int numref, int numsteps, int subdivspace, int subdivtime)
{
  double ret = 0.0;
  double tmp = 0.0;
  double minlst = 1e99;
  double maxlst = -1e99;
  try
  {
      BrickBuilderCL brick( MakePoint3D( -1.0, 0.0, 0.0),
                                   2.0*std_basis<3>(1),
                                   1.0*std_basis<3>(2),
                                   1.0*std_basis<3>(3),
                                   2*numref, numref, numref);
      MultiGridCL mg( brick);

      instat_scalar_fun_ptr sigma (0);
      SurfaceTensionCL sf( sigma, 0);
      BndCondT bc[6]= { NoBC, NoBC, NoBC, NoBC, NoBC, NoBC };
      LsetBndDataCL::bnd_val_fun bfun[6]= { 0,0,0,0,0,0};
      LsetBndDataCL lsbnd( 6, bc, bfun);

      ret = 0.0;
      for (int i = 0; i < numsteps; ++i)
      {
          const double told = (double)i/numsteps;
          const double tnew = (double)(i+1)/numsteps;

          LevelsetP2CL & lset_old( * LevelsetP2CL::Create( mg, lsbnd, sf));
          lset_old.CreateNumbering( mg.GetLastLevel());
          lset_old.Init( lset_test, told);

          LevelsetP2CL & lset_new( * LevelsetP2CL::Create( mg, lsbnd, sf));
          lset_new.CreateNumbering( mg.GetLastLevel());
          lset_new.Init( lset_test, tnew);

          TetraAccumulatorTupleCL accus ;
          // ProgressBarTetraAccumulatorCL accup(MG,"STTranspAcc");
          // accus.push_back(&accup);

          // STSurfaceTestAccumulatorCL accu (mg, NULL, &lset_old, &lset_new, told, tnew, subdivspace, subdivtime, &tmp);

          STSurfaceTestAccumulatorCL accu (mg, lset_test, NULL, NULL, told, tnew, subdivspace, subdivtime, &tmp, &minlst, &maxlst);
          MaybeAddProgressBar(mg, "Surface integration", accus, -1);
          accus.push_back(&accu);
          // std::cout << " accumulate " << std::endl;
          accumulate( accus, mg, lset_old.idx.TriangLevel(), lset_old.idx.GetBndInfo());
          // getchar();
          ret += tmp;

          delete &lset_new;
          delete &lset_old;
      }
      // ret = std::fabs(ret-std::sqrt(5.0/4.0)); //*1.05983938019);

      // ret = std::fabs(ret-0.25*M_PI); //*1.05983938019);
      ret = std::fabs(ret-0.5*M_PI);

      // ret = std::fabs(ret-0.03125*M_PI*M_PI); //*1.05983938019);
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
  std::cout << " minlst = " << minlst << std::endl;
  std::cout << " maxlst = " << maxlst << std::endl;
  return ret;
}



}//end of namespace DROPS

#define SS 3
#define ST 3

const int sr0 = 4;
const int tr0 = 2;

void output_results(DROPS::SMatrixCL<SS,ST> res)
{

    std::cout << "results:" << std::endl;

    std::cout << "|-----|";
    for (int q = 0, tr = tr0; q < ST; ++q, tr*=2)
        std::cout << "--------------|";
    std::cout << std::endl;

    std::cout << "|     |";
    for (int q = 0, tr = tr0; q < ST; ++q, tr*=2)
        std::cout << " " << std::setw(12) << tr << " |";
    std::cout << std::endl;

    std::cout << "|-----|";
    for (int q = 0, tr = tr0; q < ST; ++q, tr*=2)
        std::cout << "--------------|";
    std::cout << std::endl;

    for (int p = 0, sr = sr0; p < SS; ++p, sr*=2)
    {
        std::cout << "| " << std::setw(3) << sr << " | ";
        for (int q = 0, tr = tr0; q < ST; ++q, tr*=2)
        {
            std::cout << std::setw(12) << res(p,q) << " | ";
        }
        std::cout << std::endl;
    }

    std::cout << "|_____|";
    for (int q = 0, tr = tr0; q < ST; ++q, tr*=2)
        std::cout << "______________|";
    std::cout << std::endl;

}


int main ()
{
    DROPS::ProgressBarTetraAccumulatorCL::Activate();
    DROPS::SMatrixCL<SS,ST> res_global;

    for (int p = 0, sr = sr0; p < SS; ++p, sr*=2)
        for (int q = 0, tr = tr0; q < ST; ++q, tr*=2)
        {
            std::cout << "\n\n run with ( sr = " << sr << ", tr = " << tr << ")" << std::endl;
            res_global(p,q) = DROPS::Run(  sr, tr, 1, 1);
            std::cout << "result ( sr = " << sr << ", tr = " << tr << ") : " << res_global(p,q) << std::endl;
            //output_results(res_global);
        }

    output_results(res_global);


    DROPS::SMatrixCL<SS,ST> res_local;

    for (int p = 0, sr = sr0; p < SS; ++p, sr*=2)
        for (int q = 0, tr = tr0; q < ST; ++q, tr*=2)
        {
            std::cout << "\n\n run with ( sr = " << sr << ", tr = " << tr << ")" << std::endl;
            res_local(p,q) = DROPS::Run(  1, 1, sr, tr);
            std::cout << "result ( sr = " << sr << ", tr = " << tr << ") : " << res_local(p,q) << std::endl;
            //output_results(res_local);
        }

    output_results(res_local);

    return 0;
}
