#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "num/discretize.h"
#include "num/fe.h"
#include "misc/problem.h"
#include "levelset/levelset.h"
#include "stokes/instatstokes2phase.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include <fstream>
#include <iomanip>

using namespace DROPS;

inline int iSign( double x)
{
    return x<0 ? -1 : x>0 ? 1 : 0;
}

// \Omega_1 is the domain with phasebnd < 0.
double phasebnd (const Point3DCL& p)
{
    return p[1] + p[2] - 0.05;
//    return p[1] + p[2] - 0.55;
//    return p[0] + p[1] - 0.5;
}

double f1 (const Point3DCL& p, double)
{
    return p.norm_sq();
}

double f2 (const Point3DCL& p, double)
{
    return 3.0*p[0]*p[0] + p[1]*p[1] + 2.0*p[2]*p[2] + 2.0;
}

// double f1 (const Point3DCL& p, double)
// {
//     return p[0] + 0.5*p[1]  - 1.0;
// }
// 
// double f2 (const Point3DCL& p, double)
// {
//     return 0.7*p[0] + p[1] + 0.3*p[2];
// }

double f1sq (const Point3DCL& p, double)
{
    return std::pow( f1( p, 0.), 2);
}

double f2sq (const Point3DCL& p, double)
{
    return std::pow( f2( p, 0.), 2);
}

void InitPiecewiseP2 (instat_scalar_fun_ptr fneg, VecDescCL& un,
    instat_scalar_fun_ptr fpos, VecDescCL& up, MultiGridCL& MG, VecDescCL& Phi)
{
    const Uint idx= Phi.RowIdx->GetIdx(),
               unidx= un.RowIdx->GetIdx(),
               upidx= up.RowIdx->GetIdx();


    DROPS_FOR_TRIANG_VERTEX( MG, MG.GetLastLevel(), it) {
        switch (iSign( Phi.Data[it->Unknowns(idx)])) {
          case 1:
            un.Data[it->Unknowns( unidx)]= 0.;
            up.Data[it->Unknowns( upidx)]= fpos( it->GetCoord(), 0.);
            break;
          case -1:
            un.Data[it->Unknowns( unidx)]= fneg( it->GetCoord(), 0.);
            up.Data[it->Unknowns( upidx)]= 0.;
            break;
          default:
            un.Data[it->Unknowns( unidx)]= fneg( it->GetCoord(), 0.);
            up.Data[it->Unknowns( upidx)]= fpos( it->GetCoord(), 0.);
        }
    }
    DROPS_FOR_TRIANG_EDGE( MG, MG.GetLastLevel(), it) {
        switch (iSign( Phi.Data[it->Unknowns(idx)])) {
          case 1:
            un.Data[it->Unknowns( unidx)]= 0.;
            up.Data[it->Unknowns( upidx)]= fpos( GetBaryCenter( *it), 0.);
            break;
          case -1:
            un.Data[it->Unknowns( unidx)]= fneg( GetBaryCenter( *it), 0.);
            up.Data[it->Unknowns( upidx)]= 0.;
            break;
          default:
            un.Data[it->Unknowns( unidx)]= fneg( GetBaryCenter( *it), 0.);
            up.Data[it->Unknowns( upidx)]= fpos( GetBaryCenter( *it), 0.);
        }
    }
    InterfacePatchCL cut;
    DROPS_FOR_TRIANG_TETRA( MG, MG.GetLastLevel(), it) {
        cut.Init( *it, Phi);
        if (cut.Intersects()) {
            for (Uint i= 0; i < 4; ++i) {
                const VertexCL& v= *it->GetVertex( i);
                un.Data[v.Unknowns( unidx)]= fneg( v.GetCoord(), 0.);
                up.Data[v.Unknowns( upidx)]= fpos( v.GetCoord(), 0.);
            }
            for (Uint i= 0; i < 6; ++i) {
                const EdgeCL& e= *it->GetEdge( i);
                un.Data[e.Unknowns( unidx)]= fneg( GetBaryCenter( e), 0.);
                up.Data[e.Unknowns( upidx)]= fpos( GetBaryCenter( e), 0.);
            }
        }
    }
}

void P1XOnPart (const VecDescCL& p1x, const ExtIdxDescCL& Xidx, VecDescCL& p_part,
    const LevelsetP2CL& lset, bool posPart)
{
    const Uint lvl= p1x.RowIdx->TriangLevel,
        idxnum= p1x.RowIdx->GetIdx();
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const MultiGridCL& mg= lset.GetMG();

    p_part.SetIdx( p1x.RowIdx);
    VectorCL& pp= p_part.Data;
    pp= p1x.Data; // Our Ensight-writer ignores the extra p1x values at the end.

    // add extended pressure
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        const IdxT nr= it->Unknowns( idxnum);
        if (Xidx[nr] == NoIdx) continue;

        const bool is_pos= InterfacePatchCL::Sign( ls.val( *it)) == 1;
        if (posPart == is_pos) continue; // extended hat function == 0 on this part
        if (posPart)
            pp[nr]+= p1x.Data[Xidx[nr]];
        else
            pp[nr]-= p1x.Data[Xidx[nr]];
    }
}

class CoeffCL
{
  public:
    double mu (double) const { return 1.0; }
};

int main (int argc, char** argv)
{
  try {
    int numref=10;
    double xfemstab=0.;
    if (argc==3) {
        numref=atoi(argv[1]);
        xfemstab=atof(argv[2]);
    }

    std::cerr << "numref: " << numref << "\txfemstab: " << xfemstab <<'\n';

    DROPS::BrickBuilderCL brick( Point3DCL( -1.0),
                                  2.0*std_basis<3>(1),
                                  2.0*std_basis<3>(2),
                                  2.0*std_basis<3>(3),
                                  numref, numref, numref);
    MultiGridCL mg( brick);

//    TetraBuilderCL builder( 0);//,   std_basis<3>( 0), 2*std_basis<3>( 1),
                               //2*std_basis<3>( 2), 2*std_basis<3>( 3));
 //   MultiGridCL mg( builder);


    LevelsetP2CL lset( mg);
    CreateNumb( lset.idx.GetIdx(), lset.idx, mg, NoBndDataCL<>());
    lset.Phi.SetIdx( &lset.idx);
    lset.Init( &phasebnd);

    IdxDescCL idx( P1X_FE);
    CreateNumb( /*Level*/ mg.GetLastLevel(), idx, mg, NoBndDataCL<>());
    ExtIdxDescCL extidx( &idx, /*omit_bound=*/ xfemstab);
    extidx.UpdateXNumbering( &idx, lset, true);
    std::cerr << "P1-Unbekannte: " << extidx.GetNumUnknownsP1()
              << " P1X-Unbekannte: " << idx.NumUnknowns << '\n';

    VecDescCL beta( &idx), b( &idx);

    // For the ensight visualisation of the piecewise quadratic function f1-f2
    IdxDescCL p2idx( P2_FE);
    CreateNumb( /*Level*/ mg.GetLastLevel(), p2idx, mg, NoBndDataCL<>());
    VecDescCL uneg( &p2idx), upos( &p2idx);
    InitPiecewiseP2( f1, uneg, f2, upos, mg, lset.Phi);

    // Setup the mass matrix
    MatDescCL M( &idx, &idx);
    SetupPrMass_P1X( mg, CoeffCL(), &M, lset, extidx);

    // Setup the right hand side
    IdxT Numb[4];
    double absdet;
    InterfacePatchCL cut;
    Quad5CL<> qf;

    DROPS_FOR_TRIANG_TETRA( mg, mg.GetLastLevel(), sit) {
        GetLocalNumbP1NoBnd( Numb, *sit, idx);
        absdet= sit->GetVolume()*6.0;
        cut.Init( *sit, lset.Phi);

        cut.ComputeSubTets();
        BaryCoordCL* nodes;
        for (Uint i= 0; i < 4; ++i) {
            const IdxT xidx= extidx[Numb[i]];
            const bool have_xidx( xidx != NoIdx);
            double intp1= 0.0, intp1x= 0.0, intpos= 0.0, intneg= 0.0;
            LocalP1CL<> phip1;
            phip1[i]= 1.0;
            for (Uint k=0; k< cut.GetNumTetra(); ++k)
            {
                nodes = qf.TransformNodes(cut.GetTetra(k));
                qf.assign(*sit, cut.GetNumNegTetra()>k ? &f1 : &f2, 0.0, nodes);
                Quad5CL<> qphi(phip1, nodes);
                if (cut.GetNumNegTetra()>k)
                    intneg += Quad5CL<>(qphi*qf).quad(absdet*VolFrac(cut.GetTetra(k)));
                else
                    intpos += Quad5CL<>(qphi*qf).quad(absdet*VolFrac(cut.GetTetra(k)));
                delete nodes;
            }
            intp1+= intpos + intneg;
            if (have_xidx)
                intp1x+= cut.GetSign( i) == 1 ? -intneg : intpos;
            b.Data[Numb[i]]+= intp1;
            if (have_xidx) b.Data[xidx]+= intp1x;
        }
    }

    double intval=0.;
    DROPS_FOR_TRIANG_TETRA( mg, mg.GetLastLevel(), sit) {
        absdet= sit->GetVolume()*6.0;
        cut.Init( *sit, lset.Phi);
        cut.ComputeSubTets();
        BaryCoordCL* nodes;
        for (Uint i=0; i<cut.GetNumTetra(); i++)
        {
            nodes = qf.TransformNodes(cut.GetTetra(i));
            qf.assign( *sit, cut.GetNumNegTetra()>i ? &f1sq : &f2sq, 0.0, nodes);
            intval += qf.quad(absdet*std::fabs( VolFrac(cut.GetTetra(i))));
            delete nodes;
        }
    }

    // Solve the linear system...
    int max_iter= 200;
    double tol= 1e-16;
    PCG( M.Data, beta.Data, b.Data, JACPcCL( 1.0), max_iter, tol, /*measure_relative_tol*/ true);
    std::cerr <<  "iter: " << max_iter << "\ttol: " << tol << '\n';

    //Ensight output
    typedef NoBndDataCL<> BndDataT;
    BndDataT ubnd;
    typedef P2EvalCL<double, BndDataT, VecDescCL> DiscP2CL;
    DiscP2CL upd( &upos, &ubnd, &mg);
    DiscP2CL und( &uneg, &ubnd, &mg);
    EnsightP2SolOutCL ensight( mg, &lset.idx);
    const std::string filename= "ensight/xfem";
    const std::string datgeo= filename+".geo",
    datup = filename+".up" ,
    datun = filename+".un" ,
    datulp = filename+".ulp" ,
    datuln = filename+".uln" ,
    datscl= filename+".scl";

    ensight.CaseBegin( std::string("xfem.case").c_str());
    ensight.DescribeGeom( "Cube", datgeo);
    ensight.DescribeScalar( "Levelset", datscl);
    ensight.DescribeScalar( "up", datup);
    ensight.DescribeScalar( "un", datun);
    ensight.DescribeScalar( "ulp", datulp);
    ensight.DescribeScalar( "uln", datuln);
    ensight.putGeom( datgeo);
    ensight.putScalar( datscl, lset.GetSolution());
    ensight.putScalar( datup, upd);
    ensight.putScalar( datun, und);

    // Output the L_2-projection
    VecDescCL ulneg, ulpos;
    P1XOnPart( beta, extidx, ulpos, lset, true);
    P1XOnPart( beta, extidx, ulneg, lset, false);
    typedef P1EvalCL<double, BndDataT, VecDescCL> DiscP1CL;
    DiscP1CL ulpd( &ulpos, &ubnd, &mg);
    DiscP1CL ulnd( &ulneg, &ubnd, &mg);
    ensight.putScalar( datulp, ulpd);
    ensight.putScalar( datuln, ulnd);
    ensight.Commit();
    ensight.CaseEnd();

    //1D-Plots
    std::ofstream out ("u.txt");
    std::ofstream outpr ("up.txt");
    Point3DCL p;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    DROPS_FOR_TRIANG_VERTEX( mg, mg.GetLastLevel(), it) {
        p= it->GetCoord();
        if (p[0]==0 && p[1]==0) {
            const double u=  lset.Phi.Data[it->Unknowns( lset.idx.GetIdx())] > 0. ? upos.Data[it->Unknowns( p2idx.GetIdx())] : uneg.Data[it->Unknowns( p2idx.GetIdx())];
            const IdxT nr= it->Unknowns( idx.GetIdx());
            const double up= lset.Phi.Data[it->Unknowns( lset.idx.GetIdx())] > 0. ? ulpos.Data[nr] : ulneg.Data[nr];
            out << p[2] << " " << u << '\n';
            outpr << p[2] << " " << up << '\n';
        }
    }

    std::cerr << std::setprecision(20);
    std::cerr << "||u_l||_0^2 = " << dot (M.Data*beta.Data, beta.Data) <<'\n';
    std::cerr << "|| u ||_0^2 = " << intval << '\n';
    const double err = intval - dot (M.Data*beta.Data, beta.Data);
    std::cerr << "|| u - u_l ||_0 = " << std::scientific << std::sqrt(std::abs(err)) << '\n';
    if (err<0) std::cerr << "Fehler: Norm der Projektion > Norm der Funktion\n";
  } catch (DROPSErrCL d) {
        d.handle();
    }

    return 0;
}