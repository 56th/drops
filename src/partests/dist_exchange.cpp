/// \file dist_exchange.cpp
/// \brief test of exchange for numerical data
/// \author LNM RWTH Aachen: Patrick Esser; SC RWTH Aachen: Oliver Fortmeier

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

#include "parallel/parallel.h"
#include "misc/container.h"
#include "DiST/DiST.h"
#include "parallel/exchange.h"
#include "parallel/partime.h"
#include "parallel/loadbal.h"
#include "geom/simplex.h"
#include "geom/builder.h"
#include "geom/multigrid.h"
#include "misc/problem.h"
#include "num/discretize.h"
#include <iostream>
#include <sstream>

using namespace std;

double dist_func( const DROPS::Point3DCL& p, double) { return norm(p) - 0.5; }
double zero_func( const DROPS::Point3DCL&,   double) { return 0; }

namespace DROPS{

void InitP2( const MultiGridCL& mg, VecDescCL& vec, instat_scalar_fun_ptr func)
{
    Uint lvl= vec.GetLevel(),
         idx= vec.RowIdx->GetIdx();


    for (MultiGridCL::const_TriangVertexIteratorCL sit= mg.GetTriangVertexBegin(lvl),
        send= mg.GetTriangVertexEnd(lvl); sit != send; ++sit)  //Vertices
    {
        if (sit->Unknowns.Exist(idx))
        {
            vec.Data[sit->Unknowns(idx)]= func( sit->GetCoord(), 0);
        }
    }

    for (MultiGridCL::const_TriangEdgeIteratorCL sit= mg.GetTriangEdgeBegin(lvl),
        send= mg.GetTriangEdgeEnd(lvl); sit!=send; ++sit)      //Edges
    {
        if ( sit->Unknowns.Exist(idx))
            vec.Data[sit->Unknowns(idx)]= func( GetBaryCenter( *sit), 0);
    }

}


/// \brief Build a brick and tell parallel info class about the multigrid
void BuildBrick( MultiGridCL*& mg)
{
    Point3DCL origin, e1, e2, e3;
    e1[0]= e2[1]= e3[2]= 1.;
    Uint ref[3]= { 2, 2, 2};
    BrickBuilderCL builder( origin, e1, e2, e3, ref[0], ref[1], ref[2]);
    mg = new MultiGridCL( builder);
}

/// \brief Check, if accumulation is correct
/** Use 1 as value for all dof. Hence, the accumulated value must equal the
    number of copies of the simplex.
*/
bool CheckAccumulation( MultiGridCL& mg)
{
    const size_t num_idx=4;
    std::vector<IdxDescCL*> idx_desc(num_idx);
    std::vector<VecDescCL*> vec_desc(num_idx);
    idx_desc[0]= new IdxDescCL( P1_FE);
    idx_desc[1]= new IdxDescCL( P2_FE);
    idx_desc[2]= new IdxDescCL( vecP2_FE);
    idx_desc[3]= new IdxDescCL( P1X_FE);

    IdxDescCL lidx( P2_FE);
    lidx.CreateNumbering( mg.GetLastLevel(), mg);
    VecDescCL lset( &lidx);
    InitP2( mg, lset, dist_func);
    BndDataCL<> lbnd( mg.GetBnd().GetNumBndSeg());

    for ( size_t i=0; i<num_idx; ++i){
        if (idx_desc[i]->IsExtended())
            idx_desc[i]->CreateNumbering( mg.GetLastLevel(), mg, &lset, &lbnd);
        else
            idx_desc[i]->CreateNumbering( mg.GetLastLevel(), mg);
        vec_desc[i]= new VecDescCL( idx_desc[i]);
        vec_desc[i]->Data= 1.0;
    }

    for (size_t i = 0; i< num_idx; ++i) {
        idx_desc[i]->GetEx().Accumulate( vec_desc[i]->Data);
    }

    std::cout << " Number of neighbors of process " << ProcCL::Master()
              << ": " << idx_desc[0]->GetEx().GetNumNeighs() << std::endl;
    bool correct=true;
    // Check on vertices
    DROPS_FOR_TRIANG_VERTEX( mg, mg.GetLastLevel(), it){
        if ( it->Unknowns.Exist()){
            for ( size_t i=0; i<num_idx; ++i){
                const Uint idx=idx_desc[i]->GetIdx(),
                        lvl= idx_desc[i]->TriangLevel();
                if (it->Unknowns.Exist( idx) && it->Unknowns.InTriangLevel(lvl)){
                    const IdxT dof= it->Unknowns( idx);
                    const double numDist= idx_desc[i]->GetEx().GetNumProcs(dof)+1;
                    for ( Uint j=0; j<idx_desc[i]->NumUnknownsVertex(); ++j){
                        if ( numDist!=vec_desc[i]->Data[dof+j]){
                            correct= false;
                            printf("Proc nr %d, testcase %ld: %f should be %f, Prio is %i on level %i with bary: %f %f %f\n",
                            		ProcCL::MyRank(), i, vec_desc[i]->Data[dof+j], numDist, it->GetPrio(),
                            		it->GetLevel(), it->GetBary()[0], it->GetBary()[1], it->GetBary()[2]);
                            it->DebugInfo(cdebug);
                        }
                    }
                    if (idx_desc[i]->IsExtended() && idx_desc[i]->IsExtended(dof)) {
                        const IdxT xdof= idx_desc[i]->GetXidx()[dof];
                        const double numDistX= idx_desc[i]->GetEx().GetNumProcs(xdof)+1;
                        for ( Uint j=0; j<idx_desc[i]->NumUnknownsVertex(); ++j){
                            if ( numDistX!=vec_desc[i]->Data[xdof+j]){
                                correct= false;
                                printf("Proc nr %d, testcase %ld (XFEM): %f should be %f, Prio is %i on level %i with bary: %f %f %f\n",
                                        ProcCL::MyRank(), i, vec_desc[i]->Data[xdof+j], numDistX, it->GetPrio(),
                                        it->GetLevel(), it->GetBary()[0], it->GetBary()[1], it->GetBary()[2]);
                                it->DebugInfo(cdebug);
                            }
                        }
                    }
                }
            }
        }
    }
    // Check on edges
    DROPS_FOR_TRIANG_EDGE( mg, mg.GetLastLevel(), it){
        if ( it->Unknowns.Exist()){
            for ( size_t i=0; i<num_idx; ++i){
                const Uint idx=idx_desc[i]->GetIdx();
                if (it->Unknowns.Exist( idx) && it->Unknowns.InTriangLevel(idx_desc[i]->TriangLevel())){
                    const IdxT dof= it->Unknowns( idx);
                    const double numDist= idx_desc[i]->GetEx().GetNumProcs(dof)+1;
                    for ( Uint j=0; j<idx_desc[i]->NumUnknownsEdge(); ++j)
                        if ( numDist!=vec_desc[i]->Data[dof+j]){
                            correct= false;
                            printf("Proc nr %d, testcase %ld: %f should be %f, Prio is %i on level %i with bary: %f %f %f\n",
                            		ProcCL::MyRank(), i, vec_desc[i]->Data[dof+j], numDist, it->GetPrio(),
                            		it->GetLevel(), it->GetBary()[0], it->GetBary()[1], it->GetBary()[2]);
                            it->DebugInfo(cdebug);
                    }
                }
            }
        }
    }

    // Check ExchangeBlockCL
    ExchangeBlockCL exblock;
    exblock.AttachTo( idx_desc[2]->GetEx());
    exblock.AttachTo( idx_desc[1]->GetEx());
    exblock.AttachTo( idx_desc[0]->GetEx());
    VectorCL x( idx_desc[2]->NumUnknowns() + idx_desc[1]->NumUnknowns() + idx_desc[0]->NumUnknowns());
    x= 1.0;
    exblock.Accumulate( x);
    for ( size_t i=0; i<idx_desc[2]->NumUnknowns(); ++i)
        if ( vec_desc[2]->Data[i]!=x[i])
            correct=false;
    for ( size_t i=0; i<idx_desc[1]->NumUnknowns(); ++i)
        if ( vec_desc[1]->Data[i]!=x[i+idx_desc[2]->NumUnknowns()])
            correct=false;
    for ( size_t i=0; i<idx_desc[0]->NumUnknowns(); ++i)
        if ( vec_desc[0]->Data[i]!=x[i+idx_desc[2]->NumUnknowns()+idx_desc[1]->NumUnknowns()])
            correct=false;

    for ( size_t i=0; i<num_idx; ++i){
        delete idx_desc[i];
        delete vec_desc[i];
    }
    return ProcCL::Check(correct);
}

std::valarray<double> getReferenceValue( const MultiGridCL& mg, IdxDescCL* idxDesc)
{
    std::valarray<double> result(3);
    result = 0.0;
    const Uint idx= idxDesc->GetIdx(),
            numUnksVert= idxDesc->NumUnknownsVertex(),
            numUnksEdge= idxDesc->NumUnknownsEdge(),
            lvl= idxDesc->TriangLevel();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, mg.GetLastLevel(), it){
        if ( it->Unknowns.Exist() && it->Unknowns.Exist(idx) && it->Unknowns.InTriangLevel(lvl)){
            const bool withXdof= idxDesc->IsExtended() && idxDesc->IsExtended(it->Unknowns(idx));
            if ( it->IsLocal()){
                const int xmul= withXdof ? 2 : 1;
                result[0] += numUnksVert*xmul;
                result[1] += numUnksVert*xmul;
                result[2] += numUnksVert*xmul;
            }
            else if ( idxDesc->GetEx().AmIOwner(it->Unknowns(idx))){
                const double numDist= idxDesc->GetEx().GetNumProcs(it->Unknowns(idx))+1;
                result[0] += numUnksVert;
                result[1] += numDist*numUnksVert;
                result[2] += numDist*numDist*numUnksVert;
                if (withXdof) {
                    const double numDistX= idxDesc->GetEx().GetNumProcs(idxDesc->GetXidx()[it->Unknowns(idx)])+1;
                    result[0] += numUnksVert;
                    result[1] += numDistX*numUnksVert;
                    result[2] += numDistX*numDistX*numUnksVert;
                }
            }
        }
    }
    // Check on edges
    DROPS_FOR_TRIANG_CONST_EDGE( mg, mg.GetLastLevel(), it){
        if ( it->Unknowns.Exist() && it->Unknowns.Exist(idx) && it->Unknowns.InTriangLevel(lvl)){
            if ( it->IsLocal()){
                result[0] += numUnksEdge;
                result[1] += numUnksEdge;
                result[2] += numUnksEdge;
            }
            else if ( idxDesc->GetEx().AmIOwner(it->Unknowns(idx))){
                const double numDist= idxDesc->GetEx().GetNumProcs(it->Unknowns(idx))+1;
                result[0] += numUnksEdge;
                result[1] += numDist*numUnksEdge;
                result[2] += numDist*numDist*numUnksEdge;
            }
        }
    }
    return ProcCL::GlobalSum(result);
}

bool CheckDot( const ExchangeCL& Ex, const VectorCL& x, bool isXacc,
    const VectorCL& y, bool isYacc, double result)
{
    VectorCL x_acc(x), y_acc(y), x_acc_ref(x), y_acc_ref(y), x_acc_norm(x);
    if ( !isXacc)
        Ex.Accumulate(x_acc_ref);
    if ( !isYacc)
        Ex.Accumulate(y_acc_ref);

    const double dot1= Ex.ParDot( x, isXacc, y, isYacc);
    const double dot2= Ex.ParDot( x, isXacc, y, isYacc, &x_acc, &y_acc);
    const double normx= Ex.Norm( x, isXacc, &x_acc_norm);
    const double normx_dot= std::sqrt(Ex.ParDot( x, isXacc, x, isXacc));

    const bool same_result=   ProcCL::Check((bool)(dot1==dot2));
    const bool x_acc_correct= Ex.Norm( VectorCL(x_acc-x_acc_ref), true)<1e-10;
    const bool y_acc_correct= Ex.Norm( VectorCL(y_acc-y_acc_ref), true)<1e-10;
    const bool x_acc_norm_correct= Ex.Norm( VectorCL(x_acc_norm-x_acc_ref), true)<1e-10;
    const bool result_correct= dot1==result;
    const bool norm_correct= normx==normx_dot;

    if ( !same_result){
        std::cout << " - Call with demand for accumulated values gives different results" << std::endl;
    }
    if ( !x_acc_correct){
        std::cout << " - x is not correct accumulated" << std::endl;
    }
    if ( !y_acc_correct){
        std::cout << " - y is not correct accumulated" << std::endl;
    }
    if ( !result_correct){
        std::cout << " - Inner product is wrong. Computed " << dot1 << ", but should be " << result << std::endl;
    }
    if ( !x_acc_norm_correct){
        std::cout << " - x is not correct accumulated in norm computation" << std::endl;
    }
    if ( !norm_correct){
        std::cout << " - Norm is not correct where isXacc=" << isXacc << std::endl;
    }

    return same_result && x_acc_correct && y_acc_correct && result_correct && x_acc_norm_correct && norm_correct;
}

/** \brief Check if ExchangeBlockCL delivers the same result as all its ExchangeCLs*/
bool CheckDot( const ExchangeBlockCL& Ex, const VectorCL& x, bool isXacc,
    const VectorCL& y, bool isYacc)
{
    VectorCL x_acc(x), y_acc(y), x_acc_ref(x), y_acc_ref(y), x_acc_norm(x);
    if ( !isXacc){
        for ( size_t j=0; j<Ex.GetNumBlocks(); ++j){
            const size_t offset=Ex.GetOffset(j);
            VectorCL v( Ex.GetEx(j).GetNum());
            std::copy( Addr(x)+offset, Addr(x)+v.size()+offset, Addr(v));
            Ex.GetEx(j).Accumulate(v);
            std::copy( Addr(v), Addr(v)+v.size(), Addr(x_acc_ref)+offset);
        }
    }
    if ( !isYacc){
        for ( size_t j=0; j<Ex.GetNumBlocks(); ++j){
            const size_t offset=Ex.GetOffset(j);
            VectorCL v( Ex.GetEx(j).GetNum());
            std::copy( Addr(y)+offset, Addr(y)+v.size()+offset, Addr(v));
            Ex.GetEx(j).Accumulate(v);
            std::copy( Addr(v), Addr(v)+v.size(), Addr(y_acc_ref)+offset);
        }
    }
    double result=0;
    for ( size_t j=0; j<Ex.GetNumBlocks(); ++j){
        const size_t offset=Ex.GetOffset(j);
        VectorCL v( Ex.GetEx(j).GetNum());
        VectorCL u( Ex.GetEx(j).GetNum());
        std::copy( Addr(x)+offset, Addr(x)+v.size()+offset, Addr(v));
        std::copy( Addr(y)+offset, Addr(y)+u.size()+offset, Addr(u));
        result += Ex.GetEx(j).ParDot( v, isXacc, u, isYacc);
    }

    const double dot1= Ex.ParDot( x, isXacc, y, isYacc);
    const double dot2= Ex.ParDot( x, isXacc, y, isYacc, &x_acc, &y_acc);
    const double normx= Ex.Norm( x, isXacc, &x_acc_norm);
    const double normx_dot= std::sqrt(Ex.ParDot( x, isXacc, x, isXacc));

    const bool same_result=   ProcCL::Check((bool)(dot1==dot2));
    const bool x_acc_correct= Ex.Norm( VectorCL(x_acc-x_acc_ref), true)<1e-10;
    const bool y_acc_correct= Ex.Norm( VectorCL(y_acc-y_acc_ref), true)<1e-10;
    const bool x_acc_norm_correct= Ex.Norm( VectorCL(x_acc_norm-x_acc_ref), true)<1e-10;
    const bool result_correct= dot1==result;
    const bool norm_correct= normx==normx_dot;

    if ( !same_result){
        std::cout << " - Call with demand for accumulated values gives different results" << std::endl;
    }
    if ( !x_acc_correct){
        std::cout << " - x is not correct accumulated" << std::endl;
    }
    if ( !y_acc_correct){
        std::cout << " - y is not correct accumulated" << std::endl;
    }
    if ( !result_correct){
        std::cout << " - Inner product is wrong. Computed " << dot1 << ", but should be " << result << std::endl;
    }
    if ( !x_acc_norm_correct){
        std::cout << " - x is not correct accumulated in norm computation" << std::endl;
    }
    if ( !norm_correct){
        std::cout << " - Norm is not correct where isXacc=" << isXacc << std::endl;
    }

    return same_result && x_acc_correct && y_acc_correct && result_correct && x_acc_norm_correct && norm_correct;
}


/// \brief Check, if inner product computation is correct
/** Use 1 as value for all dof. Hence, the accumulated value must equal the
    number of copies of the simplex.
*/
bool CheckInnerProducts( MultiGridCL& mg)
{
    const size_t num_idx=4;
    bool correct= true;
    std::vector<IdxDescCL*> idx_desc(num_idx);
    std::vector<VecDescCL*> vec_desc(num_idx);
    idx_desc[0]= new IdxDescCL( P1_FE);
    idx_desc[1]= new IdxDescCL( P2_FE);
    idx_desc[2]= new IdxDescCL( vecP2_FE);
    idx_desc[3]= new IdxDescCL( P1X_FE);

    IdxDescCL lidx( P2_FE);
    lidx.CreateNumbering( mg.GetLastLevel(), mg);
    VecDescCL lset( &lidx);
    InitP2( mg, lset, dist_func);
    BndDataCL<> lbnd( mg.GetBnd().GetNumBndSeg());

    for ( size_t i=0; i<num_idx; ++i){
        std::cout << " - Check case " << (i+1) << std::endl;

        if (idx_desc[i]->IsExtended())
            idx_desc[i]->CreateNumbering( mg.GetLastLevel(), mg, &lset, &lbnd);
        else
            idx_desc[i]->CreateNumbering( mg.GetLastLevel(), mg);
        vec_desc[i]= new VecDescCL( idx_desc[i]);
        vec_desc[i]->Data= 1.0;
        VectorCL x( vec_desc[i]->Data), y( x);
        std::valarray<double> ref_result= getReferenceValue( mg, idx_desc[i]);

        bool case_correct= true;
        if ( !CheckDot( idx_desc[i]->GetEx(), x, true, y, true, ref_result[0])){
            case_correct= false;
        }
        if ( !CheckDot( idx_desc[i]->GetEx(), x, true, y, false, ref_result[1])){
            case_correct= false;
        }
        if ( !CheckDot( idx_desc[i]->GetEx(), x, false, y, true, ref_result[1])){
            case_correct= false;
        }
        if ( !CheckDot( idx_desc[i]->GetEx(), x, false, y, false, ref_result[2])){
            case_correct= false;
        }

        if ( !ProcCL::Check(case_correct)) correct=false;
        else std::cout << "   => Case " << (i+1) << " is correct" << std::endl;
    }

    ExchangeBlockCL exBlock;
    exBlock.AttachTo( idx_desc[2]->GetEx());
    exBlock.AttachTo( idx_desc[1]->GetEx());
    exBlock.AttachTo( idx_desc[0]->GetEx());

    std::cout << " - Check case 5" << std::endl;
    bool case_correct=true;
    VectorCL x( 1.0, exBlock.GetNum());
    VectorCL y( 1.0, exBlock.GetNum());
    if ( !CheckDot(exBlock, x, true, y, true)){
        case_correct= false;
    }
    if ( !CheckDot(exBlock, x, true, y, false)){
        case_correct= false;
    }
    if ( !CheckDot(exBlock, x, false, y, true)){
        case_correct= false;
    }
    if ( !CheckDot(exBlock, x, false, y, false)){
        case_correct= false;
    }
    if ( !ProcCL::Check(case_correct)) correct=false;
    else std::cout << "   => Case 5 is correct" << std::endl;

    for ( size_t i=0; i<num_idx; ++i){
        delete idx_desc[i];
        delete vec_desc[i];
    }

    return correct;
}


void MakeTimeMeasurements( MultiGridCL& mg, const size_t num_test)
{
    IdxDescCL idx( vecP2_FE);
    idx.CreateNumbering( mg.GetLastLevel(), mg);
    ParTimerCL timer;
    VectorCL x( 1., idx.NumUnknowns()), x_acc(x), y(1., idx.NumUnknowns()), y_acc(y);
    const size_t size_acc   = ProcCL::GlobalSum( x.size()),
                 size_global= idx.GetGlobalNumUnknowns();

    if ( idx.GetEx().CommViaOwner())
        std::cout << "    o using communication via owner\n";
    else
        std::cout << "    o using direct communication via neighbors\n";
    std::cout << "    o vector length: global " << size_global
              << ", accumulated " << size_acc << std::endl;

    double dummy=0;
    std::cout << " - on two distributed vectors ... ";
    timer.Reset();
    for ( size_t i=0; i<num_test; ++i){
        dummy+= idx.GetEx().LocalDot( x, false, y, false, &x_acc, &y_acc);
    }
    timer.Stop();
    std::cout << "took " << timer.GetTime() << " sec\n";

    std::cout << " - on a distributed and an accumulated vector ... ";
    timer.Reset();
    for ( size_t i=0; i<num_test; ++i){
        dummy+= idx.GetEx().LocalDot( x, true, y, false, &x_acc, &y_acc);
    }
    timer.Stop();
    std::cout << "took " << timer.GetTime() << " sec\n";

    std::cout << " - on two accumulated vectors ... ";
    timer.Reset();
    for ( size_t i=0; i<num_test; ++i){
        dummy+= idx.GetEx().LocalDot( x, true, y, true);
    }
    timer.Stop();
    std::cout << "took " << timer.GetTime() << " sec\n";
}

} // end of namespace DROPS

int main( int argc, char **argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
#endif
    try {
        DROPS::MultiGridCL* mg= 0;
        DROPS::BuildBrick( mg);
        DROPS::LoadBalCL lb( *mg);      // loadbalancing
        lb.DoMigration();        // distribute initial grid
        for (int i=0; i<2; ++i) {
            MarkAll( *mg);
            mg->Refine();
            mg->SizeInfo( std::cout);
            lb.DoMigration();
            mg->SizeInfo( std::cout);
        }

/*        for (DROPS::MultiGridCL::VertexIterator it = mg->GetAllVertexBegin(); it != mg->GetAllVertexEnd(); ++it){
      		if ( it->GetBary()[0] == 0.5 && it->GetBary()[1] == 0.0 && it->GetBary()[2] == 0.5){
  				if (it->AmIOwner())
  					printf("vertex is on proc nr %d, level: %d, prio: %d and AmIOwner is true with coord %f %f %f\n",
  							DROPS::ProcCL::MyRank(), it->GetLevel(), it->GetPrio(),
  							it->GetBary()[0], it->GetBary()[1], it->GetBary()[2]);
  				else
  					printf("vertex is on proc nr %d, level: %d, prio: %d and AmIOwner is false with coord %f %f %f\n",
  							DROPS::ProcCL::MyRank(), it->GetLevel(), it->GetPrio(),
  							it->GetBary()[0], it->GetBary()[1], it->GetBary()[2]);
        	}
        }
*/
/*        for (DROPS::MultiGridCL::TetraIterator it = mg->GetAllTetraBegin(); it != mg->GetAllTetraEnd(); ++it){
        	for (DROPS::Uint i=0; i<4; ++i) {
        		if ( it->GetVertex(i)->GetBary()[0] == 0.5 && it->GetVertex(i)->GetBary()[1] == 0.0 && it->GetVertex(i)->GetBary()[2] == 0.5){
        			if (it->HasGhost())
        				printf("proc nr %d, has ghost is true, level: %d\n", DROPS::ProcCL::MyRank(), it->GetLevel());
        			else
        				printf("proc nr %d, has ghost is false, level: %d\n", DROPS::ProcCL::MyRank(), it->GetLevel());
        		}
        	}
        }
*/


        std::cout << "Check accumulation and inner products for four cases:\n"
                  << " 1) P1_FE\n"
                  << " 2) P2_FE\n"
                  << " 3) vecP2_FE\n"
                  << " 4) P1X_FE\n"
                  << " 5) Blocked 3), 2), and 1)"
                  << std::endl;

        if ( !DROPS::CheckAccumulation( *mg)){
            std::cout << " Accumulation is broken" << std::endl;
        }
        else{
            std::cout << " Accumulation seems to be alright" << std::endl;
        }
        if ( !DROPS::CheckInnerProducts( *mg)){
            std::cout << " Computing inner products is broken" << std::endl;
        }
        else{
            std::cout << " Computing inner products seems to be alright" << std::endl;
        }

        const size_t num_test=1000;
        std::cout << "Performing time measurements for " << num_test
                  << " parallel dots ..." << std::endl;
        MakeTimeMeasurements( *mg, num_test);

        std::string filename("sane.chk");
        DROPS::ProcCL::AppendProcNum(filename);
        std::ofstream sanityfile( filename.c_str());
        const DROPS::DiST::InfoCL& info= DROPS::DiST::InfoCL::Instance();
        info.SizeInfo(sanityfile);
        info.Instance().DebugInfo(sanityfile);
        sanityfile << "\n\n===========================\n\n";

        if ( DROPS::ProcCL::Check( DROPS::DiST::InfoCL::Instance().IsSane( sanityfile))){
            std::cout << " DiST-module seems to be alright!" << std::endl;
        }
        else{
            std::cout << " DiST-module seems to be broken!" << std::endl;
        }
        delete mg; mg=0;

    }
    catch (DROPS::DROPSErrCL err) {err.handle();}
    return 0;
}
