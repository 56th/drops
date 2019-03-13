/// \file problem.cpp
/// \brief Classes that constitute a problem.
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#include "misc/problem.h"
#include "num/interfacePatch.h"
#include "num/quadrature.h"
#include "num/lattice-eval.h"
#include "parallel/exchange.h"


namespace DROPS
{

const Uint        IdxDescCL::InvalidIdx = std::numeric_limits<Uint>::max();
std::vector<bool> IdxDescCL::IdxFree;

IdxDescCL::IdxDescCL( FiniteElementT fe, const BndCondCL& bnd, double omit_bound)
    : FE_InfoCL( fe), Idx_( GetFreeIdx()), TriangLevel_( 0), NumUnknowns_( 0), Bnd_(bnd), Bnd_aux_(bnd),
      extIdx_( omit_bound != -99 ? omit_bound : IsExtended() ? 1./32. : -1.) 
// default value for omit_bound is 1./32. for XFEM and -1 otherwise
{
#ifdef _PAR
    ex_= new ExchangeCL();
#else
    ex_= new DummyExchangeCL();
#endif
}

IdxDescCL::IdxDescCL( FiniteElementT fe, const BndCondCL& bnd1, const BndCondCL& bnd2, double omit_bound)
    : FE_InfoCL( fe), Idx_( GetFreeIdx()), TriangLevel_( 0), NumUnknowns_( 0), Bnd_(bnd1), Bnd_aux_(bnd2),
      extIdx_( omit_bound != -99 ? omit_bound : IsExtended() ? 1./32. : -1.) 
// default value for omit_bound is 1./32. for XFEM and -1 otherwise
{
#ifdef _PAR
    ex_= new ExchangeCL();
#else
    ex_= new DummyExchangeCL();
#endif
}

IdxDescCL::~IdxDescCL()
{
    if (Idx_!=InvalidIdx)
        IdxFree[Idx_]= true;
    delete ex_;
}

Uint IdxDescCL::GetFreeIdx()
{
    size_t sysnum= 0;
    for (; sysnum<IdxFree.size(); ++sysnum)
        if (IdxFree[sysnum]) break;
    if (sysnum>=IdxFree.size())
        IdxFree.push_back( false);
    else
        IdxFree[sysnum]= false;
    return sysnum;
}

IdxDescCL::IdxDescCL( const IdxDescCL& orig)
 : FE_InfoCL(orig), Idx_(orig.Idx_), TriangLevel_(orig.TriangLevel_), NumUnknowns_(orig.NumUnknowns_),
   Bnd_(orig.Bnd_), Bnd_aux_(orig.Bnd_aux_), extIdx_(orig.extIdx_)
{
    // invalidate orig
    const_cast<IdxDescCL&>(orig).Idx_= InvalidIdx;
#ifdef _PAR
    ex_= new ExchangeCL(*orig.ex_);
#else
    ex_= new DummyExchangeCL(*orig.ex_);
#endif

}

void IdxDescCL::swap( IdxDescCL& obj)
/// Note, that std::swap cannot be used for IdxDescCL-objects as the
/// assignment operator is not implemented.
{
    Assert( GetFE()==obj.GetFE(), DROPSErrCL("IdxDescCL::swap: FE-types differ"), ~0);
        std::swap( Idx_,         obj.Idx_);
    std::swap( TriangLevel_, obj.TriangLevel_);
    std::swap( NumUnknowns_, obj.NumUnknowns_);
    std::swap( Bnd_,         obj.Bnd_);
    std::swap( extIdx_,      obj.extIdx_);
    std::swap( ex_,          obj.ex_);
}

bool IdxDescCL::Equal(IdxDescCL& i, IdxDescCL& j, const MultiGridCL* mg)
/// \param i The left IdxDescCL.
/// \param j The right IdxDescCL.
/// \param mg Optional pointer to a multigrid. If it is given the numbers
///     on the simplices are compared, too. This is rather expensive and
///     only needed for some correctness tests.
/// \return If mg==0: True, iff all members of i and j have the same value.
///     If mg!=0: True, iff all members of i and j have the same value and
///     all numbers on the simplices of the given triangulation are equal.
{
    const Uint lvl= i.TriangLevel();
    if (lvl != j.TriangLevel()) {
        std::cout << "Compare_Indices: Indices on different levels.\n";
        return false;
    }
    if (i.NumUnknowns() != j.NumUnknowns()) {
        std::cout << "Compare_Indices: NumUnknowns different.\n";
        return false;
    }
    if (i.GetFE() !=  j.GetFE()) {
        std::cout << "Compare_Indices: FE types different.\n";
        return false;
    }
    if (i.NumUnknownsVertex_ != j.NumUnknownsVertex_) {
        std::cout << "Compare_Indices: NumUnknownsVertex different.\n";
        return false;
    }
    if (i.NumUnknownsEdge_ != j.NumUnknownsEdge_) {
        std::cout << "Compare_Indices: NumUnknownsEdge different.\n";
        return false;
    }
    if (i.NumUnknownsFace_ != j.NumUnknownsFace_) {
        std::cout << "Compare_Indices: NumUnknownsFace different.\n";
        return false;
    }
    if (i.NumUnknownsTetra_ != j.NumUnknownsTetra_) {
        std::cout << "Compare_Indices: NumUnknownsTetra different.\n";
        return false;
    }
    if (!mg)
        return true;

    const Uint iidx= i.GetIdx(),
               jidx= j.GetIdx();
    if (iidx == jidx)
        return true;
    if ( i.NumUnknownsVertex_ != 0)
        for (MultiGridCL::const_TriangVertexIteratorCL it= mg->GetTriangVertexBegin( lvl),
             theend= mg->GetTriangVertexEnd( lvl); it != theend; ++it)
            if (it->Unknowns.Exist())
                if ( (it->Unknowns.Exist( iidx) != it->Unknowns.Exist( jidx)) ) {
                    std::cout << "Compare_Indices: Vertex difference.\n";
                    return false;
                }
    if (i.NumUnknownsEdge_ != 0)
        for (MultiGridCL::const_TriangEdgeIteratorCL it= mg->GetTriangEdgeBegin( lvl),
             theend= mg->GetTriangEdgeEnd( lvl); it != theend; ++it)
            if (it->Unknowns.Exist())
                if ( (it->Unknowns.Exist( iidx) != it->Unknowns.Exist( jidx)) ) {
                    std::cout << "Compare_Indices: Edge difference.\n";
                    return false;
                }
    if (i.NumUnknownsFace_ != 0)
        for (MultiGridCL::const_TriangFaceIteratorCL it= mg->GetTriangFaceBegin( lvl),
             theend= mg->GetTriangFaceEnd( lvl); it != theend; ++it)
            if (it->Unknowns.Exist())
                if ( (it->Unknowns.Exist( iidx) != it->Unknowns.Exist( jidx)) ) {
                    std::cout << "Compare_Indices: Face difference.\n";
                    return false;
                }
    if ( i.NumUnknownsTetra_ != 0)
        for (MultiGridCL::const_TriangTetraIteratorCL it= mg->GetTriangTetraBegin( lvl),
             theend= mg->GetTriangTetraEnd( lvl); it != theend; ++it)
            if (it->Unknowns.Exist())
                if ( (it->Unknowns.Exist( iidx) != it->Unknowns.Exist( jidx)) ) {
                    std::cout << "Compare_Indices: Tetra difference.\n";
                    return false;
                }
    return true;
}


#ifdef _PAR
/// \brief Count number of owned unknowns on this proc
/** Get number of unknowns that are local or owned by this proc*/
IdxT IdxDescCL::GetNumOwnedUnknowns() const
{
    return GetEx().GetNumOwned();
}

/// \brief Count global number of unknowns
/** Get global number over all procs of unknowns. Each unknown is just count one time*/
IdxT IdxDescCL::GetGlobalNumUnknowns() const
{
    return ProcCL::GlobalSum(GetNumOwnedUnknowns());
}
#endif

void P1XtoP1 (const IdxDescCL& xidx, const VectorCL& p1x, const IdxDescCL& idx, VectorCL& posPart, VectorCL& negPart, const VecDescCL& lset, const MultiGridCL& mg)
{
    const Uint lvl= idx.TriangLevel(),
               p1idxnum= idx.GetIdx(),
                 idxnum= xidx.GetIdx(),
               lsidxnum= lset.RowIdx->GetIdx();
    const ExtIdxDescCL& extIdx= xidx.GetXidx();
    const size_t p1unknowns = extIdx.GetNumUnknownsStdFE();

    negPart.resize(p1unknowns);
    posPart.resize(p1unknowns);
    posPart = negPart = p1x[std::slice(0, p1unknowns, 1)];

    // add extended pressure
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        const IdxT   nr= it->Unknowns( idxnum);
        if (!it->Unknowns.Exist( idxnum)) continue;
        const IdxT p1nr= it->Unknowns( p1idxnum);
        if (extIdx[nr]==NoIdx) continue;

        if (InterfacePatchCL::Sign( lset.Data[it->Unknowns(lsidxnum)]) == 1)
            negPart[p1nr]= p1x[nr] - p1x[extIdx[nr]];
        else
            posPart[p1nr]= p1x[nr] + p1x[extIdx[nr]];
    }
}

void P1toP1X (const IdxDescCL& xidx, VectorCL& p1x, const IdxDescCL& idx, const VectorCL& posPart, const VectorCL& negPart, const VecDescCL& lset, const MultiGridCL& mg)
{
    const Uint lvl= idx.TriangLevel(),
                idxnum= xidx.GetIdx(),
                p1idxnum= idx.GetIdx(),
                lsidxnum= lset.RowIdx->GetIdx();

    p1x.resize(xidx.NumUnknowns());
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
    {
        const IdxT nr= it->Unknowns(idxnum);
        if (!it->Unknowns.Exist( idxnum)) continue;
        const IdxT p1nr= it->Unknowns(p1idxnum);
        const bool is_pos= InterfacePatchCL::Sign( lset.Data[it->Unknowns(lsidxnum)])==1;
        if (is_pos)
            p1x[nr]= posPart[p1nr];
        else
            p1x[nr]= negPart[p1nr];
    }
    // add extended pressure
    const ExtIdxDescCL& extIdx= xidx.GetXidx();
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
    {
        const IdxT nr= it->Unknowns(idxnum);
        if (!it->Unknowns.Exist( idxnum)) continue;
        const IdxT p1nr= it->Unknowns(p1idxnum);
        if (extIdx[nr]==NoIdx) continue;
        p1x[extIdx[nr]]= (posPart[p1nr] - negPart[p1nr]);
    }
}

void ExtractComponent( const VectorCL& vecFE, VectorCL& scalarFE, Uint comp, Uint stride)
{
    Assert( vecFE.size()==scalarFE.size()*stride, DROPSErrCL("ExtractComponent: vector sizes do not match"), DebugNumericC);
    for (size_t i=0, s=scalarFE.size(); i<s; ++i)
        scalarFE[i]= vecFE[i*stride+comp];
}

void CreateNumbOnTetra( const Uint idx, IdxT& counter, Uint stride,
                        const MultiGridCL::TriangTetraIteratorCL& begin,
                        const MultiGridCL::TriangTetraIteratorCL& end, const Uint level)
{
    if (stride == 0) return;
    for (MultiGridCL::TriangTetraIteratorCL it=begin; it!=end; ++it)
    {
    	if (!it->Unknowns.InTriangLevel(level)) continue;
        it->Unknowns.Prepare( idx);
        it->Unknowns(idx)= counter;
        counter+= stride;
    }
}

/// \brief Routine to number unknowns on the vertices surrounding an
/// interface.
///
/// This function allocates memory for the Unknown-indices in system
/// idx on all vertices belonging to tetras between begin and end which
/// are cut by the zero level of lset. If InterfacePatchCL::IntersectsInterior() is true,
/// all vertices are numbered, if only InterfacePatchCL::Intersects() is true, only the
/// vertices with InterfacePatchCL::GetSign(vertex) == 0 are numbered. The other vertices in
/// such tetrahedra obtain NoIdx as number, but they are not counted as unknowns.
///
/// The first number used is the initial value of counter, the next
/// numbers are counter+stride, counter+2*stride, and so on.
/// Upon return, counter contains the first number, that was not used,
/// that is \# Unknowns+stride.
/// A more user friendly interface is provided by IdxDescCL::CreateNumbOnInterface.
void CreateNumbOnInterfaceVertex (const Uint idx, IdxT& counter, Uint stride,
        const MultiGridCL::TriangVertexIteratorCL& vbegin,
        const MultiGridCL::TriangVertexIteratorCL& vend,
        const MultiGridCL::TriangTetraIteratorCL& begin,
        const MultiGridCL::TriangTetraIteratorCL& end,
    const VecDescCL& ls, const BndDataCL<>& lsetbnd, double omit_bound= -1./*default to using all dof*/)
{
    if (stride == 0) return;

    LocalP2CL<> hat_sq[4]; // values of phi_i*phi_i
    for (int i= 0; i < 4; ++i)  {
        hat_sq[i][i]= 1.;
        for (int j = 0; j < 4; ++j)
            if (i != j) hat_sq[i][EdgeByVert( i, j) + 4]= 0.25;
    }
    // first set NoIdx in all vertices
    for (MultiGridCL::TriangVertexIteratorCL vit= vbegin; vit != vend; ++vit) {
        vit->Unknowns.Prepare(idx);
        vit->Unknowns.Invalidate(idx);
    }
    // then create numbering of vertices at the interface
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 2);
    std::valarray<double> ls_loc( lat.vertex_size());
    LocalP2CL<> locp2_ls;
    SPatchCL<3> patch;
    QuadDomainCodim1CL<3> qdom;
    std::valarray<double> shape_sq; // square of a P1-shape-function as integrand
    for (MultiGridCL::TriangTetraIteratorCL it= begin; it != end; ++it) {
        locp2_ls.assign( *it, ls, lsetbnd);
        evaluate_on_vertexes( locp2_ls, lat, Addr( ls_loc));
        if (equal_signs( ls_loc))
            continue;

        const double limit= omit_bound < 0. ? 0. : omit_bound*std::pow( it->GetVolume()*6, 4./3.);
        patch.make_patch<MergeCutPolicyCL>( lat, ls_loc);
        make_CompositeQuad2Domain2D( qdom, patch, *it);
        shape_sq.resize( qdom.vertex_size());
        for (Uint i= 0; i < NumVertsC; ++i) {
            UnknownHandleCL& unknowns= const_cast<VertexCL*>( it->GetVertex( i))->Unknowns;
            if (unknowns.Exist( idx))
                continue;

            evaluate_on_vertexes( hat_sq[i], qdom, Addr( shape_sq));
            if (quad_codim1( shape_sq, qdom) > limit) {
                unknowns( idx)= counter;
                counter+= stride;
            }
        }
    }
}

/// \brief Routine to number unknowns on the vertices in a strip surrounding an
/// interface(distance less than a parameter dist).
///
/// This function allocates memory for the Unknown-indices in system
/// idx on all vertices belonging to tetras between begin and end which
/// are in a narrow band of the zero level of lset.
///
/// The first number used is the initial value of counter, the next
/// numbers are counter+stride, counter+2*stride, and so on.
/// Upon return, counter contains the first number, that was not used,
/// that is \# Unknowns+stride.
/// A more user friendly interface is provided by IdxDescCL::CreateNumbNearInterface.
    void CreateNumbNearInterfaceVertex (const Uint idx, IdxT& counter, Uint stride,
                                        const MultiGridCL::TriangVertexIteratorCL& vbegin,
                                        const MultiGridCL::TriangVertexIteratorCL& vend,
                                        const MultiGridCL::TriangTetraIteratorCL& begin,
                                        const MultiGridCL::TriangTetraIteratorCL& end,
                                        const VecDescCL& ls, const BndDataCL<>& lsetbnd, double dist,double omit_bound= -1./*default to using all dof*/)
    {
        if (stride == 0) return;
        std::cout<<"Notice-----the level-set function is assume to be a signed distance function!!\n";
        std::cout<<"Test "<<dist<<std::endl;
        LocalP2CL<> hat_sq[4]; // values of phi_i*phi_i
        for (int i= 0; i < 4; ++i)  {
            hat_sq[i][i]= 1.;
            for (int j = 0; j < 4; ++j)
                if (i != j) hat_sq[i][EdgeByVert( i, j) + 4]= 0.25;
        }
        // first set NoIdx in all vertices
        for (MultiGridCL::TriangVertexIteratorCL vit= vbegin; vit != vend; ++vit) {
            vit->Unknowns.Prepare(idx);
            vit->Unknowns.Invalidate(idx);
        }
        // then create numbering of vertices at the interface
        const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 2);
        std::valarray<double> ls_loc( lat.vertex_size());
        LocalP2CL<> locp2_ls;
        //  SPatchCL<3> patch;
        // QuadDomainCodim1CL<3> qdom;
        //  std::valarray<double> shape_sq; // square of a P1-shape-function as integrand
        for (MultiGridCL::TriangTetraIteratorCL it= begin; it != end; ++it) {
            locp2_ls.assign( *it, ls, lsetbnd);
            evaluate_on_vertexes( locp2_ls, lat, Addr( ls_loc));
            if (distance( ls_loc)>dist)
                continue;

            //   patch.make_patch<MergeCutPolicyCL>( lat, ls_loc);
            //   make_CompositeQuad2Domain2D( qdom, patch, *it);
            //   shape_sq.resize( qdom.vertex_size());
            for (Uint i= 0; i < NumVertsC; ++i) {
                UnknownHandleCL& unknowns= const_cast<VertexCL*>( it->GetVertex( i))->Unknowns;
                if (unknowns.Exist( idx))
                    continue;

                //      evaluate_on_vertexes( hat_sq[i], qdom, Addr( shape_sq));
                //     if (quad_codim1( shape_sq, qdom) > 0) {
                unknowns( idx)= counter;
                counter+= stride;
                //    }
            }
        }
    }
/// \brief Routine to number P2-unknowns on the vertices and edges surrounding an
/// interface.
///
/// This function allocates memory for the Unknown-indices in system
/// idx on all vertices belonging to tetras between begin and end which
/// are cut by the zero level of lset.
///
/// The first number used is the initial value of counter, the next
/// numbers are counter+stride, counter+2*stride, and so on.
/// Upon return, counter contains the first number, that was not used,
/// that is \# Unknowns+stride.
/// A more user friendly interface is provided by IdxDescCL::CreateNumbOnInterface.
void CreateNumbOnInterfaceP2 (const Uint idx, IdxT& counter, Uint stride,
        const MultiGridCL::TriangVertexIteratorCL& vbegin,
        const MultiGridCL::TriangVertexIteratorCL& vend,
        const MultiGridCL::TriangEdgeIteratorCL& ebegin,
        const MultiGridCL::TriangEdgeIteratorCL& eend,
        const MultiGridCL::TriangTetraIteratorCL& begin,
        const MultiGridCL::TriangTetraIteratorCL& end,
        const VecDescCL& ls, const BndDataCL<>& lsetbnd, double omit_bound= -1./*default to using all dof*/)
{
    //const size_t stride= 1;

    LocalP2CL<> p2[10];
    for (int i= 0; i < 10; ++i)
        p2[i][i]= 1.;
    // first set NoIdx in all vertices
    for (MultiGridCL::TriangVertexIteratorCL vit= vbegin; vit != vend; ++vit) {
        vit->Unknowns.Prepare(idx);
        vit->Unknowns.Invalidate(idx);
    }
    for (MultiGridCL::TriangEdgeIteratorCL   eit= ebegin; eit != eend; ++eit) {
        eit->Unknowns.Prepare(idx);
        eit->Unknowns.Invalidate(idx);
    }
    // then create numbering of vertices at the interface
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( 1);
    std::valarray<double> ls_loc( lat.vertex_size());
    LocalP2CL<> locp2_ls;
    SPatchCL<3> patch;
    QuadDomainCodim1CL<3> qdom;
    std::valarray<double> qp2; // P2-shape-function as integrand
    for (MultiGridCL::TriangTetraIteratorCL it= begin; it != end; ++it) {
        locp2_ls.assign( *it, ls, lsetbnd);
        evaluate_on_vertexes( locp2_ls, lat, Addr( ls_loc));
        if (equal_signs( ls_loc))
            continue;

        const double limit= omit_bound < 0. ? 0. : omit_bound*std::pow( it->GetVolume()*6, 4./3.);
        patch.make_patch<MergeCutPolicyCL>( lat, ls_loc);
        if (patch.empty())
            continue;

        make_CompositeQuad5Domain2D( qdom, patch, *it);
        qp2.resize( qdom.vertex_size());
        for (Uint i= 0; i < 10; ++i) {
            UnknownHandleCL& unknowns= i < 4 ? const_cast<VertexCL*>( it->GetVertex( i))->Unknowns
                                             : const_cast<EdgeCL*>( it->GetEdge( i - 4))->Unknowns;
            if (unknowns.Exist( idx))
                continue;

            evaluate_on_vertexes( p2[i], qdom, Addr( qp2));
            if (quad_codim1( qp2*qp2, qdom) > limit) {
                unknowns( idx)= counter;
                counter+= stride;
            }
        }
    }
}

void IdxDescCL::CreateNumbOnInterface(Uint level, MultiGridCL& mg, const VecDescCL& ls,
                                      const BndDataCL<>& lsetbnd, double omit_bound)
/// Uses CreateNumbOnInterfaceVertex on the triangulation with level \p level on the multigrid \p mg.
/// One can create P1-and P2-elements.
{
    // set up the index description
    const Uint idxnum= GetIdx();
    TriangLevel_= level;
    NumUnknowns_= 0;

    // allocate space for indices; number unknowns in TriangLevel level
    if ((GetFE() == P1IF_FE) || (GetFE() == vecP1IF_FE))
        CreateNumbOnInterfaceVertex( idxnum, NumUnknowns_, NumUnknownsVertex(),
            mg.GetTriangVertexBegin(level), mg.GetTriangVertexEnd(level),
            mg.GetTriangTetraBegin( level), mg.GetTriangTetraEnd( level), ls, lsetbnd, omit_bound);
    else if ((GetFE() == P2IF_FE) || (GetFE() == vecP2IF_FE))
        CreateNumbOnInterfaceP2( idxnum, NumUnknowns_, NumUnknownsVertex(),
            mg.GetTriangVertexBegin(level), mg.GetTriangVertexEnd(level),
            mg.GetTriangEdgeBegin(level),   mg.GetTriangEdgeEnd(level),
            mg.GetTriangTetraBegin( level), mg.GetTriangTetraEnd( level),
            ls, lsetbnd, omit_bound);

    if (NumUnknownsFace() != 0 || NumUnknownsTetra() != 0)
        throw DROPSErrCL( "CreateNumbOnInterface: Only vertex and edge unknowns are implemented.\n" );
}

    void IdxDescCL::CreateNumbNearInterface(Uint level, MultiGridCL& mg, const VecDescCL& ls,
                                            const BndDataCL<>& lsetbnd, double dist,double omit_bound)
/// Uses CreateNumbNearInterfaceVertex on the triangulation with level \p level on the multigrid \p mg.
/// One can only create P1-elements.
    {
        // set up the index description
        const Uint idxnum= GetIdx();
        TriangLevel_= level;
        NumUnknowns_= 0;

        // allocate space for indices; number unknowns in TriangLevel level
        if (NumUnknownsVertex() != 0)
            CreateNumbNearInterfaceVertex( idxnum, NumUnknowns_, NumUnknownsVertex(),
                                           mg.GetTriangVertexBegin(level), mg.GetTriangVertexEnd(level),
                                           mg.GetTriangTetraBegin( level), mg.GetTriangTetraEnd( level), ls, lsetbnd, dist,omit_bound);

        if (NumUnknownsEdge() != 0 || NumUnknownsFace() != 0 || NumUnknownsTetra() != 0)
            throw DROPSErrCL( "CreateNumbOnInterface: Only vertex unknowns are implemented\n" );
    }

void IdxDescCL::CreateNumbStdFE( Uint level, MultiGridCL& mg)
// numbering of standard FE
{
    const Uint idxnum= GetIdx();
    TriangLevel_= level;
    NumUnknowns_ = 0;
    match_fun match= mg.GetBnd().GetMatchFun();

    // allocate space for indices; number unknowns in TriangLevel level
    if (match)
    {
        if (NumUnknownsVertex())
            CreatePeriodicNumbOnSimplex( idxnum, NumUnknowns_, NumUnknownsVertex(),
                mg.GetTriangVertexBegin(level), mg.GetTriangVertexEnd(level), Bnd_, level, match);
        if (NumUnknownsEdge())
            CreatePeriodicNumbOnSimplex( idxnum, NumUnknowns_, NumUnknownsEdge(),
                mg.GetTriangEdgeBegin(level), mg.GetTriangEdgeEnd(level), Bnd_, level, match);
        if (NumUnknownsFace())
            CreatePeriodicNumbOnSimplex( idxnum, NumUnknowns_, NumUnknownsFace(),
                mg.GetTriangFaceBegin(level), mg.GetTriangFaceEnd(level), Bnd_, level, match);
        if (NumUnknownsTetra())
            CreateNumbOnTetra( idxnum, NumUnknowns_, NumUnknownsTetra(),
                mg.GetTriangTetraBegin(level), mg.GetTriangTetraEnd(level), level);
    }
    else
    {
        if (NumUnknownsVertex())
            CreateNumbOnSimplex( idxnum, NumUnknowns_, NumUnknownsVertex(),
                mg.GetTriangVertexBegin(level), mg.GetTriangVertexEnd(level), Bnd_, level);
        if (NumUnknownsEdge())
            CreateNumbOnSimplex( idxnum, NumUnknowns_, NumUnknownsEdge(),
                mg.GetTriangEdgeBegin(level), mg.GetTriangEdgeEnd(level), Bnd_, level);
        if (NumUnknownsFace())
            CreateNumbOnSimplex( idxnum, NumUnknowns_, NumUnknownsFace(),
                mg.GetTriangFaceBegin(level), mg.GetTriangFaceEnd(level), Bnd_, level);
        if (NumUnknownsTetra())
            CreateNumbOnTetra( idxnum, NumUnknowns_, NumUnknownsTetra(),
                mg.GetTriangTetraBegin(level), mg.GetTriangTetraEnd(level), level);
    }
}

void IdxDescCL::CreateNumbering( Uint level, MultiGridCL& mg, const VecDescCL* lsetp, const BndDataCL<>* lsetbnd)
/// Memory for the Unknown-Indices on TriangLevel level is allocated
/// and the unknowns are numbered.
/// If a matching function is specified, numbering on periodic boundaries
/// is performed, too.
/// After that the extended DoFs are numbered for extended FE.
{
    if (IsOnInterface())
    {
#ifdef _PAR
        throw DROPSErrCL("IdxDescCL::CreateNumbering: Check first, if numbering on interface works in parDROPS.");
#endif
        if (lsetp == 0) throw DROPSErrCL("IdxDescCL::CreateNumbering: no level set function for interface numbering given");
        CreateNumbOnInterface( level, mg, *lsetp, *lsetbnd, GetXidx().GetBound());
    }
    else {
        CreateNumbStdFE( level, mg);
        if (IsExtended() && !IsExtendedSpaceTime()) {
            if (lsetp == 0){
                std::cout << "I am extended" << std::endl;
                throw DROPSErrCL("IdxDescCL::CreateNumbering: no level set function for XFEM numbering given");
            }
            NumUnknowns_= extIdx_.UpdateXNumbering( this, mg, *lsetp, *lsetbnd, true);
        }
        else if (IsExtendedSpaceTime()){
            if (lsetp == 0) throw DROPSErrCL("IdxDescCL::CreateNumbering: no level set function for XFEM numbering given");
            //_NumUnknowns_= extIdx_.UpdateSTXNumbering( this, mg, *lsetp, *lsetp_new, *lsetbnd, true);
            // throw DROPSErrCL("Not yet... ");
            std::cout << " there should be a STXNumbering coming ... ";
        }
    }
#ifdef _PAR
    if (!IsExtendedSpaceTime())
        ex_->CreateList(mg, this, true, true);
#endif
}


void IdxDescCL::UpdateXNumbering( MultiGridCL& mg, const VecDescCL& lset, const BndDataCL<>& lsetbnd)
{
    if (IsExtended() && !IsExtendedSpaceTime()) {
        NumUnknowns_= extIdx_.UpdateXNumbering( this, mg, lset, lsetbnd, false);
#ifdef _PAR
        ex_->CreateList(mg, this, true, true);
#endif
    }
}

    void IdxDescCL::CreateNumbering( Uint level, MultiGridCL& mg, const VecDescCL* lsetp, const BndDataCL<>* lsetbnd, double width)
/// Memory for the Unknown-Indices on TriangLevel level is allocated
/// and the unknowns are numbered.
/// If a matching function is specified, numbering on periodic boundaries
/// is performed, too.width
/// After that the extended DoFs are numbered for extended FE.
    {

        if (width>0) {
            if (lsetp == 0)
                throw DROPSErrCL("IdxDescCL::CreateNumbering: no level set function for interface numbering given");

            CreateNumbNearInterface(level, mg, *lsetp, *lsetbnd, width, GetXidx().GetBound());
        }
        else {

            if (IsNearInterface()) {
#ifdef _PAR
                throw DROPSErrCL("IdxDescCL::CreateNumbering: Check first, if numbering on interface works in parDROPS.");
#endif
                if (lsetp == 0)
                    throw DROPSErrCL("IdxDescCL::CreateNumbering: no level set function for interface numbering given");

                CreateNumbNearInterface(level, mg, *lsetp, *lsetbnd, width, GetXidx().GetBound());
            } else if (IsOnInterface()) {
#ifdef _PAR
                throw DROPSErrCL("IdxDescCL::CreateNumbering: Check first, if numbering on interface works in parDROPS.");
#endif
                if (lsetp == 0)
                    throw DROPSErrCL("IdxDescCL::CreateNumbering: no level set function for interface numbering given");
                CreateNumbOnInterface(level, mg, *lsetp, *lsetbnd, GetXidx().GetBound());
            } else {
                CreateNumbStdFE(level, mg);
                if (IsExtended()) {
                    if (lsetp == 0) {
                        std::cout << "I am extended" << std::endl;
                        throw DROPSErrCL("IdxDescCL::CreateNumbering: no level set function for XFEM numbering given");
                    }
                    NumUnknowns_ = extIdx_.UpdateXNumbering(this, mg, *lsetp, *lsetbnd, true);
                }
            }
        }
#ifdef _PAR
        ex_->CreateList(mg, this, true, true);
#endif
    }




void IdxDescCL::CreateNumbering( Uint level, MultiGridCL& mg, const BndCondCL& Bnd,
    const VecDescCL* lsetp, const BndDataCL<>* lsetbnd)
{
    CreateNumbering( level, mg, Bnd, Bnd, lsetp, lsetbnd);
}

void IdxDescCL::CreateNumbering( Uint level, MultiGridCL& mg, const BndCondCL& Bnd, const BndCondCL& Bnd_aux,
    const VecDescCL* lsetp, const BndDataCL<>* lsetbnd)
{
    Bnd_= Bnd;
    Bnd_aux_= Bnd_aux;
    CreateNumbering( level, mg, lsetp, lsetbnd);
}


void IdxDescCL::CreateNumbering( Uint level, MultiGridCL& mg, const IdxDescCL& baseIdx,
    const VecDescCL* lsetp, const BndDataCL<>* lsetbnd)
{
    Bnd_= baseIdx.Bnd_;
    CreateNumbering( level, mg, lsetp, lsetbnd);
}

void IdxDescCL::DeleteNumbering(MultiGridCL& MG)
/// This routine writes NoIdx as unknown-index for all indices of the
/// given index-description. NumUnknowns will be set to zero.
{
    const Uint idxnum = GetIdx();    // idx is the index in UnknownIdxCL
    const Uint level  = -1;          // last level
    NumUnknowns_ = 0;

    // delete memory allocated for indices
    if (NumUnknownsVertex())
        DeleteNumbOnSimplex( idxnum, MG.GetAllVertexBegin(level), MG.GetAllVertexEnd(level) );
    if (NumUnknownsEdge())
        DeleteNumbOnSimplex( idxnum, MG.GetAllEdgeBegin(level), MG.GetAllEdgeEnd(level) );
    if (NumUnknownsFace())
        DeleteNumbOnSimplex( idxnum, MG.GetAllFaceBegin(level), MG.GetAllFaceEnd(level) );
    if (NumUnknownsTetra())
        DeleteNumbOnSimplex( idxnum, MG.GetAllTetraBegin(level), MG.GetAllTetraEnd(level) );
    extIdx_.DeleteXNumbering();
#ifdef _PAR
    ex_->clear();
#endif
}

IdxT ExtIdxDescCL::UpdateXNumbering( IdxDescCL* Idx, const MultiGridCL& mg, const VecDescCL& lset, const BndDataCL<>& lsetbnd, bool NumberingChanged)
{
    const Uint sysnum= Idx->GetIdx(),
        level= Idx->TriangLevel(),
        stride= Idx->IsScalar() ? 1: 3; // scalar or vector-valued FE
    IdxT extIdx= NumberingChanged ? Idx->NumUnknowns() : Xidx_.size();
    Xidx_old_.assign( extIdx, NoIdx);
    Xidx_.swap( Xidx_old_);
    lset_= &lset;
    lsetbnd_ = &lsetbnd;
    InterfaceTetraCL cut;
    LocalP2CL<> hat_sq[4]; // values of phi_i*phi_i

    for (int i=0; i<4; ++i) //initialize hat_ii
    {
        hat_sq[i][i]=1.;
        for (int j=0; j<4; ++j)
            if (i!=j)
                hat_sq[i][EdgeByVert(i,j)+4]=0.25;
    }
    LocalP2CL<> locPhi;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, level, it)
    {
        const double h3= it->GetVolume()*6,
            h= cbrt( h3), h5= h*h*h3, // h^5
            limit= h5*omit_bound_;
        locPhi.assign( *it, lset, lsetbnd);
        cut.Init( *it, locPhi);
        SVectorCL<4> loc_int; // stores integrals \int_T p^2 dx, where p = p_i^\Gamma. Extended DoF is omitted
                              // iff this integral falls below a certain bound (omit_bound*h^5) for all tetrahedra T.
        if (cut.Intersects() )
        {

            for (int ch=0; ch<8; ++ch)
            {
                cut.ComputeCutForChild(ch);
                for (Uint i=0; i<4; ++i)
                    loc_int[i]+= cut.quad( hat_sq[i], h3, cut.GetSign(i) != 1); // integrate on other part; for sign == 0, the support is in the positive part.
            }
            if (cut.IntersectsInterior())
                for (Uint i=0; i<4; ++i)
                { // extend all DoFs
                    if (loc_int[i] < limit) continue; // omit DoFs of minor importance, which lead to unstable solutions
                    if (!it->GetVertex(i)->Unknowns.Exist( sysnum)) continue;
                    const IdxT nr= it->GetVertex(i)->Unknowns(sysnum);
                    if (Xidx_[nr]==NoIdx) // not extended yet
                        for (Uint k=0; k<stride; ++k)
                            Xidx_[nr+k]= extIdx++;
                }
            else
                for (Uint i=0; i<4; ++i)
                { // extend only DoFs on interface
                    if (cut.GetSign(i)!=0) continue;
                    if (loc_int[i] < limit) continue; // omit DoFs of minor importance, which lead to unstable solutions
                    if (!it->GetVertex(i)->Unknowns.Exist( sysnum)) continue;
                    const IdxT nr= it->GetVertex(i)->Unknowns(sysnum);
                    if (Xidx_[nr]==NoIdx) // not extended yet
                        for (Uint k=0; k<stride; ++k)
                            Xidx_[nr+k]= extIdx++;
                }
        }
    }
#ifdef _PAR
    // communicate extended dofs on vertices
    CommunicateXFEMNumbCL comm( Idx);
    comm.Call();

    // number all extended dofs from other procs (where extended dof is flagged by NoIdx-1)
    for (size_t i=0; i<Xidx_.size(); ++i)
        if (Xidx_[i] == NoIdx-1)
            Xidx_[i]= extIdx++;
#endif
    return extIdx;
}

#ifdef _PAR
bool ExtIdxDescCL::CommunicateXFEMNumbCL::Gather( const DiST::TransferableCL& t, DiST::SendStreamCL& s)
{
    VertexCL* sp = 0;
    simplex_cast( t, sp);

    if (sp->Unknowns.Exist(current_Idx_->GetIdx()))
        s << current_Idx_->IsExtended( sp->Unknowns(current_Idx_->GetIdx()));
    else
        s << false;
    return true;
}

bool ExtIdxDescCL::CommunicateXFEMNumbCL::Scatter( DiST::TransferableCL& t, const size_t numData, DiST::MPIistreamCL& r)
{
    VertexCL* sp= 0;
    simplex_cast( t, sp);

    bool RemoteExtended = false, tmp;
    for (size_t i = 0; i< numData; ++i) {
        r >> tmp;
        RemoteExtended = RemoteExtended || tmp;
    }
    if (!sp->Unknowns.Exist(current_Idx_->GetIdx()))
        return true;
    const IdxT dof= sp->Unknowns(current_Idx_->GetIdx());

    if (!current_Idx_->IsExtended( dof) && RemoteExtended)
        for (Uint i = 0; i < current_Idx_->NumUnknownsVertex(); ++i)
            current_Idx_->GetXidx()[dof+i]= NoIdx-1;
    return true;
}

void ExtIdxDescCL::CommunicateXFEMNumbCL::Call()
{
    DiST::InterfaceCL::DimListT dimlist; dimlist.push_back( 0);
    DiST::PrioListT Prios; Prios.push_back(PrioMaster);
    const Uint max_lvl= current_Idx_->TriangLevel();
    DiST::LevelListCL Levels(max_lvl);

    DiST::InterfaceCL comm( Levels, Prios, Prios, dimlist);
    comm.Communicate( *this);
}
#endif

void ExtIdxDescCL::Old2New(VecDescCL* v)
{
    VectorCL tmp( v->Data);
    // set v to zero vector with appropriate length
    v->SetIdx( v->RowIdx);
    // copy standard FE part
    const IdxT nStd= GetNumUnknownsStdFE();
#if defined(__SUNPRO_CC) || defined(DROPS_WIN)
    for (size_t i=0; i<nStd; ++i)
        v->Data[i] = tmp[i];
#else
    v->Data[std::slice( 0, nStd, 1)]= tmp[std::slice( 0, nStd, 1)];
#endif

    if (Xidx_.size() != Xidx_old_.size()) { // standard FE index changed (e.g., grid changed)
#ifndef _PAR
          std::cout << "ExtIdxDescCL::Old2New: Xidx: " << Xidx_.size()
                    << "\tXidx_old: " << Xidx_old_.size()
                    << "Extended Unknowns set to 0.\n";
#endif
        return;
    }

    // treat extended part
    IdxT ni= 0, di=0, ri= 0, ci= 0;
    for (size_t i= 0; i < Xidx_.size(); ++i) {
        if ( Xidx_old_[i] == NoIdx && Xidx_[i] != NoIdx) ++ni;
        if ( Xidx_old_[i] != NoIdx && Xidx_[i] == NoIdx) ++di;
        if ( Xidx_old_[i] != NoIdx && Xidx_[i] != NoIdx) { // extended dof was also extended before
            if ( Xidx_old_[i] != Xidx_[i])
                ++ri;
            v->Data[Xidx_[i]]= tmp[Xidx_old_[i]];
            ++ci;
        }
    }

#ifndef _PAR
    std::cout << "ExtIdxDescCL::Old2New: #P1-unknowns: " <<Xidx_.size()
              << "\t#new dof: " << ni
              << "\t#deleted dof: " << di
              << "\t#renumbered dof: " << ri
              << "\t#copied extended-dof: " << ci
              << '\n';
#endif
}

void permute_fe_basis (MultiGridCL& mg, IdxDescCL& idx, const PermutationT& p)
{
    const Uint sys= idx.GetIdx();
    const Uint lvl= idx.TriangLevel();
    const Uint num_components= idx.NumUnknownsVertex();

   if (idx.IsExtended())
        permute_fe_basis_extended_part( idx.GetXidx(), p, num_components);

    switch (idx.GetFE()) {
      case P0_FE:
        DROPS_FOR_TRIANG_TETRA( mg, lvl, it)
            if (it->Unknowns.Exist( sys))
                it->Unknowns( sys)= num_components*p[it->Unknowns( sys)/num_components];
        break;

      case P1_FE:  // fall through
      case P1X_FE: // fall through
      case P1IF_FE:
        DROPS_FOR_TRIANG_VERTEX( mg, lvl, it)
            if (it->Unknowns.Exist( sys) && it->Unknowns( sys) != NoIdx)
                it->Unknowns( sys)= num_components*p[it->Unknowns( sys)/num_components];
        break;

      case vecP1Bubble_FE: // fall through
      case P1Bubble_FE:
        DROPS_FOR_TRIANG_VERTEX( mg, lvl, it)
            if (it->Unknowns.Exist( sys))
                it->Unknowns( sys)= num_components*p[it->Unknowns( sys)/num_components];
        DROPS_FOR_TRIANG_TETRA( mg, lvl, it)
            if (it->Unknowns.Exist( sys))
                it->Unknowns( sys)= num_components*p[it->Unknowns( sys)/num_components];
        break;

      case P1D_FE:
        DROPS_FOR_TRIANG_FACE( mg, lvl, it)
            if (it->Unknowns.Exist( sys))
                it->Unknowns( sys)= num_components*p[it->Unknowns( sys)/num_components];
        break;

      case vecP2R_FE: // fall through
      case vecP2_FE:  // fall through
      case P2R_FE:    // fall through
      case P2_FE:
        DROPS_FOR_TRIANG_VERTEX( mg, lvl, it)
            if (it->Unknowns.Exist( sys))
                it->Unknowns( sys)= num_components*p[it->Unknowns( sys)/num_components];
        DROPS_FOR_TRIANG_EDGE( mg, lvl, it)
            if (it->Unknowns.Exist( sys))
                it->Unknowns( sys)= num_components*p[it->Unknowns( sys)/num_components];
        break;
      default: throw DROPSErrCL("permute_fe_basis: unknown FE type\n");
    }
}

void
LocalNumbP1CL::assign_indices_only (const TetraCL& s, const IdxDescCL& idx)
{
    const Uint sys= idx.GetIdx();
    for (Uint i= 0; i < 4; ++i)
        num[i]= s.GetVertex( i)->Unknowns.Exist( sys) ? s.GetVertex( i)->Unknowns( sys) : NoIdx;
}

LocalNumbP1CL::LocalNumbP1CL(const TetraCL& s, const IdxDescCL& idx)
/// \param s The tet, from which index-numbers are read.
/// \param idx The IdxDescCL-object to be used.
{
    this->assign_indices_only( s, idx);
}


void
LocalNumbP2CL::assign_indices_only (const TetraCL& s, const IdxDescCL& idx)
{
    const Uint sys= idx.GetIdx();
    if (!idx.IsDG())
    {
        for (Uint i= 0; i < 4; ++i)
            num[i]= s.GetVertex( i)->Unknowns.Exist( sys) ? s.GetVertex( i)->Unknowns( sys) : NoIdx;
        for(Uint i= 0; i < 6; ++i)
            num[i+4]= s.GetEdge( i)->Unknowns.Exist( sys) ? s.GetEdge( i)->Unknowns( sys)   : NoIdx;
    }
    else
    {
        Uint first = s.Unknowns(sys);
        for (int i = 0; i < 10; ++i)
            num[i] = first++;
    }
}

LocalNumbP2CL::LocalNumbP2CL(const TetraCL& s, const IdxDescCL& idx)
/// \param s The tet, from which index-numbers are read.
/// \param idx The IdxDescCL-object to be used.
{
    this->assign_indices_only( s, idx);
}

} // end of namespace DROPS
