/// \file volume_adjustment.cpp
/// \brief Maintain and adjust the volume of the subdomains with positive/negative level set value.
/// \author LNM RWTH Aachen: Joerg Grande

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
 * Copyright 2017 LNM RWTH Aachen, Germany
*/

#include "levelset/volume_adjustment.h"
#include "levelset/levelset.h"
#include "misc/params.h"
#include "num/discretize.h"
#include "num/lattice-eval.h"
#include "num/quadrature.h"

namespace DROPS
{

std::unique_ptr<VolumeAdjustmentCL> VolumeAdjustmentCL::Create (LevelsetP2CL* lset, const ParamCL& P)
{
    switch (P.get<int>("VolCorrection"))
    {
        case 0: return std::unique_ptr<VolumeAdjustmentCL> (new VolumeAdjustmentCL (lset));
        case 1: return std::unique_ptr<VolumeAdjustmentCL> (new GlobalVolumeAdjustmentCL (lset));
        case 2: return std::unique_ptr<VolumeAdjustmentCL> (new ComponentBasedVolumeAdjustmentCL (lset));
        default: throw DROPSErrCL("VolumeAdjustmentCL::Create: This case is not covered.\n");
    }
}

void VolumeAdjustmentCL::DebugOutput (std::ostream& os) const
{
    os << "VolumeAdjustmentCL: No volume adjustment.\n" << std::endl;
}

template <class VolumeFunT>
double compute_volume_correction (VolumeFunT f, double refvol, double reltol)
{
    const double tol= reltol*refvol;

    double v0=f(0.) - refvol;
    if (std::abs(v0)<=tol)
        return 0.;

    // hint: surf(Kugel) = [3/4/pi*vol(Kugel)]^(2/3) * 4pi
    double d0=0,
           d1=v0*0.23/std::pow (refvol, 2./3.);
    double v1=f(d1) - refvol;
    if (std::abs(v1)<=tol)
        return d1;

    // secant method for initial value
    while (v1*v0 > 0) // same sign
    {
        const double d2=d1-1.2*v1*(d1-d0)/(v1-v0);
        d0=d1; d1=d2; v0=v1; v1=f(d1) - refvol;
        if (std::abs(v1)<=tol)
            return d1;
    }

    // Anderson-Bjoerk for exact value
    while (true)
    {
        const double d2=(v1*d0-v0*d1)/(v1-v0),
                     v2=f(d2) - refvol;
        if (std::abs(v2)<=tol)
            return d2;

        if (v2*v1 < 0) // different signs
          { d0=d1; d1=d2; v0=v1; v1=v2; }
        else
          { const double c=1.0-v2/v1; d1=d2; v1=v2; v0*= c>0 ? c : 0.5; }
    }
    return 0.; // unreachable
}

void GlobalVolumeAdjustmentCL::AdjustVolume ()
{
    if (global_reference_volume_ == -1.0)
        throw DROPSErrCL ("GlobalVolumeAdjustmentCL::AdjustVolume: The global_reference_volume_ has not been set, but it is required for this method (maybe you forgot SetGlobalReferenceVolume somewhere).\n");

    dphi_= compute_volume_correction (
        [this](double x)->double { return lset_->GetVolume(x, num_subdivision_); },
        global_reference_volume_,
        tol_
    );
    lset_->Phi.Data+= dphi_;
}

void GlobalVolumeAdjustmentCL::DebugOutput (std::ostream& os) const
{
    os << "GlobalVolumeAdjustmentCL: shift for level set function dphi_: " << dphi_ << ", "
       "new rel. volume: " << lset_->GetVolume()/global_reference_volume_ << "." << std::endl;
}


/// \brief Compute the strongly connected components of the directed graph
//of a matrix.
/// The square matrix M is interpreted as graph with vertices
//0..M.nom_cols()-1.
/// There is a directed edge (i,j), iff M_ij != 0. (The entries in row i
//are considered as successors of i).
///
/// The strongly connected components of the graph are computed with
//breadth-first search.
class GraphComponentsCL
{
  private:
   std::vector<size_t> component_;      ///< component[i] is the number of the component of i.
   std::vector<size_t> component_size_; ///< number of elements in component i.

 /// \brief Perform a breadth-first search to discover the connected component of v.
   template < class T>
     void visit_connected_component_of (size_t v, const SparseMatBaseCL<T>& M);

  public:
    /// \brief Enumerate the connected components of the graph of M.
    template <class T>
    void number_connected_components (const SparseMatBaseCL<T>& M);
    /// \brief Number of strongly connected components.
    size_t                     num_components () const { return component_size_.size(); }
    /// \brief component_map()[i] is the number of the component, to which vertex i belongs.
    const std::vector<size_t>& component_map  () const { return component_; }
    void clear() { component_.clear(); component_size_.clear(); }

    /// \brief component(c) contains all vertices in component c. The vertices are in ascending order.
    std::vector<size_t> component  (size_t c) const;
    /// \brief the number of vertices in component i.
    const std::vector<size_t>& component_size () const { return component_size_; }
    /// \brief High-level overview of the numbering process.
    inline void stats (std::ostream& os) const;
};

std::vector<size_t> GraphComponentsCL::component (size_t c) const
{
    std::vector<size_t> ret;
    ret.reserve (component_size()[c]);
    for (size_t i= 0; i < component_.size(); ++i)
        if (component_[i] == c)
            ret.push_back( i);

    return std::move (ret);
}

template <typename T>
  void
  GraphComponentsCL::number_connected_components (const SparseMatBaseCL<T>& M)
{
    const size_t num_verts= M.num_cols();
    component_.clear();
    component_.resize( num_verts, NoVert);
    component_size_.clear();

    for (size_t v= 0; v < num_verts; ++v)
        if (component_[v] == NoVert) {
            component_size_.push_back( 0); ///< Initialize a new component.
            visit_connected_component_of( v, M);
        }
}

template <typename T>
  void
  GraphComponentsCL::visit_connected_component_of (size_t v0, const SparseMatBaseCL<T>& M)
{
    const size_t cur_component= component_size_.size() - 1;
    size_t* cur_component_size= &component_size_.back();

    std::list<size_t> deque;
    deque.push_back( v0);
    while (!deque.empty()) {
        const size_t v= deque.front();
        deque.pop_front();
        component_[v]= cur_component;
        ++*cur_component_size;

        const double* succ_val= M.GetFirstVal( v);
        for (const size_t* succ= M.GetFirstCol( v), * rowend= M.GetFirstCol( v + 1); succ != rowend; ++succ, ++succ_val) {
            if (*succ_val == T()) // Edges in the graph correspond to non-zero entries. Row i contains the successors of vertex i.
                continue;
            if (component_[*succ] == NoVert) {
                component_[*succ]= cur_component;
                deque.push_back( *succ);
            }
        }
    }
}

//*****************************************************************************
//                               AdjacencyAccuCL
//*****************************************************************************

class AdjacencyAccuCL : public TetraAccumulatorCL
{
  private:
    const LevelsetP2CL& lset;
    MatrixCL& A;
    MatrixCL& B;
    SparseMatBuilderCL<double>* mA_;
    SparseMatBuilderCL<double>* mB_;

    LocalNumbP2CL n; ///< global numbering of the P2-unknowns
    SMatrixCL<10,10> LocAdjaMat;
    SMatrixCL<10,10> LocAdjaMatB;
    LocalP2CL<> ls_loc;

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and, if required, the Dirichlet-values needed to eliminate the boundary-dof from the global system.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system.
    void update_global_system ();

  public:
    AdjacencyAccuCL (const LevelsetP2CL& ls, MatrixCL& Au, MatrixCL& Bu): lset(ls), A(Au), B(Bu){};

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new AdjacencyAccuCL ( *this); };
};

void AdjacencyAccuCL::begin_accumulation ()
{
    std::cout << "entering AdjacencyAccuCL: ";
    const size_t num_unks= lset.idx.NumUnknowns();
    mA_= new SparseMatBuilderCL<double>( &A, num_unks, num_unks);
    mB_= new SparseMatBuilderCL<double>( &B, num_unks, num_unks);
}

void AdjacencyAccuCL::finalize_accumulation ()
{
    mA_->Build();
    delete mA_;
    mB_->Build();
    delete mB_;
#ifndef _PAR
    std::cout << A.num_nonzeros() << " nonzeros in A, ";
    std::cout << B.num_nonzeros() << " nonzeros in B, ";
#endif
    std::cout << '\n';
}

void AdjacencyAccuCL::visit (const TetraCL& tet)
{
    local_setup( tet);
    update_global_system();
}

void AdjacencyAccuCL::local_setup (const TetraCL& tet)
{
    ls_loc.assign( tet, lset.Phi, lset.GetBndData());
    for (Uint e= 6; e<18; ++e) {
        const Ubyte v0= VertOfEdge( e, 0),
                    v1= VertOfEdge( e, 1);
        const bool edge= (sign(ls_loc[v0]) + sign(ls_loc[v1]) >= 1) || (sign(ls_loc[v0])+sign(ls_loc[v1]) == -2);
        LocAdjaMat( v0, v1)= LocAdjaMat( v1, v0)= edge;
        LocAdjaMatB( v0, v1)= LocAdjaMatB( v1, v0)= 1;
    }
    for (Uint e= 30; e<43; ++e) {
        const Ubyte v0= VertOfEdge( e, 0),
                    v1= VertOfEdge( e, 1);
        const bool edge= (sign(ls_loc[v0])+sign(ls_loc[v1]) >= 1) || (sign(ls_loc[v0])+sign(ls_loc[v1]) == -2);
        LocAdjaMat( v0, v1)= LocAdjaMat( v1, v0)= edge;
        LocAdjaMatB( v0, v1)= LocAdjaMatB( v1, v0)=1;
    }
    n.assign( tet, lset.idx.GetFinest(), lset.GetBndData());
}

void AdjacencyAccuCL::update_global_system ()
{
    SparseMatBuilderCL<double>& mA= *mA_;
    SparseMatBuilderCL<double>& mB= *mB_;

    for(int i= 0; i < 10; ++i){    // assemble row Numb[i]
        if (n.WithUnknowns( i)) { // dof i is not on a Dirichlet boundary
            for(int j= 0; j < 10; ++j) {
                if (n.WithUnknowns( j)) { // dof j is not on a Dirichlet boundary
                    if (LocAdjaMat(i,j) != 0.)
                        mA( n.num[i], n.num[j])= LocAdjaMat(i,j);
                    if (LocAdjaMatB(i,j) != 0.)
                        mB( n.num[i], n.num[j])= LocAdjaMatB(i,j);
                }
            }
        }
    }
}

void SetupAdjacency (MatrixCL& A, MatrixCL& B, const LevelsetP2CL& lset)
/// Set up matrix A
{
    // TimerCL time;
    // time.Start();

    AdjacencyAccuCL accu(lset, A, B);
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    accumulate( accus, lset.GetMG(), lset.idx.TriangLevel(), lset.idx.GetMatchingFunction(), lset.idx.GetBndInfo());
    // time.Stop();
    // std::cout << "setup: " << time.GetTime() << " seconds" << std::endl;
}


//*****************************************************************************
//                               ComponentBasedVolumeAdjustmentCL
//*****************************************************************************

ComponentBasedVolumeAdjustmentCL::ComponentBasedVolumeAdjustmentCL(LevelsetP2CL* lset)
    : VolumeAdjustmentCL(lset) 
{}

ComponentBasedVolumeAdjustmentCL::~ComponentBasedVolumeAdjustmentCL()
{}


void ComponentBasedVolumeAdjustmentCL::compute_indicator_functions (const MatrixCL& A)
{
    doCorrection_= std::vector<bool> (num_components(), true);

    const component_vector tmp (ExtendOneStep (A, component_of_dof_, doCorrection_));
    const component_vector tmp2 (ExtendOneStep (A, tmp, doCorrection_));

    indicator_functions_= std::vector<std::valarray<double>> (num_components(), std::valarray<double>(component_of_dof_.size()));
    // For every point in tmp2, check which component is present and write into the corresponding characteristic function.
    for(Uint j= 0; j < component_of_dof_.size(); ++j)
        indicator_functions_[tmp2[j]][j]= 1.;
}

auto ComponentBasedVolumeAdjustmentCL::ExtendOneStep(const MatrixCL& A,
    const component_vector& cp,
    std::vector<bool>& doCorrection) const -> component_vector
{
    component_vector ret (cp);
    for (size_t v0= 0; v0 < A.num_rows(); ++v0) {
        if (sign_of_component_[cp[v0]] == 1)
            continue;
        for (size_t j= A.row_beg (v0); j != A.row_beg(v0 + 1); ++j) {
            const size_t v1= A.col_ind (j);
            if (sign_of_component_[cp[v1]] == 1)
                ret[v1]= cp[v0];
            if (sign_of_component_[cp[v1]] == -1 && cp[v1] != cp[v0]) {
               doCorrection[cp[v0]]= false;
               doCorrection[cp[v1]]= false;
            }
        }
    }
    return std::move (ret);
}

void ComponentBasedVolumeAdjustmentCL::FindComponents ()
{
    MatrixCL CompAdja;
    MatrixCL MeshAdja;
    SetupAdjacency (CompAdja, MeshAdja, *lset_);
    GraphComponentsCL Split;
    Split.number_connected_components(CompAdja);

    component_of_dof_= Split.component_map();
    Volumes.resize (Split.num_components()); // neccessary to make num_components() return the current number of components.
    for (Uint c= 0; c < num_components(); ++c)
        Volumes[c]= CalculateVolume(c, 0.);

    sign_of_component_.resize (num_components());
    ComputeReferencePoints();
    compute_indicator_functions (MeshAdja);
}

void ComponentBasedVolumeAdjustmentCL::init_coord_of_dof ()
{
    coord_of_dof_.resize (lset_->Phi.Data.size());

    const Uint lvl= lset_->Phi.GetLevel(),
               idx= lset_->Phi.RowIdx->GetIdx();
    const auto& MG= lset_->GetMG();

    DROPS_FOR_TRIANG_CONST_VERTEX (MG, lvl, it)
        if (it->Unknowns.Exist (idx))
            coord_of_dof_[it->Unknowns (idx)]= it->GetCoord();

    DROPS_FOR_TRIANG_CONST_EDGE (MG, lvl, it)
        if (it->Unknowns.Exist (idx))
            coord_of_dof_[it->Unknowns (idx)]= GetBaryCenter( *it);
}


void ComponentBasedVolumeAdjustmentCL::InitVolume_impl ()
{
    init_coord_of_dof();
    coord_of_dof_backup_= coord_of_dof_;
    FindComponents();

    // Set initial target volumes.
    targetVolumes= Volumes;

    make_backup();

    DebugOutput (std::cout);
}

void ComponentBasedVolumeAdjustmentCL::Repair()
{
    coord_of_dof_backup_= std::move (coord_of_dof_);
    init_coord_of_dof();
    sign_of_component_backup_= std::move (sign_of_component_);
    FindComponents();
    MatchComponents ();
    coord_of_dof_backup_= coord_of_dof_;
    make_backup();

    DebugOutput (std::cout);
}

void ComponentBasedVolumeAdjustmentCL::DebugOutput (std::ostream& os) const
{
    os << "ComponentBasedVolumeAdjustmentCL::Number of components " << num_components() << "\nVolumes ";
    seq_out(std::begin(Volumes),std::end(Volumes),os,", ");
    os << std::endl << "ReferencePoints\n";
    seq_out(std::begin(ReferencePoints),std::end(ReferencePoints),os,", ");
    os << std::endl;
}


void ComponentBasedVolumeAdjustmentCL::ComputeReferencePoints()
{
    ReferencePoints= std::vector<Point3DCL> (num_components());

    std::vector<double> tmp (num_components(), std::numeric_limits<double>::min());
    const VectorCL& ls= lset_->Phi.Data;
    for (size_t i= 0; i < component_of_dof_.size(); ++i) {
        const Uint c= component_of_dof_[i];
        if (std::abs (ls[i]) > tmp[c]) {
            ReferencePoints[c] = coord_of_dof_[i];
            sign_of_component_[c] = ls[i] > 0. ? 1 : -1;
            tmp[c] = std::abs (ls[i]);
        }
    }
}

double ComponentBasedVolumeAdjustmentCL::CalculateVolume(Uint c, double shift) const
{
    double ret= 0.;
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance (2);
    std::valarray<double> ls_values (lat.vertex_size());
    QuadDomainCL qdom;
    LocalP2CL<> loc_phi;
    TetraPartitionCL partition;
    LocalNumbP2CL n;

    VecDescCL Copy(&lset_->idx);
    Copy.Data= lset_->Phi.Data;
    if (c > 0 && shift != 0.)
        Copy.Data+= shift*indicator_functions_[c];

    DROPS_FOR_TRIANG_TETRA( lset_->GetMG(), lset_->idx.TriangLevel(), it) {
        n.assign_indices_only(*it, lset_->idx.GetFinest());
        loc_phi.assign(*it,Copy,lset_->GetBndData());

        bool comp_exists = false;
        for (Uint a= 0; a < 10; a++)
            if (n.WithUnknowns(a)) {
                if (component_of_dof_[n.num[a]] == c)
                    comp_exists = true;
                else if (c > 0 && component_of_dof_[n.num[a]] > 0)   // by definition, component 0 is the surrounding liquid
                    loc_phi[a] = 1.0;                   // remove negative components which are not the considered component c (this change does not alter the position of the boundary of component c
        }
        if (!comp_exists)
            continue;
        evaluate_on_vertexes (loc_phi, lat, Addr(ls_values));
        partition.make_partition< SortedVertexPolicyCL,MergeCutPolicyCL> (lat, ls_values);
        make_CompositeQuad3Domain (qdom, partition);
        DROPS::GridFunctionCL<> integrand (1., qdom.vertex_size());
        ret+= quad( integrand, it->GetVolume()*6., qdom, c > 0 ? NegTetraC : PosTetraC);
    }
    return ret;
}

auto ComponentBasedVolumeAdjustmentCL::component_of_point (const std::vector<Point3DCL>& refpts,
    const std::vector<int>& sign_of_component, const component_vector& component_of_dof, const std::vector<Point3DCL>& coord_of_dof, const std::vector<int>& sign_of_component2) const -> component_vector
{
    // Temporary arrays to store the minimal distances to the points in pts and the new component of the minimizers.
    std::vector<double> distances (refpts.size(), std::numeric_limits<double>::max());
    component_vector cnew (refpts.size(), -1); // cnew[i] is the component number of refpts[i] with respect to component_of_dof.

    for (size_t i= 0; i < component_of_dof.size(); ++i) {
        for (size_t c= 0; c < refpts.size(); ++c) {
            const double tmpdistance= (coord_of_dof[i] - refpts[c]).norm();
            const bool components_have_same_sign= sign_of_component[c] == sign_of_component2[component_of_dof[i]];
            if (components_have_same_sign && tmpdistance < distances[c]) {
                distances[c]= tmpdistance;
                cnew[c]= component_of_dof[i];
            }
        }
    }
    return std::move (cnew);
}

void ComponentBasedVolumeAdjustmentCL::MatchComponents ()
{
    // cold represents a permutation which maps the new component number (from ReferencePoints_) to the old.
    const component_vector cold (component_of_point (ReferencePoints, sign_of_component_, component_of_dof_backup_, coord_of_dof_backup_, sign_of_component_backup_));
    // cnew represents a permutation which maps the old component number (from ReferencePoints_backup) to the new.
    const component_vector cnew (component_of_point (ReferencePoints_backup, sign_of_component_backup_, component_of_dof_, coord_of_dof_, sign_of_component_));


    // M will be the adjacency matrix of the undirected graph G: The nodes of G are the (old and new) components. The edges are given by the components of the reference points (in both directions).
    MatrixCL M;
    const size_t nnew= num_components(),
                 nold= ReferencePoints_backup.size(),
                 n= nold + nnew;
    SparseMatBuilderCL<double> Mb (&M, n, n);
    for (Uint i= 0; i < nold; ++i)
        Mb (cnew[i] + nold, i)= Mb (i, cnew[i] + nold)= 1.;
    for (Uint i= 0; i < nnew; ++i)
        Mb (cold[i], i + nold)= Mb (i + nold, cold[i])= 1.;
    Mb.Build();
    std::cout << "ComponentBasedVolumeAdjustmentCL::MatchComponents: M: " << M << std::endl;

    GraphComponentsCL G;
    G.number_connected_components (M);

    // The new target volumes (for the new components) are computed on the connected components of G.
    std::vector<double> newtargetVolumes (num_components());
    for (size_t i= 0; i < G.num_components(); ++i) {
        const component_vector& component (G.component (i));
        double targetvolsum= 0.,
               volsum= 0;
        for (auto c: component) // All target volumes (of the old components) are summed up and and so are all current volumes (of the new components).
            if (c < nold)
                targetvolsum+= targetVolumes[c];
            else
                volsum+= Volumes[c - nold];
        for (auto c: component) // Each new component gets its fraction (within the new components) of the total (old) targetVolume.
            if (c >= nold)
                newtargetVolumes[c - nold]= Volumes[c - nold]/volsum * targetvolsum;
    }
    targetVolumes= newtargetVolumes;
}

void ComponentBasedVolumeAdjustmentCL::AdjustVolume()
{
    FindComponents();
    MatchComponents();

    // adapt Level Set
    for (Uint i= 0; i < num_components(); ++i) {
        std::cout << "Adjustment for component " << i << ": ";
        if (sign_of_component_[i] == 1) {
            std::cout << "Skipping positive component.\n";
            continue;
        }
        if (doCorrection_[i]) {
            const double s= compute_volume_correction ([this,i](double x)->double { return CalculateVolume(i, x); },
                                                       targetVolumes[i],
                                                       tol_);
            lset_->Phi.Data+=indicator_functions_[i]*s;
            std::cout << s << "\n";
        }
        else
            std::cout << "No adjustment, components are too close.\n";
    }
    ComputeReferencePoints();
    make_backup();
    FindComponents();
    MatchComponents();
    make_backup();
}

void ComponentBasedVolumeAdjustmentCL::make_backup()
{
    component_of_dof_backup_= component_of_dof_;
    sign_of_component_backup_= sign_of_component_;
    ReferencePoints_backup= ReferencePoints;
}

} // end of namespace DROPS
