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

class VTKComponentMapCL: public VTKVariableCL
{
  public:
    using component_vector = ComponentBasedVolumeAdjustmentCL::component_vector;

  private:
    const LevelsetP2CL* lset_ = nullptr;
    const component_vector* component_of_dof_ = nullptr;

  public:
    VTKComponentMapCL (std::string var_name)
        : VTKVariableCL (var_name) {}
    VTKComponentMapCL (const LevelsetP2CL* lset, const component_vector* component_of_dof, std::string var_name)
        : VTKVariableCL (var_name), lset_(lset), component_of_dof_(component_of_dof) {}

   /// \brief Called by VTKOutCL::Write().
    void put( VTKOutCL&) const;
    Uint GetDim() const { return 1; }
};

void VTKComponentMapCL::put( VTKOutCL& cf) const
{
    if (lset_ == nullptr)
        return;
    VecDescCL vd(const_cast<MLIdxDescCL*>(&lset_->idx));
    for (size_t i = 0; i < component_of_dof_->size(); ++i)
        vd.Data[i] = (*component_of_dof_)[i];
    cf.PutScalar( make_P2Eval( lset_->GetMG(), lset_->GetBndData(), vd), varName());
}


std::unique_ptr<VolumeAdjustmentCL> VolumeAdjustmentCL::Create (LevelsetP2CL* lset, const ParamCL& P)
{
    const char* methods[] = { "None", "Global", "ComponentBased" };
    std::string string_method;

    try {
        string_method = P.get<std::string>("VolCorrection");
    }
    catch (DROPSParamErrCL& ) {
        std::cerr << "VolumeAdjustmentCL::Create: Warning: The VolCorrection parameter expects a string (\"None\", \"Global\" or \"ComponentBased\"). Using no volume correction (\"None\") for this simulation.\n";
        return std::unique_ptr<VolumeAdjustmentCL> (new VolumeAdjustmentCL (lset));
    }
    if (string_method == "0" || string_method == "1" || string_method == "2")
        std::cerr << "VolumeAdjustmentCL::Create: Deprecated use of the VolCorrection parameter. The parameter expects a (non-empty) string (\"None\", \"Global\" or \"ComponentBased\").\n";
    if (string_method == "" || string_method == "0") {
        if (string_method == "")
            std::cerr << "Using no volume correction (\"None\") for this simulation.\n";
        return std::unique_ptr<VolumeAdjustmentCL> (new VolumeAdjustmentCL (lset));
    }
    else if (string_method == methods[0])
        return std::unique_ptr<VolumeAdjustmentCL> (new VolumeAdjustmentCL (lset));
    else if (string_method == "1" || string_method == methods[1])
        return std::unique_ptr<VolumeAdjustmentCL> (new GlobalVolumeAdjustmentCL (lset));
    else if (string_method == "2" || string_method == methods[2]) {
#ifdef _PAR
        throw DROPSErrCL("ComponentBasedVolumeAdjustmentCL is not yet implemented in parallel.\n");
#endif
        return std::unique_ptr<VolumeAdjustmentCL> (new ComponentBasedVolumeAdjustmentCL (lset));
    }
    else
        throw DROPSErrCL("VolumeAdjustmentCL::Create: Please specify which volume correction method you want to use. The VolCorrection parameter expects one of the following strings: \"None\", \"Global\" or \"ComponentBased\".\n");
}

void VolumeAdjustmentCL::DebugOutput (std::ostream& os) const
{
    os << "VolumeAdjustmentCL: No volume adjustment.\n" << std::endl;
}

VTKVariableCL& VolumeAdjustmentCL::make_VTKComponentMap(std::string var_name)
{
    return *new VTKComponentMapCL(var_name);
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

/// \brief Compute the strongly connected components of a directed graph.
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
   template <class GraphT>
     void visit_connected_component_of (size_t v, const GraphT& M);

  public:
    /// \brief Enumerate the connected components of the graph M.
    template <class GraphT>
    void number_connected_components (const GraphT& M);
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

template <typename GraphT>
  void
  GraphComponentsCL::number_connected_components (const GraphT& M)
{
    const size_t num_verts= M.num_verts();
    component_.clear();
    component_.resize( num_verts, NoVert);
    component_size_.clear();

    for (size_t v= 0; v < num_verts; ++v)
        if (component_[v] == NoVert) {
            component_size_.push_back( 0); ///< Initialize a new component.
            visit_connected_component_of( v, M);
        }
}

template <typename GraphT>
  void
  GraphComponentsCL::visit_connected_component_of (size_t v0, const GraphT& M)
{
    const size_t cur_component= component_size_.size() - 1;
    size_t* cur_component_size= &component_size_.back();

    std::vector<size_t> successors; // all succsessors of v, see below.
    std::list<size_t> deque;
    deque.push_back( v0);
    while (!deque.empty()) {
        const size_t v= deque.front();
        deque.pop_front();
        component_[v]= cur_component;
        ++*cur_component_size;
        successors= M.get_all_succsessors (v);
        for (size_t w: successors)
            if (component_[w] == NoVert) {
                component_[w]= cur_component;
                deque.push_back (w);
            }
    }
}

//*****************************************************************************
//                               SparseMatrixGraphCL
//*****************************************************************************
class SparseMatrixGraphCL
{
  public:
    using MatrixT= SparseMatBaseCL<unsigned char>;

  private:
    const MatrixT& M_;
    unsigned char edge_limit= 0; // M_.val() > edge_val indicates an edge between the vertices given by row and column index.

  public:
    SparseMatrixGraphCL (const MatrixT& M) : M_ (M) {}

    SparseMatrixGraphCL& set_edge_limit (unsigned char el) {
        edge_limit= el;
        return *this;
    }

    size_t num_verts () const  { return M_.num_cols(); }

    std::vector<size_t> get_all_succsessors (size_t v) const {
        std::vector<size_t> successors;
        successors.reserve (M_.row_beg (v + 1) - M_.row_beg (v));
        size_t w;
        for (size_t e= M_.row_beg (v); e < M_.row_beg (v + 1); ++e)
            if (M_.val( e) > edge_limit && (w= M_.col_ind( e)) != v)
                successors.push_back( w);
        return successors;
    }
};

//*****************************************************************************
//                               AdjacencyAccuCL
//*****************************************************************************
class AdjacencyAccuCL : public TetraAccumulatorCL
{
  private:
    static constexpr Uint edges_[25]= { // All edges of the regular children of the reference tetra, see VertOfEdgeAr.
         6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
        30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
        42
    };

    const LevelsetP2CL& lset;
    // M encodes the adjacency matrix of two graphs G_comp and G_mesh: a value > 0 is an edge in G_mesh, a value > 1 is an edge within a G_comp (i.e. the vertices of the edge have the same sign of the level set function).
    SparseMatBaseCL<unsigned char>& M;
    SparseMatBuilderCL<unsigned char>* mA_;

    LocalNumbP2CL n; ///< global numbering of the P2-unknowns
    unsigned char LocAdjaMat[10][10];
    LocalP2CL<> ls_loc;

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and, if required, the Dirichlet-values needed to eliminate the boundary-dof from the global system.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system.
    void update_global_system ();

  public:
    AdjacencyAccuCL (const LevelsetP2CL& ls, SparseMatBaseCL<unsigned char>& Marg) : lset (ls), M (Marg) {};

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new AdjacencyAccuCL (*this); };
};

constexpr Uint AdjacencyAccuCL::edges_[25];

void AdjacencyAccuCL::begin_accumulation ()
{
    std::cout << "entering AdjacencyAccuCL: ";
    const size_t num_unks= lset.idx.NumUnknowns();
    mA_= new SparseMatBuilderCL<unsigned char> (&M, num_unks, num_unks);
}

void AdjacencyAccuCL::finalize_accumulation ()
{
    mA_->Build();
    delete mA_;
#ifndef _PAR
//    std::cout << M.num_nonzeros() << " nonzeros in M, ";
#endif
//    std::cout << '\n';
}

void AdjacencyAccuCL::visit (const TetraCL& tet)
{
    local_setup( tet);
    update_global_system();
}

void AdjacencyAccuCL::local_setup (const TetraCL& tet)
{
    ls_loc.assign (tet, lset.Phi, lset.GetBndData());
    std::memset (LocAdjaMat, 0, sizeof (LocAdjaMat));
    n.assign( tet, lset.idx.GetFinest(), lset.GetBndData());
    for (Uint e: edges_) {
        const Ubyte v0= VertOfEdge (e, 0),
                    v1= VertOfEdge (e, 1);
        // The domain is partitioned into (\Omega_pos+\Gamma)  \disjoint_union  \Omega_neg
        const bool edge= (sign(ls_loc[v0]) + sign(ls_loc[v1]) >= 1) || (sign(ls_loc[v0])+sign(ls_loc[v1]) == -2) || (sign(ls_loc[v0])==0 && sign(ls_loc[v1])==0);
        LocAdjaMat[v0][v1]= LocAdjaMat[v1][v0]= 1 + edge;
    }
}

void AdjacencyAccuCL::update_global_system ()
{
    SparseMatBuilderCL<unsigned char>& mA= *mA_;

    for(int i= 0; i < 10; ++i) // assemble row Numb[i]
        if (n.WithUnknowns (i)) // dof i is not on a Dirichlet boundary
            for(int j= 0; j < 10; ++j)
                if (n.WithUnknowns (j) && LocAdjaMat[i][j] != 0)
                    mA (n.num[i], n.num[j])= LocAdjaMat[i][j];
}


//*****************************************************************************
//                               VolumeAccuCL
//*****************************************************************************

class VolumeAccuCL : public TetraAccumulatorCL
{
  private:
    const LevelsetP2CL& lset_;
    const VecDescCL& ls_vd_;
    size_t c_;
    std::vector<int> sign_of_component_;
    const ComponentBasedVolumeAdjustmentCL::component_vector* component_of_dof_ = nullptr;
    const std::vector<std::valarray<double>>* indicator_functions_ = nullptr;

    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance (2);
    std::valarray<double> ls_values= std::valarray<double> (lat.vertex_size());
    QuadDomainCL qdom;
    LocalP2CL<> loc_phi;
    TetraPartitionCL partition;
    LocalNumbP2CL n;

    std::vector<double> volumes_omp;
    double* thread_volume= 0;
    double volume;

  public:
    VolumeAccuCL (const LevelsetP2CL& lset, const VecDescCL& ls_vd, size_t c, const std::vector<int>& sign_of_component, const ComponentBasedVolumeAdjustmentCL::component_vector& component_of_dof)
        : lset_ (lset), ls_vd_ (ls_vd), c_ (c), sign_of_component_ (sign_of_component), component_of_dof_ (&component_of_dof) {}

    void use_indicator_functions(const std::vector<std::valarray<double>>& indicator_functions) {
        indicator_functions_ = &indicator_functions;
    }

    void begin_accumulation () {
//         std::cout << "entering VolumeAccuCL: ";
        volume= 0.;
        volumes_omp= std::vector<double> (omp_get_max_threads());
        thread_volume= &volumes_omp[0];
    }
    void finalize_accumulation() {
        volume= std::accumulate (volumes_omp.begin(), volumes_omp.end(), 0.);
//         std::cout << volume << ".\n";
    }
    void visit (const TetraCL& sit);
    TetraAccumulatorCL* clone (int tid) {
        VolumeAccuCL* cl= new VolumeAccuCL ( *this);
        cl->thread_volume= &volumes_omp[tid];
        return cl;
    }

    double get_volume () const { return volume; }
};

void VolumeAccuCL::visit (const TetraCL& tet)
{
    n.assign_indices_only (tet, lset_.idx.GetFinest());
    loc_phi.assign (tet, ls_vd_, lset_.GetBndData());
    // Shift in this case. Components have to be at least 4 P2 DOFs apart. Integrate over all tetrahedra in the two-neighbourhood (without any modifications).
    if (indicator_functions_ != nullptr) {
        double is= 0.;
        for (Uint i= 0; i < 10; ++i)
            if (n.WithUnknowns(i))
                is+= (*indicator_functions_)[c_][n.num[i]];
        if (is == 0.)
            return;
    }
    // No shift in this case. Components may be arbitrarily close, therefore one has to handle the case with multiple negative components on a single tetrahedron.
    else {
        bool comp_exists= false;
        for (Uint i= 0; i < 10; i++)
            if (n.WithUnknowns (i)) {
                if ((*component_of_dof_)[n.num[i]] == c_)
                    comp_exists= true;
                // Remove different components (!= c_) of the same sign (this change does not alter
                // the position of the boundary of component c_.
                else if (sign_of_component_[c_] == sign_of_component_[(*component_of_dof_)[n.num[i]]]) {
                    loc_phi[i]= -sign_of_component_[c_];
                }
            }
        if (!comp_exists)
            return;
    }

    evaluate_on_vertexes (loc_phi, lat, Addr (ls_values));
    partition.make_partition< SortedVertexPolicyCL,MergeCutPolicyCL> (lat, ls_values);
    make_CompositeQuad3Domain (qdom, partition);
    GridFunctionCL<> integrand (1., qdom.vertex_size());
    *thread_volume+= quad( integrand, 6.*tet.GetVolume(), qdom, sign_of_component_[c_] == -1 ? NegTetraC : PosTetraC);
}


//*****************************************************************************
//                               ComponentBasedVolumeAdjustmentCL
//*****************************************************************************

ComponentBasedVolumeAdjustmentCL::ComponentBasedVolumeAdjustmentCL(LevelsetP2CL* lset)
    : VolumeAdjustmentCL(lset)
{}

ComponentBasedVolumeAdjustmentCL::~ComponentBasedVolumeAdjustmentCL()
{}


void ComponentBasedVolumeAdjustmentCL::compute_indicator_functions (const SparseMatBaseCL<unsigned char>& A)
{
    doCorrection_= std::vector<bool> (num_components(), true);

    if (num_components() == 2) {
        indicator_functions_.clear();
        for (Uint c=0; c<2; ++c)
            indicator_functions_.push_back(sign_of_component_[c] == 1 ? std::valarray<double>() : std::valarray<double>(1.0, component_of_dof_.size()));
        return;
    }

    const component_vector tmp (ExtendOneStep (A, component_of_dof_, doCorrection_));
    const component_vector tmp2 (ExtendOneStep (A, tmp, doCorrection_));

    indicator_functions_= std::vector<std::valarray<double>> (num_components(), std::valarray<double>(component_of_dof_.size()));
    // For every point in tmp2, check which component is present and write into the corresponding characteristic function.
    for(Uint j= 0; j < component_of_dof_.size(); ++j)
        indicator_functions_[tmp2[j]][j]= 1.;
}

auto ComponentBasedVolumeAdjustmentCL::ExtendOneStep(const SparseMatBaseCL<unsigned char>& A,
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
    return ret;
}

void ComponentBasedVolumeAdjustmentCL::FindComponents ()
{
    // TimerCL time;
    // time.Start();
    // M encodes the adjacency matrix of two graphs G_comp and G_mesh: a value > 0 is an edge in G_mesh, a value > 1 is an edge within a G_comp (i.e. the vertices of the edge have the same sign of the level set function).
    SparseMatBaseCL<unsigned char> M;
    AdjacencyAccuCL accu (*lset_, M);
    TetraAccumulatorTupleCL accus;
    accus.push_back (&accu);
    accumulate (accus, lset_->GetMG(), lset_->idx.TriangLevel(), lset_->idx.GetBndInfo());
    // time.Stop();
    // std::cout << "setup: " << time.GetTime() << " seconds" << std::endl;

    GraphComponentsCL Split;
    Split.number_connected_components (SparseMatrixGraphCL (M).set_edge_limit (1)); // Operate on the graph representing the components.

    component_of_dof_= Split.component_map();
    Volumes.resize (Split.num_components()); // neccessary to make num_components() return the current number of components.

    sign_of_component_.resize (num_components());
    ComputeReferencePoints();
    for (Uint c= 0; c < num_components(); ++c)
        Volumes[c]= CalculateVolume(c, 0.);
    compute_indicator_functions (M); // Uses the graph representing the whole mesh.
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

    // Compatibility with global volume correction (in case the intial volume is know analytically, it is enforced prior to the simulation).
    if (num_components() == 2){
        const size_t cneg = sign_of_component_[0] == -1 ? 0 : 1;
        if(global_reference_volume_ != -1.0)
            targetVolumes[cneg] = global_reference_volume_;
    }

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
    os << "ComponentBasedVolumeAdjustmentCL: Number of components " << num_components() << "\nVolumes ";
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
    VecDescCL Copy (&lset_->idx);
    Copy.Data= lset_->Phi.Data;

    VolumeAccuCL accu(*lset_, Copy, c, sign_of_component_, component_of_dof_);
    if (shift != 0.) {
        Copy.Data+= shift*indicator_functions_[c];
        accu.use_indicator_functions(indicator_functions_);
    }
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    accumulate( accus, lset_->GetMG(), lset_->idx.TriangLevel(), lset_->idx.GetBndInfo());
    return accu.get_volume();
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

    if (std::count(cnew.begin(), cnew.end(),-1) > 0)
        throw DROPSErrCL("ComponentBasedVolumeAdjustmentCL::component_of_point: Your problem changed from two-phase to one-phase or vice versa... this is not yet supported (and we encourage you to think about what you are doing!).\n");

    return cnew;
}

void ComponentBasedVolumeAdjustmentCL::MatchComponents ()
{
    // cold represents a permutation which maps the new component number (from ReferencePoints_) to the old.
    const component_vector cold (component_of_point (ReferencePoints, sign_of_component_, component_of_dof_backup_, coord_of_dof_backup_, sign_of_component_backup_));
    // cnew represents a permutation which maps the old component number (from ReferencePoints_backup) to the new.
    const component_vector cnew (component_of_point (ReferencePoints_backup, sign_of_component_backup_, component_of_dof_, coord_of_dof_, sign_of_component_));


    // M will be the adjacency matrix of the undirected graph G: The nodes of G are the (old and new) components. The edges are given by the components of the reference points (in both directions).
    SparseMatBaseCL<unsigned char> M;
    const size_t nnew= num_components(),
                 nold= ReferencePoints_backup.size(),
                 n= nold + nnew;
    SparseMatBuilderCL<unsigned char> Mb (&M, n, n);
    for (Uint i= 0; i < nold; ++i)
        Mb (cnew[i] + nold, i)= Mb (i, cnew[i] + nold)= 1;
    for (Uint i= 0; i < nnew; ++i)
        Mb (cold[i], i + nold)= Mb (i + nold, cold[i])= 1;
    Mb.Build();
    // std::cout << "ComponentBasedVolumeAdjustmentCL::MatchComponents: M: " << M << std::endl;

    GraphComponentsCL G;
    G.number_connected_components (SparseMatrixGraphCL (M));

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

VTKVariableCL& ComponentBasedVolumeAdjustmentCL::make_VTKComponentMap(std::string var_name)
{
    return *new VTKComponentMapCL( lset_, &component_of_dof_, var_name);
}

} // end of namespace DROPS
