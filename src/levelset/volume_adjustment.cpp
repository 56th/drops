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
 std::vector<size_t> component_backup_;

 std::vector<size_t> extendedcomponent_;      ///< extendedcomponent[i] is the number of the extended component of i.
 std::vector<bool> DoCorrection_;
 std::vector<std::valarray<double> > comp_char_functions;


 /// \brief Extend all components but the one with number 0 by one level
 void ExtendOneStep( const MatrixCL& A, const std::vector<size_t>& helper);

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
 const std::vector<size_t>& component_map_backup () const { return component_backup_; }
 const std::vector<size_t>& extended_component_map  () const { return extendedcomponent_; }
 bool DoCorrection (size_t x) const { return DoCorrection_[x]; }
 void make_backup() {component_backup_=component_;}
 void clear() {component_.clear(); component_backup_.clear(); extendedcomponent_.clear(); DoCorrection_.clear(); comp_char_functions.clear(); }

 std::vector<size_t>& write_component_map  () { return component_; }

 /// \brief component(c) contains all vertices in component c. The vertices are in ascending order.
 std::vector<size_t> component  (size_t c) const;
 /// \brief the number of vertices in component i.
 const std::vector<size_t>& component_size () const { return component_size_; }
 /// \brief High-level overview of the numbering process.
 inline void stats (std::ostream& os) const;
 /// \brief Renumbers the components in such a way, that the outer phase is always component 0
 void renumber_components(const LevelsetP2CL& lset);
 /// \brief Extends the connectivity components by a neighborhood of 2 DOFs
 void ExtendComponents(const MatrixCL& A);
 /// \brief This routine is called later in the process to generate characteristic functions for all the components
 void GenerateCharFunctions();

 const std::valarray<double>& GetCharFunc (Uint Komp) const {return comp_char_functions[Komp];}
};

void GraphComponentsCL::GenerateCharFunctions(){
    // initialize functions... give them the right dimension first
    comp_char_functions.clear();
    comp_char_functions.resize(num_components(),std::valarray<double>(extendedcomponent_.size()));
    // check every point in the extendedcomponent_-function for which component is present and write into the corresponding characteristic function
    for(Uint j=0; j<extendedcomponent_.size(); ++j)
        comp_char_functions[extendedcomponent_[j]][j]=1;
}

void GraphComponentsCL::ExtendComponents(const MatrixCL& A){
    ExtendOneStep(A,component_);
    std::vector<size_t> temp(extendedcomponent_);
    ExtendOneStep(A,temp);
    GenerateCharFunctions();
}

void GraphComponentsCL::ExtendOneStep(const MatrixCL& A, const std::vector<size_t>& helper){
    extendedcomponent_=helper;
    for (size_t i=0; i<A.num_rows(); ++i)
    {
        const size_t v0=i;
        if (helper[v0]==0)
            continue;
        for(size_t j=A.row_beg(i); j!=A.row_beg(i+1); j++)
        {
            const size_t v1=A.col_ind(j);
            if (helper[v1]==0)
                extendedcomponent_[v1]=helper[v0];
            if(helper[v1]!=0 && helper[v1]!=helper[v0])
            {
               DoCorrection_[helper[v0]]=false;
               DoCorrection_[helper[v1]]=false;
            }
        }
    }
}

void GraphComponentsCL::renumber_components(const LevelsetP2CL& lset)
{
    // This routine ensures that component 0 is always the _outer_ component, means the surrounding fluid, where the levelset function takes positiv values
    Uint komponente0=0;
    Uint ticker=0;
    for(; ticker<lset.Phi.Data.size(); ++ticker)
        if(lset.Phi.Data[ticker]>0)
        {
            komponente0=component_map()[ticker];
            break;
        }
    if(ticker==lset.Phi.Data.size())
        throw DROPSErrCL("GraphComponentsCL::renumber_components: komponente0 not found");
    if(komponente0==0) return;
    std::swap(component_size_[0], component_size_[komponente0]);
    for (uint a=0; a<component_.size(); ++a)
        if(component_[a]==0)
            component_[a]=komponente0;
        else if(component_[a]==komponente0)
            component_[a]=0;
}

std::vector<size_t> GraphComponentsCL::component (size_t c) const
{
    std::vector<size_t> ret;
    ret.reserve( component_size()[c]);

    for (size_t i= 0; i < component_.size(); ++i)
        if (component_[i] == c)
            ret.push_back( i);

    return ret;
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
    DoCorrection_.resize(0);
    DoCorrection_.resize(num_components(),true);
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

void SetupAdjacency( MatrixCL& A, MatrixCL& B, const LevelsetP2CL& lset)
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


void ComponentCL::DebugOutput (std::ostream& os) const
{
     os << "c_: " << c_ << ", volume_: " << volume_ << ", reference_volume_: " << reference_volume_
        << ", refPoint_: " << refPoint_ << ", doCorrection_: " << doCorrection_
        << ", char_function_.sum(): " << char_function_.sum() << '\n';
}

ComponentBasedVolumeAdjustmentCL::ComponentBasedVolumeAdjustmentCL(LevelsetP2CL* lset)
    : VolumeAdjustmentCL(lset), Split(*new GraphComponentsCL()) 
{}

ComponentBasedVolumeAdjustmentCL::~ComponentBasedVolumeAdjustmentCL()
{
    delete &Split;
}


int ComponentBasedVolumeAdjustmentCL::GetNumberOfComponents() const
{
    return Split.num_components();
}

GraphComponentsCL& ComponentBasedVolumeAdjustmentCL::GetSplit()
{
    return Split;
}

void ComponentBasedVolumeAdjustmentCL::InitVolume_impl()
{              
    MatrixCL CompAdja;
    MatrixCL MeshAdja;
    SetupAdjacency (CompAdja,MeshAdja, *lset_);
    Split.number_connected_components(CompAdja);
    Split.renumber_components(*lset_);
    CalculateInitialVolumes();
    FindReferencePoints();
    make_backup(true);
    Split.ExtendComponents(MeshAdja);
    for (int c=0; c < (int)Split.num_components(); ++c) {
        ComponentCL cp{c, Volumes[c], Volumes[c], ReferencePoints[c], Split.DoCorrection(c), Split.GetCharFunc (c)};
        components_.push_back (cp);
        components_backup_.push_back (cp);
         
    }
    component_of_dof_= Split.component_map();
    component_of_dof_backup_= Split.component_map();

    DebugOutput (std::cout);
}
      
void ComponentBasedVolumeAdjustmentCL::Repair()
{
    MatrixCL CompAdja;
    MatrixCL MeshAdja;
    SetupAdjacency (CompAdja,MeshAdja, *lset_);
    Split.clear();
    Split.number_connected_components(CompAdja);
    MatchComponents();
    FindReferencePoints();
    make_backup();
    Split.ExtendComponents(MeshAdja);

    components_.clear();
    components_backup_.clear();
    for (int c=0; c < (int)Split.num_components(); ++c) {
        ComponentCL cp{c, Volumes[c], Volumes[c], ReferencePoints[c], Split.DoCorrection(c), Split.GetCharFunc (c)};
        components_.push_back (cp);
        components_backup_.push_back (cp);
         
    }
    component_of_dof_= Split.component_map();
    component_of_dof_backup_= Split.component_map();
    DebugOutput (std::cout);
}

void ComponentBasedVolumeAdjustmentCL::DebugOutput (std::ostream& os) const
{
    os << "ComponentBasedVolumeAdjustmentCL::Number of components " << GetNumberOfComponents() << "\nVolumes ";
    seq_out(std::begin(Volumes),std::end(Volumes),os,", ");
    os << std::endl << "ReferencePoints\n";
    seq_out(std::begin(ReferencePoints),std::end(ReferencePoints),os,", ");
    os << std::endl << "components_:\n";
    for (auto cp: components_)
        cp.DebugOutput (os);
    os << "components_backup_:\n";
    for (auto cp: components_backup_)
        cp.DebugOutput (os);
    os << std::endl;
}


//*****************************************************************************
//                               ComponentBasedVolumeAdjustmentCL
//*****************************************************************************

void ComponentBasedVolumeAdjustmentCL::FindReferencePoints() {
    std::vector<double> CurrentAbsMax(GetNumberOfComponents(),std::numeric_limits<double>::min());
    ReferencePoints.clear();
    ReferencePoints.resize(GetNumberOfComponents());
    LocalNumbP2CL n;
    LocalP2CL<> loc_phi;
    DROPS_FOR_TRIANG_TETRA( lset_->GetMG(), lset_->idx.TriangLevel(), it) {
        loc_phi.assign(*it,lset_->Phi,lset_->GetBndData());
        n.assign_indices_only(*it, lset_->idx.GetFinest());
        for (Uint a=0;a<10;a++) {
            if (n.WithUnknowns(a)) {
                const Uint CurrentComponent = Split.component_map()[n.num[a]];
                if (CurrentAbsMax[CurrentComponent] < std::abs(loc_phi[a])) {
                    ReferencePoints[CurrentComponent] = a<4? it->GetVertex(a)->GetCoord() : GetBaryCenter(*it->GetEdge(a-4));
                    CurrentAbsMax[CurrentComponent] = std::abs(loc_phi[a]);
                }
            }
        }
    }
}

void ComponentBasedVolumeAdjustmentCL::CalculateInitialVolumes() {
    CalculateVolumes(Volumes, 0);
    std::cout << "Initial Volumes" << std::endl;
    for (Uint i=0; i<Volumes.size(); ++i)
        std::cout << "Component "<< i << ": " << Volumes[i] << std::endl;
}

void ComponentBasedVolumeAdjustmentCL::CalculateVolumes(std::valarray<double>& volumes, std::valarray<double>* epsilons) const {
    volumes.resize(0);
    volumes.resize(Split.num_components());
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance ( 2);
    std::valarray<double> ls_values (lat.vertex_size());
    QuadDomainCL qdom;
    LocalP2CL<> loc_phi;
    TetraPartitionCL partition;
    LocalNumbP2CL n;

    VecDescCL Copy(&lset_->idx);
    Copy.Data=lset_->Phi.Data;

    if (epsilons!=0)
        for (Uint i=0; i<epsilons->size(); ++i)
            Copy.Data+=(*epsilons)[i]*Split.GetCharFunc(i);

//    DROPS_FOR_TRIANG_TETRA( lset_->GetMG(), lset_->idx.TriangLevel(), it) {
//        loc_phi.assign(*it,Copy,lset_->GetBndData());
//        evaluate_on_vertexes (loc_phi, lat, Addr(ls_values));
//        partition.make_partition< SortedVertexPolicyCL,MergeCutPolicyCL>(lat, ls_values);
//        make_CompositeQuad3Domain( qdom, partition);
//        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
//        n.assign_indices_only(*it, lset_->idx);
//        for(uint a=0;a<10;a++)
//            if(n.WithUnknowns(a) && Split.component_map()[n.num[a]]){
//                // run over all DOF within the tetrahedron and check the values of component_map()[i] to decide, to which component the volume should be added
//                volumes[Split.component_map()[n.num[a]]]+= quad( integrand, it->GetVolume()*6., qdom, NegTetraC);
//                break;         // component 0 is the surrounding fluid
//            }
//        volumes[0]+= quad( integrand, it->GetVolume()*6., qdom, PosTetraC); // Volume of the positive part is always added to the surrounding fluid
//    }

    LocalP2CL<double> changed_loc_phi;
    DROPS_FOR_TRIANG_TETRA( lset_->GetMG(), lset_->idx.TriangLevel(), it) {
        n.assign_indices_only(*it, lset_->idx.GetFinest());
        loc_phi.assign(*it,Copy,lset_->GetBndData());

        for (uint k=1;k<Split.num_components() ; k++) {
            changed_loc_phi = loc_phi;
            bool comp_exists = false;
            for (uint a=0;a<10;a++)
                 if (n.WithUnknowns(a)) {
                     if (Split.component_map()[n.num[a]] == k)
                         comp_exists = true;
                     else if (Split.component_map()[n.num[a]] > 0)   // by definition, component 0 is the surrounding liquid
                         changed_loc_phi[a] = 1.0;                   // remove negative components which are not the considered component k (this change does not alter the position of the boundary of component k
                 }
            if (!comp_exists)
                continue;
            evaluate_on_vertexes (changed_loc_phi, lat, Addr(ls_values));
            partition.make_partition< SortedVertexPolicyCL,MergeCutPolicyCL>(lat, ls_values);
            make_CompositeQuad3Domain( qdom, partition);
            DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
            volumes[k]+= quad( integrand, it->GetVolume()*6., qdom, NegTetraC);   // all remaining negative parts are of component k
        }

        evaluate_on_vertexes (loc_phi, lat, Addr(ls_values));
        partition.make_partition< SortedVertexPolicyCL,MergeCutPolicyCL>(lat, ls_values);
        make_CompositeQuad3Domain( qdom, partition);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        volumes[0]+= quad( integrand, it->GetVolume()*6., qdom, PosTetraC); // Volume of the positive part is always added to the surrounding fluid
    }
}

double ComponentBasedVolumeAdjustmentCL::ComputeComponentAdjustment (int compnumber)
{
    return compute_volume_correction (
        [this,compnumber](double x)->double {
            std::valarray<double> shift(Split.num_components());
            shift[compnumber]=x;
            CalculateVolumes(Volumes, &shift);
            return Volumes[compnumber];
        },
        Volumes_backup[compnumber],
        tol_
    );
}

void ComponentBasedVolumeAdjustmentCL::MatchComponents() {
    Uint RPBS = ReferencePoints_backup.size();

    // Create and initialize temporary arrays to store the minimal distances to the reference points and the corresponding minimizer
    std::vector<double> distances(RPBS, std::numeric_limits<double>::max());
    std::vector<int> globalIDs(RPBS, -1);

    double tempdistance = 0.0;
    LocalNumbP2CL n;
    DROPS_FOR_TRIANG_TETRA( lset_->GetMG(), lset_->idx.TriangLevel(), it) {
        n.assign_indices_only(*it, lset_->idx.GetFinest());
        for (Uint a=0;a<10;a++) {
            // if there is a DOF
            if (n.WithUnknowns(a)) {
                // Calculate the distance to every reference point... if it's smaller: store the new distance and global index
                for (Uint ap=0; ap<RPBS; ++ap) {
                    tempdistance = ((a<4 ? it->GetVertex(a)->GetCoord() : GetBaryCenter(*it->GetEdge(a-4)))-ReferencePoints_backup[ap]).norm();
                    if (tempdistance < distances[ap]) {
                        distances[ap] = tempdistance;
                        globalIDs[ap] = n.num[a];
                    }
                }
            }
        }
    }

    // Create a temporary function to store the component in each DOF
    std::vector<size_t> temp(Split.component_map().size());
    for (Uint old_component_number=0; old_component_number<RPBS; ++old_component_number) {
        const Uint new_component_number = Split.component_map()[globalIDs[old_component_number]];
        const std::vector<size_t>& componentI = Split.component(new_component_number); // gather all points that belong to the component with the new component number
        for (Uint it=0; it<componentI.size(); ++it) {
            temp[componentI[it]] = old_component_number;                              // "color" them with respect to the old component number, i.e. relabel them
        }
    }
    Split.write_component_map() = temp;                                // fill the new information into the private member "component_" of Split
}

void ComponentBasedVolumeAdjustmentCL::AdjustVolume() {
    MatrixCL ComponentAdja;
    MatrixCL FullAdja;
    SetupAdjacency(ComponentAdja, FullAdja, *lset_);
    Split.number_connected_components(ComponentAdja);
    Split.renumber_components(*lset_); // after this step component 0 is the surrounding liquid
    FindReferencePoints();
    CalculateVolumes(Volumes, 0);
    if (!Handle_topo_change())
        MatchComponents();

    Split.ExtendComponents(FullAdja);
    // adapt Level Set
    for (Uint i= 1; i<Split.num_components(); ++i) {
        if (Split.DoCorrection(i)) {
            const double s= ComputeComponentAdjustment(i);
            lset_->Phi.Data+=Split.GetCharFunc(i)*s;
            std::cout << "Adjustment for component " << i << ": " << s << "\n";
        }
        else {
            std::cout << "Adjustment for component " << i << ": There hasn't been made any adjustment, components to close \n";
        }
    }
    FindReferencePoints();
    make_backup();
    SetupAdjacency(ComponentAdja, FullAdja, *lset_);
    Split.number_connected_components(ComponentAdja);
    MatchComponents();
    FindReferencePoints();
    make_backup();
}

void ComponentBasedVolumeAdjustmentCL::make_backup(bool complete) {
    Split.make_backup();
    ReferencePoints_backup=ReferencePoints;
    if(complete)
        Volumes_backup=Volumes;
}

bool ComponentBasedVolumeAdjustmentCL::Handle_topo_change(){
    bool change=false;
    Uint RPS =ReferencePoints.size();
    Uint RPBS=ReferencePoints_backup.size();
    if (RPS!=RPBS) {
        change=true;
        if (RPS==RPBS+1) { // New component: find out which one of the old components split up
            // look up the old component numbers of the new reference points in the old map and store the results temporarily
            std::vector<size_t> temp(RPS); // global index of the DOF closest to the referencepoint
            // std::vector<size_t> tempOldGrid(RPS); not used yet
            std::vector<size_t> ORPN(RPBS); // ORPN = Old reference point number ... global index of the DOF closest to the referencepoint_backup
            std::vector<double> distances(std::numeric_limits<double>::max(),RPS);
            // std::vector<double> DistanzenOldGrid(std::numeric_limits<double>::max(),RPS); not used yet
            std::vector<double> distancesORPN(std::numeric_limits<double>::max(),RPBS);
            double tempdistance;
            // run over the up2date grid. for each old and new reference point one gets the closest DOF, including the distance of the DOF to the considered reference point
            LocalNumbP2CL n;
            DROPS_FOR_TRIANG_TETRA( lset_->GetMG(), lset_->idx.TriangLevel(), it) {
                n.assign_indices_only(*it, lset_->idx.GetFinest());
                for (Uint a=0; a<10; a++) {
                    // if there is a DOF
                    if (n.WithUnknowns(a)) {
                        // Calculate the distance to every reference point... if it's smaller: store the new distance and global index
                        for (Uint ap=0; ap<RPS; ++ap) {
                            tempdistance = ((a<4 ? it->GetVertex(a)->GetCoord() : GetBaryCenter(*it->GetEdge(a-4)))-ReferencePoints[ap]).norm();
                            if (tempdistance < distances[ap]) {
                                distances[ap] = tempdistance;
                                temp[ap] = n.num[a];
                            }
                        }
                        for (Uint ap2=0; ap2<RPBS; ++ap2) {
                            tempdistance = ((a<4 ? it->GetVertex(a)->GetCoord() : GetBaryCenter(*it->GetEdge(a-4)))-ReferencePoints_backup[ap2]).norm();
                            if (tempdistance < distancesORPN[ap2]) {
                                distancesORPN[ap2] = tempdistance;
                                ORPN[ap2] = n.num[a];
                            }
                        }
                    }
                }
            }
            std::vector<size_t> OldAffiliation(RPS);
            for (Uint i=1; i<RPS; ++i)
                OldAffiliation[i] = Split.component_map_backup()[temp[i]];///!!!
            // in OldAffiliation one finds the hypothetical affiliation of the new reference points with respect to the old components

            std::vector<size_t> SortingCopy(OldAffiliation);
            std::vector<size_t>::iterator finder;
            std::vector<size_t>::iterator finder_twofold;

            // compare pairwise and check which one appears twice
            std::sort(SortingCopy.begin(),SortingCopy.end());
            finder = std::adjacent_find(SortingCopy.begin(),SortingCopy.end());
            if (finder == SortingCopy.end())
                DROPSErrCL("ComponentBasedVolumeAdjustmentCL::Handle_topo_change(): New Component, but no old component appeared twice");
            // in finder there is the number of the old volume, which split up

            // get the old volume and redistribute it on a percentage basis
            double VolumeToSplit = Volumes_backup[*finder];

            // find the first component that appears twice according to the old numbering
            finder_twofold = std::find(OldAffiliation.begin(),OldAffiliation.end(),*finder);
            size_t Comp1 = *finder_twofold;
            // find the second component that appears twice according to the old numbering
            // (there has to be a second one, otherwise there would have been an error message earlier)
            size_t Comp2 = *(std::find(finder_twofold+1,OldAffiliation.end(),*finder));


            double VolumeSum=Volumes[Comp1]+Volumes[Comp2];
            // calculate the volumes as to how they should have been, in order to avoid any mass loss
            Volumes[Comp1]=Volumes[Comp1]/VolumeSum*VolumeToSplit;
            Volumes[Comp2]=Volumes[Comp2]/VolumeSum*VolumeToSplit;

            for (Uint a=0; a<RPBS; ++a) {
                if (a == *finder)
                    continue;
                Volumes[Split.component_map()[ORPN[a]]]=Volumes_backup[Split.component_map_backup()[ORPN[a]]];
            }
            // new status for the further simulation
            Volumes_backup.resize(Volumes.size());
            Volumes_backup=Volumes;
        }
        else if (RPS==RPBS-1) { // component vanished
            // find out which of the old components coalesced

            std::vector<size_t> temp(RPBS);
            std::vector<double> distances(std::numeric_limits<double>::max(),RPBS);

            // run through grid and find the distance minimizers to all old reference points
            double tempdistance;

            LocalNumbP2CL n;
            DROPS_FOR_TRIANG_TETRA( lset_->GetMG(), lset_->idx.TriangLevel(), it) {
                n.assign_indices_only(*it, lset_->idx.GetFinest());
                for (Uint a=0; a<10; a++) {
                    // if there is a DOF
                    if (n.WithUnknowns(a)) {
                        // Calculate the distance to every reference point... if it's smaller: store the new distance and global index
                        for (Uint ap=0; ap<RPBS; ++ap) {
                            tempdistance = ((a<4 ? it->GetVertex(a)->GetCoord() : GetBaryCenter(*it->GetEdge(a-4)))-ReferencePoints_backup[ap]).norm();
                            if (tempdistance < distances[ap]) {
                                distances[ap] = tempdistance;
                                temp[ap] = n.num[a];
                            }
                        }
                    }
                }
            }
            // initialize new volumes
            for (Uint a=0; a<RPS; ++a) Volumes[a]=0;
            // run through all old reference points, check which old volume was present and add it to the new volume at the same place
            for (Uint b=0; b<RPBS; ++b)
                Volumes[Split.component_map()[temp[b]]]+=Volumes_backup[Split.component_map_backup()[temp[b]]];
            // from here on forward the simulation is newly set up... the component numbers are not conserved in this method
            // new status for the further simulation
            Volumes_backup.resize(Volumes.size());
            Volumes_backup=Volumes;
        }
        else {
            DROPSErrCL("ComponentBasedVolumeAdjustmentCL::Handle_topo_change() : Numbers of components differ by a number bigger than one... case not considered yet ... I give up...");
        }
    }
    return change;
}





} // end of namespace DROPS
