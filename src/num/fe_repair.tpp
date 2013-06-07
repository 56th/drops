/// \file fe_repair.tpp
/// \brief Repair-classes with respect to multigrid-changes for finite-elements.
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande; SC RWTH Aachen:

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

namespace DROPS
{

/// RepairFEDataCL

template <class LocalFEDataT>
void RepairFEDataCL<LocalFEDataT>::repair (AugmentedDofVecT& dof, VectorCL& newdata) const
{
    AugmentedDofVecT::iterator d;
    BaryCoordCL tmp( Uninitialized);
    for (typename ChildVecT::const_iterator old_ch= data.begin(); !(dof.empty() || old_ch == data.end()); ++old_ch) {
        const SMatrixCL<4,4>& to_old_child= parent_to_child_bary( old_ch->first);
        d= dof.begin();
        while (d != dof.end()) {
            tmp= to_old_child*d->second;
            if (contained_in_reference_tetra( tmp, 8*std::numeric_limits<double>::epsilon())) {
                DoFHelperCL<typename LocalFEDataT::value_type, VectorCL>::set( newdata, d->first, old_ch->second( tmp));
                d= dof.erase( d);
            }
            else
                ++d;
        }
    }
    if (!dof.empty())
        throw DROPSErrCL("RepairFEDataCL::repair: Could not locate all new dof.\n");
}

template <class LocalFEDataT>
bool RepairFEDataCL<LocalFEDataT>::has_child (Ubyte ch) const
{
    for (size_t i= 0; i < data.size(); ++i)
        if (data[i].first == ch)
            return true;
    return false;
}


/// RepairFECL

template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
  RepairFECL<ValueT, LocalFEDataT, BndDataT>::RepairFECL (const MultiGridCL& mg, const VecDescCL& old, const BndDataT<value_type>& bnd)
        : mg_( mg), old_vd_ ( old), bnd_( bnd)
{
    localfedata_::Init();
    pre_refine();
}

template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
  void
  RepairFECL<ValueT, LocalFEDataT, BndDataT>::pre_refine ()
{
    parent_data_.clear();
    level0_leaves_.clear();
    repair_needed_.clear();

    Uint lvl= old_vd_.RowIdx->TriangLevel();
    LocalFECL lp2;
    DROPS_FOR_TRIANG_CONST_TETRA( mg_, lvl, it) {
        if (!it->IsUnrefined())
            continue;
        if (it->GetLevel() > 0) {
            // These could be deleted by the refinement algo
            /// \todo To store less, one can generally add "if ( it->IsMarkedForRemovement() || !it->IsRegular())" here. However, the TetrabuilderCL::BogoReMark-algorithm does not obey the user-marks, thus we postpone this optimization.
            const TetraCL* p= it->GetParent();
            const Ubyte ch= std::find( p->GetChildBegin(), p->GetChildEnd(), &*it) - p->GetChildBegin();
            const RefRuleCL& rule= p->GetRefData();
            lp2.assign_on_tetra( *it, old_vd_, bnd_);
            parent_data_[p].data.push_back( std::make_pair( rule.Children[ch], lp2));
        }
        else {
            // Leaves in level 0 can give birth to children, which cannot easily be distinguished from tetras, which just remained over one refinement step, cf. case (2) and (4b) in repair(). We memoize the leaves in level 0, as they tend to be few.
            level0_leaves_.insert( &*it);
        }
    }
#ifdef _PAR
    HandlerParentDataCL handler( parent_data_);
    DiST::InterfaceCL* interf;
    // init the interface
    DiST::PrioListT from, to;
    from.push_back( PrioGhost);
    to.push_back( PrioMaster);
    DiST::InterfaceCL::DimListT dimlist;
    dimlist.push_back( DiST::GetDim<TetraCL>());
    DiST::LevelListCL lvls( lvl); // levels 0,..,TriangLevel

    interf = new DiST::InterfaceCL( lvls, from, to, dimlist, true);
    interf->Communicate( handler);
    delete interf;
#endif
}

template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
  AugmentedDofVecT
  RepairFECL<ValueT, LocalFEDataT, BndDataT>::collect_unrepaired_dofs (const TetraCL& t)
{
    LocalNumbCL n_new( t, *new_vd_->RowIdx);
    AugmentedDofVecT dof;
    dof.reserve( localfedata_::NumUnknownsOnTetra);
    for (Uint i= 0; i < localfedata_::NumUnknownsOnTetra; ++i)
        if (n_new.WithUnknowns( i) && repair_needed( n_new.num[i])) {
            dof.push_back( std::make_pair( n_new.num[i], localfedata_::p_dof_[i]));
            mark_as_repaired( n_new.num[i]); // All callers will repair all dofs or else throw an exception.
        }
    return dof;
}

template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
void RepairFECL<ValueT, LocalFEDataT, BndDataT>::unchanged_refinement (const TetraCL& t)
{
    const VectorCL& olddata= old_vd_.Data;
          VectorCL& newdata= new_vd_->Data;
    LocalNumbCL n_old( t, *old_vd_.RowIdx);
    LocalNumbCL n_new( t, *new_vd_->RowIdx);
    for (Uint i= 0; i < localfedata_::NumUnknownsOnTetra; ++i)
        if (n_new.WithUnknowns( i) && repair_needed( n_new.num[i])) {
            Assert( n_old.WithUnknowns( i), DROPSErrCL( "TetraRepairFECL::unchanged_refinement: "
                "Old and new function must use the same boundary-data-types.\n"), DebugNumericC);
            const value_type& tmp= DoFHelperCL<value_type, VectorCL>::get( olddata, n_old.num[i]);
            DoFHelperCL<value_type, VectorCL>::set( newdata, n_new.num[i], tmp);
            mark_as_repaired( n_new.num[i]);
        }
}

template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
void RepairFECL<ValueT, LocalFEDataT, BndDataT>::regular_leaf_refinement (const TetraCL& t)
{
    const AugmentedDofVecT& dof= collect_unrepaired_dofs( t);
    if (dof.empty())
        return;

    const TetraCL* p= t.GetParent();
    const Ubyte ch= std::find( p->GetChildBegin(), p->GetChildEnd(), &t) - p->GetChildBegin();
    const SMatrixCL<4,4>& T= child_to_parent_bary( p->GetRefData().Children[ch]);

    LocalFECL oldp2;
    oldp2.assign_on_tetra( *p, old_vd_, bnd_);
    for (AugmentedDofVecT::const_iterator d= dof.begin(); d != dof.end(); ++d)
        DoFHelperCL<value_type, VectorCL>::set( new_vd_->Data, d->first, oldp2( T*d->second));
}

template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
void RepairFECL<ValueT, LocalFEDataT, BndDataT>::genuine_refinement (const TetraCL& t, const RepairFEDataCL<LocalFECL>& repairdata)
{
    AugmentedDofVecT dof= collect_unrepaired_dofs( t);
    if (dof.empty())
        return;

    const TetraCL* p= t.GetParent();
    const Ubyte ch= std::find( p->GetChildBegin(), p->GetChildEnd(), &t) - p->GetChildBegin();
    const SMatrixCL<4,4>& T= child_to_parent_bary( p->GetRefData().Children[ch]);
    const TetraCL* gp= p->GetParent();
    const Ubyte gpch= std::find( gp->GetChildBegin(), gp->GetChildEnd(), p) - gp->GetChildBegin();
    const SMatrixCL<4,4>& S= child_to_parent_bary( gp->GetRefData().Children[gpch]);
    for (AugmentedDofVecT::iterator d= dof.begin(); d != dof.end(); ++d)
        d->second= S*(T*d->second);
    repairdata.repair( dof, new_vd_->Data);
}

template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
void RepairFECL<ValueT, LocalFEDataT, BndDataT>::unrefinement (const TetraCL& t, const RepairFEDataCL<LocalFECL>& repairdata)
{
    AugmentedDofVecT dof= collect_unrepaired_dofs( t);
    if (dof.empty())
        return;

    QRDecompCL<4> T;
    BaryCoordCL tmp;
    repairdata.repair( dof, new_vd_->Data);
}

template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
void RepairFECL<ValueT, LocalFEDataT, BndDataT>::changed_refinement (const TetraCL& t, Ubyte ch, const RepairFEDataCL<LocalFECL>& repairdata)
{
    AugmentedDofVecT dof= collect_unrepaired_dofs( t);
    if (dof.empty())
        return;

    const SMatrixCL<4,4>& to_parent= child_to_parent_bary( ch);
    for (AugmentedDofVecT::iterator d= dof.begin(); d != dof.end(); ++d)
        d->second= to_parent*d->second;
    repairdata.repair( dof, new_vd_->Data);
}

template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
  void
  RepairFECL<ValueT, LocalFEDataT, BndDataT>::repair (VecDescCL& new_vd)
{
    new_vd_= &new_vd;

    const Uint lvl= new_vd_->RowIdx->TriangLevel();
    Assert( lvl == old_vd_.RowIdx->TriangLevel() || lvl ==  old_vd_.RowIdx->TriangLevel() - 1,
        DROPSErrCL( "RepairFECL<ValueT>::repair: Different levels\n"), DebugNumericC);
    if (lvl == old_vd_.RowIdx->TriangLevel() - 1)
        std::cout << "old level: " << old_vd_.RowIdx->TriangLevel() << " mg_.GetLastLevel(): " << mg_.GetLastLevel() << '\n';

    VectorCL& newdata= new_vd_->Data;
    repair_needed_.resize( newdata.size(), true);

    DROPS_FOR_TRIANG_CONST_TETRA( mg_, lvl, t) {
        if (parent_data_.count( &*t) == 1) // Case 1
            unrefinement( *t, parent_data_[&*t]);
        // From here on, t has no parent-data itself.
        else if (t->GetLevel() == 0) // Case 2
            // If t had arrived in the triangulation via coarsening, it would have had parent-data.
            unchanged_refinement( *t);
        else { // t has no parent-data and t->GetLevel() > 0
            const TetraCL* p= t->GetParent();
            if (parent_data_.count( p) == 1) { // Case 3
                const Ubyte chpos= std::find( p->GetChildBegin(), p->GetChildEnd(), &*t) - p->GetChildBegin(),
                            ch= p->GetRefData().Children[chpos]; // Child-number of t per topo.h.
                if (parent_data_[p].has_child( ch)) // t was a child of p before refinement, handle as case (2b).
                    unchanged_refinement( *t);
                else
                    changed_refinement( *t, ch, parent_data_[p]);
            }
            else { // t has no repair-data, t->GetLevel() > 0, and p has no repair-data
                if ((p->GetLevel() > 0 && parent_data_.count( p->GetParent()) == 1))
                    genuine_refinement( *t, parent_data_[p->GetParent()]); // Case (4a).
                else if (level0_leaves_.count( p) == 1)
                    regular_leaf_refinement( *t); // Case (4b).
                else // Case (2b)
                    unchanged_refinement( *t);
            }
        }
    }
}

/******************************************************************************
* R E P A I R  M A P   CL                                                     *
******************************************************************************/

/// \brief class for storing repair data on master and ghost tetrahedra
template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
class RepairFECL<ValueT, LocalFEDataT, BndDataT>::RepairMapCL {
  private:
    typedef typename localfedata_::LocalFECL   LocalFECL;
    typedef DROPS_STD_UNORDERED_MAP<const TetraCL*, RepairFEDataCL<LocalFECL> > RepairMapT;
    RepairMapT master_data_;
    RepairMapT ghost_data_;

  public:
    void clear() { master_data_.clear(); ghost_data_.clear(); }
    void push_back(const std::pair<const TetraCL*, RepairFEDataCL<LocalFECL> >& e){
        if (e.first.IsGhost())
            ghost_data_.push_back(e);
        else
            master_data_.push_back(e);
    }
    size_t count(const TetraCL* t){
        if (t->IsGhost())
            return ghost_data_.count(t);
        else
            return master_data_.count(t);
    }
    RepairFEDataCL<LocalFECL> & operator[](const TetraCL* t){
        if (t->IsGhost())
            return ghost_data_[t];
        else
            return master_data_[t];
    }

};


/******************************************************************************
* H A N D L E R  P A R E N T D A T A   C L                                    *
******************************************************************************/
#ifdef _PAR
/// \brief Handler for generating information about parent data on other processes
template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
class RepairFECL<ValueT, LocalFEDataT, BndDataT>::HandlerParentDataCL
{
  private:
    RepairMapCL& parent_data_;

  public:
    HandlerParentDataCL( RepairMapCL& data)
        : parent_data_(data) {}
    /// \brief Gather information about distributed dof on sender proc
    bool Gather( DiST::TransferableCL&, DiST::SendStreamCL&);
    /// \brief Scatter information about distributed dof on sender proc
    bool Scatter( DiST::TransferableCL&, const size_t&, DiST::MPIistreamCL&);
};

template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
bool RepairFECL<ValueT, LocalFEDataT, BndDataT>::HandlerParentDataCL::Gather(
    DiST::TransferableCL& t, DiST::SendStreamCL& send)
{
    TetraCL* tetra;
    simplex_cast(t, tetra);
    if (parent_data_.count( tetra) == 1){
        const typename RepairFEDataCL<LocalFECL>::ChildVecT& data = parent_data_[tetra].data;
        send << data.size();
        for (size_t i = 0; i<data.size(); ++i){
            send << data[i].first;
            for (size_t k=0; k< localfedata_::NumUnknownsOnTetra; ++k) // write LocalPiCL
                send << data[i].second[k];
        }
        return true;
    }
    return false;
}

template <class ValueT, template<class> class LocalFEDataT, template<class> class BndDataT>
bool RepairFECL<ValueT, LocalFEDataT, BndDataT>::HandlerParentDataCL::Scatter(
    DiST::TransferableCL& t, __UNUSED__ const size_t& numData,
    DiST::MPIistreamCL& recv)
{
    Assert(numData == 1, DROPSErrCL("HandlerParentDataCL::Scatter: more than one ghost exists"), DebugDiSTC);

    TetraCL* tetra;
    simplex_cast(t, tetra);
    Ubyte child;
    LocalFECL lp2;
    size_t num_childs;
    recv >> num_childs;
    for (size_t ch =0; ch< num_childs; ++ch){
        recv >> child;
        for (size_t k=0; k< localfedata_::NumUnknownsOnTetra; ++k) // read LocalPiCL
            recv >> lp2[k];
        parent_data_[tetra].data.push_back( std::make_pair( child, lp2));
    }
    return true;
}
#endif
} // end of namespace DROPS
