/// \file prolongation.tpp
/// \brief prolongation and restriction routines for (parallel) multigrid solver
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

namespace DROPS{

template <class ValueT>
void ProlongationCL<ValueT>::Create(IdxDescCL* coarse, IdxDescCL* fine)
{
    coarse_ = coarse;
    fine_   = fine;
    const FiniteElementT fe= coarse->GetFE();
    if (fe == P1_FE){
#ifdef _PAR
        throw DROPSErrCL("P1-prolongation not yet implemented in parallel");
#endif
        BuildP1ProlongationMatrix( *coarse_, *fine);
        return;
    }

    if (fe != vecP2_FE && fe != P2_FE)
        throw DROPSErrCL("ProlongationCL: FE type not supported, yet");

    scale_.resize(fine_->NumUnknowns());

    const IdxT f_idx = fine_->GetIdx();
    const Uint f_level = fine_->TriangLevel();
    const Uint c_level = coarse_->TriangLevel();
    const Uint ndofs= fine_->NumUnknownsVertex();

    // calculate weights for prolongation/restriction
    for (MultiGridCL::const_TetraIterator sit= mg_.GetAllTetraBegin( c_level),
         theend= mg_.GetAllTetraEnd( c_level); sit != theend; ++sit){
        if (!sit->IsInTriang( c_level)) continue;
#ifdef _PAR
        if (sit->HasGhost()) continue;
#endif
        if (!sit->IsInTriang( f_level)) {

            const std::vector<IdxT> tmp ( CollectChildUnknownsP2( *sit, f_idx));
            for (size_t i =0; i< tmp.size(); ++i)
                if (tmp[i] != NoIdx)
                    for (Uint k=0; k<ndofs; ++k)
                        scale_[tmp[i]+k] += 1.0;
        }
        else { // coarse and fine tetra are identical; the restriction is trivial.
            for (Uint i=0; i<4; ++i)
                if (sit->GetVertex( i)->Unknowns.Exist()
                    && sit->GetVertex( i)->Unknowns.Exist( f_idx))
                    for (Uint k=0; k<ndofs; k++)
                        scale_[sit->GetVertex( i)->Unknowns( f_idx)+k] +=1.0;
            for (Uint i=0; i<6; ++i)
                if (sit->GetEdge( i)->Unknowns.Exist()
                    && sit->GetEdge( i)->Unknowns.Exist( f_idx))
                    for (Uint k=0; k<ndofs; k++)
                        scale_[sit->GetEdge( i)->Unknowns( f_idx)+k] +=1.0;
            }
    }

#ifdef _PAR
    fine_->GetEx().Accumulate(scale_);
#endif

    scale_ = 1.0/scale_;
    BuildP2ProlongationMatrix( *coarse_, *fine);
}

template <class ValueT>
std::vector<IdxT> ProlongationCL<ValueT>::CollectChildUnknownsP2(const TetraCL& t, const Uint f_idx) const
{
    const RefRuleCL& rule= t.GetRefData();
    typedef std::map<Ubyte, IdxT>::value_type map_entryT;
    std::map<Ubyte, IdxT> childvertex;
    std::map<Ubyte, IdxT> childedge;

    for (Uint j= 0; j < rule.ChildNum; ++j) {
        const ChildDataCL& child= GetChildData( rule.Children[j]);
        const TetraCL* const childp= t.GetChild( j);
        for (Ubyte k= 0; k<4; ++k) {
            const VertexCL* const p= childp->GetVertex( k);
            childvertex[child.Vertices[k]]= p->Unknowns.Exist()
                                            && p->Unknowns.Exist( f_idx)
                                            ? p->Unknowns( f_idx) : NoIdx;
        }
        for (Ubyte k= 0; k<6; ++k) {
            const EdgeCL* const p= childp->GetEdge( k);
            childedge[child.Edges[k]]= p->Unknowns.Exist()
                                       && p->Unknowns.Exist( f_idx)
                                       ? p->Unknowns( f_idx) : NoIdx;
        }
    }
    std::vector<IdxT> ret( childvertex.size() + childedge.size());
    std::transform( childvertex.begin(), childvertex.end(), ret.begin(),
                    select2nd<map_entryT>());
    std::transform( childedge.begin(), childedge.end(), ret.begin() + childvertex.size(),
                    select2nd<map_entryT>());
    return ret;
}

/// \name Coefficient-tables for P2-prolongation.
/// For each refinement rule the local P2 prolongation matrix is stored
/// as a sequence of matrix rows in P2_local_prolongation_row_idx.
/// Beginning and one-past-the-end of each sequence are stored in
/// P2_local_prolongation_mat_beg.
/// Conceptually, rows of the local prolongation matrices are represented
/// by an index into an array of the 35 possible different rows.
/// Considered as a matrix, this array is sparse. Consequently, it is
/// stored in compressed row storage format.
/// \{
const unsigned int P2_local_prolongation_mat_beg[65]= {
    0, 10, 24, 38, 56, 70, 88, 106,
    128, 142, 160, 178, 200, 219, 241, 263,
    289, 303, 321, 340, 362, 380, 402, 425,
    451, 469, 491, 514, 541, 564, 591, 618,
    649, 663, 682, 700, 723, 741, 764, 786,
    812, 830, 853, 875, 902, 925, 952, 979,
    1010, 1028, 1051, 1074, 1101, 1123, 1150, 1177,
    1208, 1230, 1256, 1283, 1314, 1341, 1372, 1403,
    1438
};

const unsigned char P2_local_prolongation_rows[1438]= {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
    12, 13, 0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 14, 15, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7,
    8, 9, 10, 11, 14, 15, 16, 12, 13, 18, 0, 1, 2, 3, 6, 4, 5, 7, 8, 9, 19, 20, 21,
    18, 0, 1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 11, 19, 20, 21, 12, 13, 17, 0, 1, 2, 3,
    5, 6, 4, 7, 8, 9, 14, 15, 19, 20, 21, 16, 17, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8,
    9, 10, 11, 14, 15, 19, 20, 21, 16, 12, 18, 17, 13, 0, 1, 2, 3, 7, 4, 5, 6, 8,
    9, 22, 23, 24, 25, 0, 1, 2, 3, 4, 7, 5, 6, 8, 9, 10, 11, 22, 23, 24, 12, 13,
    26, 0, 1, 2, 3, 5, 7, 4, 6, 8, 9, 14, 15, 22, 23, 16, 25, 17, 27, 0, 1, 2, 3,
    4, 5, 7, 6, 8, 9, 10, 11, 14, 15, 22, 23, 16, 12, 13, 27, 26, 18, 0, 1, 2, 3,
    6, 7, 4, 5, 8, 9, 19, 20, 22, 23, 21, 24, 25, 18, 28, 0, 1, 2, 3, 4, 6, 7, 5,
    8, 9, 10, 11, 19, 20, 22, 23, 21, 24, 12, 13, 26, 17, 0, 1, 2, 3, 5, 6, 7, 4,
    8, 9, 14, 15, 19, 20, 22, 23, 21, 16, 25, 17, 27, 13, 0, 1, 2, 3, 4, 5, 6, 7,
    8, 9, 10, 11, 14, 15, 19, 20, 22, 23, 21, 16, 12, 27, 26, 18, 17, 13, 0, 1, 2,
    3, 8, 4, 5, 6, 7, 9, 29, 30, 31, 26, 0, 1, 2, 3, 4, 8, 5, 6, 7, 9, 10, 11, 29,
    30, 31, 12, 13, 25, 0, 1, 2, 3, 5, 8, 4, 6, 7, 9, 14, 15, 29, 30, 31, 16, 26,
    17, 28, 0, 1, 2, 3, 4, 5, 8, 6, 7, 9, 10, 11, 14, 15, 29, 30, 31, 16, 12, 13,
    25, 18, 0, 1, 2, 3, 6, 8, 4, 5, 7, 9, 19, 20, 29, 30, 21, 26, 18, 32, 0, 1, 2,
    3, 4, 6, 8, 5, 7, 9, 10, 11, 19, 20, 29, 30, 21, 12, 13, 32, 25, 17, 0, 1, 2,
    3, 5, 6, 8, 4, 7, 9, 14, 15, 19, 20, 29, 30, 21, 16, 26, 17, 32, 13, 28, 0, 1,
    2, 3, 4, 5, 6, 8, 7, 9, 10, 11, 14, 15, 19, 20, 29, 30, 21, 16, 12, 32, 25, 18,
    17, 13, 0, 1, 2, 3, 7, 8, 4, 5, 6, 9, 22, 23, 29, 30, 31, 24, 25, 12, 0, 1,
    2, 3, 4, 7, 8, 5, 6, 9, 10, 11, 22, 23, 29, 30, 31, 24, 13, 26, 25, 12, 0, 1,
    2, 3, 5, 7, 8, 4, 6, 9, 14, 15, 22, 23, 29, 30, 31, 16, 25, 17, 27, 12, 28, 0,
    1, 2, 3, 4, 5, 7, 8, 6, 9, 10, 11, 14, 15, 22, 23, 29, 30, 31, 16, 13, 27, 26,
    25, 12, 18, 28, 0, 1, 2, 3, 6, 7, 8, 4, 5, 9, 19, 20, 22, 23, 29, 30, 21, 24,
    25, 18, 32, 12, 28, 0, 1, 2, 3, 4, 6, 7, 8, 5, 9, 10, 11, 19, 20, 22, 23, 29,
    30, 21, 24, 13, 32, 26, 25, 12, 17, 28, 0, 1, 2, 3, 5, 6, 7, 8, 4, 9, 14, 15,
    19, 20, 22, 23, 29, 30, 21, 16, 25, 17, 32, 27, 12, 13, 28, 0, 1, 2, 3, 4, 5,
    6, 7, 8, 9, 10, 11, 14, 15, 19, 20, 22, 23, 29, 30, 21, 16, 32, 27, 26, 25, 12,
    18, 17, 13, 28, 0, 1, 2, 3, 9, 4, 5, 6, 7, 8, 33, 34, 32, 27, 0, 1, 2, 3, 4,
    9, 5, 6, 7, 8, 10, 11, 33, 34, 32, 27, 12, 13, 28, 0, 1, 2, 3, 5, 9, 4, 6, 7,
    8, 14, 15, 33, 34, 32, 16, 17, 24, 0, 1, 2, 3, 4, 5, 9, 6, 7, 8, 10, 11, 14,
    15, 33, 34, 32, 16, 12, 13, 24, 18, 28, 0, 1, 2, 3, 6, 9, 4, 5, 7, 8, 19, 20,
    33, 34, 21, 27, 18, 31, 0, 1, 2, 3, 4, 6, 9, 5, 7, 8, 10, 11, 19, 20, 33, 34,
    21, 27, 12, 13, 31, 17, 28, 0, 1, 2, 3, 5, 6, 9, 4, 7, 8, 14, 15, 19, 20, 33,
    34, 21, 16, 17, 31, 24, 13, 0, 1, 2, 3, 4, 5, 6, 9, 7, 8, 10, 11, 14, 15, 19,
    20, 33, 34, 21, 16, 12, 31, 24, 18, 17, 13, 0, 1, 2, 3, 7, 9, 4, 5, 6, 8, 22,
    23, 33, 34, 32, 24, 25, 16, 0, 1, 2, 3, 4, 7, 9, 5, 6, 8, 10, 11, 22, 23, 33,
    34, 32, 24, 12, 13, 16, 26, 28, 0, 1, 2, 3, 5, 7, 9, 4, 6, 8, 14, 15, 22, 23,
    33, 34, 32, 25, 17, 27, 24, 16, 0, 1, 2, 3, 4, 5, 7, 9, 6, 8, 10, 11, 14, 15,
    22, 23, 33, 34, 32, 12, 13, 27, 24, 16, 26, 18, 28, 0, 1, 2, 3, 6, 7, 9, 4, 5,
    8, 19, 20, 22, 23, 33, 34, 21, 24, 25, 18, 31, 16, 28, 0, 1, 2, 3, 4, 6, 7, 9,
    5, 8, 10, 11, 19, 20, 22, 23, 33, 34, 21, 24, 12, 13, 31, 16, 26, 17, 28, 0, 1,
    2, 3, 5, 6, 7, 9, 4, 8, 14, 15, 19, 20, 22, 23, 33, 34, 21, 25, 17, 31, 27,
    24, 16, 13, 28, 0, 1, 2, 3, 4, 5, 6, 7, 9, 8, 10, 11, 14, 15, 19, 20, 22, 23,
    33, 34, 21, 12, 31, 27, 24, 16, 26, 18, 17, 13, 28, 0, 1, 2, 3, 8, 9, 4, 5, 6,
    7, 29, 30, 33, 34, 31, 27, 26, 21, 0, 1, 2, 3, 4, 8, 9, 5, 6, 7, 10, 11, 29,
    30, 33, 34, 31, 27, 12, 13, 21, 25, 28, 0, 1, 2, 3, 5, 8, 9, 4, 6, 7, 14, 15,
    29, 30, 33, 34, 31, 16, 26, 17, 21, 24, 28, 0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10,
    11, 14, 15, 29, 30, 33, 34, 31, 16, 12, 13, 21, 24, 25, 18, 28, 0, 1, 2, 3, 6,
    8, 9, 4, 5, 7, 19, 20, 29, 30, 33, 34, 27, 26, 18, 32, 31, 21, 0, 1, 2, 3, 4,
    6, 8, 9, 5, 7, 10, 11, 19, 20, 29, 30, 33, 34, 27, 12, 13, 32, 31, 21, 25, 17,
    28, 0, 1, 2, 3, 5, 6, 8, 9, 4, 7, 14, 15, 19, 20, 29, 30, 33, 34, 16, 26, 17,
    32, 31, 21, 24, 13, 28, 0, 1, 2, 3, 4, 5, 6, 8, 9, 7, 10, 11, 14, 15, 19, 20,
    29, 30, 33, 34, 16, 12, 32, 31, 21, 24, 25, 18, 17, 13, 28, 0, 1, 2, 3, 7, 8,
    9, 4, 5, 6, 22, 23, 29, 30, 33, 34, 31, 24, 25, 21, 16, 12, 0, 1, 2, 3, 4, 7,
    8, 9, 5, 6, 10, 11, 22, 23, 29, 30, 33, 34, 31, 24, 13, 21, 16, 26, 25, 12, 0,
    1, 2, 3, 5, 7, 8, 9, 4, 6, 14, 15, 22, 23, 29, 30, 33, 34, 31, 25, 17, 21, 27,
    24, 16, 12, 28, 0, 1, 2, 3, 4, 5, 7, 8, 9, 6, 10, 11, 14, 15, 22, 23, 29, 30,
    33, 34, 31, 13, 21, 27, 24, 16, 26, 25, 12, 18, 28, 0, 1, 2, 3, 6, 7, 8, 9, 4,
    5, 19, 20, 22, 23, 29, 30, 33, 34, 24, 25, 18, 32, 31, 21, 16, 12, 28, 0, 1, 2,
    3, 4, 6, 7, 8, 9, 5, 10, 11, 19, 20, 22, 23, 29, 30, 33, 34, 24, 13, 32, 31,
    21, 16, 26, 25, 12, 17, 28, 0, 1, 2, 3, 5, 6, 7, 8, 9, 4, 14, 15, 19, 20, 22,
    23, 29, 30, 33, 34, 25, 17, 32, 31, 21, 27, 24, 16, 12, 13, 28, 0, 1, 2, 3, 4,
    5, 6, 7, 8, 9, 10, 11, 14, 15, 19, 20, 22, 23, 29, 30, 33, 34, 32, 31, 21, 27,
    24, 16, 26, 25, 12, 18, 17, 13, 28
};

/// Column index for compressed row format storage.
const unsigned char P2_prolongation_col_ind[116]= {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 4, 0, 1, 4, 0, 1, 4, 7,
    8, 0, 1, 4, 5, 6, 0, 2, 5, 0, 2, 5, 0, 2, 5, 7, 9, 0, 2, 4,
    5, 6, 1, 2, 4, 5, 6, 1, 2, 6, 1, 2, 6, 1, 2, 6, 8, 9, 0, 3,
    7, 0, 3, 7, 0, 3, 5, 7, 9, 0, 3, 4, 7, 8, 1, 3, 4, 7, 8, 2,
    3, 5, 7, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 3, 8, 1, 3, 8,
    1, 3, 6, 8, 9, 2, 3, 6, 8, 9, 2, 3, 9, 2, 3, 9
};

/// Beginning of rows for compressed row format storage.
const unsigned char P2_prolongation_row_beg[36]= {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 16, 21, 26, 29, 32, 37, 42,
    47, 50, 53, 58, 61, 64, 69, 74, 79, 84, 94, 97, 100, 105, 110, 113, 116
};

/// Coefficients of the prolongation matrices as index into coeff.
const unsigned char P2_prolongation_coeff_idx[116]= {
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 0, 4, 0, 2, 4, 0, 0, 1, 3, 3,
    0, 0, 1, 3, 3, 2, 0, 4, 0, 2, 4, 0, 0, 1, 3, 3, 0, 0, 3, 1, 3, 0, 0, 3, 3, 1,
    2, 0, 4, 0, 2, 4, 0, 0, 1, 3, 3, 2, 0, 4, 0, 2, 4, 0, 0, 3, 1, 3, 0, 0, 3, 1,
    3, 0, 0, 3, 3, 1, 0, 0, 3, 3, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 0, 4, 0, 2,
    4, 0, 0, 3, 1, 3, 0, 0, 3, 3, 1, 2, 0, 4, 0, 2, 4
};

/// The coefficients in the prolongation matrices.
const double P2_prolongation_coeff[6]= {
    -0.125, 0.25, 0.375, 0.5, 0.75, 1.0
};
/// \}

#ifdef _PAR
template <class ValueT>
void ProlongationCL<ValueT>::Restrict(const TetraCL* sit, const VecDescCL& fine, const bool doinjection) const
{
    LocalP2CL<ValueT> lp2; // for ghost - master communication

    if (sit->HasGhost()) return;
    if (!sit->IsGhost()) return;

    const IdxDescCL& fIdx = *fine.RowIdx;

    const Uint f_idx= fIdx.GetIdx();

    // sit->IsInTriang( f_level) != sit->IsInTriang( c_level)
    // Collect the fine indices.
    const std::vector<IdxT> fUnknowns( CollectChildUnknownsP2( *sit, f_idx));

    // Here, the green rule for regular refinement and the
    // regular rule are treated in the same way: & 63 cuts off
    // the bit that distinguishes them.
    const Uint rule= sit->GetRefRule() & 63;
    for (Uint i= P2_local_prolongation_mat_beg[rule]; // Construct the rows into mat.
         i != P2_local_prolongation_mat_beg[rule+1]; ++i) {
        const IdxT thefUnknown= fUnknowns[i-P2_local_prolongation_mat_beg[rule]];
        if (thefUnknown != NoIdx) {
            const Uint row= P2_local_prolongation_rows[i];
            for (Uint j= P2_prolongation_row_beg[row]; // Construct a single row into mat.
                 j != P2_prolongation_row_beg[row+1]; ++j) {
                const value_type f_val = DoFHelperCL<value_type, VectorCL>::get( fine.Data, thefUnknown);
                const value_type sca_val = DoFHelperCL<value_type, VectorCL>::get( scale_, thefUnknown);
                const double wheight_val = P2_prolongation_coeff[P2_prolongation_coeff_idx[j]];

                if (doinjection && wheight_val == 1.0)
                    lp2[P2_prolongation_col_ind[j]] += wheight_val * sca_val * f_val;
                if (!doinjection)
                    lp2[P2_prolongation_col_ind[j]] += wheight_val * sca_val * f_val;

            }
        }
    }
    data_[&*sit] = lp2;
}

template <class ValueT>
void ProlongationCL<ValueT>::Prolong(const TetraCL* sit, const LocalP2CL<ValueT>& coarse, VecDescCL& fine) const
{
    const IdxDescCL& fIdx = *fine.RowIdx;
    const Uint f_idx= fIdx.GetIdx();

    // Collect the fine indices.
    const std::vector<IdxT> fUnknowns( CollectChildUnknownsP2( *sit, f_idx));

    // Here, the green rule for regular refinement and the
    // regular rule are treated in the same way: & 63 cuts off
    // the bit that distinguishes them.
    const Uint rule= sit->GetRefRule() & 63;
    for (Uint i= P2_local_prolongation_mat_beg[rule]; // Construct the rows into mat.
         i != P2_local_prolongation_mat_beg[rule+1]; ++i) {
        const IdxT thefUnknown= fUnknowns[i-P2_local_prolongation_mat_beg[rule]];
        if (thefUnknown != NoIdx) {
            const Uint row= P2_local_prolongation_rows[i];
            for (Uint j= P2_prolongation_row_beg[row]; // Construct a single row into mat.
                 j != P2_prolongation_row_beg[row+1]; ++j) {
                const IdxT thecUnknown= P2_prolongation_col_ind[j];
                const value_type f_val = DoFHelperCL<value_type, VectorCL>::get( fine.Data, thefUnknown);
                const value_type sca_val = DoFHelperCL<value_type, VectorCL>::get( scale_, thefUnknown);
                const double wheight_val = P2_prolongation_coeff[P2_prolongation_coeff_idx[j]];
                const value_type tmp = f_val + wheight_val * sca_val * coarse[thecUnknown];
                DoFHelperCL<value_type, VectorCL>::set( fine.Data, thefUnknown, tmp);
            }
        }
    }
}
#endif

template <class ValueT>
void ProlongationCL<ValueT>::BuildP1ProlongationMatrix( const IdxDescCL& cIdx, const IdxDescCL& fIdx)
{
    const Uint c_level= cIdx.TriangLevel();
    const Uint c_idx= cIdx.GetIdx();
    const Uint f_idx= fIdx.GetIdx();
    MatrixBuilderCL mat( &prolongation_, fIdx.NumUnknowns(), cIdx.NumUnknowns());
    IdxT i;

    // setup index part of matrix
    // Iterate over all edges, interpolate values on new mid vertices
    for (MultiGridCL::const_EdgeIterator sit= mg_.GetAllEdgeBegin( c_level),
         theend= mg_.GetAllEdgeEnd( c_level); sit!=theend; ++sit)
        if ( sit->IsRefined() && sit->GetMidVertex()->Unknowns.Exist()
             && sit->GetMidVertex()->Unknowns.Exist(f_idx)
             && !sit->GetMidVertex()->Unknowns.Exist(c_idx) )  {
            // only new non-boundary vertices are interpolated
// if(!(*sit->GetMidVertex()->Unknowns.Exist(idx)) std::cout << "no unknown index in mid vertex" << std::endl;
            i= sit->GetMidVertex()->Unknowns(f_idx);
            if (sit->GetVertex(0)->Unknowns.Exist() && sit->GetVertex(0)->Unknowns.Exist(c_idx))
                mat(i,sit->GetVertex(0)->Unknowns(c_idx))= 0.5;
            if (sit->GetVertex(1)->Unknowns.Exist() && sit->GetVertex(1)->Unknowns.Exist(c_idx))
                mat(i,sit->GetVertex(1)->Unknowns(c_idx))= 0.5;
        }
    // Iterate over the vertices of the coarse triangulation and copy values
    for (MultiGridCL::const_TriangVertexIteratorCL sit= mg_.GetTriangVertexBegin( c_level),
         theend= mg_.GetTriangVertexEnd( c_level); sit!=theend; ++sit)
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(c_idx) ) {
            mat( sit->Unknowns(f_idx), sit->Unknowns(c_idx) )= 1.0;
        }

    mat.Build();
}

template <class ValueT>
void ProlongationCL<ValueT>::BuildP2ProlongationMatrix( const IdxDescCL& cIdx, const IdxDescCL& fIdx)
{
    const Uint c_level= cIdx.TriangLevel();
    const Uint f_level= fIdx.TriangLevel();
    const Uint c_idx= cIdx.GetIdx();
    const Uint f_idx= fIdx.GetIdx();
// It is only a copy of the 1D case 'ndofs' times
    const Uint ndofs= fIdx.NumUnknownsVertex();

    MatrixBuilderCL mat( &prolongation_, fIdx.NumUnknowns(), cIdx.NumUnknowns());
    MatrixBuilderCL mat_inj( &injection_, fIdx.NumUnknowns(), cIdx.NumUnknowns());

    DROPS_FOR_TRIANG_CONST_TETRA(mg_, c_level, sit){

#ifdef _PAR
        if (sit->HasGhost()) continue;
#endif

        if (!sit->IsInTriang( f_level)) {
            IdxT cUnknowns[10];
            // Collect the coarse indices.
            for (Uint i=0; i<4; ++i)
                cUnknowns[i]= (sit->GetVertex( i)->Unknowns.Exist()
                               && sit->GetVertex( i)->Unknowns.Exist( c_idx))
                              ? sit->GetVertex( i)->Unknowns( c_idx) : NoIdx;
            for (Uint i=0; i<6; ++i)
                cUnknowns[i+4]= (sit->GetEdge( i)->Unknowns.Exist()
                                 && sit->GetEdge( i)->Unknowns.Exist( c_idx))
                                ? sit->GetEdge( i)->Unknowns( c_idx) : NoIdx;
            // Collect the fine indices.
            const std::vector<IdxT> fUnknowns( CollectChildUnknownsP2( *sit, f_idx));

            // Here, the green rule for regular refinement and the
            // regular rule are treated in the same way: & 63 cuts off
            // the bit that distinguishes them.
            const Uint rule= sit->GetRefRule() & 63;
            for (Uint i= P2_local_prolongation_mat_beg[rule]; // Construct the rows into mat.
                 i != P2_local_prolongation_mat_beg[rule+1]; ++i) {
                const IdxT thefUnknown= fUnknowns[i-P2_local_prolongation_mat_beg[rule]];
                if (thefUnknown != NoIdx) {
                    const Uint row= P2_local_prolongation_rows[i];
                    for (Uint j= P2_prolongation_row_beg[row]; // Construct a single row into mat.
                         j != P2_prolongation_row_beg[row+1]; ++j) {

                        const double wheight_val = P2_prolongation_coeff[P2_prolongation_coeff_idx[j]];

                        const IdxT thecUnknown= cUnknowns[P2_prolongation_col_ind[j]];
                        if (thecUnknown != NoIdx){
                            if (wheight_val ==1)
                                for (Uint k=0; k<ndofs; k++)
                                    mat_inj(thefUnknown+k, thecUnknown+k) += wheight_val * scale_[thefUnknown+k];
                            for (Uint k=0; k<ndofs; k++)
                                mat(thefUnknown+k, thecUnknown+k) += wheight_val * scale_[thefUnknown+k];
                        }
                    }
                }
            }
        }
        else { // coarse and fine tetra are identical; the restriction is trivial.
            for (Uint i=0; i<4; ++i)
                if (sit->GetVertex( i)->Unknowns.Exist()
                    && sit->GetVertex( i)->Unknowns.Exist( f_idx))
                    for (Uint k=0; k<ndofs; k++){
                        mat( sit->GetVertex( i)->Unknowns( f_idx)+k, sit->GetVertex( i)->Unknowns( c_idx)+k) += scale_[sit->GetVertex( i)->Unknowns( f_idx)+k];
                        mat_inj( sit->GetVertex( i)->Unknowns( f_idx)+k, sit->GetVertex( i)->Unknowns( c_idx)+k) += scale_[sit->GetVertex( i)->Unknowns( f_idx)+k];
                    }
            for (Uint i=0; i<6; ++i)
                if (sit->GetEdge( i)->Unknowns.Exist()
                    && sit->GetEdge( i)->Unknowns.Exist( f_idx))
                    for (Uint k=0; k<ndofs; k++){
                        mat( sit->GetEdge( i)->Unknowns( f_idx)+k, sit->GetEdge( i)->Unknowns( c_idx)+k) += scale_[sit->GetEdge( i)->Unknowns( f_idx)+k];
                        mat_inj( sit->GetEdge( i)->Unknowns( f_idx)+k, sit->GetEdge( i)->Unknowns( c_idx)+k) += scale_[sit->GetEdge( i)->Unknowns( f_idx)+k];
                    }
        }
    }
    mat.Build();
    mat_inj.Build();
}



/******************************************************************************
* H A N D L E R  P A R E N T D A T A   C L                                    *
******************************************************************************/
#ifdef _PAR

/// \brief Handler for generating information about parent data on other processes
template <class ValueT>
class ProlongationCL<ValueT>::HandlerParentDataCL
{
  private:
    DataT& parent_data_;

  public:
    HandlerParentDataCL( DataT& data)
        : parent_data_(data) {}
    /// \brief Gather information about distributed dof on sender proc
    bool Gather( DiST::TransferableCL&, DiST::SendStreamCL&);
    /// \brief Scatter information about distributed dof on sender proc
    bool Scatter( DiST::TransferableCL&, const size_t&, DiST::MPIistreamCL&);
};

template <class ValueT>
bool ProlongationCL<ValueT>::HandlerParentDataCL::Gather(
    DiST::TransferableCL& t, DiST::SendStreamCL& send)
{
    TetraCL* tetra;
    simplex_cast(t, tetra);
    if (parent_data_.count( tetra) == 1){
        const LocalP2CL<value_type> & tmp = parent_data_[&*tetra];
        for (size_t k=0; k< 10; ++k) // write LocalPiCL
            send << tmp[k];
        parent_data_.erase(tetra);
        return true;
    }
    return false;
}

template <class ValueT>
bool ProlongationCL<ValueT>::HandlerParentDataCL::Scatter(
    DiST::TransferableCL& t, __UNUSED__ const size_t& numData,
    DiST::MPIistreamCL& recv)
{
    Assert(numData == 1, DROPSErrCL("HandlerParentDataCL::Scatter: more than one ghost exists"), DebugDiSTC);

    TetraCL* tetra;
    simplex_cast(t, tetra);
    LocalP2CL<value_type> lp2;
    for (size_t k=0; k< 10; ++k) // read LocalPiCL
        recv >> lp2[k];
    parent_data_[&*tetra]= lp2;
    return true;
}

#endif

#ifdef _PAR
template <class ValueT>
VectorCL ProlongationCL<ValueT>::operator* (const VectorCL& vec) const {
    VecDescCL old(coarse_);
    old.Data = vec;
    VecDescCL neu(fine_);
    neu.Data = 0.0;
    const Uint lvl= coarse_->TriangLevel();
    const Uint idx= coarse_->GetIdx();

    typedef VecDescCL::DataType VecT;
    typedef DoFHelperCL<value_type, VecT> DoFT;

    DROPS_FOR_TRIANG_CONST_TETRA( mg_, lvl, it) {
        if (!it->HasGhost()) continue;
        LocalP2CL<ValueT> lp2; // for master - ghost communication

        for (Uint i= 0; i< NumVertsC; ++i)
            if (it->GetVertex( i)->Unknowns.Exist() && it->GetVertex( i)->Unknowns.Exist(idx))
                lp2[i]= DoFT::get( vec, it->GetVertex( i)->Unknowns( idx));
            else
                lp2[i]= value_type(0.0);
        for (Uint i= 0; i< NumEdgesC; ++i)
            if ( it->GetEdge( i)->Unknowns.Exist() &&  it->GetEdge( i)->Unknowns.Exist(idx))
                lp2[i+NumVertsC]= DoFT::get( vec, it->GetEdge( i)->Unknowns( idx));
            else
                lp2[i+NumVertsC]=  value_type(0.0);
        data_[&*it] = lp2;
    }

    if (prolongation_.num_nonzeros() != 0)
        neu.Data = prolongation_* old.Data;

    HandlerParentDataCL handler( data_);
    DiST::InterfaceCL* interf;
    // init the interface
    DiST::PrioListT from, to;
    from.push_back( PrioMaster);
    to.push_back( PrioGhost);
    DiST::InterfaceCL::DimListT dimlist;
    dimlist.push_back( DiST::GetDim<TetraCL>());
    DiST::LevelListCL lvls( lvl); // levels 0,..,TriangLevel

    interf = new DiST::InterfaceCL( lvls, from, to, dimlist, true);
    interf->Communicate( handler);
    delete interf;

    for (typename DataT::const_iterator it = data_.begin(); it != data_.end(); ++it) {
        Prolong(it->first, it->second, neu);
    }

    data_.clear();
    fine_->GetEx().Accumulate( neu.Data);
    return neu.Data;
}
#else
template <class ValueT>
VectorCL ProlongationCL<ValueT>::operator* (const VectorCL& vec) const {
    return prolongation_*vec;
}
#endif

#ifdef _PAR
template <class ValueT>
VectorCL ProlongationCL<ValueT>::mytransp_mul(const VectorCL& vec, const bool doinjection) const
{

    VecDescCL old(fine_);
    old.Data = vec;
    VecDescCL neu(coarse_);

    Uint lvl= coarse_->TriangLevel();

    // restriction in MG-cycle: fine data is distributed, result is distributed
    fine_->GetEx().Accumulate(old.Data);

    for (MultiGridCL::const_TetraIterator sit= mg_.GetAllTetraBegin( lvl),
         theend= mg_.GetAllTetraEnd( lvl); sit != theend; ++sit){
        if (!sit->IsInTriang(lvl)) continue;
        Restrict(&*sit, old, doinjection);
    }

    if (doinjection){
        if (injection_.num_nonzeros() != 0)
            neu.Data = transp_mul(injection_, old.Data);
    }
    else{
        if (prolongation_.num_nonzeros() != 0)
            neu.Data = transp_mul(prolongation_, old.Data);
    }

    HandlerParentDataCL handler( data_);
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

    for (typename DataT::const_iterator it = data_.begin(); it != data_.end(); ++it) {
        const TetraCL* tet = it->first;
        for (Uint v=0; v<4; ++v)
            if (tet->GetVertex(v)->Unknowns.Exist() && tet->GetVertex(v)->Unknowns.Exist(coarse_->GetIdx())){
                const IdxT pos = tet->GetVertex(v)->Unknowns(coarse_->GetIdx());
                value_type tmp= DoFHelperCL<value_type, VectorCL>::get( neu.Data, pos);
                tmp += data_[tet][v];
                DoFHelperCL<value_type, VectorCL>::set( neu.Data, pos, tmp);
            }
        for (Uint e=0; e<6; ++e)
            if (tet->GetEdge(e)->Unknowns.Exist() && tet->GetEdge(e)->Unknowns.Exist(coarse_->GetIdx())){
                const IdxT pos = tet->GetEdge(e)->Unknowns(coarse_->GetIdx());
                value_type tmp= DoFHelperCL<value_type, VectorCL>::get( neu.Data, pos);
                tmp += data_[tet][e+4];
                DoFHelperCL<value_type, VectorCL>::set( neu.Data, pos, tmp);
            }
    }

    data_.clear();
    return neu.Data;
}
#else
template <class ValueT>
VectorCL ProlongationCL<ValueT>::mytransp_mul(const VectorCL& vec, const bool doinjection) const
{
    if (doinjection)
        return transp_mul(injection_, vec);
    else
        return transp_mul(prolongation_,vec);
}
#endif

} // end of namespace DROPS
