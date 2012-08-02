/// \file combinatorialcut.cpp
/// \brief tests the PrincipalLattice-class
/// \author LNM RWTH Aachen: Joerg Grande; SC RWTH Aachen:

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
 * Copyright 2011 LNM/SC RWTH Aachen, Germany
*/

#include "geom/principallattice.h"
#include "geom/subtriangulation.h"
#include "num/quadrature.h"
#include "misc/container.h"
#include "num/discretize.h"
#include "num/lattice-eval.h"
#include "geom/multigrid.h"

#include <iostream>
#include <sstream>
#include <tr1/unordered_map>

namespace DROPS {

namespace RefEdge {

const Uint NumVertsC= 2;

/// Vertices of a given edge
inline Ubyte VertOfEdge (Ubyte, Ubyte num) { return num; }

/// Vertex opposing a given vertex
/// Note, that this is different from the n-dimensional cases, n > 2, as vertexes and facets are the same in 1d. Therefore, they should not be numbered differently. This is a natural ordering, which also could have been used for tetras. However, changing this now would be painful.
inline Ubyte OppVert (Ubyte vert) { return 1 - vert; }

} // end of namespace DROPS::RefEdge

namespace RefTri {

const Uint NumVertsC= 3;
const Uint NumEdgesC= 3;

/// Vertices of a given edge
const Ubyte vert_of_edge_ar[NumEdgesC][2]= {
    {0,1}, {0,2}, {1,2}
};
inline Ubyte VertOfEdge (Ubyte edge, Ubyte num) { return vert_of_edge_ar[edge][num]; }

/// Edge opposing a given vertex
/// Note, that this is different from the higher-dimensional cases as edges and facets are the same in 2d. Therefore, they should not be numbered differently. This is a natural ordering, which also could have been used for tetras. However, changing this now would be painful.
inline Ubyte OppEdge (Ubyte vert) { return 2 - vert; }

} // end of namespace DROPS::RefTri

namespace RefPenta {

const Uint NumVertsC= 5;
const Uint NumEdgesC= 10;
const Uint NumTrisC=  10;
const Uint NumTetrasC= 5;

/// Vertices of a given edge
const Ubyte vert_of_edge_ar[NumEdgesC][2]= {
    {0,1}, {0,2}, {1,2}, {0,3}, {1,3}, {2,3}, // The first six edges are from the reference-tetra
    {0,4}, {1,4}, {2,4}, {3,4}
};
inline Ubyte VertOfEdge (Ubyte edge, Ubyte num) { return vert_of_edge_ar[edge][num]; }

/// Vertices of a given tetra
const Ubyte vert_of_tetra_ar[NumVertsC][4]= {
    {1,2,3,4}, {0,2,3,4}, {0,1,3,4}, {0,1,2,4}, {0,1,2,3}
};
inline Ubyte VertOfTetra (Ubyte tetra, Ubyte num) { return vert_of_tetra_ar[tetra][num]; }

/// Tetra opposing a given vertex
inline Ubyte OppTetra (Ubyte vert) { return vert; }

} // end of namespace DROPS::RefPenta

template <Uint Dim>
inline Ubyte VertOfEdge (Ubyte edge, Ubyte num); // not defined, but full specializations for Dim=1..4 are.

template <>
inline Ubyte VertOfEdge<1> (Ubyte, Ubyte num)
{ return RefEdge::VertOfEdge( 0, num); }

template <>
inline Ubyte VertOfEdge<2> (Ubyte edge, Ubyte num)
{ return RefTri::VertOfEdge( edge, num); }

template <>
inline Ubyte VertOfEdge<3> (Ubyte edge, Ubyte num)
{ return VertOfEdge( edge, num); }

template <>
inline Ubyte VertOfEdge<4> (Ubyte edge, Ubyte num)
{ return RefPenta::VertOfEdge( edge, num); }


template <Uint Dim= 3>
class SignTraitCL;

template <Uint Dim>
std::ostream&
operator<< (std::ostream& out, const SignTraitCL<Dim>& c);

template <Uint Dim>
class SignTraitCL
{
  private:
    enum { NumVerts= Dim + 1,
           NumEdges= (Dim*(Dim + 1))/2,
// dimension vs. max. number of roots:
// 1 2, the 0-pattern
// 2 3, the 0-pattern
// 3 4, 2 vertexes on each side (and 0-pattern)
// 4 6, 2 verts on one side, three on the other: 10 - 1 - 3
// in general:
// d < 3: see above,
// else: distribute vertexes as evenly as possible on both sides:
//       d odd:  (d+1 over 2) - 2*((d+1)/2 over 2)              = ((d+1)/2)^2
//       d even: (d+1 over 2) - (d/2 over 2) - (d/2 + 1 over 2) = d*(d+2)/4
           max_num_root= Dim < 4 ? Dim + 1 : (Dim & 1u ? NumVerts*NumVerts : Dim*(Dim + 2))/4
    };

    Ubyte num_root_vert_; ///< number of vertices, where the level set function is zero.
    Ubyte num_root_;      ///< number of roots of the level set function; invariant: num_root_vert <= num_root <= max_num_root.
    byte sign_[NumVerts]; ///< Sign of the level set function of the vertices; \f$\in\{-1,0,1\}\f$

    ///\brief local number with respect to the reference simplex of the object on the cut: [0,num_root_vert): vertex numbers in (0..NumVerts-1); [num_root_vert, num_root): edge numbers in (0..NumEdges-1). Both parts are sorted in increasing order.
    Ubyte cut_simplex_    [max_num_root];
    Ubyte cut_simplex_rep_[max_num_root]; ///< local number of the object on the cut: (0..NumVerts+NumEdges-1)

    void compute_cuts ();

  public:
    SignTraitCL () : num_root_vert_( NumVerts), num_root_( NumVerts) {} ///< Uninitialized default state
    SignTraitCL (const byte   ls[NumVerts]) { assign( ls); }            ///< Assign the sign pattern on the vertices.
    SignTraitCL (const double ls[NumVerts]) { assign( ls); }            ///< Assign the sign pattern on the vertices.
    void assign (const byte   ls[NumVerts]); ///< Assign a sign pattern on the vertices.
    void assign (const double ls[NumVerts]); ///< Assign a sign pattern on the vertices.

    bool empty () const { return num_root_ == 0; } ///< True, iff there is no intersection.

    bool has_codim_le_1 () const { return num_root_ >= Dim; }        ///< True, iff the intersection has non-zero (Dim-1)-dimensional measure.
    bool has_codim_0 () const { return num_root_vert_ == NumVerts; } ///< True, exactly for the zero-pattern.
    bool no_zero_vertex () const { return num_root_vert_ == 0; }     ///< True, iff there is no vertex, in which ls vanishes.

    Ubyte num_cut_simplexes () const { return num_root_; }      ///< Number of edges and vertices with a root of ls.
    Ubyte num_zero_vertexes () const { return num_root_vert_; } ///< Number of vertices of the simplex that are roots of ls.

    byte sign (int i) const { return sign_[i]; } ///< -1,0,1; sign of vertex i.

    /// Return local number of edges/verts with a root of ls. For edges, [] returns a edge number in 0..NumEdges-1, and () returns an extended vertex number in NumVerts..NumVerts+NumEdges-1.
    ///@{
    Ubyte operator[] (int i) const { return cut_simplex_[i]; }
    Ubyte operator() (int i) const { return cut_simplex_rep_[i]; }
    ///@}
    friend std::ostream& operator<< <Dim>(std::ostream&, const SignTraitCL<Dim>&); ///< Debug-output to a stream (dumps all members)
};

template <Uint Dim>
void
SignTraitCL<Dim>::compute_cuts ()
{
    for (Ubyte i= 0; i < NumVerts; ++i)
        if (sign( i) == 0)
            cut_simplex_[num_root_vert_++]= i;
    num_root_= num_root_vert_;
    for (Ubyte i= 0; i < NumEdges; ++i)
        if (sign( VertOfEdge<Dim>( i, 0))*sign( VertOfEdge<Dim>( i, 1)) == -1)
            cut_simplex_[num_root_++]= i;
    std::memcpy( cut_simplex_rep_, cut_simplex_, num_root_*sizeof(byte));
    for (int i= num_root_vert_; i < num_root_; ++i)
        cut_simplex_rep_[i]+= NumVerts;
}

template <Uint Dim>
void
SignTraitCL<Dim>::assign (const byte ls[NumVerts])
{
    num_root_vert_= num_root_= 0;

    byte sum= 0;
    for (Ubyte i= 0; i < NumVerts; ++i)
        sum+= (sign_[i]= ls[i]);
    if (std::abs( sum) == static_cast<int>( NumVerts)) // optimize the case of uncut tetras
        return;

    compute_cuts ();
}

template <Uint Dim>
inline void
SignTraitCL<Dim>::assign (const double ls[NumVerts])
{
    byte ls_byte[NumVerts];
    std::transform( ls + 0, ls + NumVerts, ls_byte + 0, DROPS::sign);
    this->assign( ls_byte);

}

template <Uint Dim>
std::ostream&
operator<< (std::ostream& out, const SignTraitCL<Dim>& c)
{
    const int NumVerts= Dim + 1;
    out << static_cast<int>( c.num_root_vert_) << ' ' << static_cast<int>( c.num_root_) << '\n';
    for (int i= 0; i < NumVerts; ++i)
        out << static_cast<int>( c.sign_[i]) << ' ';
    out << '\n';
    for (int i= 0; i < NumVerts; ++i)
        out << static_cast<int>( c.cut_simplex_[i]) << ' ';
    out << '\n';
    for (int i= 0; i < NumVerts; ++i)
        out << static_cast<int>( c.cut_simplex_rep_[i]) << ' ';
    return out << '\n';
}

typedef SignTraitCL<3> TetraSignTraitCL;

} // end of namespace DROPS


template <DROPS::Uint Dim>
void
write_sign_trait (const DROPS::SignTraitCL<Dim>& c, std::ostream& os)
{
    os << "signs: ";
    for (DROPS::Uint i= 0; i < Dim+1; ++i)
        os << (int) c.sign( i) << ' ';
    os << '\n';
    os << (int) c.empty() << ' ' << (int) c.has_codim_le_1() << ' ' << (int) c.has_codim_0() << ' '
       << (int) c.no_zero_vertex() << ' ' <<  (int) c.num_cut_simplexes()  << ' ' << (int) c.num_zero_vertexes () << '\n';
    os << "cut_simplexes: ";
    for (DROPS::Uint i= 0; i < c.num_cut_simplexes(); ++i)
        os << (int) c[i] << ' ' << (int) c(i) << "  ";
    os << '\n';
}

void
write_sign_pattern_trait (const DROPS::SignPatternTraitCL& c, std::ostream& os)
{
    os << "signs: ";
    for (DROPS::Uint i= 0; i < 4; ++i)
        os << (int) c.sign( i) << ' ';
    os << '\n';
    os << (int) c.empty() << ' ' << (int) c.is_2d() << ' ' << (int) c.is_3d() << ' '
       << (int) c.no_zero_vertex() << ' ' <<  (int) c.num_cut_simplexes()  << ' ' << (int) c.num_zero_vertexes () << '\n';
    os << "cut_simplexes: ";
    for (DROPS::Uint i= 0; i < c.num_cut_simplexes(); ++i)
        os << (int) c[i] << ' ' << (int) c(i) << "  ";
    os << '\n';
}

void write_sign_traits_1_2_3_4 ()
{
    std::ofstream f1( "sign_trait1.txt");
    std::ofstream f2( "sign_trait2.txt");
    std::ofstream f3( "sign_trait3.txt");
    std::ofstream f4( "sign_trait4.txt");
    DROPS::byte ls[5];
    int c= 0;
    for (ls[0]= -1; ls[0] <= 1; ++ls[0]) {
      for (ls[1]= -1; ls[1] <= 1; ++ls[1]) {
        write_sign_trait( DROPS::SignTraitCL<1>( ls), f1);
        for (ls[2]= -1; ls[2] <= 1; ++ls[2]) {
          write_sign_trait( DROPS::SignTraitCL<2>( ls), f2);
          for (ls[3]= -1; ls[3] <= 1; ++ls[3]) {
            write_sign_trait( DROPS::SignTraitCL<3>( ls), f3);
            for (ls[4]= -1; ls[4] <= 1; ++ls[4], c++) {
                std::cout << "c: " << c << " ls: " << (int) ls[0] << ' ' << (int) ls[1]
                          << ' ' << (int) ls[2] << ' ' << (int) ls[3]  << ' ' << (int) ls[4]<< std::endl;
                write_sign_trait( DROPS::SignTraitCL<4>( ls), f4);
            }
          }
        }
      }
    }

}
void test_sign_traits_cut ()
{
    std::ofstream f0( "sign_pattern_trait.txt");
    std::ofstream f1( "sign_trait_3.txt");
    DROPS::byte ls[4];
    int c= 0;
    for (int i= -1; i <= 1; ++i)
      for (int j= -1; j <= 1; ++j)
        for (int k= -1; k <= 1; ++k)
          for (int l= -1; l <= 1; ++l, ++c) {
              ls[0]= i; ls[1]= j; ls[2]= k; ls[3]= l;
              std::cout << "c: " << c << " ls: " << (int) ls[0] << ' ' << (int) ls[1] << ' ' << (int) ls[2] << ' ' << (int) ls[3] << std::endl;
              write_sign_pattern_trait( DROPS::SignPatternTraitCL( ls), f0);
              write_sign_trait(         DROPS::SignTraitCL<3>(     ls), f1);
          }
}

void test_tetra_cut ()
{
    DROPS::GridFunctionCL<> ls( 4);
    ls[0]= -1.; ls[1]= 0.; ls[2]= 0.; ls[3]= 0.;
    DROPS::TetraPartitionCL tet;
    // tet.partition_principal_lattice<DROPS::SortedVertexPolicyCL, DROPS::MergeCutPolicyCL> ( 1, ls);
    // std::cerr << tet;
    int c= 0;
    for (int i= -1; i <= 1; ++i)
      for (int j= -1; j <= 1; ++j)
        for (int k= -1; k <= 1; ++k)
          for (int l= -1; l <= 1; ++l, ++c) {
              if (i == 0 && j == 0 && k == 0 && l == 0) continue;
              ls[0]= i; ls[1]= j; ls[2]= k; ls[3]= l;
              std::cout << "c: " << c << " ls: " << ls[0] << ' ' << ls[1] << ' ' << ls[2] << ' ' << ls[3] << std::endl;
              DROPS::RefTetraPartitionCL cut( static_cast<double*>(&ls[0]));
              DROPS::SignPatternTraitCL comb_cut( static_cast<double*>(&ls[0]));
              tet.make_partition<DROPS::SortedVertexPolicyCL, DROPS::MergeCutPolicyCL> ( DROPS::PrincipalLatticeCL::instance( 1), ls);
//              if (c == 5) {
//                  std::cerr << comb_cut << std::endl;
//                  std:: cerr << cut << std::endl;
//                  std::cerr << tet << std::endl << std::endl;
//              }
              std::ostringstream name;
              name << "hallo" << c << ".vtu";
              std::ofstream file( name.str().c_str());
              DROPS::write_paraview_vtu( file, tet);
          }
}

void test_cut_surface ()
{
    DROPS::GridFunctionCL<> ls( 4);
    ls[0]= -1.; ls[1]= 0.; ls[2]= 0.; ls[3]= 0.;
    DROPS::SurfacePatchCL tet;
    // tet.partition_principal_lattice ( 1, ls);
    // std::cerr << tet;
    int c= 0;
    for (int i= -1; i <= 1; ++i)
      for (int j= -1; j <= 1; ++j)
        for (int k= -1; k <= 1; ++k)
          for (int l= -1; l <= 1; ++l, ++c) {
              if (i == 0 && j == 0 && k == 0 && l == 0) continue;
              ls[0]= i; ls[1]= j; ls[2]= k; ls[3]= l;
              std::cout << "c: " << c << " ls: " << ls[0] << ' ' << ls[1] << ' ' << ls[2] << ' ' << ls[3] << std::endl;
              DROPS::RefTetraPartitionCL cut( static_cast<double*>(&ls[0]));
              DROPS::SignPatternTraitCL comb_cut( static_cast<double*>(&ls[0]));
              tet.make_patch<DROPS::MergeCutPolicyCL>( DROPS::PrincipalLatticeCL::instance( 1), ls);
              std::ostringstream name;
              name << "hallo_surf" << c << ".vtu";
              std::ofstream file( name.str().c_str());
              DROPS::write_paraview_vtu( file, tet);
          }
}

void test_principal_lattice ()
{
    for (int i= 1; i <= 4; ++i) {
        const DROPS::PrincipalLatticeCL& lat= DROPS::PrincipalLatticeCL::instance( i);
        std::cout << "=======================================" << lat.num_intervals() << ' ' << lat.vertex_size() << " " << lat.tetra_size() << std::endl;
        for (DROPS::PrincipalLatticeCL::const_vertex_iterator v= lat.vertex_begin(), end= lat.vertex_end(); v != end; ++v) {
            std::cout << lat.num_intervals()*(*v) << std::endl;
        }
        std:: cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
        for (DROPS::PrincipalLatticeCL::const_tetra_iterator v= lat.tetra_begin(), end= lat.tetra_end(); v != end; ++v) {
            std::cout << (*v)[0] << ' '  << (*v)[1] << ' ' << (*v)[2] << ' ' << (*v)[3] << ' ' << std::endl;
        }
    }
}

inline double sphere (const DROPS::Point3DCL& p)
{
    return p.norm() - 0.5;
}

inline double sphere_instat (const DROPS::Point3DCL& p, double)
{
    return sphere( p);
}

void test_sphere_cut ()
{
    DROPS::TetraBuilderCL tetrabuilder( 0);
    DROPS::MultiGridCL mg( tetrabuilder);

    const DROPS::PrincipalLatticeCL& lat= DROPS::PrincipalLatticeCL::instance( 10);
    DROPS::GridFunctionCL<> ls( lat.vertex_size());
    evaluate_on_vertexes( &sphere_instat, *mg.GetAllTetraBegin(), lat, 0., Addr( ls));
    DROPS::TetraPartitionCL tet;
    tet.make_partition<DROPS::SortedVertexPolicyCL, DROPS::MergeCutPolicyCL>( lat, ls);
    std::ostringstream name;
    name << "sphere.vtu";
    std::ofstream file( name.str().c_str());
    DROPS::write_paraview_vtu( file, tet, DROPS::NegTetraC);
    file.close();

    DROPS::SurfacePatchCL surf;
    surf.make_patch<DROPS::MergeCutPolicyCL>( lat, ls);
    name.str( "");
    name << "sphere_surf.vtu";
    file.open( name.str().c_str());
    DROPS::write_paraview_vtu( file, surf);
}


void test_sphere_integral ()
{
    std::cout << "Enter the number of subdivisions of the cube: ";
    DROPS::Uint num_sub;
    std::cin >> num_sub;
    std::cout << "Enter the number of subdivisions of the principal lattice: ";
    DROPS::Uint num_sub_lattice;
    std::cin >> num_sub_lattice;
    DROPS::BrickBuilderCL brick(DROPS::Point3DCL( -1.),
                                2.*DROPS::std_basis<3>(1),
                                2.*DROPS::std_basis<3>(2),
                                2.*DROPS::std_basis<3>(3),
                                num_sub, num_sub, num_sub);
    DROPS::MultiGridCL mg( brick);
    const DROPS::PrincipalLatticeCL& lat= DROPS::PrincipalLatticeCL::instance( num_sub_lattice);
    DROPS::GridFunctionCL<> ls( lat.vertex_size());
    DROPS::TetraPartitionCL tet;
    // DROPS::SurfacePatchCL patch;
    double vol_neg= 0., vol_pos= 0.;
    // double surf= 0.;
    DROPS::QuadDomainCL qdom;

    DROPS_FOR_TRIANG_TETRA( mg, 0, it) {
        evaluate_on_vertexes( sphere_instat, *it, lat, 0., Addr( ls));
        // tet.make_partition<DROPS::UnorderedVertexPolicyCL, DROPS::MergeCutPolicyCL>( lat, ls);
        // tet.make_partition<DROPS::SortedVertexPolicyCL, DROPS::MergeCutPolicyCL>( lat, ls);
        tet.make_partition<DROPS::PartitionedVertexPolicyCL, DROPS::MergeCutPolicyCL>( lat, ls);
        // patch.make_partition<DROPS::SortedVertexPolicyCL, DROPS::MergeCutPolicyCL>( lat, ls);
        DROPS::make_CompositeQuad5Domain( qdom, tet);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        double tmp_neg, tmp_pos;
        quad( integrand, it->GetVolume()*6., qdom, tmp_neg, tmp_pos);
        vol_neg+= tmp_neg; vol_pos+= tmp_pos;
        // q5_2d.assign( patch);
        // DROPS::GridFunctionCL<> surf_integrand( 1., q5_2d.size());
        // surf+= q5_2d.quad( surf_integrand);

    }
    std::cout << "Volume of the negative part: " << vol_neg << ", volume of the positive part: " << vol_pos << std::endl;
}

void test_extrapolated_sphere_integral ()
{
    std::cout << "Enter the number of subdivisions of the cube: ";
    DROPS::Uint num_sub;
    std::cin >> num_sub;
    std::cout << "Enter the number of extrapolation levels: ";
    DROPS::Uint num_level;
    std::cin >> num_level;
    DROPS::BrickBuilderCL brick(DROPS::Point3DCL( -1.),
                                2.*DROPS::std_basis<3>(1),
                                2.*DROPS::std_basis<3>(2),
                                2.*DROPS::std_basis<3>(3),
                                num_sub, num_sub, num_sub);
    DROPS::MultiGridCL mg( brick);
    double vol_neg= 0., vol_pos= 0.;
    DROPS::QuadDomainCL qdom;
    DROPS::ExtrapolationToZeroCL extra( num_level, DROPS::RombergSubdivisionCL());
    // DROPS::ExtrapolationToZeroCL extra( num_level, DROPS::HarmonicSubdivisionCL());

    DROPS_FOR_TRIANG_TETRA( mg, 0, it) {
        DROPS::LocalP2CL<> ls_loc( *it, &sphere_instat);
        make_ExtrapolatedQuad5Domain( qdom, ls_loc, extra);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        double tmp_neg, tmp_pos;
        quad( integrand, it->GetVolume()*6., qdom, tmp_neg, tmp_pos);
        vol_neg+= tmp_neg; vol_pos+= tmp_pos;
    }
    std::cout << "Volume of the negative part: " << vol_neg << ", volume of the positive part: " << vol_pos << std::endl;
}

void test_sphere_surface_integral ()
{
    std::cout << "Enter the number of subdivisions of the cube: ";
    DROPS::Uint num_sub;
    std::cin >> num_sub;
    std::cout << "Enter the number of subdivisions of the principal lattice: ";
    DROPS::Uint num_sub_lattice;
    std::cin >> num_sub_lattice;
    DROPS::BrickBuilderCL brick(DROPS::Point3DCL( -1.),
                                2.*DROPS::std_basis<3>(1),
                                2.*DROPS::std_basis<3>(2),
                                2.*DROPS::std_basis<3>(3),
                                num_sub, num_sub, num_sub);
    DROPS::MultiGridCL mg( brick);
    const DROPS::PrincipalLatticeCL& lat= DROPS::PrincipalLatticeCL::instance( num_sub_lattice);
    DROPS::GridFunctionCL<> ls( lat.vertex_size());
    DROPS::SurfacePatchCL patch;
    double surf= 0.;
    DROPS::QuadDomain2DCL qdom;

    DROPS_FOR_TRIANG_TETRA( mg, 0, it) {
        evaluate_on_vertexes( &sphere_instat, *it, lat, 0., Addr( ls));
        patch.make_patch<DROPS::MergeCutPolicyCL>( lat, ls);
        DROPS::make_CompositeQuad5Domain2D( qdom, patch, *it);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        surf+= quad_2D( integrand, qdom);
    }
    std::cout << "Surface: " <<surf << std::endl;
}

void test_extrapolated_sphere_surface_integral ()
{
    std::cout << "Enter the number of subdivisions of the cube: ";
    DROPS::Uint num_sub;
    std::cin >> num_sub;
    std::cout << "Enter the number of extrapolation levels: ";
    DROPS::Uint num_level;
    std::cin >> num_level;
    DROPS::BrickBuilderCL brick(DROPS::Point3DCL( -1.),
                                2.*DROPS::std_basis<3>(1),
                                2.*DROPS::std_basis<3>(2),
                                2.*DROPS::std_basis<3>(3),
                                num_sub, num_sub, num_sub);
    DROPS::MultiGridCL mg( brick);
    double surf= 0.;
    DROPS::QuadDomain2DCL qdom;
    DROPS::ExtrapolationToZeroCL extra( num_level, DROPS::RombergSubdivisionCL());
    // DROPS::ExtrapolationToZeroCL extra( num_level, DROPS::HarmonicSubdivisionCL());

    DROPS_FOR_TRIANG_TETRA( mg, 0, it) {
        DROPS::LocalP2CL<> ls_loc( *it, &sphere_instat);
        make_ExtrapolatedQuad5Domain2D( qdom, ls_loc, *it, extra);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        surf+= quad_2D( integrand, qdom);
    }
    std::cout << "Surface: " << surf << std::endl;
}

int main()
{
    try {
        // test_tetra_cut();
        // test_cut_surface();
        // test_principal_lattice();
        // test_sphere_cut();
        // test_sphere_integral();
        // test_extrapolated_sphere_integral();
        // test_sphere_surface_integral();
        // test_extrapolated_sphere_surface_integral();
        test_sign_traits_cut();
        write_sign_traits_1_2_3_4();
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
    return 0;
}
