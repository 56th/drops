/// \file spacetime.cpp
/// \brief classes that constitute space time geometries - 
/// Note that in contrast to VertexCL, ..., SimplexCL, ... these geometries only appear element-local 
/// and their purpose is the decomposition of a cut space-time primitive into several uncut ones
/// main components:
///  * Point4DContainer (gathered collection of unique Points)
///  * geometric primitive classes:
///    * Tetra4DCL
///    * HyperTrigCL
///    * GeneralizedPrism4CL
///    * PentatopeCL
/// PentatopeCL basically does all the work. With the help of levelset values a cut Pentatope
/// can be decomposed into uncut Pentatopes, HyperTrigs and GeneralizedPrisms. Each class can
/// decompose itself into Pentatopes. The Pentatope class also decomposes the interface into 
/// Tetra4Ds. 
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

#include "num/spacetime_geom.h"

namespace DROPS
{

inline void PentatopeCL::calc_abs_determinant() const{
    if (!absdet_of_trafo_initialized)
    {
        for (Uint i = 0; i < 5; ++i)
            for (Uint j = i+1; j < 5; ++j)
            {
                if (i==j) continue;
                Point4DCL diff = *p[j]-*p[i];
                if (diff.norm() < 1e-15)
                {
                    absdet_of_trafo = 0.0;
                    absdet_of_trafo_initialized = true;
                    return;
                }
            }
        QRDecompCL<4> Tmat;
        SMatrixCL<4,4> & trafomat = Tmat.GetMatrix();
        for (Uint i = 0; i < 4; ++i)
            for (Uint j = 0; j < 4; ++j)
            {
                trafomat(i,j) = (*p[j+1])[i]-(*p[0])[i];
            }
        if (Tmat.prepare_solve(false))
            absdet_of_trafo = 0.0;
        else
            absdet_of_trafo = std::fabs(Tmat.Determinant_R());
        absdet_of_trafo_initialized = true;
    }
}

double PentatopeCL::Measure() const{
    calc_abs_determinant();
    return absdet_of_trafo*1.0/24.0;
}

double PentatopeCL::AbsDet() const{
    calc_abs_determinant();
    return absdet_of_trafo;
}

/**
   a cut exists if the vertex values have a different sign. Here we count 0 as
   a value with a positive sign, s.t. for vertex value -1 and 0 we have a
   (degenerate) cut. As the resulting pentas have 0-measure they will be neglected 
   anyway...
 */
bool cutting_criteria_fulfilled(const double a, const double b)
{
    return (((a >= 0.0 && b < 0.0) ||(b >= 0.0 && a < 0.0))
            ||((a < 0.0 && b > 0.0)||(b < 0.0 && a > 0.0)));
}

// decomposition into positive and negative pentatopes. The level set values already 
// define a P1 representation. They could have been obtained via interpolation 
// or L2-Projection or whatever... 
void PentatopeCL::decompose_add_signed_pentas_4dtets(SArrayCL<double,5> & lsetvals, 
                                                     std::vector<PentatopeCL> & negpentas, 
                                                     std::vector<PentatopeCL> & pospentas,
                                                     std::vector<Tetra4DCL> & iftets) 
{
    const double eps = 1e-15;
    Uint ncuts = 0; //number of cutted edges

    Uint cutcounter[] = {0,0,0,0,0};
    for (Uint i = 0; i < 5; ++i)
    {
        for (Uint j = i+1; j < 5; ++j)
        {
            if (cutting_criteria_fulfilled(lsetvals[i],lsetvals[j]))
            {
                cutcounter[i] ++;
                cutcounter[j] ++;
                ncuts++;
            }
        }
    }
        
    // so far we can only deal with zero levels which are not crossing a vertex:
    if (ncuts!=0 && ncuts!=4 && ncuts!=6){
        std::cout << " ERROR-output: \n";
        std::cout << " *this = \n" << *this << std::endl;
        std::cout << " lsetvals = " << lsetvals << std::endl;
        std::cout << " ncuts = " << ncuts << std::endl;
        throw DROPSErrCL("PentatopeCL::decompose_add_signed_pentas_4dtets : number of edge cuts not in {0,4,6}");
    }

    //first (trivial) case:
    if (ncuts==0){
        if (lsetvals[0] > 0)
            pospentas.push_back(PentatopeCL(*this));
        else
            negpentas.push_back(PentatopeCL(*this));
        return;
    }

    //sorting of vertices - the latest one (or two) is separated from the other.
    //this facilitates the construction of the help-geometries prism4 and hypertrig.
    std::vector<std::pair<Uint,Uint> > verts(5);
    for (Uint i = 0; i < 5; ++i)
        verts[i] = std::pair<Uint,Uint>(cutcounter[i],Uint(i));
    std::sort(verts.begin(),verts.end());

    if (ncuts==4) // all edges leaving verts[4].second are cutted
    {
        Uint singlevertex = verts[4].second;            // the number of the vertex that is separated
        const Point4DCL & x4(*(p[singlevertex]));       // the vertex that is separated
        bool sign_of_sv = (lsetvals[singlevertex] > 0); // is the separated vertex positive?

        SArrayCL<const Point4DCL*,4> b( Uninitialized); // cutpoints (new points)
        SArrayCL<const Point4DCL*,4> x( Uninitialized); // basepoints (belonging to pentatope)

        Uint cuts = 0;
        for (Uint i = 0; i < 5; ++i)
        {
            if (i==singlevertex) continue;
            double alpha = - lsetvals[i]/(lsetvals[singlevertex]-lsetvals[i]);
            x[cuts] = p[i];
            b[cuts] = pcont((1-alpha) * (*x[cuts]) + alpha * x4);
            cuts++;
        }
            
        GeneralizedPrism4CL genprism4(pcont,x,b);
        PentatopeCL pp(pcont,*(b[0]),*(b[1]),*(b[2]),*(b[3]),x4, /* already_in_pcont = */ true);
        
        if (sign_of_sv){
            if (pp.Measure() > eps)
                pospentas.push_back(pp);
            genprism4.decompose_add_to_pentas(negpentas);
        }
        else
        {
            if (pp.Measure() > eps)
                negpentas.push_back(pp);
            genprism4.decompose_add_to_pentas(pospentas);
        }
            
        // the interface
        Tetra4DCL tet(pcont,*(b[0]),*(b[1]),*(b[2]),*(b[3]), /*already_in_pcont = */ true);
        tet.SetHelpPoint(x4,sign_of_sv,true);
        Point4DCL n = tet.Normal();
        if (tet.Measure() > eps)
            if (std::abs(std::abs(n[3])-1.0) < 1e-6){
                std::cout << " n = " << n << std::endl;
                std::cout << " strange tet found , tet = " << tet << std::endl;
                getchar();
            }
        if (tet.Measure() > eps)
            iftets.push_back(tet);
        
        return;
    }
    else if (ncuts==6) // all edges leaving verts[3].second or verts[4].second 
                       // and not connecting those two vertices are cutted
    {
        Uint sepvertex1 = verts[3].second; // the number of one of the vertices that are separated
        Uint sepvertex2 = verts[4].second; // the number of one of the vertices that are separated
        const Point4DCL & x4(*p[sepvertex1]);     // one of the vertices that is separated
        const Point4DCL & x5(*p[sepvertex2]);     // one of the vertices that is separated

        bool sign_of_sv = (lsetvals[sepvertex1] > 0); // are the separated vertex positive?

        SArrayCL<const Point4DCL*,3> c( Uninitialized); // first group of cutpoints (new points) (connected to sepvertex1)
        SArrayCL<const Point4DCL*,3> d( Uninitialized); // second group of cutpoints (new points) (connected to sepvertex2)
        SArrayCL<const Point4DCL*,3> x( Uninitialized); // basepoints (belonging to pentatope)

        Uint dcuts = 0; // "double cuts" - those two cuts which are on the two edges 
                        // connecting one vertex with the vertices sepvertex1/2 
        for (Uint i = 0; i < 5; ++i)
        {
            if (i==sepvertex1 || i==sepvertex2) continue;
            double alpha4 = - lsetvals[i]/(lsetvals[sepvertex1]-lsetvals[i]);
            double alpha5 = - lsetvals[i]/(lsetvals[sepvertex2]-lsetvals[i]);
            x[dcuts] = p[i];
            c[dcuts] = pcont((1-alpha4) * *(x[dcuts]) + alpha4 * x4);
            d[dcuts] = pcont((1-alpha5) * *(x[dcuts]) + alpha5 * x5);
            dcuts++;
        }

        SArrayCL<const Point4DCL*,4> extc( Uninitialized); // first group of cutpoints together with sepvertex1
        SArrayCL<const Point4DCL*,4> extd( Uninitialized); // second group of cutpoints together with sepvertex1
        for (Uint i = 0; i < 3; ++i)
        {
            extc[i] = c[i];
            extd[i] = d[i];
        }
        extc[3] = &x4;
        extd[3] = &x5;


        GeneralizedPrism4CL genprism4(pcont,extc,extd);
        HyperTrigCL hypertrig(pcont,x,c,d);

        if (sign_of_sv){
            genprism4.decompose_add_to_pentas(pospentas);
            hypertrig.decompose_add_to_pentas(negpentas);
        }
        else
        {
            genprism4.decompose_add_to_pentas(negpentas);
            hypertrig.decompose_add_to_pentas(pospentas);
        }

        // the interface ( a prism-3 in 4D)
        Tetra4DCL tet1(pcont,*(c[0]),*(c[1]),*(c[2]),*(d[2]), /* already_in_pcont = */ true);
        tet1.SetHelpPoint(x4,sign_of_sv, /* already_in_pcont = */ true);
        Point4DCL n = tet1.Normal();
        if (tet1.Measure() > eps)
            if (std::abs(std::abs(n[3])-1.0) < 1e-6){
                std::cout << " n = " << n << std::endl;
                std::cout << " strange tet found , tet = " << tet1 << std::endl;
                getchar();
            }

        Tetra4DCL tet2(pcont,*(c[0]),*(c[1]),*(d[1]),*(d[2]), /* already_in_pcont = */ true);
        tet2.SetHelpPoint(x4,sign_of_sv, /* already_in_pcont = */ true);
        n = tet2.Normal();
        if (tet2.Measure() > eps)
            if (std::abs(std::abs(n[3])-1.0) < 1e-6){
                std::cout << " n = " << n << std::endl;
                std::cout << " strange tet found , tet = " << tet2 << std::endl;
                getchar();
            }
        Tetra4DCL tet3(pcont,*(c[0]),*(d[0]),*(d[1]),*(d[2]), /* already_in_pcont = */ true);
        tet3.SetHelpPoint(x4,sign_of_sv, /* already_in_pcont = */ true);
        n = tet3.Normal();
        if (tet3.Measure() > eps)
            if (std::abs(std::abs(n[3])-1.0) < 1e-6){
                std::cout << " n = " << n << std::endl;
                std::cout << " strange tet found , tet = " << tet3 << std::endl;
                getchar();
            }

        if (tet1.Measure() > eps)
            iftets.push_back(tet1);
        if (tet2.Measure() > eps)
            iftets.push_back(tet2);
        if (tet3.Measure() > eps)
            iftets.push_back(tet3);

        return;
    }

}

std::ostream& operator<<(std::ostream& os, PentatopeCL penta)
{ 
    os << " Output of PentatopeCL:\n vertices:\n  ";
    for (Uint i = 0; i < 5; ++i)
        os << *(penta.p[i]) << '\n';
    os << "measure = " << penta.Measure() << "\n";
    os << "absdet = " << penta.AbsDet() << "\n";
    return os;
}
    
/**
   Decompose the generalized prism4 into pentatopes. Depending on plattice_1d_els divide the 
   "limit tetrahedra" x[0]-x[3] and y[0]-y[3] according to a principal lattice structure first,
   i.e. do a uniform spatial refinement.
*/
void GeneralizedPrism4CL::decompose_add_to_pentas (std::vector<PentatopeCL> & pentas, int plattice_1d_els, int time_1d_els) const 
{
    const double eps = 1e-15; // for now

    if (plattice_1d_els==1 && time_1d_els==1)
    { // the simple case
        PentatopeCL P1(pcont,*(x[0]),*(x[1]),*(x[2]),*(x[3]),*(y[3]), /*already_in_pcont = */ true);
        PentatopeCL P2(pcont,*(x[0]),*(x[1]),*(x[2]),*(y[2]),*(y[3]), /*already_in_pcont = */ true);
        PentatopeCL P3(pcont,*(x[0]),*(x[1]),*(y[1]),*(y[2]),*(y[3]), /*already_in_pcont = */ true);
        PentatopeCL P4(pcont,*(x[0]),*(y[0]),*(y[1]),*(y[2]),*(y[3]), /*already_in_pcont = */ true);

        if (P1.Measure() > eps)
            pentas.push_back(P1);
        if (P2.Measure() > eps)
            pentas.push_back(P2);
        if (P3.Measure() > eps)
            pentas.push_back(P3);
        if (P4.Measure() > eps)
            pentas.push_back(P4);
    }
    else
    { // the general case
        SMatrixCL<4,3> xtrafo, ytrafo;
        for (Uint i = 0; i < 4; ++i)
            for (Uint j = 0; j < 3; ++j)
            {
                xtrafo(i,j) = ( (*x[j+1])[i]- (*x[0])[i]);
                ytrafo(i,j) = ( (*y[j+1])[i]- (*y[0])[i]);
            }

        for (int ti = 0; ti < time_1d_els; ++ti)
        {
            Point4DCL z_x((1.0-((double)ti)/time_1d_els)* (*x[0])+(((double)ti)/time_1d_els)* (*y[0]));
            Point4DCL z_y((1.0-(ti+1.0)/time_1d_els)* (*x[0])+((ti+1.0)/time_1d_els)* (*y[0]));
            

            const PrincipalLatticeCL& pl( PrincipalLatticeCL::instance( plattice_1d_els));
            const std::vector<BaryCoordCL>& plvertices (pl.vertices());
            Point3DCL refcoord;
            SArrayCL<const Point4DCL*,4> tx,ty;

            for (PrincipalLatticeCL::const_tetra_iterator it= pl.tetra_begin(); it != pl.tetra_end(); ++it) {
                for ( Uint i = 0; i < 4; ++i)
                {
                    refcoord = BaryToRefCoord(plvertices[(*it)[i]]);
                    Point4DCL tmpx = z_x + xtrafo * refcoord;
                    tx[i] = pcont(tmpx);
                    Point4DCL tmpy = z_y + ytrafo * refcoord;
                    ty[i] = pcont(tmpy);
                }

                PentatopeCL P1(pcont, *(tx[0]), *(tx[1]), *(tx[2]), *(tx[3]), *(ty[3]), /*already_in_pcont = */ true);
                PentatopeCL P2(pcont, *(tx[0]), *(tx[1]), *(tx[2]), *(ty[2]), *(ty[3]), /*already_in_pcont = */ true);
                PentatopeCL P3(pcont, *(tx[0]), *(tx[1]), *(ty[1]), *(ty[2]), *(ty[3]), /*already_in_pcont = */ true);
                PentatopeCL P4(pcont, *(tx[0]), *(ty[0]), *(ty[1]), *(ty[2]), *(ty[3]), /*already_in_pcont = */ true);
                const double p1m = P1.Measure(); P1.Measure(); /// <- the second one makes things work with optimization... without it things get really strange...
                const double p2m = P2.Measure(); // P2.Measure();
                const double p3m = P3.Measure(); // P3.Measure();
                const double p4m = P4.Measure(); // P4.Measure();
                if (p1m > eps)
                    pentas.push_back(P1);
                if (p2m > eps)
                    pentas.push_back(P2);
                if (p3m > eps)
                    pentas.push_back(P3);
                if (p4m > eps)
                    pentas.push_back(P4);
            }
        }
    }
}


std::ostream& operator<<(std::ostream& os, GeneralizedPrism4CL genprism4)
{
    os << "Output of GeneralizedPrism4:\n x-points:\n";
    for (Uint i = 0; i < 4; ++i)
        os <<  *(genprism4.x[i]) << '\n';
    os << " --- \n y-points:\n";
    for (Uint i = 0; i < 4; ++i)
        os <<  *(genprism4.y[i]) << "\n\n";
    return os;
}

/**
   Decomposition of the hyper triangle into 6 pentatopes. The structure here is simpel:
   - the "diagonal triangle" u[0],v[1],w[2] is part of each pentatope. 
   - additionally the points which not already part of the "diagonal triangle" and 
     belong to one of 6  the "main triangles {u[0],u[1],u[2]},..,{w[0],w[1],w[2]},
     {u[0],v[0],v[0]},..,{u[2],v[2],v[2]} are added.
*/
void HyperTrigCL::decompose_add_to_pentas (std::vector<PentatopeCL> & pentas)
{
    const double eps = 1e-15; // for now

    PentatopeCL Du( pcont, *(u[0]), *(u[1]), *(u[2]), *(v[1]), *(w[2]), /*already_in_pcont = */ true);
    PentatopeCL Dv( pcont, *(u[0]), *(v[0]), *(v[1]), *(v[2]), *(w[2]), /*already_in_pcont = */ true);
    PentatopeCL Dw( pcont, *(u[0]), *(v[1]), *(w[0]), *(w[1]), *(w[2]), /*already_in_pcont = */ true);
    PentatopeCL D1( pcont, *(u[0]), *(v[0]), *(v[1]), *(w[0]), *(w[2]), /*already_in_pcont = */ true);
    PentatopeCL D2( pcont, *(u[0]), *(u[1]), *(v[1]), *(w[1]), *(w[2]), /*already_in_pcont = */ true);
    PentatopeCL D3( pcont, *(u[0]), *(u[2]), *(v[1]), *(v[2]), *(w[2]), /*already_in_pcont = */ true);
    if (Du.Measure() > eps)
        pentas.push_back(Du);
    if (Dv.Measure() > eps)
        pentas.push_back(Dv);
    if (Dw.Measure() > eps)
        pentas.push_back(Dw);
    if (D1.Measure() > eps)
        pentas.push_back(D1);
    if (D2.Measure() > eps)
        pentas.push_back(D2);
    if (D3.Measure() > eps)
        pentas.push_back(D3);
}

std::ostream& operator<<(std::ostream& os, HyperTrigCL hypert)
{
    os << "Output of HyperTrig:\n u-points:\n";
    for (Uint i = 0; i < 3; ++i)
        os <<  *(hypert.u[i]) << '\n';
    os << " --- \n v-points:\n";
    for (Uint i = 0; i < 3; ++i)
        os <<  *(hypert.v[i]) << '\n';
    os << " --- \n w-points:\n";
    for (Uint i = 0; i < 3; ++i)
        os <<  *(hypert.w[i]) << "\n\n";
    return os;
}
    

inline void Tetra4DCL::calc_abs_determinant() const
{
    if (!absdet_of_trafo_initialized)
    {

        for (Uint i = 0; i < 4; ++i)
            for (Uint j = i+1; j < 4; ++j)
            {
                if (i==j) continue;
                Point4DCL diff = *(p[j])-*(p[i]);
                if (diff.norm() < 1e-15)
                {
                    absdet_of_trafo = 0.0;
                    absdet_of_trafo_initialized = true;
                    return;
                }
            }
        
        normal = cross_product(*(p[1])-*(p[0]),*(p[2])-*(p[0]),*(p[3])-*(p[0]));
        absdet_of_trafo = normal.norm();
        normal /= absdet_of_trafo;
        

        // orientation - Normal should point from neg to pos
        // (works only if helppoint is prescribed)
        if (has_helppoint)
        {
            Point4DCL d = *helppoint-*(p[0]);
            const double sdir = inner_prod(d,normal);
            if (!(( sdir > 0 ) == (helpsign))) // normal is not pointing towards pos. region
                normal *= -1.0;
        }
        else
        {
            std::cout << " WARNING: No helppoint... " << std::endl;
        }

        // old version (before the 4D-cross-product):
        // QRDecompCL<3> Tmat;
        // SMatrixCL<3,3> & T = Tmat.GetMatrix();
        // SMatrixCL<4,3> trafomat;
        // for (Uint i = 0; i < 4; ++i)
        //     for (Uint j = 0; j < 3; ++j)
        //     {
        //         trafomat(i,j) = p[j+1][i]-p[0][i];
        //     }
        // T = GramMatrix(trafomat);
        // if (Tmat.prepare_solve(false))
        //     absdet_of_trafo = 0.0;
        // else
        //     absdet_of_trafo = std::sqrt(std::fabs(Tmat.Determinant_R()));

        absdet_of_trafo_initialized = true;
        normal_initialized = true;
    }
}

double Tetra4DCL::Measure() const{
    calc_abs_determinant();
    return absdet_of_trafo*1.0/6.0;
}

double Tetra4DCL::AbsDet() const{
    calc_abs_determinant();
    return absdet_of_trafo;
}

// return space time normal
Point4DCL Tetra4DCL::Normal() const{
    calc_normal();
    return normal;
}

// return 1.0/(sqrt(1+w_n^2)) = ||n_s|| where n_s is the spatial part of the space time normal 
double Tetra4DCL::Nu() const{
    calc_normal();
    return std::sqrt( std::pow(normal[0],2) + std::pow(normal[1],2) + std::pow(normal[2],2));
}

inline void Tetra4DCL::calc_normal() const{
    if (!normal_initialized)
    {
        // old version (before the 4D-cross-product):
        // There holds A^T n = 0 by definition
        // After QR-decomposition there further holds R^T Q^T n = 0
        // As a has full rank there holds R^T y = 0 => y = c * (0,0,0,1)
        // We consider y = (0,0,0,1) as it is normalized to ||y||_2 = 1
        // Thus there holds n = Q * y with ||n|| = 1
        
        // QRDecompCL<4,3> Tmat;
        // SMatrixCL<4,3> & trafomat = Tmat.GetMatrix();
        // for (Uint i = 0; i < 4; ++i)
        //     for (Uint j = 0; j < 3; ++j)
        //         trafomat(i,j) = p[j+1][i]-p[0][i];
        // Tmat.prepare_solve(false);
        // normal[0]=normal[1]=normal[2]=0.0; normal[3] = 1.0; // y = (0,0,0,1)
        // Tmat.apply_Q(normal);

        normal = cross_product(*(p[1])-*(p[0]),*(p[2])-*(p[0]),*(p[3])-*(p[0]));
        absdet_of_trafo = normal.norm();
        normal /= absdet_of_trafo;

        // orientation - Normal should point from neg to pos
        // (works only if helppoint is prescribed)
        if (has_helppoint)
        {
            Point4DCL d = *helppoint-*(p[0]);
            const double sdir = inner_prod(d,normal);
            if (!(( sdir > 0 ) == (helpsign))) // normal is not pointing towards pos. region
                normal *= -1.0;
        }
        else
        {
            std::cout << " WARNING: No helppoint... " << std::endl;
        }

        normal_initialized = true;
        absdet_of_trafo_initialized = true;

    }
}

std::ostream& operator<<(std::ostream& os, const Tetra4DCL & tet4d)
{ 
    os << " Output of Tetra4D:\n vertices:\n  ";
    for (Uint i = 0; i < 4; ++i)
        os << tet4d.p[i] << '\n';
    os << "measure = " << tet4d.Measure() << "\n";
    os << "absdet = " << tet4d.AbsDet() << "\n";
    return os;
}



} // end of namespace DROPS
