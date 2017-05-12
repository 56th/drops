/// \file discretize.cpp
/// \brief discretizations for several PDEs and FE types
/// \author LNM RWTH Aachen: Patrick Esser, Sven Gross, Trung Hieu Nguyen, Joerg Peters, Volker Reichelt; SC RWTH Aachen:

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

#include "num/discretize.h"
#include "num/interfacePatch.h"
#include "num/fe.h"

namespace DROPS
{

static inline BaryCoordCL*
TransformNodes (const SArrayCL<BaryCoordCL,4>& M, BaryCoordCL* p, Uint NumNodes, const BaryCoordCL* Node)
{
    if (!p) p= new BaryCoordCL[NumNodes];
    for (Uint i=0; i < NumNodes; ++i)
        //p[i]=M*Node[i]; M (als Matrix) ist spaltenweise gespeichert!
        for (Uint k= 0; k < 4; ++k)
            p[i][k]= M[0][k]*Node[i][0] + M[1][k]*Node[i][1]
                   + M[2][k]*Node[i][2] + M[3][k]*Node[i][3];
    return p;
}

BaryCoordCL Quad2DataCL::Node[NumNodesC];

const double Quad2DataCL::Wght[2]= {
    1./120., /* Node[0] bis Node[3]*/
    2./15., /* 2./15., Node[4]*/
};

const double Quad2DataCL::Weight[NumNodesC]= {
    1./120.,
    1./120.,
    1./120.,
    1./120.,

    2./15.
};

Quad2DataCL::Quad2DataCL()
{
    Node[0]= MakeBaryCoord( 1., 0., 0., 0.);
    Node[1]= MakeBaryCoord( 0., 1., 0., 0.);
    Node[2]= MakeBaryCoord( 0., 0., 1., 0.);
    Node[3]= MakeBaryCoord( 0., 0., 0., 1.);
    Node[4]= BaryCoordCL( 0.25);
}

const double Quad2Data_Mul_P2_CL::weights_[10][NumNodesC]= {
    { 1./360., 0., 0., 0., -1./90. },
    { 0., 1./360., 0., 0., -1./90. },
    { 0., 0., 1./360., 0., -1./90. },
    { 0., 0., 0., 1./360., -1./90. },

    { 1./180., 1./180., 0., 0., 1./45. },
    { 1./180., 0., 1./180., 0., 1./45. },
    { 0., 1./180., 1./180., 0., 1./45. },
    { 1./180., 0., 0., 1./180., 1./45. },
    { 0., 1./180., 0., 1./180., 1./45. },
    { 0., 0., 1./180., 1./180., 1./45. }
};

const double Quad2Data_Mul_P1_CL::weights_[4][NumNodesC]= {
    { 1./120., 0., 0., 0., 1./30. },
    { 0., 1./120., 0., 0., 1./30. },
    { 0., 0., 1./120., 0., 1./30. },
    { 0., 0., 0., 1./120., 1./30. }
};

///\brief special implementation to copy the first 4 nodes.
BaryCoordCL*
Quad2DataCL::TransformNodes (const SArrayCL<BaryCoordCL,4>& M, BaryCoordCL* p)
{
    if (!p) p= new BaryCoordCL[NumNodesC];
    //tN[i]=M*Node[i]; M (als Matrix) ist spaltenweise gespeichert!
    std::memcpy( p, M.begin(), 4*sizeof( BaryCoordCL));
    for (Uint k= 0; k < 4; ++k)
            p[4][k]= 0.25*( M[0][k] + M[1][k] + M[2][k] + M[3][k]);
    return p;
}

namespace {
    Quad2DataCL theQuad2DataInitializer_; // The constructor sets up the static arrays
} // end of anonymous namespace

//**************************************************************************
// Class: Quad3DataCL                                                      *
//**************************************************************************

BaryCoordCL Quad3DataCL::Node[NumNodesC];

const double Quad3DataCL::Wght[2]= {
    -2./15., /* -(n+1)^2/[4(n+2)] /6       Node[0]*/
    3./40.  /* (n+3)^2/[4(n+1)(n+2)] /6 , Node[1] bis Node[4]*/
};

const double Quad3DataCL::Weight[5]= {
    -2./15., /* -(n+1)^2/[4(n+2)] /6       Node[0]*/
    3./40.,  /* (n+3)^2/[4(n+1)(n+2)] /6 , Node[1] bis Node[4]*/
    3./40.,
    3./40.,
    3./40.
};

std::valarray<double> Quad3DataCL::P2_Val[10]; // P2_Val[i] contains FE_P2CL::H_i( Node).

Quad3DataCL::Quad3DataCL()
{
    Node[0]= BaryCoordCL( 0.25);
    const double A= 1./6.,
                 B= 0.5;
    Node[1]= MakeBaryCoord( A,A,A,B);
    Node[2]= MakeBaryCoord( A,A,B,A);
    Node[3]= MakeBaryCoord( A,B,A,A);
    Node[4]= MakeBaryCoord( B,A,A,A);

    FE_P2CL::ApplyAll( NumNodesC, Node, P2_Val);
}

BaryCoordCL* Quad3DataCL::TransformNodes (const SArrayCL<BaryCoordCL,4>& M, BaryCoordCL* p)
{
    return DROPS::TransformNodes( M, p, NumNodesC, Node);
}

namespace {
    Quad3DataCL theQuad3DataInitializer_; // The constructor sets up the static arrays
} // end of anonymous namespace

//**************************************************************************
// Class: Quad5DataCL                                                      *
//**************************************************************************

BaryCoordCL Quad5DataCL::Node[NumNodesC];

const double Quad5DataCL::Wght[4]= {
    0.01975308641975308641975308641975308641975, /* 8./405.,                                   Node[0]*/
    0.01198951396316977000173064248499538621936, /* (2665.0 + 14.0*std::sqrt( 15.0))/226800.0, Node[1] bis Node[4]*/
    0.01151136787104539754677023934921978132914, /* (2665.0 - 14.0*std::sqrt( 15.0))/226800.0, Node[5] bis Node[8]*/
    0.008818342151675485008818342151675485008818 /* 5./567.,                                   Node[9] bis Node[14]*/
};

const double Quad5DataCL::Weight[NumNodesC]= {
    0.01975308641975308641975308641975308641975,

    0.01198951396316977000173064248499538621936,
    0.01198951396316977000173064248499538621936,
    0.01198951396316977000173064248499538621936,
    0.01198951396316977000173064248499538621936,

    0.01151136787104539754677023934921978132914,
    0.01151136787104539754677023934921978132914,
    0.01151136787104539754677023934921978132914,
    0.01151136787104539754677023934921978132914,

    0.008818342151675485008818342151675485008818,
    0.008818342151675485008818342151675485008818,
    0.008818342151675485008818342151675485008818,
    0.008818342151675485008818342151675485008818,
    0.008818342151675485008818342151675485008818,
    0.008818342151675485008818342151675485008818
};

std::valarray<double> Quad5DataCL::P2_Val[10]; // P2_Val[i] contains FE_P2CL::H_i( Node).

Quad5DataCL::Quad5DataCL()
{
    Node[0]= BaryCoordCL( 0.25);
    const double A1= (7.0 - std::sqrt( 15.0))/34.0,
                 B1= (13.0 + 3.0*std::sqrt( 15.0))/34.0;
    Node[1]= MakeBaryCoord( A1,A1,A1,B1);
    Node[2]= MakeBaryCoord( A1,A1,B1,A1);
    Node[3]= MakeBaryCoord( A1,B1,A1,A1);
    Node[4]= MakeBaryCoord( B1,A1,A1,A1);
    const double A2= (7.0 + std::sqrt( 15.0))/34.0,
                 B2= (13.0 - 3.0*std::sqrt( 15.0))/34.0;
    Node[5]= MakeBaryCoord( A2,A2,A2,B2);
    Node[6]= MakeBaryCoord( A2,A2,B2,A2);
    Node[7]= MakeBaryCoord( A2,B2,A2,A2);
    Node[8]= MakeBaryCoord( B2,A2,A2,A2);
    const double A3= (10.0 - 2.0*std::sqrt( 15.0))/40.0,
                 B3= (10.0 + 2.0*std::sqrt( 15.0))/40.0;
    Node[9] = MakeBaryCoord( A3,A3,B3,B3);
    Node[10]= MakeBaryCoord( A3,B3,A3,B3);
    Node[11]= MakeBaryCoord( A3,B3,B3,A3);
    Node[12]= MakeBaryCoord( B3,A3,A3,B3);
    Node[13]= MakeBaryCoord( B3,A3,B3,A3);
    Node[14]= MakeBaryCoord( B3,B3,A3,A3);

    FE_P2CL::ApplyAll( NumNodesC, Node, P2_Val);
}

BaryCoordCL*
Quad5DataCL::TransformNodes (const SArrayCL<BaryCoordCL,4>& M, BaryCoordCL* p)
{
    return DROPS::TransformNodes( M, p, NumNodesC, Node);
}

namespace {
    Quad5DataCL theQuad5DataInitializer_; // The constructor sets up the static arrays
} // end of anonymous namespace

//**************************************************************************
// Class: Quad9_1DDataCL                                                   *
//**************************************************************************
const double Quad9_1DDataCL::Node[NumNodesC] {
     0,
    -0.5384693101,
     0.5384693101,
    -0.9061798459,
     0.9061798459
};

const double Quad9_1DDataCL::Weight[NumNodesC] {
    0.5688888889, 
    0.4786286705,
    0.4786286705,
    0.2369268851,
    0.2369268851
};

//**************************************************************************
// Class: Quad5_2DDataCL                                                   *
//**************************************************************************
Point3DCL Quad5_2DDataCL::Node[NumNodesC];  // Barycentric coord for 2D

const double Quad5_2DDataCL::Wght[3]= {
      9./80.,                          /*Node[0]*/
      (155. - std::sqrt( 15.0))/2400., /*Node[1] to Node[3]*/
      (155. + std::sqrt( 15.0))/2400.  /*Node[4] to Node[6]*/
};

const double Quad5_2DDataCL::Weight[NumNodesC]= {
      9./80.,                          /*Node[0]*/
      (155. - std::sqrt( 15.0))/2400., /*Node[1] to Node[3]*/
      (155. - std::sqrt( 15.0))/2400.,
      (155. - std::sqrt( 15.0))/2400.,
      (155. + std::sqrt( 15.0))/2400., /*Node[4] to Node[6]*/
      (155. + std::sqrt( 15.0))/2400.,
      (155. + std::sqrt( 15.0))/2400.
};

Quad5_2DDataCL::Quad5_2DDataCL ()
{
    Node[0]= Point3DCL( 1./3.);
    const double A1= (6.0 - std::sqrt( 15.0))/21.0,
    B1= (9.0 + 2.0*std::sqrt( 15.0))/21.0;
    Node[1]= MakePoint3D( A1,A1,B1);
    Node[2]= MakePoint3D( A1,B1,A1);
    Node[3]= MakePoint3D( B1,A1,A1);
    const double A2= (6.0 + std::sqrt( 15.0))/21.0,
    B2= (9.0 - 2.0*std::sqrt( 15.0))/21.0;
    Node[4]= MakePoint3D( A2,A2,B2);
    Node[5]= MakePoint3D( A2,B2,A2);
    Node[6]= MakePoint3D( B2,A2,A2);
}

namespace {
    Quad5_2DDataCL theQuad52DDataInitializer_; // The constructor sets up the static arrays
} // end of anonymous namespace

//**************************************************************************
// Class: Quad3_4DDataCL                                                   *
//**************************************************************************
STBaryCoordCL Quad3_4DDataCL::Node[NumNodesC]; 

const double Quad3_4DDataCL::Weight[NumNodesC]= {
      1./120.,                          /*Node[0]*/
      1./120.,                          /*Node[1]*/
      1./120.,                          /*Node[2]*/
      1./120.,                          /*Node[3]*/
      1./120.                           /*Node[4]*/
};

Quad3_4DDataCL::Quad3_4DDataCL ()
{
    const double A = 0.118350341907227374; 
    const double B = 0.526598632371090503;
    Node[0]= MakeSTBaryCoord( A, A, A, A, 1-4.0*A);
    Node[1]= MakeSTBaryCoord( B, A, A, A, 1-3.0*A-B);
    Node[2]= MakeSTBaryCoord( A, B, A, A, 1-3.0*A-B);
    Node[3]= MakeSTBaryCoord( A, A, B, A, 1-3.0*A-B);
    Node[4]= MakeSTBaryCoord( A, A, A, B, 1-3.0*A-B);
}

namespace {
    Quad3_4DDataCL theQuad34DDataInitializer_; // The constructor sets up the static arrays
} // end of anonymous namespace


//**************************************************************************
// Class: Quad5_4DDataCL                                                   *
//**************************************************************************
STBaryCoordCL Quad5_4DDataCL::Node[NumNodesC]; 

const double Quad5_4DDataCL::Weight[NumNodesC]= {
    1.852751404612798e-03, 1.124562934816807e-03, 1.124562934816807e-03, 1.124562934816807e-03, 1.124562934816807e-03, 
    1.079714963974758e-03, 1.079714963974758e-03, 1.079714963974758e-03, 1.079714963974758e-03, 8.271211627735704e-04, 
    8.271211627735704e-04, 8.271211627735704e-04, 8.271211627735704e-04, 8.271211627735704e-04, 8.271211627735704e-04, 
    2.153172719775622e-03, 1.306909403436344e-03, 1.306909403436344e-03, 1.306909403436344e-03, 1.306909403436344e-03, 
    1.254789390403849e-03, 1.254789390403849e-03, 1.254789390403849e-03, 1.254789390403849e-03, 9.612378213284026e-04, 
    9.612378213284026e-04, 9.612378213284026e-04, 9.612378213284026e-04, 9.612378213284026e-04, 9.612378213284026e-04, 
    8.403303534373656e-04, 5.100545956270171e-04, 5.100545956270171e-04, 5.100545956270171e-04, 5.100545956270171e-04, 
    4.897134364759201e-04, 4.897134364759201e-04, 4.897134364759201e-04, 4.897134364759201e-04, 3.751474792131096e-04, 
    3.751474792131096e-04, 3.751474792131096e-04, 3.751474792131096e-04, 3.751474792131096e-04, 3.751474792131096e-04, 
    9.201712711248374e-05, 5.585155691227325e-05, 5.585155691227325e-05, 5.585155691227325e-05, 5.585155691227325e-05, 
    5.362417690682013e-05, 5.362417690682013e-05, 5.362417690682013e-05, 5.362417690682013e-05, 4.107907460378738e-05, 
    4.107907460378738e-05, 4.107907460378738e-05, 4.107907460378738e-05, 4.107907460378738e-05, 4.107907460378738e-05
};


Quad5_4DDataCL::Quad5_4DDataCL ()
{
   Node[ 0] = MakeSTBaryCoord( 2.394617701415296e-01, 2.394617701415297e-01, 2.394617701415297e-01, 2.394617701415297e-01, 4.215291943388128e-02);
   Node[ 1] = MakeSTBaryCoord( 8.809422860931937e-02, 8.809422860931940e-02, 8.809422860931940e-02, 6.935643947381607e-01, 4.215291943388128e-02);
   Node[ 2] = MakeSTBaryCoord( 8.809422860931937e-02, 8.809422860931940e-02, 6.935643947381607e-01, 8.809422860931940e-02, 4.215291943388128e-02);
   Node[ 3] = MakeSTBaryCoord( 8.809422860931937e-02, 6.935643947381607e-01, 8.809422860931940e-02, 8.809422860931940e-02, 4.215291943388128e-02);
   Node[ 4] = MakeSTBaryCoord( 6.935643947381605e-01, 8.809422860931940e-02, 8.809422860931940e-02, 8.809422860931940e-02, 4.215291943388128e-02);
   Node[ 5] = MakeSTBaryCoord( 3.063133928002588e-01, 3.063133928002589e-01, 3.063133928002589e-01, 3.890690216534197e-02, 4.215291943388128e-02);
   Node[ 6] = MakeSTBaryCoord( 3.063133928002588e-01, 3.063133928002589e-01, 3.890690216534197e-02, 3.063133928002589e-01, 4.215291943388128e-02);
   Node[ 7] = MakeSTBaryCoord( 3.063133928002588e-01, 3.890690216534197e-02, 3.063133928002589e-01, 3.063133928002589e-01, 4.215291943388128e-02);
   Node[ 8] = MakeSTBaryCoord( 3.890690216534187e-02, 3.063133928002589e-01, 3.063133928002589e-01, 3.063133928002589e-01, 4.215291943388128e-02);
   Node[ 9] = MakeSTBaryCoord( 5.397548057923118e-02, 5.397548057923109e-02, 4.249480597038283e-01, 4.249480597038283e-01, 4.215291943388128e-02);
   Node[10] = MakeSTBaryCoord( 5.397548057923118e-02, 4.249480597038283e-01, 5.397548057923109e-02, 4.249480597038283e-01, 4.215291943388128e-02);
   Node[11] = MakeSTBaryCoord( 5.397548057923118e-02, 4.249480597038283e-01, 4.249480597038283e-01, 5.397548057923109e-02, 4.215291943388128e-02);
   Node[12] = MakeSTBaryCoord( 4.249480597038282e-01, 4.249480597038283e-01, 5.397548057923109e-02, 5.397548057923109e-02, 4.215291943388128e-02);
   Node[13] = MakeSTBaryCoord( 4.249480597038282e-01, 5.397548057923109e-02, 4.249480597038283e-01, 5.397548057923109e-02, 4.215291943388128e-02);
   Node[14] = MakeSTBaryCoord( 4.249480597038282e-01, 5.397548057923109e-02, 5.397548057923109e-02, 4.249480597038283e-01, 4.215291943388128e-02);
   Node[15] = MakeSTBaryCoord( 1.975708074923217e-01, 1.975708074923217e-01, 1.975708074923217e-01, 1.975708074923217e-01, 2.097167700307132e-01);
   Node[16] = MakeSTBaryCoord( 7.268320062726330e-02, 7.268320062726334e-02, 7.268320062726334e-02, 5.722336280874968e-01, 2.097167700307132e-01);
   Node[17] = MakeSTBaryCoord( 7.268320062726341e-02, 7.268320062726334e-02, 5.722336280874968e-01, 7.268320062726334e-02, 2.097167700307132e-01);
   Node[18] = MakeSTBaryCoord( 7.268320062726341e-02, 5.722336280874968e-01, 7.268320062726334e-02, 7.268320062726334e-02, 2.097167700307132e-01);
   Node[19] = MakeSTBaryCoord( 5.722336280874968e-01, 7.268320062726334e-02, 7.268320062726334e-02, 7.268320062726334e-02, 2.097167700307132e-01);
   Node[20] = MakeSTBaryCoord( 2.527275411247960e-01, 2.527275411247959e-01, 2.527275411247959e-01, 3.210060659489898e-02, 2.097167700307132e-01);
   Node[21] = MakeSTBaryCoord( 2.527275411247960e-01, 2.527275411247959e-01, 3.210060659489898e-02, 2.527275411247959e-01, 2.097167700307132e-01);
   Node[22] = MakeSTBaryCoord( 2.527275411247960e-01, 3.210060659489898e-02, 2.527275411247959e-01, 2.527275411247959e-01, 2.097167700307132e-01);
   Node[23] = MakeSTBaryCoord( 3.210060659489911e-02, 2.527275411247959e-01, 2.527275411247959e-01, 2.527275411247959e-01, 2.097167700307132e-01);
   Node[24] = MakeSTBaryCoord( 4.453311806941906e-02, 4.453311806941899e-02, 3.506084969152243e-01, 3.506084969152243e-01, 2.097167700307132e-01);
   Node[25] = MakeSTBaryCoord( 4.453311806941906e-02, 3.506084969152243e-01, 4.453311806941899e-02, 3.506084969152243e-01, 2.097167700307132e-01);
   Node[26] = MakeSTBaryCoord( 4.453311806941918e-02, 3.506084969152243e-01, 3.506084969152243e-01, 4.453311806941899e-02, 2.097167700307132e-01);
   Node[27] = MakeSTBaryCoord( 3.506084969152244e-01, 3.506084969152243e-01, 4.453311806941899e-02, 4.453311806941899e-02, 2.097167700307132e-01);
   Node[28] = MakeSTBaryCoord( 3.506084969152244e-01, 4.453311806941899e-02, 3.506084969152243e-01, 4.453311806941899e-02, 2.097167700307132e-01);
   Node[29] = MakeSTBaryCoord( 3.506084969152244e-01, 4.453311806941899e-02, 4.453311806941899e-02, 3.506084969152243e-01, 2.097167700307132e-01);
   Node[30] = MakeSTBaryCoord( 1.339616115220625e-01, 1.339616115220626e-01, 1.339616115220626e-01, 1.339616115220626e-01, 4.641535539117497e-01);
   Node[31] = MakeSTBaryCoord( 4.928237531745672e-02, 4.928237531745671e-02, 4.928237531745671e-02, 3.879993201358802e-01, 4.641535539117497e-01);
   Node[32] = MakeSTBaryCoord( 4.928237531745672e-02, 4.928237531745671e-02, 3.879993201358802e-01, 4.928237531745671e-02, 4.641535539117497e-01);
   Node[33] = MakeSTBaryCoord( 4.928237531745672e-02, 3.879993201358802e-01, 4.928237531745671e-02, 4.928237531745671e-02, 4.641535539117497e-01);
   Node[34] = MakeSTBaryCoord( 3.879993201358801e-01, 4.928237531745671e-02, 4.928237531745671e-02, 4.928237531745671e-02, 4.641535539117497e-01);
   Node[35] = MakeSTBaryCoord( 1.713602789541758e-01, 1.713602789541757e-01, 1.713602789541757e-01, 2.176560922572299e-02, 4.641535539117497e-01);
   Node[36] = MakeSTBaryCoord( 1.713602789541758e-01, 1.713602789541757e-01, 2.176560922572299e-02, 1.713602789541757e-01, 4.641535539117497e-01);
   Node[37] = MakeSTBaryCoord( 1.713602789541758e-01, 2.176560922572299e-02, 1.713602789541757e-01, 1.713602789541757e-01, 4.641535539117497e-01);
   Node[38] = MakeSTBaryCoord( 2.176560922572301e-02, 1.713602789541757e-01, 1.713602789541757e-01, 1.713602789541757e-01, 4.641535539117497e-01);
   Node[39] = MakeSTBaryCoord( 3.019539343085142e-02, 3.019539343085138e-02, 2.377278296132737e-01, 2.377278296132737e-01, 4.641535539117497e-01);
   Node[40] = MakeSTBaryCoord( 3.019539343085142e-02, 2.377278296132737e-01, 3.019539343085138e-02, 2.377278296132737e-01, 4.641535539117497e-01);
   Node[41] = MakeSTBaryCoord( 3.019539343085142e-02, 2.377278296132737e-01, 2.377278296132737e-01, 3.019539343085138e-02, 4.641535539117497e-01);
   Node[42] = MakeSTBaryCoord( 2.377278296132738e-01, 2.377278296132737e-01, 3.019539343085138e-02, 3.019539343085138e-02, 4.641535539117497e-01);
   Node[43] = MakeSTBaryCoord( 2.377278296132738e-01, 3.019539343085138e-02, 2.377278296132737e-01, 3.019539343085138e-02, 4.641535539117497e-01);
   Node[44] = MakeSTBaryCoord( 2.377278296132738e-01, 3.019539343085138e-02, 3.019539343085138e-02, 2.377278296132737e-01, 4.641535539117497e-01);
   Node[45] = MakeSTBaryCoord( 6.536944720772242e-02, 6.536944720772242e-02, 6.536944720772242e-02, 6.536944720772242e-02, 7.385222111691103e-01);
   Node[46] = MakeSTBaryCoord( 2.404839412561921e-02, 2.404839412561919e-02, 2.404839412561919e-02, 1.893326064540321e-01, 7.385222111691103e-01);
   Node[47] = MakeSTBaryCoord( 2.404839412561921e-02, 2.404839412561919e-02, 1.893326064540321e-01, 2.404839412561919e-02, 7.385222111691103e-01);
   Node[48] = MakeSTBaryCoord( 2.404839412561921e-02, 1.893326064540321e-01, 2.404839412561919e-02, 2.404839412561919e-02, 7.385222111691103e-01);
   Node[49] = MakeSTBaryCoord( 1.893326064540322e-01, 2.404839412561919e-02, 2.404839412561919e-02, 2.404839412561919e-02, 7.385222111691103e-01);
   Node[50] = MakeSTBaryCoord( 8.361893068710013e-02, 8.361893068710009e-02, 8.361893068710009e-02, 1.062099676958939e-02, 7.385222111691103e-01);
   Node[51] = MakeSTBaryCoord( 8.361893068710013e-02, 8.361893068710009e-02, 1.062099676958939e-02, 8.361893068710009e-02, 7.385222111691103e-01);
   Node[52] = MakeSTBaryCoord( 8.361893068710013e-02, 1.062099676958939e-02, 8.361893068710009e-02, 8.361893068710009e-02, 7.385222111691103e-01);
   Node[53] = MakeSTBaryCoord( 1.062099676958939e-02, 8.361893068710009e-02, 8.361893068710009e-02, 8.361893068710009e-02, 7.385222111691103e-01);
   Node[54] = MakeSTBaryCoord( 1.473449113046366e-02, 1.473449113046365e-02, 1.160044032849812e-01, 1.160044032849812e-01, 7.385222111691103e-01);
   Node[55] = MakeSTBaryCoord( 1.473449113046366e-02, 1.160044032849812e-01, 1.473449113046365e-02, 1.160044032849812e-01, 7.385222111691103e-01);
   Node[56] = MakeSTBaryCoord( 1.473449113046366e-02, 1.160044032849812e-01, 1.160044032849812e-01, 1.473449113046365e-02, 7.385222111691103e-01);
   Node[57] = MakeSTBaryCoord( 1.160044032849812e-01, 1.160044032849812e-01, 1.473449113046365e-02, 1.473449113046365e-02, 7.385222111691103e-01);
   Node[58] = MakeSTBaryCoord( 1.160044032849812e-01, 1.473449113046365e-02, 1.160044032849812e-01, 1.473449113046365e-02, 7.385222111691103e-01);
   Node[59] = MakeSTBaryCoord( 1.160044032849812e-01, 1.473449113046365e-02, 1.473449113046365e-02, 1.160044032849812e-01, 7.385222111691103e-01);
}

namespace {
    Quad5_4DDataCL theQuad54DDataInitializer_; // The constructor sets up the static arrays
} // end of anonymous namespace






const double Quad3PosWeightsCL::_points[8][4]= {
    {1.,0.,0.,0.}, {0.,1.,0.,0.}, {0.,0.,1.,0.}, {0.,0.,0.,1.},
    {1./3.,1./3.,1./3.,0.}, {1./3.,1./3.,0.,1./3.}, {1./3.,0.,1./3.,1./3.},
    {0.,1./3.,1./3.,1./3.}
    };

const double FaceQuad2CL::_points[4][2]= {
    {0., 0.}, {1., 0.}, {0., 1.}, {1./3., 1./3.} };

const double P1BubbleDiscCL::_points1[26][4]= {
    {1.,0.,0.,0.}, {1./3.,1./3.,1./3.,0.}, {1./3.,1./3.,0.,1./3.}, {1./3.,0.,1./3.,1./3.},
    {0.,1.,0.,0.}, {1./3.,1./3.,1./3.,0.}, {1./3.,1./3.,0.,1./3.}, {0.,1./3.,1./3.,1./3.},
    {0.,0.,1.,0.}, {1./3.,1./3.,1./3.,0.}, {1./3.,0.,1./3.,1./3.}, {0.,1./3.,1./3.,1./3.},
    {0.,0.,0.,1.}, {1./3.,1./3.,0.,1./3.}, {1./3.,0.,1./3.,1./3.}, {0.,1./3.,1./3.,1./3.},
    {1.,0.,0.,0.}, {0.,1.,0.,0.}, {0.,0.,1.,0.}, {0.,0.,0.,1.}, {.5,.5,0.,0.}, {.5,0.,.5,0.}, {.5,0.,0.,.5}, {0.,.5,.5,0.}, {0.,.5,0.,.5}, {0.,0.,.5,.5}
    };

void P2DiscCL::GetGradientsOnRef( LocalP1CL<Point3DCL> GRef[10])
{
    for (int i= 0; i < 10; ++i)
    {
        GRef[i][0]= FE_P2CL::DHRef( i, 0,0,0);
        GRef[i][1]= FE_P2CL::DHRef( i, 1,0,0);
        GRef[i][2]= FE_P2CL::DHRef( i, 0,1,0);
        GRef[i][3]= FE_P2CL::DHRef( i, 0,0,1);
    }
}

void P2DiscCL::GetGradientsOnRef( Quad2CL<Point3DCL> GRef[10])
{
    for (int i=0; i<10; ++i)
    {
        GRef[i][0]= FE_P2CL::DHRef( i, 0,0,0);
        GRef[i][1]= FE_P2CL::DHRef( i, 1,0,0);
        GRef[i][2]= FE_P2CL::DHRef( i, 0,1,0);
        GRef[i][3]= FE_P2CL::DHRef( i, 0,0,1);
        GRef[i][4]= FE_P2CL::DHRef( i, 0.25,0.25,0.25);
    }
}

void P2DiscCL::GetGradientsOnRef( Quad5CL<Point3DCL> GRef[10])
{
    for (int i=0; i<10; ++i)
        for (int j=0; j<Quad5DataCL::NumNodesC; ++j)
        {
            const BaryCoordCL& Node= Quad5DataCL::Node[j];
            GRef[i][j]= FE_P2CL::DHRef( i, Node[1], Node[2], Node[3]);
        }
}

void P2DiscCL::GetGradientsOnRef( Quad5_2DCL<Point3DCL> GRef[10],
    const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    Quad5_2DCL<Point3DCL>::SetInterface( p, NodeInTetra);
    for (int i= 0; i < 10; ++i)
        for (int j= 0; j < Quad5_2DDataCL::NumNodesC; ++j) {
            const BaryCoordCL& Node= NodeInTetra[j];
            GRef[i][j]= FE_P2CL::DHRef( i, Node[1], Node[2], Node[3]);
        }
}

void P2DiscCL::GetP2Basis( Quad5_2DCL<> p2[10], const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    Quad5_2DCL<>::SetInterface( p, NodeInTetra);
    for (int j= 0; j < Quad5_2DDataCL::NumNodesC; ++j) {
        const BaryCoordCL& Node= NodeInTetra[j];
        p2[0][j]= FE_P2CL::H0( Node);
        p2[1][j]= FE_P2CL::H1( Node);
        p2[2][j]= FE_P2CL::H2( Node);
        p2[3][j]= FE_P2CL::H3( Node);
        p2[4][j]= FE_P2CL::H4( Node);
        p2[5][j]= FE_P2CL::H5( Node);
        p2[6][j]= FE_P2CL::H6( Node);
        p2[7][j]= FE_P2CL::H7( Node);
        p2[8][j]= FE_P2CL::H8( Node);
        p2[9][j]= FE_P2CL::H9( Node);
    }
}

void P1DiscCL::GetP1Basis( Quad5_2DCL<> p1[4], const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    Quad5_2DCL<>::SetInterface( p, NodeInTetra);
    for (int j= 0; j < Quad5_2DDataCL::NumNodesC; ++j) {
        const BaryCoordCL& Node= NodeInTetra[j];
        p1[0][j]= FE_P1CL::H0( Node);
        p1[1][j]= FE_P1CL::H1( Node);
        p1[2][j]= FE_P1CL::H2( Node);
        p1[3][j]= FE_P1CL::H3( Node);
    }
}

void P2RidgeDiscCL::GetEnrichmentFunction( LocalP2CL<>& ridgeFunc_p, LocalP2CL<>& ridgeFunc_n, const LocalP2CL<>& lset)
/// initialize ridge enrichment function (to be interpreted as isoP2 function = P1 on child)
{
    int sgn[10]; // sign pattern of lset
    for (int v=0; v<4; ++v) {
        const double absval= std::abs(lset[v]);
        sgn[v]= InterfacePatchCL::Sign(lset[v]);
        ridgeFunc_p[v]= sgn[v]== 1 ? 0 : 2*absval;
        ridgeFunc_n[v]= sgn[v]==-1 ? 0 : 2*absval;
    }

    for (int e=0; e<6; ++e) {
        sgn[e+4]= InterfacePatchCL::Sign(lset[e+4]);
        const int v1= VertOfEdge(e,0),
                  v2= VertOfEdge(e,1);
        const bool vzw= sgn[v1] != sgn[v2] || sgn[v1] != sgn[e+4]; // change of sign
        const double eVal= lset[e+4],
               absModLset= vzw ? std::max( std::abs(eVal), std::max(std::abs(lset[v1]), std::abs(lset[v2]))) : std::abs(eVal);
        ridgeFunc_p[e+4]= absModLset - eVal;
        ridgeFunc_n[e+4]= absModLset + eVal;
    }
}

void P2RidgeDiscCL::GetExtBasisOnChildren( LocalP2CL<> p1ridge_p[4][8], LocalP2CL<> p1ridge_n[4][8], const LocalP2CL<>& lset)
/// returns extended basis functions per child on pos./neg. part (P2 on child)
{
    LocalP2CL<> Fabs_p, Fabs_n;       // enrichment function
    LocalP2CL<> extFabs_p, extFabs_n; // extension of enrichment function from child to parent
    GetEnrichmentFunction( Fabs_p, Fabs_n, lset);
    for (int ch= 0; ch < 8; ++ch) {
        // extend P1 values on child (as Fabs has to be interpreted as isoP2 function) to whole parent
        ExtendP1onChild( Fabs_p, ch, extFabs_p);
        ExtendP1onChild( Fabs_n, ch, extFabs_n);
        for (int i=0; i<4; ++i) { // init extended basis functions: p1r = p1 * Fabs
            LocalP2CL<> &p1r_p= p1ridge_p[i][ch],
                        &p1r_n= p1ridge_n[i][ch];
            // p1_i[i] == 1
            p1r_p[i]= extFabs_p[i];
            p1r_n[i]= extFabs_n[i];
            // p1_i == 0 on opposite face
            const Ubyte oppFace= OppFace(i);
            for (Ubyte j=0; j<3; ++j) {
                const Ubyte vf= VertOfFace(oppFace,j),
                            ef= EdgeOfFace(oppFace,j) + 4;
                p1r_p[vf]= 0;
                p1r_p[ef]= 0;
                p1r_n[vf]= 0;
                p1r_n[ef]= 0;
                // p1_i == 0.5 on edges connecting vert i with opposite face
                const Ubyte e= EdgeByVert(i,vf) + 4;
                p1r_p[e]= 0.5*extFabs_p[e];
                p1r_n[e]= 0.5*extFabs_n[e];
            }
        }
    }
}

void P2RidgeDiscCL::GetExtBasisPointwise( LocalP2CL<> p1ridge_p[4], LocalP2CL<> p1ridge_n[4], const LocalP2CL<>& lset)
/// returns extended basis functions on pos./neg. part (to be interpreted pointwise in P2 degrees of freedom)
{
    LocalP2CL<> Fabs_p, Fabs_n;       // enrichment function
    GetEnrichmentFunction( Fabs_p, Fabs_n, lset);
    for (int i=0; i<4; ++i) { // init extended basis functions: p1r = p1 * Fabs
        LocalP2CL<> &p1r_p= p1ridge_p[i],
                    &p1r_n= p1ridge_n[i];
        // p1_i[i] == 1
        p1r_p[i]= Fabs_p[i];
        p1r_n[i]= Fabs_n[i];
        // p1_i == 0 on opposite face
        const Ubyte oppFace= OppFace(i);
        for (Ubyte j=0; j<3; ++j) {
            const Ubyte vf= VertOfFace(oppFace,j),
                        ef= EdgeOfFace(oppFace,j) + 4;
            p1r_p[vf]= 0;
            p1r_p[ef]= 0;
            p1r_n[vf]= 0;
            p1r_n[ef]= 0;
            // p1_i == 0.5 on edges connecting vert i with opposite face
            const Ubyte e= EdgeByVert(i,vf) + 4;
            p1r_p[e]= 0.5*Fabs_p[e];
            p1r_n[e]= 0.5*Fabs_n[e];
        }
    }
}

void P2RtoP2( const IdxDescCL& p2ridx, const VectorCL& p2r, const IdxDescCL& p2idx, VectorCL& posPart, VectorCL& negPart, const VecDescCL& lset, const BndDataCL<>& lsetbnd, const MultiGridCL& mg)
{
    const ExtIdxDescCL& extIdx= p2ridx.GetXidx();
    const int lvl= p2ridx.TriangLevel();
    const int stride= p2ridx.NumUnknownsVertex();
    const size_t p2unknowns = extIdx.GetNumUnknownsStdFE();
    if (p2unknowns != p2idx.NumUnknowns())
        throw DROPSErrCL( "P2RtoP2: inconsistent indices\n");

    LocalP2CL<> loc_phi;
    LocalP2CL<> Fabs_p, Fabs_n;    // enrichment function
    LocalP1CL<Point3DCL> linear;   // linear function to be mulitplied with enrichment function
    LocalP2CL<Point3DCL> ext_p, ext_n, linP2; // extended part which comes additionally to the std part
    LocalNumbP2CL numb; // std P2 numbering
    IdxT xnr[4];        // extended numbering
    InterfaceTetraCL patch;

    negPart.resize(p2unknowns);
    posPart.resize(p2unknowns);
    posPart = negPart = p2r[std::slice(0, p2unknowns, 1)]; // copy std part

    for (MultiGridCL::const_TriangTetraIteratorCL sit = mg.GetTriangTetraBegin(lvl), send=mg.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        loc_phi.assign( *sit, lset, lsetbnd);
        patch.Init( *sit, loc_phi);
        if (patch.Intersects()) { // we are at the phase boundary.
            // compute extended part
            P2RidgeDiscCL::GetEnrichmentFunction( Fabs_p, Fabs_n, loc_phi);
            numb.assign( *sit, p2ridx, p2ridx.GetBndInfo());
            for (int v=0; v<4; ++v) {
                xnr[v]= numb.num[v]!=NoIdx ? extIdx[numb.num[v]] : NoIdx;
            }
            for (int v=0; v<4; ++v)
                if (xnr[v] != NoIdx) {
                    for (int k=0; k<stride; ++k)
                        linear[v][k]= p2r[xnr[v]+k];
                }
            linP2.assign( linear);
            ext_p= Fabs_p*linP2;
            ext_n= Fabs_n*linP2;
            // write posPart/negPart
            for (int i=0; i<10; ++i) {
                const IdxT nr= numb.num[i];
                if (nr!=NoIdx)
                    for (int k=0; k<stride; ++k) {
                        posPart[nr+k]= p2r[nr+k] + ext_p[i][k];
                        negPart[nr+k]= p2r[nr+k] + ext_n[i][k];
                    }
            }
        }
    }
}


} // end of namespace DROPS
