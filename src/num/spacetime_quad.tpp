/// \file spacetime_quad.tpp
/// \brief template implementations
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
#include "num/spacetime_quad.h"

namespace DROPS
{

template <class IntRuleData, class Geom>
void Gather4DIntegrationPoints(const std::vector<Geom> & simplex4d, GridFunctionCL<Point4DCL> & points)
{
    // ScopeTimer scopetiming("Gather4DIntegrationInformation");
    points.resize(simplex4d.size()*IntRuleData::NumNodesC);
    for(Uint p = 0; p < simplex4d.size(); ++p)
    {
        const Uint offset = p*IntRuleData::NumNodesC;
        SArrayCL<Point4DCL,IntRuleData::NumNodesC> spoints 
            = Transform4DIntegrationPoints<IntRuleData>(simplex4d[p]);
        for (Uint i = 0; i < IntRuleData::NumNodesC; ++i)
            points[offset+i] = spoints[i];
    }
}

template <class IntRuleData, class Geom>
void Gather4DIntegrationWeights(const std::vector<Geom> & simplex4d, GridFunctionCL<double> & weights)
{
    // ScopeTimer scopetiming("Gather4DIntegrationInformation");
    weights.resize(simplex4d.size()*IntRuleData::NumNodesC);
    for(Uint p = 0; p < simplex4d.size(); ++p)
    {
        const Uint offset = p*IntRuleData::NumNodesC;
        SArrayCL<double,IntRuleData::NumNodesC> pentaweights 
            = Transform4DIntegrationWeights<IntRuleData>(simplex4d[p]);
        for (Uint i = 0; i < IntRuleData::NumNodesC; ++i)
            weights[offset+i] = pentaweights[i];
    }
}


template <class IntRuleData>
void Gather4DNormals(const std::vector<Tetra4DCL> & tets, 
                                GridFunctionCL<Point4DCL> & normals)
{
    normals.resize(tets.size()*IntRuleData::NumNodesC);
    for(Uint p = 0; p < tets.size(); ++p)
    {
        const Uint offset = p*IntRuleData::NumNodesC;
        Point4DCL normal = tets[p].Normal();
        for (Uint i = 0; i < IntRuleData::NumNodesC; ++i)
            normals[offset+i] = normal;
    }
}

template <class IntRuleData>
void Gather4DNu(const std::vector<Tetra4DCL> & tets, 
                GridFunctionCL<double> & nus)
{
    nus.resize(tets.size()*IntRuleData::NumNodesC);
    for(Uint p = 0; p < tets.size(); ++p)
    {
        const Uint offset = p*IntRuleData::NumNodesC;
        const double nu = tets[p].Nu();
        for (Uint i = 0; i < IntRuleData::NumNodesC; ++i)
            nus[offset+i] = nu;
    }
}


template <class IntRuleData, class Geom>
SArrayCL<Point4DCL,IntRuleData::NumNodesC> Transform4DIntegrationPoints(const Geom & penta)
{
    SArrayCL<Point4DCL,IntRuleData::NumNodesC> ret(Point4DCL(0.0,0.0,0.0,0.0));

    for (Uint i = 0; i < IntRuleData::NumNodesC; ++i)
        for (Uint j = 0; j <= IntRuleData::Dim ; ++j)
            ret[i] += IntRuleData::Node[i][j] * *(penta.GetPoints()[j]);
    return ret;
}

template <class IntRuleData, class Geom>
SArrayCL<double,IntRuleData::NumNodesC> Transform4DIntegrationWeights(const Geom & simplex4d)
{
    SArrayCL<double,IntRuleData::NumNodesC> ret(0.0);
    const double absdet = simplex4d.AbsDet();
    for (Uint i = 0; i < IntRuleData::NumNodesC; ++i)
        ret[i] = IntRuleData::Weight[i] * absdet;
    return ret;
}

template <class IntRuleData4D>
CompositeSTQuadCL<IntRuleData4D>::CompositeSTQuadCL(const TetraCL& tet, const TimeInterval& ti,
                                                    const LocalP2CL<double>& lsetold, 
                                                    const LocalP2CL<double>& lsetnew, 
                                                    Uint ints_per_space_edge,
                                                    Uint subtimeintervals)
{
    Initialize(&tet,&ti,/*prism4*/NULL,&lsetold,&lsetnew,/*f*/NULL,ints_per_space_edge,subtimeintervals);
}

template <class IntRuleData4D>
CompositeSTQuadCL<IntRuleData4D>::CompositeSTQuadCL(const TetraCL& tet, const TimeInterval& ti,
                                                    instat_scalar_fun_ptr f,
                                                    Uint ints_per_space_edge,
                                                    Uint subtimeintervals)
{
    Initialize(&tet,&ti,/*prism4*/NULL,/*lsetold*/NULL,/*lsetnew*/NULL,f,ints_per_space_edge,subtimeintervals);
}

template <class IntRuleData4D>
CompositeSTQuadCL<IntRuleData4D>::CompositeSTQuadCL(const GeneralizedPrism4CL& refprism4,
                                                    instat_scalar_fun_ptr f,
                                                    Uint ints_per_space_edge,
                                                    Uint subtimeintervals)
{
    Initialize(/*tet*/NULL,/*time int*/NULL,&refprism4,/*lsetold*/NULL,/*lsetnew*/NULL,f,ints_per_space_edge,subtimeintervals);
}

template <class IntRuleData4D>
CompositeSTQuadCL<IntRuleData4D>::CompositeSTQuadCL(const GeneralizedPrism4CL& refprism4,
                                                    const LocalP2CL<double>& lsetold, 
                                                    const LocalP2CL<double>& lsetnew, 
                                                    Uint ints_per_space_edge,
                                                    Uint subtimeintervals)
{
    Initialize(/*tet*/NULL,/*time int*/NULL,&refprism4,&lsetold,&lsetnew,/*f*/NULL,ints_per_space_edge,subtimeintervals);
}


// Initialize creates a generalized prism4, divides it into subprism4 according 
// to the subdivision rule prescribed (ints_per_space_edge, subtimeintervals).
// Those subprism4 are decomposed into pentatopes. On these the level set function
// is evaluated. The level set function is either 
// given by the space level set function at begin and end time or directly as 
// a function pointer. Once the level set function is evaluated, each cut pentatope
// is decomposed into uncut pentatopes. All resulting pentatopes in positive / negative
// phase are gather in corresponding std::vectors, the same holds for interface tet4ds.
// Then integration points are distributed across the pentatopes and correspondingly
// integration points and weights are gathered for pos / neg phase and interface.
template <class IntRuleData4D>
void CompositeSTQuadCL<IntRuleData4D>::Initialize(const TetraCL* tet, const TimeInterval* ti,
                                                  const GeneralizedPrism4CL* refprism4,
                                                  const LocalP2CL<double>* lsetold,
                                                  const LocalP2CL<double>* lsetnew,
                                                  instat_scalar_fun_ptr f,
                                                  Uint ints_per_space_edge,
                                                  Uint subtimeintervals)
{
    ScopeTimer timingInitialize("CompositeSTQuadCL<..>::Initialize");
    GeneralizedPrism4CL backup_refprism4(pcont);
    const GeneralizedPrism4CL& prefprism4 = refprism4 != NULL ? *refprism4 : backup_refprism4; 
    std::vector<PentatopeCL> mainpentas; mainpentas.clear();
    negpentas.clear(), pospentas.clear(), iftetras.clear();
    prefprism4.decompose_to_pentas(mainpentas, 
                                    /*principal lattice intervals per edge=*/ ints_per_space_edge,
                                    /*subintervals per time intervals=*/ subtimeintervals);
    SArrayCL<double,5> lsetvals;

    for(std::vector<PentatopeCL>::iterator it = mainpentas.begin(); 
        it != mainpentas.end(); 
        it++)
    {
        if (lsetold !=NULL)
            // using two LocalP2 functions
            (*it).eval_timeinterpol_func_at_verts(*lsetold,*lsetnew,lsetvals);
        else
            // using a given (space,time)-function
            if (tet != NULL && ti != NULL)
                (*it).eval_func_at_verts(f, *tet, *ti, lsetvals);
            else
                throw DROPSErrCL("Cannot evaluate arbitrary function yet ..");
        (*it).decompose_add_signed_pentas_4dtets(lsetvals,negpentas,pospentas,iftetras);
    }

    const bool hasnegvals = (negpentas.size() > 0); 
    const bool hasposvals = (pospentas.size() > 0); 

    Gather4DIntegrationPoints<typename IntRuleData4D::Volume>(negpentas,ips_neg);
    Gather4DIntegrationWeights<typename IntRuleData4D::Volume>(negpentas,ipw_neg);
    Gather4DIntegrationPoints<typename IntRuleData4D::Volume>(pospentas,ips_pos);
    Gather4DIntegrationWeights<typename IntRuleData4D::Volume>(pospentas,ipw_pos);

    p1p1signs.resize(8);

    hasinterface = false; // if not changed hereafter

    if (hasposvals && hasnegvals)
    {
        if (iftetras.size() == 0){
            // std::cout << " WARNING: no interface tetras, but interface " << std::endl;
        }
        else
        {
            hasinterface = true;
            Gather4DIntegrationPoints<typename IntRuleData4D::Surface>(iftetras,ips_if);
            Gather4DIntegrationWeights<typename IntRuleData4D::Surface>(iftetras,ipw_if);
            Gather4DNormals<typename IntRuleData4D::Surface>(iftetras,ipn_if);
            Gather4DNu<typename IntRuleData4D::Surface>(iftetras,ipnu_if);
            ipsn_if.resize(ipn_if.size());
            for (Uint i = 0; i < ipn_if.size(); ++i)
            {
                for (int j = 0; j < 3; ++j)
                    ipsn_if[i][j] = ipn_if[i][j] * 1.0/ ipnu_if[i];
            }

            if (lsetold != NULL)
            {
                for (Uint i = 0 ; i < 4; ++i)
                    p1p1signs[i] = ((*lsetold)[i] >= 0);
                for (Uint i = 0 ; i < 4; ++i)
                    p1p1signs[4+i] = ((*lsetnew)[i] >= 0);
            }
            else
            {
                for (Uint i = 0 ; i < 4; ++i)
                    p1p1signs[i] = (f(tet->GetVertex(i)->GetCoord(),ti->first) >= 0);
                for (Uint i = 0 ; i < 4; ++i)
                    p1p1signs[4+i] = (f(tet->GetVertex(i)->GetCoord(),ti->second) >= 0);
            }
        }
    }

    if (!hasinterface)
    {
        ips_if.resize(0);
        ipw_if.resize(0);
        ipn_if.resize(0);
        ipsn_if.resize(0);
        ipnu_if.resize(0);
        p1p1signs = hasposvals;
    }
    
}

template <class IntRuleData4D>
void CompositeSTQuadCL<IntRuleData4D>::Report(std::ostream & out, std::string linehead, std::string linetail) const
{
    double posmeas = 0, negmeas = 0, ifmeas = 0;
    for (Uint i = 0; i < negpentas.size(); i++)
        negmeas += negpentas[i].Measure();
    
    for (Uint i = 0; i < pospentas.size(); i++)
        posmeas += pospentas[i].Measure();

    for (Uint i = 0; i < iftetras.size(); i++)
        ifmeas += iftetras[i].Measure();

    out << linehead << "number of  if. tetras = " << iftetras.size() << linetail << std::endl;
    out << linehead << "number of pos. pentas = " << pospentas.size() << linetail << std::endl;
    out << linehead << "number of neg. pentas = " << negpentas.size() << linetail << std::endl;
    out << linehead << "total neg measure = " << negmeas << linetail << std::endl;
    out << linehead << "total pos measure = " << posmeas << linetail << std::endl;
    out << linehead << "total     measure = " << posmeas+negmeas << linetail << std::endl;
    out << linehead << "total if. measure = " << ifmeas << linetail << std::endl;
    pcont.Report(out);
}

template <class IntRuleData4D>
const GridFunctionCL<Point4DCL> & CompositeSTQuadCL<IntRuleData4D>::GetVolumeIntegrationPoints ( bool posPart) const
{
    if (posPart)
        return ips_pos;
    else
        return ips_neg;
}

template <class IntRuleData4D>
const GridFunctionCL<double> & CompositeSTQuadCL<IntRuleData4D>::GetVolumeIntegrationWeights  (bool posPart) const 
{
    if (posPart)
        return ipw_pos;
    else
        return ipw_neg;
}

template <class IntRuleData4D>
Uint CompositeSTQuadCL<IntRuleData4D>::NumberOfIntegrationPoints ( bool posPart) const
{
    if (posPart)
        return ips_pos.size();
    else
        return ips_neg.size();
}

template <class IntRuleData4D>
const GridFunctionCL<Point4DCL> & CompositeSTQuadCL<IntRuleData4D>::GetInterfaceIntegrationPoints ( ) const
{
    return ips_if;
}

template <class IntRuleData4D>
const GridFunctionCL<double> & CompositeSTQuadCL<IntRuleData4D>::GetInterfaceIntegrationWeights  ( ) const 
{
    return ipw_if;
}

template <class IntRuleData4D>
Uint CompositeSTQuadCL<IntRuleData4D>::NumberOfInterfaceIntegrationPoints ( ) const
{
    return ips_if.size();
}


// Eval functions provide the evaluation of space-discrete and linear-in-time functions or 
// functions given as function pointers on the integrations points of CompositeSTQuadCL
// (after CompositeSTQuadCL as been initialized)

template <class IntRuleData4D>
GridFunctionCL<double> CompositeSTQuadCL<IntRuleData4D>::Eval ( instat_scalar_fun_ptr f, 
                                                                const GridFunctionCL<Point4DCL> & points,
                                                                const SpaceTimeMapping * map) const
{
    ScopeTimer scopetiming("CompositeSTQuadCL<..>::Eval(f)");
    GridFunctionCL<double> vals(points.size());
    for (Uint i = 0; i < points.size(); ++i){
        Point4DCL p(map->Map(points[i]));
        vals[i] = f(MakePoint3D(p[0],p[1],p[2]),p[3]);
    }
    return vals;
}

template <class IntRuleData4D>
GridFunctionCL<Point3DCL> CompositeSTQuadCL<IntRuleData4D>::Eval ( instat_vector_fun_ptr f, 
                                                                const GridFunctionCL<Point4DCL> & points,
                                                                const SpaceTimeMapping * map) const
{
    ScopeTimer scopetiming("CompositeSTQuadCL<..>::Eval(f)");
    GridFunctionCL<Point3DCL> vals(points.size());
    for (Uint i = 0; i < points.size(); ++i){
        Point4DCL p(map->Map(points[i]));
        vals[i] = f(MakePoint3D(p[0],p[1],p[2]),p[3]);
    }
    return vals;
}

template <class IntRuleData4D> 
template <class T>
GridFunctionCL<T> CompositeSTQuadCL<IntRuleData4D>::EvalLinear ( const LocalP2CL<T>& fold, 
                                                                 const LocalP2CL<T>& fnew,
                                                                 const GridFunctionCL<Point4DCL> & points) const
{
    ScopeTimer scopetiming("CompositeSTQuadCL<..>::Eval(lp2)");
    GridFunctionCL<T> vals(points.size());
    for (Uint i = 0; i < points.size(); ++i)
    {
        const T oldval = fold(MakeBaryCoord(RestrictToPoint3D(points[i])));
        const T newval = fnew(MakeBaryCoord(RestrictToPoint3D(points[i])));
        vals[i] = (1-points[i][3])*oldval + points[i][3]*newval;
    }
    return vals;
}

template <class IntRuleData4D>
GridFunctionCL<double> CompositeSTQuadCL<IntRuleData4D>::EvalOnInterface ( instat_scalar_fun_ptr f,
                                                                           const SpaceTimeMapping * map) const
{
    return Eval(f,ips_if,map);
}

template <class IntRuleData4D>
GridFunctionCL<double> CompositeSTQuadCL<IntRuleData4D>::EvalOnPart ( instat_scalar_fun_ptr f, 
                                                                      bool posPart, 
                                                                      const SpaceTimeMapping * map) const
{
    if (posPart)
        return Eval(f,ips_pos,map);
    else
        return Eval(f,ips_neg,map);
}

template <class IntRuleData4D>
GridFunctionCL<Point3DCL> CompositeSTQuadCL<IntRuleData4D>::EvalOnPart ( instat_vector_fun_ptr f, 
                                                                      bool posPart, 
                                                                      const SpaceTimeMapping * map) const
{
    if (posPart)
        return Eval(f,ips_pos,map);
    else
        return Eval(f,ips_neg,map);
}

template <class IntRuleData4D>
template <class T>
GridFunctionCL<T> CompositeSTQuadCL<IntRuleData4D>::EvalLinearOnInterface ( const LocalP2CL<T>& fold, 
                                                                            const LocalP2CL<T>& fnew) const
{
    return EvalLinear(fold,fnew,ips_if);
}

template <class IntRuleData4D>
template <class T>
GridFunctionCL<T> CompositeSTQuadCL<IntRuleData4D>::EvalLinearOnPart ( const LocalP2CL<T>& fold, 
                                                                       const LocalP2CL<T>& fnew, 
                                                                       bool posPart) const
{
    if (posPart)
        return EvalLinear(fold,fnew,ips_pos);
    else
        return EvalLinear(fold,fnew,ips_neg);
}


template <class IntRuleData4D>
const GridFunctionCL<Point4DCL> & CompositeSTQuadCL<IntRuleData4D>::GetNormalsOnInterface ( ) const
{
    return ipn_if;
}

template <class IntRuleData4D>
const GridFunctionCL<Point3DCL> & CompositeSTQuadCL<IntRuleData4D>::GetSpaceNormalsOnInterface ( ) const
{
    return ipsn_if;
}


template <class IntRuleData4D>
const GridFunctionCL<double> & CompositeSTQuadCL<IntRuleData4D>::GetNuOnInterface ( ) const
{
    static bool first = true;
    if (first)
    {
        std::cout << " WARNING : untransformed nonsense.... " << std::endl;
        std::cout << " WARNING : if you are working on the reference element just ignore the warning.... " << std::endl;
        getchar();
        first = false;
    }
    return ipnu_if;
}

template <class IntRuleData4D>
bool CompositeSTQuadCL<IntRuleData4D>::HasInterface ( ) const
{
    return hasinterface;
}

template <class IntRuleData4D>
bool CompositeSTQuadCL<IntRuleData4D>::GetVertexSign (Uint i, bool newtime ) const
{
    return p1p1signs[(newtime?4:0)+i];
}

template <class IntRuleData4D>
const std::valarray<bool> & CompositeSTQuadCL<IntRuleData4D>::GetVertexSigns ( ) const
{
    return p1p1signs;
}

// applying quadrature simply means forming the weighted sum where the gridfunction
// describes the function values at the integration points (of one domain: pos/neg or interface)
// and weights the corresponding integration weights.
// Weights are the weights for the integration w.r.t. the reference Prism: Tx[0,t], with T the 
// reference tetrahedra. Thus all transformation weights should be included in the Gridfunction f.

template <class IntRuleData4D>
template <class T>
T CompositeSTQuadCL<IntRuleData4D>::Quad ( const GridFunctionCL<T> & f, const GridFunctionCL<double> & weights) const
{
    ScopeTimer scopetiming("CompositeSTQuadCL<..>::Quad");
    T sum(0.0);
    for (Uint i = 0; i < weights.size(); ++i)
        sum += weights[i] * f[i];
    return sum;
}

template <class IntRuleData4D>
template <class T>
T CompositeSTQuadCL<IntRuleData4D>::QuadOnPart ( const GridFunctionCL<T> & f, bool posPart) const
{
    if (posPart)
        return Quad(f,ipw_pos);
    else
        return Quad(f,ipw_neg);
}

template <class IntRuleData4D>
template <class T>
T CompositeSTQuadCL<IntRuleData4D>::QuadOnInterface ( const GridFunctionCL<T> & f) const
{
    return Quad(f,ipw_if);
}


} // end of namespace DROPS
