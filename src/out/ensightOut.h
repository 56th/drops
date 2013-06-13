/// \file ensightOut.h
/// \brief solution output in Ensight6 Case format
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

/// \file
/// \todo (merge) EnsightOutCL differs a lot!

#ifndef DROPS_ENSIGHTOUT_H
#define DROPS_ENSIGHTOUT_H

#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <map>
#include "geom/multigrid.h"
#include "misc/problem.h"
#include "num/fe.h"

#ifndef _PAR
namespace DROPS
{

/// \brief Helper union for writing integers in binary files
union showInt
{
    int i;
    char s[sizeof(int)];
};

/// \brief Helper union for writing floats in binary files
union showFloat
{
    float f;
    char s[sizeof(float)];
};


class Ensight6VariableCL; //forward declaration

/// \brief Class for writing out results of a simulation in Ensight6 Case format.
///
/// Register subclasses of Ensight6VariableCL to output the geometry and scalar/vector-valued functions
/// in Ensight6-Case format.
class Ensight6OutCL
{
  private:
    char               decDigits_; ///< Number of digits in the decimal representation of numsteps_
    int                timestep_;  ///< Current timestep
    int                numsteps_;  ///< Total number of timesteps (required for the names of transient output files)
    double             time_;      ///< Time of current timestep
    std::ostringstream geomdesc_,  ///< Geometry-section of the Case-file
                       vardesc_,   ///< Variable-section of the Case-file
                       timestr_;   ///< Time-section of the Case-file
    std::string        casefile_;  ///< Name of the Case-file
    const bool         binary_;    ///< type of output
    bool               timedep_;   ///< true, if there are time-dependent variables
    std::map<std::string, Ensight6VariableCL*> vars_;        ///< The variables and geometry stored by varName.

    /// \brief Internal helper
    ///@{
    void OpenFile       (std::ofstream& of, std::string varName); ///< Append timecode for transient output and check the stream
    bool putTime        (double t);                               ///< Advance timestep_, time_, timestr_, if t_> time_ and set time_= t_; returns true, if t_>time_.
    void CheckFile      (const std::ofstream&) const;
    void CommitCaseFile ();                                       ///< (Re)write case file
    ///@}

  public:
    Ensight6OutCL  (std::string casefileName, Uint numsteps= 0, bool binary= true);
    ~Ensight6OutCL ();

    /// \brief Register a variable or the geometry for output with Write().
    ///
    /// The class takes ownership of the objects, i. e. it destroys them with delete in its destructor.
    void Register (Ensight6VariableCL& var);
    /// \brief Write the registered Ensight6-variables
    ///
    /// For t==0, write all registered objects to their files;if t>0 and t has increased with regard
    /// to the last call, write all time-dependent objects. Only the first of multiple calls with identical t has an effect.
    void Write (double t= 0.);

    /// \brief Append the current timestep-value with the required number of leading '0' to str.
    void AppendTimecode(std::string& str) const;

    /// \brief Interface for classes that implement the Ensight6VariableCL-interface, i.e. output of specific varibles; should probably be private.
    ///@{
    /// \brief Describe a geometry model
    void DescribeGeom (std::string geoName);
    /// \brief Describe a finite element function
    void DescribeVariable (std::string varName, bool isscalar);

    /// \brief Write the geometry into a file
    void putGeom   (MultiGridCL& mg, int lvl, std::string geoName);
    /// \brief Write a scalar value finite element function into a file
    template<class DiscScalT>
    void putScalar (const DiscScalT& v, std::string varName);
    /// \brief Write a vector value finite element function into a file
    template<class DiscVecT>
    void putVector (const DiscVecT& v, std::string varName);
    ///@}
};

/// \brief Base-class for the output of a single function in Ensight6 Case format.
///
/// We employ the command pattern: 'Describe' is the interface for registration in Ensight6OutCL.
/// 'put' is called for the output of the function at time t. The command objects are stored in Ensight6OutCL.
class Ensight6VariableCL
{
  private:
    std::string varName_,
                fileName_;
    bool        timedep_;

  public:
    Ensight6VariableCL (std::string varName, std::string fileName, bool timedep)
        : varName_( varName), fileName_( fileName), timedep_( timedep) {}
    virtual ~Ensight6VariableCL () {}

    std::string varName  () const { return varName_; }  ///< Name of the variable in einsight; also used as identifier in Ensight6OutCL.
    std::string fileName () const { return fileName_; } ///< Name of the file; for time-dependent objects, the timecode is attached by Ensight6OutCL.
    bool        Timedep  () const { return timedep_; }  ///< Is the object time-dependent?

    /// \brief Called by Ensight6OutCL::Register().
    virtual void Describe (Ensight6OutCL&) const= 0;
    /// \brief Called by Ensight6OutCL::Write().
    virtual void put      (Ensight6OutCL&) const= 0;
    /// \brief In the parallel version, the geometry information must be written before other variables.
    virtual bool is_geom  ()               const { return false; }
};

///\brief Output a geometry.
///
/// This outputs a triangulation of a multigrid.
class Ensight6GeomCL : public Ensight6VariableCL
{
  private:
    MultiGridCL* mg_;  ///< The multigrid
    int          lvl_; ///< Level of the triangulation

  public:
    Ensight6GeomCL (MultiGridCL& mg, int lvl, std::string varName, std::string fileName, bool timedep= false)
        : Ensight6VariableCL( varName, fileName, timedep), mg_( &mg), lvl_( lvl) {}

    void Describe (Ensight6OutCL& cf) const { cf.DescribeGeom( this->varName());  }
    void put      (Ensight6OutCL& cf) const { cf.putGeom( *mg_, lvl_, varName()); }
    bool is_geom  ()                  const { return true; }
};

///\brief Create an Ensight6GeomCL with operator new.
///
/// This is just for uniform code; the analoguous functions for scalars and vectors are more useful because
/// they help to avoid template parameters in user code.
inline Ensight6GeomCL&
make_Ensight6Geom (MultiGridCL& mg, int lvl, std::string varName, std::string fileName, bool timedep= false)
{
    return *new Ensight6GeomCL( mg, lvl, varName, fileName, timedep);
}

///\brief Represents a scalar Drops-function (P1 or P2, given as PXEvalCL) as Ensight6 variable.
template <class DiscScalarT>
class Ensight6ScalarCL : public Ensight6VariableCL
{
  private:
    const DiscScalarT f_;

  public:
    Ensight6ScalarCL (const DiscScalarT& f, std::string varName, std::string fileName, bool timedep= false)
        : Ensight6VariableCL( varName, fileName, timedep), f_( f) {}

    void Describe (Ensight6OutCL& cf) const { cf.DescribeVariable( this->varName(), true); }
    void put      (Ensight6OutCL& cf) const { cf.putScalar( f_, varName()); }
};

///\brief Create an Ensight6ScalarCL<> with operator new.
///
/// This function does the template parameter deduction for user code.
template <class DiscScalarT>
  Ensight6ScalarCL<DiscScalarT>&
    make_Ensight6Scalar (const DiscScalarT& f, std::string varName, std::string fileName, bool timedep= false)
{
    return *new Ensight6ScalarCL<DiscScalarT>( f, varName, fileName, timedep);
}

///\brief Represents a vector Drops-function (P1 or P2, given as PXEvalCL) as Ensight6 variable.
template <class DiscVectorT>
class Ensight6VectorCL : public Ensight6VariableCL
{
  private:
    const DiscVectorT f_;

  public:
    Ensight6VectorCL (const DiscVectorT& f, std::string varName, std::string fileName, bool timedep= false)
        : Ensight6VariableCL( varName, fileName, timedep), f_( f) {}

    void Describe (Ensight6OutCL& cf) const { cf.DescribeVariable( this->varName(), false); }
    void put      (Ensight6OutCL& cf) const { cf.putVector( f_, varName()); }
};

///\brief Create an Ensight6VectorCL<> with operator new.
///
/// This function does the template parameter deduction for user code.
template <class DiscVectorT>
  Ensight6VectorCL<DiscVectorT>&
    make_Ensight6Vector (const DiscVectorT& f, std::string varName, std::string fileName, bool timedep= false)
{
    return *new Ensight6VectorCL<DiscVectorT>( f, varName, fileName, timedep);
}

///\brief Represents a scalar P1X function as two Ensight6 variables.
///
/// The class registers two P1-functions, neg and pos, which do the actual output. In order to prepare the data
/// needed by them, it registers itself with a put-routine that sets up the piecewise P1-data. Note, that we depend
/// on the order
/// in which Ensight6OutCL calls the routines: Lexicographically (a std::map<> varName --> Ensight6VariableCL*), the
/// P1-variables follow after the variable itself as we append "Neg" resp. "Pos" to its name.
///
/// The filenames are also generated by appending "Neg" and "Pos" to the provided filename.
class Ensight6P1XScalarCL : public Ensight6VariableCL
{
  private:
    const VecDescCL& v_;

    mutable IdxDescCL p1idx_;
    mutable VecDescCL vneg_,
                      vpos_;

    const VecDescCL& lset_;
    BndDataCL<> bnd_;
    MultiGridCL& mg_;

  public:
    Ensight6P1XScalarCL (MultiGridCL& mg, const VecDescCL& lset, const VecDescCL& v, const BndDataCL<>& bnd,
        std::string varName, std::string fileName, bool timedep= false)
        : Ensight6VariableCL( varName, fileName, timedep),
          v_( v), vneg_( &p1idx_), vpos_( &p1idx_), lset_( lset), bnd_( bnd), mg_( mg) {}
    ~Ensight6P1XScalarCL () { if (p1idx_.NumUnknowns() != 0) p1idx_.DeleteNumbering( mg_); }

    void Describe (Ensight6OutCL& cf) const {
        cf.Register( make_Ensight6Scalar( make_P1Eval( mg_, bnd_, vneg_), varName()+"Neg", fileName()+"Neg", Timedep()));
        cf.Register( make_Ensight6Scalar( make_P1Eval( mg_, bnd_, vpos_), varName()+"Pos", fileName()+"Pos", Timedep()));
    }
    void put      (Ensight6OutCL&)    const {
        if (p1idx_.NumUnknowns() != 0)
            p1idx_.DeleteNumbering( mg_);
        p1idx_.CreateNumbering( v_.RowIdx->TriangLevel(), mg_, *v_.RowIdx);
        P1XtoP1 ( *v_.RowIdx, v_.Data, p1idx_, vpos_.Data, vneg_.Data, lset_, mg_);
    }
};

///\brief Create an Ensight6P1XScalarCL with operator new.
///
/// This is just for uniform code; the analoguous functions for scalars and vectors are more useful because
/// they help to avoid template parameters in user code.
inline Ensight6P1XScalarCL&
make_Ensight6P1XScalar (MultiGridCL& mg, const VecDescCL& lset, const VecDescCL& v,
    std::string varName, std::string fileName, bool timedep= false)
{
    return *new Ensight6P1XScalarCL( mg, lset, v, BndDataCL<>( 0), varName, fileName, timedep);
}

inline Ensight6P1XScalarCL&
make_Ensight6P1XScalar (MultiGridCL& mg, const VecDescCL& lset, const VecDescCL& v, const BndDataCL<>& bnd,
    std::string varName, std::string fileName, bool timedep= false)
{
    return *new Ensight6P1XScalarCL( mg, lset, v, bnd, varName, fileName, timedep);
}

///\brief Represents a vector-valued P2R function as two Ensight6 variables.
///
/// The class registers two P2-functions, neg and pos, which do the actual output. In order to prepare the data
/// needed by them, it registers itself with a put-routine that sets up the piecewise P2-data. Note, that we depend
/// on the order
/// in which Ensight6OutCL calls the routines: Lexicographically (a std::map<> varName --> Ensight6VariableCL*), the
/// P2-variables follow after the variable itself as we append "Neg" resp. "Pos" to its name.
///
/// The filenames are also generated by appending "Neg" and "Pos" to the provided filename.
class Ensight6P2RVectorCL : public Ensight6VariableCL
{
  private:
    const VecDescCL& v_;

    mutable IdxDescCL p2idx_;
    mutable VecDescCL vneg_,
                      vpos_;

    const VecDescCL& lset_;
    const BndDataCL<> lsetbnd_;
    BndDataCL<Point3DCL> bnd_;
    MultiGridCL& mg_;

  public:
    Ensight6P2RVectorCL (MultiGridCL& mg, const VecDescCL& lset, const BndDataCL<>& lsetbnd, const VecDescCL& v, const BndDataCL<Point3DCL>& bnd,
        std::string varName, std::string fileName, bool timedep= false)
        : Ensight6VariableCL( varName, fileName, timedep),
          v_( v), p2idx_(vecP2_FE), vneg_( &p2idx_), vpos_( &p2idx_), lset_( lset), lsetbnd_( lsetbnd), bnd_( bnd), mg_( mg) {}
    ~Ensight6P2RVectorCL () { if (p2idx_.NumUnknowns() != 0) p2idx_.DeleteNumbering( mg_); }

    void Describe (Ensight6OutCL& cf) const {
        cf.Register( make_Ensight6Vector( make_P2Eval( mg_, bnd_, vneg_), varName()+"Neg", fileName()+"Neg", Timedep()));
        cf.Register( make_Ensight6Vector( make_P2Eval( mg_, bnd_, vpos_), varName()+"Pos", fileName()+"Pos", Timedep()));
    }
    void put      (Ensight6OutCL&)    const;
};

///\brief Create an Ensight6P2RVectorCL with operator new.
///
/// This is just for uniform code; the analoguous functions for scalars and vectors are more useful because
/// they help to avoid template parameters in user code.
inline Ensight6P2RVectorCL&
make_Ensight6P2RVector (MultiGridCL& mg, const VecDescCL& lset, const BndDataCL<>& lsetbnd, const VecDescCL& v, const BndDataCL<Point3DCL>& bnd,
    std::string varName, std::string fileName, bool timedep= false)
{
    return *new Ensight6P2RVectorCL( mg, lset, lsetbnd, v, bnd, varName, fileName, timedep);
}


/// \brief Class for writing out 2 phase flows in ensight format
template<typename StokesT, typename LevelsetT>
class Ensight2PhaseOutCL
{
  private:
    Ensight6OutCL ensight_;
    bool adaptive_;                                         // changing geometry
    bool timedep_;                                          // time-dependent problem
    const StokesT&   stokes_;                               // (Navier-)Stokes problem
    const LevelsetT& lset_;                                 // Levelset problem

  public:
    /**
       \brief Construct a class for writing out a 2 phase flow problem in ensight format

       All information to this class is given by this constructor:
       \param mg Multigrid of the problem
       \param idx P2 index class
       \param stokes Stokes problem class
       \param lset   Levelset problem class
       \param directory directory of ensight files
       \param caseName  name of the case file
       \param geomName name of the geometry
       \param numsteps number of time steps
       \param adaptive flag if the geometry changes over time
       \param binary    binary of ascii output
       \param masterout (in parallel) flag if the output should be done by master
    */

    Ensight2PhaseOutCL( MultiGridCL& mg, const IdxDescCL* idx,
                        const StokesT& stokes, const LevelsetT& lset,
                        const std::string& directory, const std::string& caseName,
                        const std::string& geomName, bool adaptive,
                        Uint numsteps=0, bool binary=false)
        : ensight_(caseName + ".case", numsteps, binary), adaptive_(adaptive),
          timedep_(numsteps>0), stokes_(stokes), lset_(lset)
    {
        std::string ensf( directory + "/" + caseName);
        ensight_.Register( make_Ensight6Geom      ( mg, idx->TriangLevel(),  geomName,        ensf + ".geo", adaptive_));
        ensight_.Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",      ensf + ".scl", timedep_));
        ensight_.Register( make_Ensight6Scalar    ( stokes_.GetPrSolution(),  "Pressure",      ensf + ".pr",  timedep_));
        ensight_.Register( make_Ensight6Vector    ( stokes_.GetVelSolution(), "Velocity",      ensf + ".vel", timedep_));
    }

    /// \brief Write ensight files for a time step
    void write(){
        ensight_.Write( timedep_  ? stokes_.t : -1);
    }
};

class ReadEnsightP2SolCL
// read solution from Ensight6 Case format
{
  private:
    const MultiGridCL* _MG;
    const bool         binary_;

    void CheckFile( const std::ifstream&) const;

  public:
    ReadEnsightP2SolCL( const MultiGridCL& mg, bool binary=true)
      : _MG(&mg), binary_(binary)
    {}

    template<class BndT>
    void ReadScalar( const std::string&, VecDescCL&, const BndT&) const;
    template<class BndT>
    void ReadVector( const std::string&, VecDescCL&, const BndT&) const;
};


//=====================================================
//              template definitions
//=====================================================

template<class DiscScalT>
void Ensight6OutCL::putScalar (const DiscScalT& v, std::string varName)
{
    const MultiGridCL& mg= v.GetMG();
    const Uint lvl= v.GetLevel();
    char buffer[80];
    std::memset(buffer,0,80);
    showFloat sFlo;

    std::ofstream os;
    OpenFile( os, varName);

    if(binary_)
    {
        std::strcpy(buffer,"DROPS data file, scalar variable:");
        os.write(buffer,80);

        DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
            sFlo.f=v.val(*it);
            os.write(sFlo.s,sizeof(float));
        }

        DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
            sFlo.f=v.val(*it,0.5);
            os.write(sFlo.s,sizeof(float));
        }
    }
    else //ASCII-Ausgabe
    {
        int cnt=0;
        os.flags(std::ios_base::scientific);
        os.precision(5);
        os.width(12);

        os << "DROPS data file, scalar variable:\n";
        DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
            os << std::setw(12) << v.val( *it);
            if ( (++cnt)==6)
            { // Ensight expects six real numbers per line
                cnt= 0;
                os << '\n';
            }
        }

        DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
            os << std::setw(12) << v.val( *it, 0.5);
            if ( (++cnt)==6)
            { // Ensight expects six real numbers per line
                cnt= 0;
                os << '\n';
            }
        }
        os << '\n';
    }
}

template<class DiscVecT>
void Ensight6OutCL::putVector (const DiscVecT& v, std::string varName)
{
    const MultiGridCL& mg= v.GetMG();
    const Uint lvl= v.GetLevel();
    char buffer[80];
    std::memset(buffer,0,80);
    showFloat sFlo;

    std::ofstream os;
    OpenFile( os, varName);

    if(binary_)
    {
        std::strcpy(buffer,"DROPS data file, vector variable:");
        os.write( buffer, 80);

        DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
            for (int i=0; i<3; ++i)
            {
                sFlo.f=v.val( *it)[i];
                os.write(sFlo.s,sizeof(float));
            }
        }
        DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
            for (int i=0; i<3; ++i)
            {
                 sFlo.f=v.val( *it)[i];
                 os.write(sFlo.s,sizeof(float));
            }
        }
    }
    else
    { // ASCII
        int cnt=0;
        os.flags(std::ios_base::scientific);
        os.precision(5);
        os.width(12);

        os << "DROPS data file, vector variable:\n";

        DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
            for (int i=0; i<3; ++i)
                os << std::setw(12) << v.val( *it)[i];
            if ( (cnt+=3)==6)
            { // Ensight expects six real numbers per line
                 cnt= 0;
                 os << '\n';
            }
        }
        DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
            for (int i=0; i<3; ++i)
                os << std::setw(12) << v.val( *it)[i];
            if ( (cnt+=3)==6)
            { // Ensight expects six real numbers per line
                 cnt= 0;
                 os << '\n';
            }
        }
        os << '\n';
    }
}


// ========== ReadEnsightP2SolCL ==========

template <class BndT>
void ReadEnsightP2SolCL::ReadScalar( const std::string& file, VecDescCL& v, const BndT& bnd) const
{
    const Uint lvl= v.GetLevel(),
               idx= v.RowIdx->GetIdx();
    std::string fileName(file);

    std::ifstream is( fileName.c_str());
    CheckFile( is);

    if (binary_)
    {
        showFloat fl;
        char buffer[80];
        is.read( buffer, 80);           //ignore first 80 characters

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
             end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            is.read( fl.s, sizeof(float));
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            v.Data[it->Unknowns(idx)]= (double)fl.f;
        }

        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            is.read( fl.s, sizeof(float));
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            v.Data[it->Unknowns(idx)]= fl.f;
        }
    }
    else
    { // ASCII
        char buf[256];
        double d= 0;

        is.getline( buf, 256); // ignore first line

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
             end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            is >> d;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            v.Data[it->Unknowns(idx)]= d;
        }
        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            is >> d;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            v.Data[it->Unknowns(idx)]= d;
        }
    }

    CheckFile( is);
}

template <class BndT>
void ReadEnsightP2SolCL::ReadVector( const std::string& file, VecDescCL& v, const BndT& bnd) const
{
    const Uint lvl= v.GetLevel(),
               idx= v.RowIdx->GetIdx();
    std::string fileName(file);

    std::ifstream is( fileName.c_str());
    CheckFile( is);

    double d0= 0, d1= 0, d2= 0;

    if(binary_)
    {
        showFloat fl;
        //std::cout<<"READVECTOR: "<<file.c_str()<<"\n";
        char buffer[80];
        is.read( buffer, 80);       //ignore first 80 characters

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
            end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            is.read( fl.s, sizeof(float));
            d0=fl.f;
            is.read( fl.s, sizeof(float));
            d1=fl.f;
            is.read( fl.s, sizeof(float));
            d2=fl.f;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            const IdxT Nr= it->Unknowns(idx);
            v.Data[Nr]= d0;    v.Data[Nr+1]= d1;    v.Data[Nr+2]= d2;
        }

        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
             end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            is.read( fl.s, sizeof(float));
            d0=fl.f;
            is.read( fl.s, sizeof(float));
            d1=fl.f;
            is.read( fl.s, sizeof(float));
            d2=fl.f;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            const IdxT Nr= it->Unknowns(idx);
            v.Data[Nr]= d0;    v.Data[Nr+1]= d1;    v.Data[Nr+2]= d2;
        }
    }
    else
    { // ASCII
        char buf[256];

        is.getline( buf, 256); // ignore first line

        int count=0;
        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
            end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            is >> d0 >> d1 >> d2;
            count +=3;

            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            const IdxT Nr= it->Unknowns(idx);
            v.Data[Nr]= d0;    v.Data[Nr+1]= d1;    v.Data[Nr+2]= d2;
        }

        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            is >> d0 >> d1 >> d2;
            count +=3;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;

            const IdxT Nr= it->Unknowns(idx);
            v.Data[Nr]= d0;    v.Data[Nr+1]= d1;    v.Data[Nr+2]= d2;
        }
    }

    CheckFile( is);
}

} // end of namespace DROPS
#endif

#endif
