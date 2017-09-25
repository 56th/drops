/// \file params.h
/// \brief read parameters from file.
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt, Thorolf Schulte; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_MISC_PARAMS_H
#define DROPS_MISC_PARAMS_H

#include "misc/container.h"
#include "misc/utils.h"
#include <string>
#include <map>
#include <fstream>

namespace DROPS
{

/// \brief Specific error class for parameter-files.
class DROPSParamErrCL : public DROPSErrCL
{
  public:
    DROPSParamErrCL() : DROPSErrCL() {}
    DROPSParamErrCL(const std::string& mesg) : DROPSErrCL( mesg) {}
};

/// \brief Parser and container for JSON parameter files
///        This class is based on BOOST Property Tree
/// Usage:
///   - read in JSON file via ParamCL P << ifstream file
///   - use get routines to access data
/// For further information see TRAC
class ParamCL
{
  public:
    typedef boost::property_tree::ptree ptree_type;
    typedef ptree_type::iterator        ptree_iterator_type;
    typedef ptree_type::const_iterator  ptree_const_iterator_type;
    typedef ptree_type::path_type       path_type;
    typedef ParamCL                     self_type;

    ParamCL ();
    ParamCL (std::string path);
    ParamCL (const ptree_type& p) : pt( p) {}

    void read_json (std::string path);

    /// \brief Adds and/or overrides parameters of this->pt as specified in src.
    void update_from (const ptree_type& src);

    const ptree_type& get_child(const std::string& path) const {
      try {
        return pt.get_child( path);
      }
      catch(boost::property_tree::ptree_error& e) {
        throw DROPSParamErrCL( "ParamCL::get_child: Trying to get '" + path + "' failed.\n");
      }
    }

    /// \brief standard get routine
    /// template parameter is the desired type of returned variable
    /// \param pathInPT Path in tree hierarchy
    ///        e.g. for SimParams { "Timestep" : 0.1 };
    ///             you use: myvar = P.get("SimParams.TimeStep");
    /// \returns value stored in node
    /// \exception DROPSParamErrCL if node does not exists
    template <typename OutType>
    OutType get(const std::string& pathInPT) const
    {
      try {
        return this->pt.get<OutType>(pathInPT);
      }
      catch(boost::property_tree::ptree_error& e) {
        throw DROPSParamErrCL( "ParamCL::get: Trying to get '" + pathInPT + "' failed.\n");
      }
    }

    /// \brief get routine with default value
    /// \param pathInPT Path in tree hierarchy
    /// \param default_value value returned if node does not exist
    ///        type of returned variable is determined by the type of default_value
    ///        if node does not exist, node will be created with default_value as value
    /// \param default_value for the node (used if node does not exist yet)
    /// \returns value stored in node
    template <typename OutType>
    OutType get(const std::string & pathInPT, OutType default_val)
    {
      //value in container?
      try {
          return get<OutType>(pathInPT);
      }
      //no? then add for next time
      catch (DROPSParamErrCL& e) {
          this->pt.put(pathInPT, default_val);
          return default_val;
      }
    }

    /// \brief routine to insert a child ptree
    /// \param pathToNode Path in tree hierarchy
    /// \param child ptree to be inserted
    void put_child(const std::string & pathToNode, ptree_type child)
    {
        this->pt.put_child(pathToNode, child);
    }

    /// \brief routine to assign a value to node/create nodes manually
    /// \param pathToNode Path in tree hierarchy
    /// \param value This value will be assigned to node
    template <typename InType>
    void put(const std::string & pathToNode, InType value)
    {
        this->pt.put(pathToNode, value);
    }

    /// \brief routine to assign a value to node/create nodes manually
    /// only if node is not set yet!
    /// \param pathToNode Path in tree hierarchy
    /// \param value This value will be assigned to node
    /// \returns if value has been set already
    template <typename InType>
    bool put_if_unset(const std::string & pathToNode, InType value)
    {
      //value in container?
      try {
          get<InType>(pathToNode);
          return true;
      }
      //no? then add for next time
      catch (DROPSParamErrCL& e) {
          this->pt.put(pathToNode, value);
          return false;
      }
    }

    bool exists(const std::string &pathToNode, const boost::property_tree::ptree *ptTmp = 0 ) const;

    ptree_const_iterator_type begin () const { return pt.begin(); }
    ptree_const_iterator_type end () const { return pt.end(); }

    friend std::istream &operator>>(std::istream& stream, ParamCL& P);
    friend std::ostream &operator<<(std::ostream& stream, ParamCL& P);


  private:
    //container for data
    boost::property_tree::ptree pt;

    std::ostream& print(std::ostream& s);
    void print(boost::property_tree::ptree child, std::string level, std::ostream& s);
};

/// \brief specialisation of standard get routine for Point3DCL
template<>
DROPS::Point3DCL ParamCL::get<DROPS::Point3DCL>(const std::string & pathInPT) const;

/// \brief specialisation of standard get routine for Point2DCL
template<>
DROPS::Point2DCL ParamCL::get<DROPS::Point2DCL>(const std::string & pathInPT) const;

/// \brief specialisation of standard get routine for std::vector<std::string>
template<>
std::vector<std::string > ParamCL::get<std::vector<std::string> >(const std::string & pathInPT) const;

/// \brief specialisation of standard get routine for SArrayCL<Uint, 3>
template<>
DROPS::SArrayCL<Uint, 3> ParamCL::get<DROPS::SArrayCL<Uint, 3> >(const std::string & pathInPT) const;


/// \brief Read P from the file argv[1] or from default_file if argv[1] does not exist; else throw DROPSErrCL.
void read_parameter_file_from_cmdline (ParamCL& P, int argc, char** argv, std::string default_file= std::string());

/// \brief Adds and/or overrides parameters of P directly from the command line
void apply_parameter_modifications_from_cmdline (ParamCL& P, int argc, char **argv);

} // end of namespace DROPS

#endif




