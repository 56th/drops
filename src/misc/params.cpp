/// \file params.cpp
/// \brief read parameters from file.
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt, Thorolf Schulte ; SC RWTH Aachen: Oliver Fortmeier

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

#include "misc/params.h"

namespace DROPS
{
  // =====================================================
  //                    ParamCL
  // =====================================================

  ParamCL::ParamCL() {}
  ParamCL::ParamCL(std::string path)
  {
    this->read_json( path);
  }

  void ParamCL::read_json (std::string path)
  {
    try {
        boost::property_tree::read_json(path, this->pt);
    } catch (boost::property_tree::ptree_error& e) {
          throw DROPSParamErrCL( "ParamCL::read_json: Error while opening or reading file.\n" + std::string(e.what()) + '\n');
    }
  }

  bool ParamCL::exists(const std::string &pathToNode, const boost::property_tree::ptree *ptTmp) const
  {
      std::string::size_type n = pathToNode.find(".");
      if( !ptTmp ) ptTmp = &pt;

      if( n != std::string::npos )
      {
          std::string parent = pathToNode.substr(0,n);
          std::string child  = pathToNode.substr(n+1);
          try{
              ptTmp = &(ptTmp->get_child( parent ));
              return exists( child, ptTmp );
          }
          catch( DROPSErrCL e)
          {
              return false;
          }
      }

      //typedef boost::property_tree::ptree ptree;
      //ptree::const_assoc_iterator it = pt.find( pathToNode );
      auto it = ptTmp->find( pathToNode );
      if( it == ptTmp->not_found() )
          return false;

      return true;
  }

  std::istream &operator>>(std::istream& stream, ParamCL& P)
  {
    try {
        boost::property_tree::read_json(stream, P.pt);
    } catch (boost::property_tree::ptree_error& e) {
          throw DROPSParamErrCL( "ParamCL::operator>>: Read error.\n" + std::string(e.what()) + '\n');
    }
    return stream;
  }

  std::ostream &operator<<(std::ostream& stream, ParamCL& P)
  {
    boost::property_tree::write_json(stream, P.pt);
    return stream;
  }


template <typename T, Uint Size>
   SArrayCL<T, Size>& get_SArray (const ParamCL& p, const std::string& path, SArrayCL<T, Size>& a)
  {
    using boost::property_tree::ptree;
    Uint i= 0;
    try {
        const ptree& node= p.get_child( path);
        for (ptree::const_iterator it= node.begin(); it != node.end(); ++it)
            a[i++]= it->second.get_value<T>();
    }
    catch (boost::property_tree::ptree_error& e) {
        throw DROPSParamErrCL( std::string( "get_SArray: Trying to get '") + path + std::string( "' failed.\n"));
    }
    Assert( i == Size, DROPSParamErrCL( "get_SArray: Wrong number of components in parameter.\n"), ~0);

    return a;
  }

  template<>
  Point3DCL ParamCL::get<Point3DCL>(const std::string & pathInPT) const
  {
    Point3DCL ret( Uninitialized);
    get_SArray( *this, pathInPT, ret);
    return ret;
  }

  template<>
  Point2DCL ParamCL::get<Point2DCL>(const std::string & pathInPT) const
  {
    Point2DCL ret( Uninitialized);
    get_SArray( *this, pathInPT, ret);
    return ret;
  }

  template<>
  SArrayCL<Uint, 3> ParamCL::get<SArrayCL<Uint, 3> >(const std::string & pathInPT) const
  {
    SArrayCL<Uint, 3> ret;
    get_SArray( *this, pathInPT, ret);
    return ret;
  }


  std::ostream& ParamCL::print(std::ostream& s)
  {
    using boost::property_tree::ptree;

    for (ptree::const_iterator it = this->pt.begin(); it != this->pt.end(); ++it) {
        s << it->first << ": " << it->second.get_value<std::string>() << "\n";
        print(it->second, std::string("\t"), s);
    }
    return s;
  }

    template <>
    std::vector<std::string> ParamCL::get<std::vector<std::string> >(const
    std::string &path) const
    {
      std::vector<std::string> v;
      using boost::property_tree::ptree;
      try {
          const ptree& node= this->get_child( path);
          for (ptree::const_iterator it= node.begin(); it != node.end(); ++it)
            v.push_back(it->second.get_value<std::string>());
      }
      catch (boost::property_tree::ptree_error& e) {
          throw DROPSParamErrCL( std::string( "get: Trying to get '") + path + std::string( "' failed.\n"));
      }
      return v;
    }


  void ParamCL::print(boost::property_tree::ptree child, std::string level, std::ostream& s)
  {

    using boost::property_tree::ptree;

    for (ptree::const_iterator it = child.begin(); it != child.end(); ++it) {
        s << level << it->first << ": " << it->second.get_value<std::string>() << "\n";
        print(it->second, level+"\t", s);
    }
  }

// The function updates the ptree dst with the values in src. The variable path is only used for error messages. Defined below.
void update_parameters (const boost::property_tree::ptree& src,  boost::property_tree::ptree& dst, std::string path);

void ParamCL::update_from (const ptree_type& src)
{
    update_parameters (src, pt, "");
}


void read_parameter_file_from_cmdline (ParamCL& P, int argc, char **argv, std::string default_file)
{
    // first read default parameters from default.json, searching in 1.) local directory or 2.) path from default_file
    bool fail= true;
    std::string default_params= default_file.substr( 0, default_file.rfind( "/", default_file.size() )) + "/default.json";
    for (int i=0; i<2 && fail; ++i) {
        std::string fil= i==0 ? "default.json" : default_params;
        std::ifstream test_if_default_exists( fil.c_str());
        if (test_if_default_exists) {
            test_if_default_exists.close();
            fail= false;
            P.read_json( fil);
        }
    }
    if (fail)
        throw DROPSErrCL("read_parameter_file_from_cmdline: Unable to read default parameters from './default.json' and '"
                + default_params + "'\n");
    
    // then read parameters from JSON file specified on command line (otherwise use default_file)
    boost::property_tree::ptree params;
    if ( (argc == 1) || (argv[1][0] == '-') ) { // no file specified on command line
        if (default_file == std::string())
            throw DROPSErrCL(
                "read_parameter_file_from_cmdline: You must specify a parameter file on the command line.\n"
                "        " + std::string( argv[0]) + " <path_to_parameter_file>\n");
        std::cout << "Using fall-back parameter file '" << default_file << "'." << std::endl;
        boost::property_tree::read_json( default_file, params);
    }
    else {
        std::cout << "Using parameter file '" << argv[1] << "'." << std::endl;
        boost::property_tree::read_json( argv[1], params);
    }
    P.update_from (params);

    // finally, read parameters from command line specified by --add-param
    apply_parameter_modifications_from_cmdline( P, argc, argv);
}


void apply_parameter_modifications_from_cmdline (ParamCL& P, int argc, char **argv)
  /** Allows to apply changes to the property tree using the command line, without having to change the .json file.
      Usage: --add-param '{"path_to_node":value}'
             --add-param '...some data in JSON format...'
      Note: use enclosing single quotation marks to prevent the shell from parsing. */
{
    int param_pos = 0;

    for(int i = 1; i < argc && !param_pos; i++) {

        if (std::string (argv[i]) == "--add-param") {

            param_pos = i;
            break;
        }
    }

    if (param_pos) {

        boost::property_tree::ptree changes;
        std::stringstream input;

        for(int i = param_pos+1; i<argc; i++)
            input << argv[i];

        try {
            boost::property_tree::read_json(input,changes);
        } catch (boost::property_tree::ptree_error& e) {
              throw DROPSParamErrCL( "ParamCL::apply_parameter_modifications_from_cmdline: Error while reading from command line.\n" + std::string(e.what()) + '\n');
        }
        P.update_from (changes);
    }
}



/// A ptree is a dictionary if it does not have data, and it has at least one child. In addition, at least one child must have a nontrivial (!= "") key.
void update_parameters_dictionary (const boost::property_tree::ptree& src,  boost::property_tree::ptree& dst, std::string path)
{
    if (src.data() != boost::property_tree::ptree::data_type())
        throw DROPSErrCL("update_parameters_dictionary: The dictionary in src at `" + path + "' contains unkeyed data.\n");

    for (auto& dstch: dst) { // Update parameters which are already in dst.
        const boost::property_tree::ptree* srcch= nullptr;
        try {
            srcch= &src.get_child (dstch.first);
        } catch (boost::property_tree::ptree_error) {
            continue;
        }
        update_parameters (*srcch, dstch.second, path + "." + dstch.first);
    }
    for (auto& srcch: src) { // Add parameters from src to dst which are not already in dst.
        try {
            dst.get_child (srcch.first);
        } catch (boost::property_tree::ptree_error) {
            dst.put_child (srcch.first, srcch.second);
        }
    }
}

/// A ptree is an array if it does not have data, and it has at least one child. In addition, all keys are == "".
void update_parameters_array (const boost::property_tree::ptree& src,  boost::property_tree::ptree& dst, std::string path)
{
    if (src.data() != boost::property_tree::ptree::data_type())
        throw DROPSErrCL("update_parameters_array: The array at `" + path + "' contains unindexed data.\n");
    if (src.size() < dst.size())
        throw DROPSErrCL("update_parameters_array: Refusing to update destination array with strictly smaller array from source at `" + path + "'.\n");

    auto srcit= src.begin();
    size_t i= 0;
    std::stringstream s;
    for (auto& dstch: dst) { // Note that src.size() >= dst.size() at this point.
        if (srcit->first != "")
            throw DROPSErrCL("update_parameters_array: Expected an array in the source at `" + path + "'.\n");
        std::stringstream s;
        s << '[' << i++ << ']';
        update_parameters (srcit++->second, dstch.second, path + s.str());
    }
    dst.insert (dst.end(), srcit, src.end());
}

void update_parameters (const boost::property_tree::ptree& src,  boost::property_tree::ptree& dst, std::string path)
{
    if (src == boost::property_tree::ptree()) // nothing to do.
        return;

    if (dst == boost::property_tree::ptree()) //dst has no children and no data --> (deep) copy src to dst.
        dst= src;
    else if (dst.empty()) { // dst has no children but nontrivial data --> update scalar value.
        if (!src.empty()) // Verify that src is not a dictionary or array.
            throw DROPSErrCL("update_parameters: Expected a scalar type in src at `" + path + "'.\n");
        dst.data()= src.data();
    }
    else { // dst has children.
        if (dst.data() != boost::property_tree::ptree::data_type())
            throw DROPSErrCL("update_parameters: Expected array or dictionary at `" + path + "'.\n");
        std::string alldstkeys;
        for (auto& dstch: dst)
            alldstkeys+= dstch.first;
        if (alldstkeys == "")
            update_parameters_array (src, dst, path);
        else
            update_parameters_dictionary (src, dst, path);
    }
}

} // end of namespace DROPS


