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

  bool ParamCL::exists(const std::string &pathToNode) const
  {
      typedef boost::property_tree::ptree ptree;

      ptree::const_assoc_iterator it = pt.find( pathToNode );
      if( it == pt.not_found() )
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

  void read_parameter_file_from_cmdline (ParamCL& P, int argc, char **argv, std::string default_file)
{
    if (argc == 1) {
        if (default_file == std::string())
            throw DROPSErrCL(
                "read_parameter_file_from_cmdline: You must specify a parameter file on the command line.\n"
                "        " + std::string( argv[0]) + " <path_to_parameter_file>\n");
        std::cout << "Using default parameter file '" << default_file << "'." << std::endl;
        P.read_json( default_file);
    }
    else {
        std::cout << "Using  parameter file '" << argv[1] << "'." << std::endl;
        P.read_json( argv[1]);
        apply_parameter_modifications_from_cmdline(P,argc,argv);
    }
}

void update_parameters (const boost::property_tree::ptree& pt, std::string key, ParamCL& P)
{
    using boost::property_tree::ptree;
    std::string nkey;

    if (!key.empty()){
        nkey = key + ".";
    }

    for (ptree::const_iterator it = pt.begin(); it != pt.end(); ++it) {
        update_parameters(it->second, nkey + it->first, P);
        P.put(nkey + it->first, it->second.data());
    }
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
        update_parameters(changes,"",P);

    }

}

  // =====================================================
  //                    ReadParamsCL
  // =====================================================

  void ReadParamsCL::SetInfo( std::string& bez, char typ)
  {
      bez= group_ + bez;
      if (!IsKnown(bez))
      {
          info_[bez].first= typ;
          info_[bez].second= false;
      }
      else
          throw DROPSErrCL( "ReadParamsCL: Parameter "+bez+" already registered");
  }

  void ReadParamsCL::RegInt( int& ref, std::string bez)
  {
      SetInfo( bez, 'i');
      iparam_[bez]= &ref;
  }

  void ReadParamsCL::RegDouble( double& ref, std::string bez)
  {
      SetInfo( bez, 'd');
      dparam_[bez]= &ref;
  }

  void ReadParamsCL::RegCoord( Point3DCL& ref, std::string bez)
  {
      SetInfo( bez, 'c');
      cparam_[bez]= &ref;
  }

  void ReadParamsCL::RegString( std::string& ref, std::string bez)
  {
      SetInfo( bez, 's');
      sparam_[bez]= &ref;
  }

  void ReadParamsCL::RegInt( int& ref, std::string bez, int defaultvalue)
  {
      SetInfo( bez, 'i');
      iparam_[bez]= &ref;
      iparam_def_[bez]= defaultvalue;
  }

  void ReadParamsCL::RegDouble( double& ref, std::string bez, double defaultvalue)
  {
      SetInfo( bez, 'd');
      dparam_[bez]= &ref;
      dparam_def_[bez]= defaultvalue;
  }

  void ReadParamsCL::RegCoord( Point3DCL& ref, std::string bez, Point3DCL defaultvalue)
  {
      SetInfo( bez, 'c');
      cparam_[bez]= &ref;
      cparam_def_[bez]= defaultvalue;
  }

  void ReadParamsCL::RegString( std::string& ref, std::string bez, std::string defaultvalue)
  {
      SetInfo( bez, 's');
      sparam_[bez]= &ref;
      sparam_def_[bez]= defaultvalue;
  }

  void ReadParamsCL::BeginGroup( const std::string& s)
  {
      group_+= s;
      group_+= ":";
  }

  void ReadParamsCL::EndGroup()
  {
      if (group_.empty())
          throw DROPSErrCL("ReadParamsCL::EndGroup: missing BeginGroup!");
      int pos= group_.rfind( ':', group_.size()-2);
      group_.erase( pos+1);
  }


  inline void SkipComment( std::istream& s)
  {
    if (!s.eof()) s >> std::ws;
    while (!s.eof() && s.peek()=='#')
    {
      while (!s.eof() && s.get()!='\n') {}
      if (!s.eof()) s >> std::ws;
    }
  }

  void ReadParamsCL::ReadEntry( std::istream& is)
  {
      SkipComment( is);
      while (!is.eof())
      {
          switch(is.peek())
          {
            case '#': SkipComment( is); break;
            case '=': is.get(); SkipComment( is);
                      name_= group_ + name_;
                      ReadData( is);
                      name_.clear();
                      return;
            case '{': is.get(); if (!name_.empty())
                                 { BeginGroup( name_); name_.clear(); }
                                 else throw DROPSErrCL("ReadParamsCL: group name missing before '{'");
                                 break;
            case '}': is.get(); if (name_.empty()) EndGroup();
                                 else throw DROPSErrCL("ReadParamsCL: end of group "+group_+" before defining parameter "+name_);
                                 break;
            default: name_+= is.get(); break;
          }
          if (!is.eof()) is >> std::ws;
      }
  }

  void ReadParamsCL::ReadData( std::istream& is)
  {
      if (!IsKnown( name_))
      {
          std::cout << "Skipping unknown parameter " << name_ << std::endl;
          is.ignore( 256, '\n');
          return;
      }
      int i= 0; double d= 0; Point3DCL p(0.); std::string s;
      char typ= info_[name_].first;
      switch (typ)
      {
          case 'i': is >> i; break;
          case 'd': is >> d; break;
          case 'c': is >> p[0] >> p[1] >> p[2]; break;
          case 's': is >> s;
      }
      if (!is)
      {
          throw DROPSErrCL( "ReadParamsCL: reading of data failed for parameter "+name_);
      }
      switch (typ)
      {
          case 'i': *iparam_[name_]= i; break;
          case 'd': *dparam_[name_]= d; break;
          case 'c': *cparam_[name_]= p; break;
          case 's': *sparam_[name_]= s; break;
      }
      info_[name_].second= true;
  }

  void ReadParamsCL::Clear()
  {
      info_.clear(); iparam_.clear(); dparam_.clear(); cparam_.clear(); sparam_.clear();
  }

  void ReadParamsCL::ReadParams( std::istream& is)
  {
      if (!group_.empty())
          throw DROPSErrCL("ReadParamCL: group "+group_+" not terminated properly");

      if (!is)
          throw DROPSErrCL("ReadParamsCL: file error");

      name_.clear(); group_.clear();
      // reset info: all parameters uninitialized
      for (InfoT::iterator it= info_.begin(), end= info_.end(); it!=end; ++it)
          it->second.second= false;

      while (!is.eof())
      {
          ReadEntry( is);
          if (!name_.empty())
              throw DROPSErrCL("ReadParamsCL: error while reading parameter "+name_);
      }
      UseDefaults();
      PrintWarning();
      if (!group_.empty())
          throw DROPSErrCL("ReadParamCL: group "+group_+" not terminated properly");
  }

  /// \note This routine can be used to generate a standard parameter file
  ///
  void ReadParamsCL::WriteParams( std::ostream& os) const
  {
      os << "#=============================================\n"
         << "#    DROPS parameter file\n"
         << "#=============================================\n\n";
      for (InfoT::const_iterator it= info_.begin(), end= info_.end();
          it!=end; ++it)
      {
          std::string bez= it->first;
          char typ= it->second.first;
          os << bez << "\t=\t";
          switch(typ)
          {
            case 'i': os << *(iparam_.find(bez)->second); break;
            case 'd': os << *(dparam_.find(bez)->second); break;
            case 'c': os << *(cparam_.find(bez)->second); break;
            case 's': os << *(sparam_.find(bez)->second); break;
          }
          os << '\n';
      }
  }

  void ReadParamsCL::UseDefaults()
  {
      for (InfoT::const_iterator it= info_.begin(), end= info_.end(); it!=end; ++it)
          if(!it->second.second) // parameter uninitialized
          {
              switch(it->second.first){
                  case 'i':
                  {
                      std::map<std::string,int>::iterator f = iparam_def_.find(it->first);
                      if ( f != iparam_def_.end()){
                       *iparam_[f->first]= f->second;
                       info_[f->first].second= true;
                       std::cout << "WARNING: Parameter " << f->first << " initialized by default value " << f->second << "!\n";
                      }
                      break;
                  }
                  case 'd':
                  {
                      std::map<std::string,double>::iterator f = dparam_def_.find(it->first);
                      if ( f != dparam_def_.end()){
                       *dparam_[f->first]= f->second;
                       info_[f->first].second= true;
                       std::cout << "WARNING: Parameter " << f->first << " initialized by default value " << f->second << "!\n";
                      }
                      break;
                  }
                  case 'c':
                  {
                      std::map<std::string,Point3DCL>::iterator f = cparam_def_.find(it->first);
                      if ( f != cparam_def_.end()){
                       *cparam_[f->first]= f->second;
                       info_[f->first].second= true;
                       std::cout << "WARNING: Parameter " << f->first << " initialized by default value " << f->second << "!\n";
                      }
                      break;
                  }
                  case 's':
                  {
                      std::map<std::string,std::string>::iterator f = sparam_def_.find(it->first);
                      if ( f != sparam_def_.end()){
                       *sparam_[f->first]= f->second;
                       info_[f->first].second= true;
                       std::cout << "WARNING: Parameter " << f->first << " initialized by default value " << f->second << "!\n";
                      }
                      break;
                  }
                  default:
                      break;

              }
          }
  }

  void ReadParamsCL::PrintWarning() const
  {
      bool AllOk= true;
      for (InfoT::const_iterator it= info_.begin(), end= info_.end(); it!=end; ++it)
          if(!it->second.second) // parameter uninitialized
          {
              std::cout << "WARNING: Parameter " << it->first << " uninitialized!\n";
              AllOk= false;
          }
      if (!AllOk) throw DROPSErrCL("ReadParamsCL: Parameters above are missing in parameter file!");
  }

  } // end of namespace DROPS


