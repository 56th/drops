/// \file distributeddatatypes.h
/// \brief  Redefinintion of the types from the DDD library
/// \author LNM RWTH Aachen: ; SC RWTH Aachen: Oliver Fortmeier, Alin Bastea

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

#include <stdio.h>

#ifndef DROPS_HEADER_H
#define DROPS_HEADER_H

#define OWN_DDD

#ifdef _PAR
# ifdef USE_DDD
#   include <ddd.h>
# endif
#endif


namespace DROPS
{
#ifdef USE_DDD
    typedef DDD_HEADER HEADERT;                  ///< Object header
    typedef DDD_RET RETT;                        ///< Return types for DDD functions
    typedef DDD_IF IFT;                          ///< Interface type
    typedef ComProcPtr ComProcPtrT;              ///< Function type
    typedef ComProcXPtr ComProcXPtrT;            ///< Function type
    typedef DDD_TYPE TypeT;                      ///< Type of a distributed object
    typedef DDD_ELEM_TYPE ElemT;                 ///< Type of an element
    typedef DDD_ATTR ATTRT;                      ///< attribute of a distributed object
    typedef DDD_IF_DIR IF_DIRT;                  ///< direction of an interface
    typedef ExecProcPtr ExecProcPtrT;            ///< Function type
    typedef DDD_HDR HDRT;                        ///< each object stores a header which can be found here
    typedef DDD_GID GIDT;                        ///< global id type
    typedef DDD_PROC PROCT;                      ///< process numbers
    typedef DDD_PRIO PrioT;                      ///< priority
    typedef DDD_OBJ OBJT;                        ///< pointer to an object
    typedef HandlerCONSTRUCTOR HandlrContrctrT;  ///< Handler how to construct an object
    typedef HandlerDELETE HandlrDltT;            ///< Handler how to delete an object
    typedef HandlerXFERCOPY HandlerXfrcpT;       ///< Handler how to transfer and copy an object
    typedef HandlerXFERGATHER HandlerXfrGthrT;   ///< Handler how to transfer and gather an object
    typedef HandlerXFERSCATTER HandlerXfrSctrT;  ///< Handler how to transfer and scatter an object
    typedef HandlerUPDATE HandlerUpdtT;          ///< Handler how to update an object
    typedef HandlerOBJMKCONS HandlerObjMkConsT;  ///< Handler how to make an object consistant
    typedef HandlerSETPRIORITY HandlerSetPrioT;  ///< Handler how to set the priority
    typedef HandlerXFERDELETE HandlerXferDltT;   ///< Hanler how to transfer and delete an object
#define MarkHdrInvalid(DDD_HDR) DDD_MarkHdrInvalid(DDD_HDR); ///< Macro definition of the creator of the header
#endif

#ifdef OWN_DDD
    typedef int HEADERT;                            ///< Object header
    typedef int RETT;                               ///< Return types for DDD functions
    typedef int IFT;                                ///< Interface type
    typedef int ComProcPtrT;                        ///< Function type
    typedef int ComProcXPtrT;                       ///< Function type
    typedef int TypeT;                              ///< Type of a distributed object
    typedef int ElemT;                              ///< Type of an element
    typedef int ATTRT;                              ///< attribute of a distributed object
    typedef int IF_DIRT;                            ///< direction of an interface
    typedef int ExecProcPtrT;                       ///< Function type
    typedef int HDRT;                               ///< each object stores a header which can be found here
    typedef int GIDT;                               ///< global id type
    typedef int PROCT;                              ///< process numbers
    typedef int PrioT;                              ///< priority
    typedef int OBJT;                               ///< pointer to an object
    typedef int HandlrContrctrT;                    ///< Handler how to construct an object
    typedef int HandlrDltT;                         ///< Handler how to delete an object
    typedef int HandlerXfrcpT;                      ///< Handler how to transfer and copy an object
    typedef int HandlerXfrGthrT;                    ///< Handler how to transfer and gather an object
    typedef int HandlerXfrSctrT;                    ///< Handler how to transfer and scatter an object
    typedef int HandlerUpdtT;                       ///< Handler how to update an object
    typedef int HandlerObjMkConsT;                  ///< Handler how to make an object consistant
    typedef int HandlerSetPrioT;                    ///< Handler how to set the priority
    typedef int HandlerXferDltT;                    ///< Hanler how to transfer and delete an object
#define MarkHdrInvalid(DDD_HDR) (void(0));          ///< Macro definition of the creator of the header
#endif


/// \brief This class is the wrapper around the methods declared in geom/simplex.h
///        which because of the TypeDeclare method had to be included inside parallel/pardistributeddata.h
class DynamicDataInterfaceExtraCL
{
  public:
    static int InfoIsLocal (HDRT);

    /// \brief Return list of couplings of certain object
    static int * InfoProcList (HDRT);

    /// \brief Initiate object's \ddd{header}.
    static void HdrConstructor (HDRT, TypeT, PrioT, ATTRT);

    /// \brief Create DDD_HDR copy inside local memory simultaneously destruct original DDD_HDR
    static void HdrConstructorMove (HDRT, HDRT);

    /// \brief Transfer-command for deleting a local DDD object.
    static void XferDeleteObj (HDRT);
};
}
#endif
