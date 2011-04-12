/// \file pardistributeddatat.h
/// \brief This is the interface between DROPS and a library to handle distributed data, such as DDD
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

//How to use: 1. Read the comments before trying to edit this file.
//            2. Comments precede the data they describe.


#include <stdio.h>
#include "geom/simplex.h"
#include "misc/container.h"
#include "parallel/addeddata.h"
//#include "parallel/distributeddatatypes.h"


#ifndef DROPS_PARDISTRIBUTEDDATA_H
#define DROPS_PARDISTRIBUTEDDATA_H


namespace DROPS
{

/// \brief Interface for handling dynamic distributed data
/**    Here are declared and defined all the members of the wrapper.
       You can find the implementation of this class (wrapper) in the parddd.cpp file
       */
class DynamicDataInterfaceCL
{
    // The declaration of the class' members
  public:
    /// \name Declaring all members
    //@{
    /// \brief Declare a dynamic distributed data type
    static TypeT TypeDeclare (char *name);

    /// \brief Starts transfer phase.
    static void XferBegin();

    /// \brief End of transfer phase.
    static void XferEnd();

    /// \brief Exchanges data via an interface
    static void IFExchange(IFT, size_t, ComProcPtrT, ComProcPtrT);

    static void IFAOnewayX(IFT, ATTRT, IF_DIRT, size_t, ComProcXPtrT,ComProcXPtrT);

    /// \brief Transfer one way
    static void IFOneway(IFT,IF_DIRT,size_t, ComProcPtrT,ComProcPtrT);

    /// \brief Transfer one way on a given level (attribute)
    static void IFAOneway(IFT ,ATTRT ,IF_DIRT ,size_t , ComProcPtrT,ComProcPtrT);

    /// \brief Execute a function for all objects in the interface on a given level (attribute)
    static void IFAExecLocal (IFT, ATTRT, ExecProcPtrT);

    /// \brief Indentifying object which are created by different processes
    static void IdentifyObject (HDRT, PROCT, HDRT);

    /// \brief Execute a function for all objects in the interface
    static void IFExecLocal (IFT, ExecProcPtrT);

    /// \brief Identify an object by a number
    static void IdentifyNumber (HDRT, PROCT, int);

    static IFT IFDefine (int, TypeT *, int, PrioT *, int, PrioT *);

    /// \brief Allows to define a textual description for interfaces
    static void IFSetName (IFT, char *);

    /// \brief Transfer-command for copying a local DDD object to another processor.
    static void XferCopyObj (HDRT, PROCT, PrioT);

    /// \brief Remove object's header from DDD management
    static void HdrDestructor (HDRT);

    /// \brief Transfer array of additional data objects with a DDD local object
    static void XferAddData (int, TypeT);

    /// \brief Creates C++-Objects
    static void SetHandlerCONSTRUCTOR (TypeT, HandlrContrctrT);

    /// \brief Pass a function to the library which is been called to delete an object
    static void SetHandlerDELETE (TypeT , HandlrDltT);

    /// \brief Pass a function to the library which is been called to copy an object
    static void SetHandlerXFERCOPY (TypeT, HandlerXfrcpT);

    /// \brief Pass a function to the library which is been called to gather an object
    static void SetHandlerXFERGATHER (TypeT, HandlerXfrGthrT);

    /// \brief Pass a function to the library which is been called to scatter an object
    static void SetHandlerXFERSCATTER (TypeT, HandlerXfrSctrT);

    /// \brief Pass a function to the library which is been called to update an object
    static void SetHandlerUPDATE (TypeT, HandlerUpdtT);

    /// \brief Pass a function to the library which is been called to make an object consistent
    static void SetHandlerOBJMKCONS (TypeT, HandlerObjMkConsT);

    /// \brief Pass a function to the library which is been called to set the priority of an object
    static void SetHandlerSETPRIORITY (TypeT, HandlerSetPrioT);

    /// \brief Pass a function to the library which is been called to transfer and delete an object
    static void SetHandlerXFERDELETE (TypeT, HandlerXferDltT);

    /// \brief Searching object which are created by different processes
    static HDRT SearchHdr (GIDT);

    /// \brief Get a list of processes owning an object
    static int * InfoProcList (HDRT);

    /// \brief Show defined TypeT
    static void TypeDisplay (TypeT);

    /// \brief Display overview of all DDD interfaces.
    static void IFDisplayAll (void);

    /// \brief Check DDD runtime consistency. Returns total number of errors
    static int ConsCheck (void);

    /// \brief Display overview of single DDD interface.
    static void IFDisplay (IFT);

    /* This is an elegant implementation of a wrapper around the variable number of
     * arguments function DDD_TypeDefine. However the macros work only on the GCC compiler !!!*/

    /*#define FORCE_INLINE __attribute__((__always_inline__))

    static FORCE_INLINE void TypeDefine (TypeT t,...)
    {
        printf("%d ",__builtin_va_arg_pack_len ());
        DDD_TypeDefine(t, __builtin_va_arg_pack ());
    }*/

    /// \brief these are the callings of the variable number of arguements function DDD_TypeDefine
    static void TypeDefineAddedVec(TypeT&, DROPS::AddedVecCL*&, ElemT, DROPS::Uint*, long unsigned int, ElemT, DROPS::Point3DCL*, long unsigned int, ElemT, DROPS::AddedVecCL*);

    static void TypeDefineAddedScal(TypeT&, DROPS::AddedScalCL*&, ElemT, DROPS::Uint*, long unsigned int, ElemT, double*, long unsigned int, ElemT, DROPS::AddedScalCL*);

    static void TypeDefineChildPtrT(TypeT&, DROPS::TetraCL**&, ElemT, DROPS::TetraCL**, long unsigned int, TypeT, ElemT, DROPS::TetraCL**);

    static void TypeDefineBndPtT(TypeT&, DROPS::BndPointCL*&, ElemT, DROPS::BndIdxT*, long unsigned int, ElemT, DROPS::Point2DCL*, long unsigned int, ElemT, DROPS::BndPointCL*);

    static void TypeDefineTetra(TypeT& xt, DROPS::TetraCL*& a, ElemT b, HEADERT* c, ElemT d, DROPS::Usint* e, long unsigned int f, 
                                ElemT g, DROPS::Usint* h, long unsigned int i, ElemT j, DROPS::SArrayCL<DROPS::VertexCL*, 4u>* l, long unsigned int m, TypeT n, 
                                ElemT o, DROPS::SArrayCL<DROPS::EdgeCL*, 6u>* q, long unsigned int r, TypeT s, ElemT t, DROPS::SArrayCL<DROPS::FaceCL*, 4u>* v,
                                long unsigned int w, TypeT x, ElemT y, DROPS::TetraCL** z, long unsigned int a1, TypeT b1, ElemT c1, DROPS::UnknownHandleCL* d1, 
                                long unsigned int e1, ElemT f1, DROPS::TetraCL* g1);

    static void TypeDefineFace(TypeT&, DROPS::FaceCL*&, ElemT, HEADERT*, ElemT, const DROPS::BndIdxT*, 
                               long unsigned int, ElemT, bool*, long unsigned int, ElemT, DROPS::FaceCL*);

    static void TypeDefineEdge(TypeT&, DROPS::EdgeCL*&, ElemT, HEADERT*, ElemT, DROPS::SArrayCL<DROPS::VertexCL*, 2u>*, long unsigned int, TypeT, ElemT, DROPS::VertexCL**, 
                               long unsigned int, TypeT, ElemT, DROPS::SArrayCL<short unsigned int, 2u>*, long unsigned int, ElemT, short int*, long unsigned int, 
                               ElemT, short int*, long unsigned int, ElemT, bool*, long unsigned int, ElemT, DROPS::EdgeCL*);

    static void TypeDefineVertex(TypeT&, DROPS::VertexCL*&, ElemT, DROPS::IdCL<DROPS::VertexCL>*, long unsigned int, ElemT, DROPS::Point3DCL*, long unsigned int, ElemT, 
                                 std::vector<DROPS::BndPointCL, std::allocator<DROPS::BndPointCL> >**, long unsigned int, ElemT, bool*, long unsigned int, ElemT, 
                                 HEADERT*, ElemT, DROPS::UnknownHandleCL*, long unsigned int, ElemT, DROPS::VertexCL*);

    /// \brief Get number of used processes
    static PROCT InfoProcs (void);

    /// \brief Initialization of the library.
    static void Init (int *, char ***);

    /// \brief Get rank of the process
    static PROCT InfoMe (void);

    /// \brief Clean-up of the DDD library.
    static void Exit (void);

    /// \brief Begin identification phase.
    static void IdentifyBegin (void);

    /// \brief End of identification phase.
    static RETT IdentifyEnd (void);
    //@}

};  // end of the class (the wrapper)

}   // end of namespace
#endif
