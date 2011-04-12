/// \file parddd.cpp
/// \brief Implementing the interface DynamicDataInterfaceCL by wrapping DDD functions
/// between the DDD-library and a future, improved library
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

//  !!!After the implementation of the methods is carried out delete the call to the DDD functions and also this comment line.

//How to use: 1. Read the comments before trying to edit this file.
//            2. Comments precede the data they describe.

//Include section ----please do not edit this section because it needs to stay like this
//for the validation of the wrapper

#include "parallel/pardistributeddata.h"


//End of the include section

namespace DROPS
{

#ifdef USE_DDD
TypeT DynamicDataInterfaceCL::TypeDeclare (char *name)
{
    throw DROPSErrCL("Access to DDD function");
  //Call of the DDD function
    return DDD_TypeDeclare(name);
}

/** A call to this method establishes a global transfer operation.
    It should be issued on all processors. After this call an arbitrary
    series of {\bf Xfer}-commands may be issued. The global transfer operation
    is carried out via a \funk{XferEnd} call on each processor.
*/
void DynamicDataInterfaceCL::XferBegin()
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_XferBegin();
}

/** This method starts the object transfer process. After a call to
    this method (on all processors) all {\bf Transfer}-commands since
    the last call to \funk{XferBegin} are executed. This involves
    a set of local communications between the processors.
*/
void DynamicDataInterfaceCL::XferEnd()
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_XferEnd();
}

void DynamicDataInterfaceCL::IFExchange(IFT interf, size_t n, ComProcPtrT a, ComProcPtrT b)
{
    throw DROPSErrCL("Access to DDD function");
  //Call of the DDD function
    DDD_IFExchange(interf, n, a, b);
}

void DynamicDataInterfaceCL::IFAOnewayX(IFT interf, ATTRT attribute, IF_DIRT direction, size_t n, ComProcXPtrT a,ComProcXPtrT b)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_IFAOnewayX(interf, attribute, direction, n, a, b);
}

void DynamicDataInterfaceCL::IFOneway(IFT interf, IF_DIRT direction, size_t n, ComProcPtrT a,ComProcPtrT b)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_IFOneway(interf, direction, n, a, b);
}

void DynamicDataInterfaceCL::IFAOneway(IFT interf,ATTRT attribute,IF_DIRT direction,size_t n, ComProcPtrT a,ComProcPtrT b)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_IFAOneway(interf, attribute, direction, n, a, b);
}

void DynamicDataInterfaceCL::IFAExecLocal(IFT interf, ATTRT attribute, ExecProcPtrT exec)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_IFAExecLocal(interf, attribute, exec);
}

/**
    DDD Object identification via another DDD Object.
    After an initial call to \funk{IdentifyBegin}, this function
    identifies two object copies on separate processors. It has to be
    called on both processors with the same identification object.
    The necessary actions (e.g. message transfer) are executed via the
    final call to \funk{IdentifyEnd}; therefore a whole set of
    \funk{Identify}-operations is accumulated.

    After the identification both objects have the same library global
    object ID, which is build using the minimum of both local object IDs.

    The identification object {\em ident} must be either a distributed
    object known to both processors issueing the \funk{IdentifyObject}-command
    or a local object which is not known to these two processors, but which
    will also be identified during the current {\bf Identify}-process.

    The identification specified here may be detailed even further by
    additional calls to {\bf Identify}-operations with the same
    local object. This will construct an identification tupel from
    all {\bf Identify}-commands for this local object.
	 */
void DynamicDataInterfaceCL::IdentifyObject(HDRT ob1, PROCT process, HDRT ob2)
{
    throw DROPSErrCL("Access to DDD function");
  //Call of the DDD function
    DDD_IdentifyObject(ob1, process, ob2);
}

void DynamicDataInterfaceCL::IFExecLocal (IFT interf, ExecProcPtrT exec)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_IFExecLocal(interf, exec);
}

/**
       DDD Object identification via integer number.
       After an initial call to \funk{IdentifyBegin}, this function
       identifies two object copies on separate processors. It has to be
       called on both processors with the same identification value.
       The necessary actions (e.g. message transfer) are executed via the
       final call to \funk{IdentifyEnd}; therefore a whole set of
       \funk{Identify}-operations is accumulated.

       After the identification both objects have the same DDD global
       object ID, which is build using the minimum of both local object IDs.

       The identification specified here may be detailed even further by
       additional calls to {\bf Identify}-operations with the same
       local object. This will construct an identification tupel from
       all {\bf Identify}-commands for this local object.
       */
void DynamicDataInterfaceCL::IdentifyNumber (HDRT ob, PROCT process, int n)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_IdentifyNumber(ob, process, n);
}

/**
    This function defines a new {interface}. Its argument list contains
    three arrays: the first one specifies a subset of the global DDD object set,
    the second and third ones specify a subset of all DDD priorities.
    After the initial creation of the new interface its ID is returned.
    During all following DDD operations ({\bf Identify}- as well as
    {\bf Transfer}-operations) the interface will be kept consistent, and can
    be used for communication via \funk{IFExchange} and
    \funk{IFOneway} and analogous functions.
    */
IFT DynamicDataInterfaceCL::IFDefine (int n1, TypeT tip[], int n2, PrioT priority1[], int n3, PrioT priority2[])
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    return DDD_IFDefine(n1,tip,n2,priority1,n3,priority2);
}

void DynamicDataInterfaceCL::IFSetName (IFT interf, char * s)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_IFSetName(interf,s);
}

/** After an initial call to \funk{XferBegin}, this function
    creates a copy of one local DDD object on another processor with a certain
    priority. The necessary actions (packing/unpacking of object data, message
    transfer) are executed via the final call to \funk{XferEnd}; therefore
    a whole set of {\bf Transfer}-operations is accumulated.

    Caution: As the original object data is not copied throughout this
    call due to efficiency reasons (transferring a large number of objects
    would result in a huge amount of memory copy operations),
    the object may not be changed or deleted until the actual transfer
    has happened. Otherwise the changes will be sent, too.
    */
void DynamicDataInterfaceCL::XferCopyObj (HDRT ob, PROCT process, PrioT priority)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_XferCopyObj(ob, process, priority);
}

/** This function removes an object from DDD-management
	via destructing its \ddd{header}.
	{\em Note:} The \ddd{object} will be destructed, but its copies
	on remote processors will not be informed by \funk{HdrDestructor}.
	There are two consistent possibilities to delete \ddd{objects} which
	have copies on remote processors.
	*/
void DynamicDataInterfaceCL::HdrDestructor (HDRT ob)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_HdrDestructor(ob);
}

/** This function transfers an array of additional data objects
    corresponding to one DDD object to the same destination processor.
    Therefore the latest call to \funk{XferCopyObj} defines the
    {\em primary} object and the destination. An arbitrary number of
    \funk{XferAddData}-calls may be issued to transfer various
    data object arrays at once. This serves as a mechanism to send
    object references which cannot be handled by the native DDD
    pointer conversion technique
    */
void DynamicDataInterfaceCL::XferAddData (int n, TypeT tip)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_XferAddData(n,tip);
}

void DynamicDataInterfaceCL::SetHandlerCONSTRUCTOR (TypeT t, HandlrContrctrT ha)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_SetHandlerCONSTRUCTOR(t,ha);
}

void DynamicDataInterfaceCL::SetHandlerDELETE (TypeT t, HandlrDltT ha)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_SetHandlerDELETE(t,ha);
}

void DynamicDataInterfaceCL::SetHandlerXFERCOPY (TypeT t, HandlerXfrcpT ha)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_SetHandlerXFERCOPY(t, ha);
}

void DynamicDataInterfaceCL::SetHandlerXFERGATHER (TypeT t, HandlerXfrGthrT ha)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_SetHandlerXFERGATHER(t,ha);
}

void DynamicDataInterfaceCL::SetHandlerXFERSCATTER (TypeT t, HandlerXfrSctrT ha)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_SetHandlerXFERSCATTER(t,ha);
}

void DynamicDataInterfaceCL::SetHandlerUPDATE (TypeT t, HandlerUpdtT ha)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_SetHandlerUPDATE(t, ha);
}

void DynamicDataInterfaceCL::SetHandlerOBJMKCONS (TypeT t, HandlerObjMkConsT ha)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_SetHandlerOBJMKCONS(t, ha);
}

void DynamicDataInterfaceCL::SetHandlerSETPRIORITY (TypeT t, HandlerSetPrioT ha)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_SetHandlerSETPRIORITY(t, ha);
}

HDRT DynamicDataInterfaceCL::SearchHdr (GIDT info)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    return DDD_SearchHdr(info);
}

void DynamicDataInterfaceCL::TypeDisplay (TypeT t)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_TypeDisplay(t);
}

/**
    This function displays an overview table containing all interfaces,
    their definition parameters and the current number of constituing objects
    on the calling processor.

    For each interface and each neighbour processor corresponding
    to that interface a relation line is displayed, which contains
    the overall number of objects inside the interface, the number of
    oneway relations outwards, the number of oneway relations inwards,
    the number of exchange relations and the neighbour processor number.
	 */
void DynamicDataInterfaceCL::IFDisplayAll (void)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_IFDisplayAll();
}

/** This function performs a combined local/global consistency
    check on the object data structures and interfaces managed.
    This may be used for debugging purposes; if errors are detected,
    then some understanding of internal structures will be useful.
	 */
int DynamicDataInterfaceCL::ConsCheck (void)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    return DDD_ConsCheck();
}

/**
    This function displays an overview table for one DDD-interface,
    its definition parameters and the current number of constituting objects
    on the calling processor.

    For each neighbor processor corresponding
    to that interface a relation line is displayed, which contains
    the overall number of objects inside the interface, the number of
    oneway relations outwards, the number of oneway relations inwards,
    the number of exchange relations and the neighbor processor number.
    */
void DynamicDataInterfaceCL::IFDisplay (IFT interf)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_IFDisplay (interf);
}

/**	   These functions are the separate initialization of the DDD_TypeDefine method
       They are used to define types for transferring distributed data. Combining all eight
       functions to a single one is desirable. This is quite ugly ...
		 */
void DynamicDataInterfaceCL::TypeDefineAddedVec(DDD_TYPE& t, DROPS::AddedVecCL*& a, DDD_ELEM_TYPE b, DROPS::Uint* c, long unsigned int d, DDD_ELEM_TYPE e, DROPS::Point3DCL* f, long unsigned int g, DDD_ELEM_TYPE h, DROPS::AddedVecCL*i)
{
    throw DROPSErrCL("Access to DDD function");
    DDD_TypeDefine(t,a,b,c,d,e,f,g,h,i);
}

void DynamicDataInterfaceCL::TypeDefineAddedScal(DDD_TYPE& t, DROPS::AddedScalCL*& a, DDD_ELEM_TYPE b, DROPS::Uint* c, long unsigned int d, DDD_ELEM_TYPE e, double* f, long unsigned int g, DDD_ELEM_TYPE h, DROPS::AddedScalCL* i)
{
    throw DROPSErrCL("Access to DDD function");
    DDD_TypeDefine(t,a,b,c,d,e,f,g,h,i);
}

void DynamicDataInterfaceCL::TypeDefineChildPtrT(DDD_TYPE& t, DROPS::TetraCL**& a, DDD_ELEM_TYPE b, DROPS::TetraCL** c, long unsigned int d, DDD_TYPE e, DDD_ELEM_TYPE f, DROPS::TetraCL** g)
{
    throw DROPSErrCL("Access to DDD function");
    DDD_TypeDefine(t,a,b,c,d,e,f,g);
}

void DynamicDataInterfaceCL::TypeDefineBndPtT(DDD_TYPE& t, DROPS::BndPointCL*& a, DDD_ELEM_TYPE b, DROPS::BndIdxT* c, long unsigned int d, DDD_ELEM_TYPE e, DROPS::Point2DCL* f, long unsigned int g, DDD_ELEM_TYPE h, DROPS::BndPointCL* i)
{
    throw DROPSErrCL("Access to DDD function");
    DDD_TypeDefine(t,a,b,c,d,e,f,g,h,i);
}

void DynamicDataInterfaceCL::TypeDefineTetra(DDD_TYPE& xt, DROPS::TetraCL*& a, DDD_ELEM_TYPE b, DDD_HEADER* c, DDD_ELEM_TYPE d, DROPS::Usint* e, long unsigned int f, DDD_ELEM_TYPE g, DROPS::Usint* h, long unsigned int i, DDD_ELEM_TYPE j, DROPS::SArrayCL<DROPS::VertexCL*, 4u>* l, long unsigned int m, DDD_TYPE n, DDD_ELEM_TYPE o, DROPS::SArrayCL<DROPS::EdgeCL*, 6u>* q, long unsigned int r, DDD_TYPE s, DDD_ELEM_TYPE t, DROPS::SArrayCL<DROPS::FaceCL*, 4u>* v, long unsigned int w, DDD_TYPE x, DDD_ELEM_TYPE y, DROPS::TetraCL** z, long unsigned int a1, DDD_TYPE b1, DDD_ELEM_TYPE c1, DROPS::UnknownHandleCL* d1, long unsigned int e1, DDD_ELEM_TYPE f1, DROPS::TetraCL* g1)
{
    throw DROPSErrCL("Access to DDD function");
    DDD_TypeDefine(xt,a,b,c,d,e,f,g,h,i,j,l,m,n,o,q,r,s,t,v,w,x,y,z,a1,b1,c1,d1,e1,f1,g1);
}

void DynamicDataInterfaceCL::TypeDefineFace(DDD_TYPE& t, DROPS::FaceCL*& a, DDD_ELEM_TYPE b, DDD_HEADER* c, DDD_ELEM_TYPE d, const DROPS::BndIdxT* e, long unsigned int f, DDD_ELEM_TYPE g, bool* h, long unsigned int i, DDD_ELEM_TYPE j, DROPS::FaceCL* k)
{
    throw DROPSErrCL("Access to DDD function");
    DDD_TypeDefine(t,a,b,c,d,e,f,g,h,i,j,k);
}

void DynamicDataInterfaceCL::TypeDefineEdge(DDD_TYPE& xt, DROPS::EdgeCL*& a, DDD_ELEM_TYPE b, DDD_HEADER* c, DDD_ELEM_TYPE d, DROPS::SArrayCL<DROPS::VertexCL*, 2u>* e, long unsigned int f, DDD_TYPE g, DDD_ELEM_TYPE h, DROPS::VertexCL** i, long unsigned int j, DDD_TYPE k, DDD_ELEM_TYPE l, DROPS::SArrayCL<short unsigned int, 2u>* m, long unsigned int n, DDD_ELEM_TYPE o, short int* p, long unsigned int q, DDD_ELEM_TYPE r, short int* s, long unsigned int t, DDD_ELEM_TYPE u, bool* v, long unsigned int w, DDD_ELEM_TYPE x, DROPS::EdgeCL* y)
{
    throw DROPSErrCL("Access to DDD function");
    DDD_TypeDefine(xt,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y);
}

void DynamicDataInterfaceCL::TypeDefineVertex(DDD_TYPE& xt, DROPS::VertexCL*& a, DDD_ELEM_TYPE b, DROPS::IdCL<DROPS::VertexCL>* c, long unsigned int d, DDD_ELEM_TYPE e, DROPS::Point3DCL* f, long unsigned int g, DDD_ELEM_TYPE h, std::vector<DROPS::BndPointCL, std::allocator<DROPS::BndPointCL> >** j, long unsigned int k, DDD_ELEM_TYPE l, bool* m, long unsigned int n, DDD_ELEM_TYPE o, DDD_HEADER* p, DDD_ELEM_TYPE q, DROPS::UnknownHandleCL* r, long unsigned int s, DDD_ELEM_TYPE t, DROPS::VertexCL* u)
{
    throw DROPSErrCL("Access to DDD function");
    DDD_TypeDefine(xt,a,b,c,d,e,f,g,h,j,k,l,m,n,o,p,q,r,s,t,u);
}

void DynamicDataInterfaceCL::SetHandlerXFERDELETE (TypeT t, HandlerXferDltT ha)
{
    throw DROPSErrCL("Access to DDD function");
    //Call of the DDD function
    DDD_SetHandlerXFERDELETE(t,ha);
}

PROCT DynamicDataInterfaceCL::InfoProcs (void)
{
    throw DROPSErrCL("Access to DDD function");
    return DDD_InfoProcs();
}

/**
    This function has to be called before any other function
    of the DDD library is called. It initializes the underlying
    PPIF-library, sets all DDD options to their default values
    and initiates all DDD subsystems.

    As some of the memory handler calls will be initiated during
    the execution of this function, the memory manager has to be
    initialized before calling \funk{Init}.

		@param  a      pointer to argc (the application's parameter count)
		@param  b      pointer to argv (the application's parameter list)
	 */
void DynamicDataInterfaceCL::Init (int * a, char *** b)
{
    throw DROPSErrCL("Access to DDD function");
    /// Call of the DDD_function
    DDD_Init(a,b);
}

PROCT DynamicDataInterfaceCL::InfoMe (void)
{
    throw DROPSErrCL("Access to DDD function");
    return DDD_InfoMe();
}

/**
    This function frees memory previously allocated by DDD and finally
    finishes up the PPIF library. After the call to \funk{Exit}
    further usage of the DDD library is no longer possible during this program
    run.
    The clean-up of the memory manager should happen afterwards and is left
    to the DDD application programmer.
	 */
void DynamicDataInterfaceCL::Exit (void)
{
    throw DROPSErrCL("Access to DDD function");
    // Call of the DDD function
    DDD_Exit();
}

/**
    A call to this function establishes a global identification operation.
    It should be issued on all processors. After this call an arbitrary
    series of {\bf Identify}-commands may be issued. The global
    identification operation is carried out via a \funk{IdentifyEnd}
    call on each processor.
	 */
void DynamicDataInterfaceCL::IdentifyBegin (void)
{
    throw DROPSErrCL("Access to DDD function");
    // Call of the DDD function
    DDD_IdentifyBegin();
}

/**
    This function starts the object identification process. After a call to
    this function (on all processors) all {\bf Identify}-commands since the
    last call to \funk{IdentifyBegin} are executed. This involves a set
    of local communications between the processors.
	 */
RETT DynamicDataInterfaceCL::IdentifyEnd (void)
{
    throw DROPSErrCL("Access to DDD function");
    return DDD_IdentifyEnd();
}

int DynamicDataInterfaceExtraCL::InfoIsLocal (HDRT a)
{
    throw DROPSErrCL("Access to DDD function");
    // Call of the DDD function
    return DDD_InfoIsLocal(a);
}

int * DynamicDataInterfaceExtraCL::InfoProcList (DDD_HDR a)
{
    throw DROPSErrCL("Access to DDD function");
    // Call of the DDD function
    return DDD_InfoProcList(a);
}

/**
    This function registers a \ddd{object} via constructing
    its DDD-header. Each \ddd{object} is given a unique {\em global ID},
    which is stored in the DDD-header together with the object's
    properties (type\_id, prio, attr) and additional data used by DDD.
	 */
void DynamicDataInterfaceExtraCL::HdrConstructor (DDD_HDR a, DDD_TYPE b, DDD_PRIO c, DDD_ATTR d)
{
    throw DROPSErrCL("Access to DDD function");
    // Call of the DDD function
    DDD_HdrConstructor(a,b,c,d);
}

void DynamicDataInterfaceExtraCL::HdrConstructorMove (DDD_HDR a, DDD_HDR b)
{
    throw DROPSErrCL("Access to DDD function");
    // Call of the DDD function
    return DDD_HdrConstructorMove(a,b);
}


/**
    This function is regarded as a {\bf Transfer}-operation due
    to its influence on DDD management information on neighboring
    processors. Therefore the function has to be issued between
    a starting \funk{XferBegin} and a final \funk{XferEnd} call.
	 */
void DynamicDataInterfaceExtraCL::XferDeleteObj (DDD_HDR a)
{
    throw DROPSErrCL("Access to DDD function");
    // Call of the DDD function
    return DDD_XferDeleteObj(a);
}

#endif //USE_DDD


#ifdef OWN_DDD
TypeT DynamicDataInterfaceCL::TypeDeclare (char *)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::XferBegin()
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::XferEnd()
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::IFExchange(IFT , size_t , ComProcPtrT, ComProcPtrT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::IFAOnewayX(IFT, ATTRT, IF_DIRT, size_t, ComProcXPtrT, ComProcXPtrT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::IFOneway(IFT, IF_DIRT, size_t, ComProcPtrT, ComProcPtrT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::IFAOneway(IFT, ATTRT, IF_DIRT, size_t, ComProcPtrT, ComProcPtrT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::IFAExecLocal(IFT, ATTRT, ExecProcPtrT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::IdentifyObject(HDRT, PROCT, HDRT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::IFExecLocal (IFT, ExecProcPtrT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::IdentifyNumber (HDRT, PROCT, int)
{
    throw DROPSErrCL("Access to DDD function");
}

IFT DynamicDataInterfaceCL::IFDefine (int, TypeT*, int, PrioT*, int n3, PrioT*)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::IFSetName (IFT, char *)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::XferCopyObj (HDRT, PROCT, PrioT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::HdrDestructor (HDRT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::XferAddData (int, TypeT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::SetHandlerCONSTRUCTOR (TypeT, HandlrContrctrT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::SetHandlerDELETE (TypeT, HandlrDltT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::SetHandlerXFERCOPY (TypeT, HandlerXfrcpT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::SetHandlerXFERGATHER (TypeT, HandlerXfrGthrT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::SetHandlerXFERSCATTER (TypeT t, HandlerXfrSctrT ha)
{
    //Call of the DDD function
    DDD_SetHandlerXFERSCATTER(t,ha);
}

void DynamicDataInterfaceCL::SetHandlerUPDATE (TypeT, HandlerUpdtT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::SetHandlerOBJMKCONS (TypeT, HandlerObjMkConsT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::SetHandlerSETPRIORITY (TypeT, HandlerSetPrioT)
{
    throw DROPSErrCL("Access to DDD function");
}

HDRT DynamicDataInterfaceCL::SearchHdr (GIDT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::TypeDisplay (TypeT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::IFDisplayAll (void)
{
    throw DROPSErrCL("Access to DDD function");
}

int DynamicDataInterfaceCL::ConsCheck (void)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::IFDisplay (IFT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::TypeDefineAddedVec(TypeT&, DROPS::AddedVecCL*&, ElemT, DROPS::Uint*, long unsigned int, ElemT, DROPS::Point3DCL*, long unsigned int, ElemT, DROPS::AddedVecCL*)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::TypeDefineAddedScal(TypeT&, DROPS::AddedScalCL*&, ElemT, DROPS::Uint*, long unsigned int, ElemT, double*, long unsigned int, ElemT, DROPS::AddedScalCL*)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::TypeDefineChildPtrT(TypeT&, DROPS::TetraCL**&, ElemT, DROPS::TetraCL**, long unsigned int, TypeT, ElemT, DROPS::TetraCL**)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::TypeDefineBndPtT(TypeT&, DROPS::BndPointCL*&, ElemT, DROPS::BndIdxT*, long unsigned int, ElemT, DROPS::Point2DCL*, long unsigned int, ElemT, DROPS::BndPointCL*)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::TypeDefineTetra(TypeT&, DROPS::TetraCL*&, ElemT, HEADERT*, ElemT, DROPS::Usint*, long unsigned int, ElemT, DROPS::Usint*, long unsigned int, ElemT, 
                                DROPS::SArrayCL<DROPS::VertexCL*, 4u>*, long unsigned int, TypeT, HEADERT, DROPS::SArrayCL<DROPS::EdgeCL*, 6u>*, 
                                long unsigned int, TypeT, HEADERT, DROPS::SArrayCL<DROPS::FaceCL*, 4u>*, long unsigned int, TypeT, ElemT, DROPS::TetraCL**, long unsigned int, 
                                TypeT, ElemT, DROPS::UnknownHandleCL*, long unsigned int, ElemT, DROPS::TetraCL*)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::TypeDefineFace(TypeT&, DROPS::FaceCL*&, ElemT, HEADERT*, ElemT, const DROPS::BndIdxT*, long unsigned int, ElemT, bool*, long unsigned int, ElemT, DROPS::FaceCL*)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::TypeDefineEdge(TypeT&, DROPS::EdgeCL*&, ElemT, HEADERT*, ElemT, DROPS::SArrayCL<DROPS::VertexCL*, 2u>*, long unsigned int, TypeT, ElemT, DROPS::VertexCL**, 
                               long unsigned int, TypeT, ElemT, DROPS::SArrayCL<short unsigned int, 2u>*, long unsigned int, ElemT, short int*, long unsigned int, 
                               ElemT, short int*, long unsigned int, ElemT, bool*, long unsigned int, ElemT, DROPS::EdgeCL*)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::TypeDefineVertex( TypeT&, DROPS::VertexCL*&, ElemT, DROPS::IdCL<DROPS::VertexCL>*, long unsigned int, ElemT, DROPS::Point3DCL*, long unsigned int, ElemT, 
                                 std::vector<DROPS::BndPointCL, std::allocator<DROPS::BndPointCL> >**, long unsigned int, ElemT, bool*, long unsigned int, ElemT, 
                                 HEADERT*, ElemT, DROPS::UnknownHandleCL*, long unsigned int, ElemT, DROPS::VertexCL*)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::SetHandlerXFERDELETE (TypeT, HandlerXferDltT)
{
    throw DROPSErrCL("Access to DDD function");
}

PROCT DynamicDataInterfaceCL::InfoProcs (void)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::Init (int*, char***)
{
    throw DROPSErrCL("Access to DDD function");
}

PROCT DynamicDataInterfaceCL::InfoMe (void)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::Exit (void)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceCL::IdentifyBegin (void)
{
    throw DROPSErrCL("Access to DDD function");
}

RETT DynamicDataInterfaceCL::IdentifyEnd (void)
{
    throw DROPSErrCL("Access to DDD function");
}

int DynamicDataInterfaceExtraCL::InfoIsLocal (HDRT)
{
    throw DROPSErrCL("Access to DDD function");
}

int * DynamicDataInterfaceExtraCL::InfoProcList (HDRT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceExtraCL::HdrConstructor (HDRT, TypeT, PrioT, ATTRT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceExtraCL::HdrConstructorMove (HDRT, HDRT)
{
    throw DROPSErrCL("Access to DDD function");
}

void DynamicDataInterfaceExtraCL::XferDeleteObj (HDRT)
{
    throw DROPSErrCL("Access to DDD function");
}

#endif // OWN_DDD


}   // end of namespace
