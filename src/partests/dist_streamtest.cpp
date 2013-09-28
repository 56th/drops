/// \file dist_streamtest.cpp
/// \brief test streaming for DiST
/// \author LNM RWTH Aachen: ; SC RWTH Aachen: Oliver Fortmeier, Daniel Medina Cardona

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

#include "parallel/parallel.h"
#include "DiST/DiST.h"
#include <iostream>
#include <sstream>

using namespace std;


int main( int argc, char **argv)
{
#ifdef _PAR
    DROPS::ProcCL::Instance(&argc, &argv);
#endif
  try {
        int myrank = DROPS::ProcCL::MyRank();  // init with rank of process
        int size   = DROPS::ProcCL::Size();    // init with number of processes
        bool binary= true;
        if (argc>1)
            binary= atoi( argv[1]);
        if(myrank==0) cerr << "**************************************************************************\n";
        if(myrank==0) cerr << "***  This program tests the classes DROPS::DiST::SendStreamCL  ***\n";
        if(myrank==0) cerr << "***  and DROPS::DiST::RecvStreamCL which represent a           ***\n";
        if(myrank==0) cerr << "***  data-sending stream and a data-receiving stream respectively.     ***\n";
        if(myrank==0) cerr << "**************************************************************************\n";

        DROPS::DiST::SendStreamCL tosend(binary); //My sending stream
        DROPS::DiST::RecvStreamCL torecv(binary); //My receiving stream


//******************** I N S E R T I O N  P H A S E **********************
        if(myrank==0) cerr << "\n********************* I N S E R T I O N  P H A S E ***********************\n";

        stringstream debugOutputInsert;
        debugOutputInsert << "\nI'm PROCESS " << myrank << " and I'm sending to PROCESS " << (myrank+1)%size <<":\n";

        int intValue=24453+myrank;
        debugOutputInsert << "+ an integer (" << intValue << ")\n";
        tosend << intValue;

        float floatValue=244.53+myrank;
        debugOutputInsert << "+ a float (" << floatValue << ")\n";
        tosend << floatValue;

        double doubleValue=0.53121+myrank;
        debugOutputInsert << "+ a double (" << doubleValue << ")\n";
        tosend << doubleValue;

        DROPS::Uint UintValue=244+myrank;
        debugOutputInsert << "+ an unsigned integer (" << UintValue << ")\n";
        tosend << UintValue;

        bool boolValue=true;
        debugOutputInsert << "+ a boolean (" << boolValue << ")\n";
        tosend << boolValue;

        DROPS::Ulint UlintValue=2442342+myrank;
        debugOutputInsert << "+ an unsigned long integer (" << UlintValue << ")\n";
        tosend << UlintValue;

        DROPS::Usint UsintValue=242+myrank;
        debugOutputInsert << "+ an unsigned short integer (" << UsintValue << ")\n";
        tosend << UsintValue;

        DROPS::byte byteValue=-24+myrank;
        debugOutputInsert << "+ a signed char (" << byteValue << ")\n";
        tosend << byteValue;

        char sepValue= DROPS::DiST::SendRecvStreamAsciiTerminatorC;
        debugOutputInsert << "+ the ascii item-terminator (" << sepValue << ")\n";
        tosend << sepValue;

        DROPS::Ubyte UbyteValue='x';
        debugOutputInsert << "+ an unsigned char (" << UbyteValue << ")\n";
        tosend << UbyteValue;

        size_t SizetValue=8+myrank;
        debugOutputInsert << "+ a size_t value (" << SizetValue << ")\n";
        tosend << SizetValue;

        if (!binary && myrank == 0) {
            cerr << " Orig is: " << tosend.str()     << " position: " << tosend.tellp()     << " state: " << tosend.rdstate()     << '\n';
            // tosend.setstate( ios_base::failbit);
            DROPS::DiST::SendStreamCL tosendcopy= tosend;
            cerr << " Copy is: " << tosendcopy.str() << " position: " << tosendcopy.tellp() << " state: " << tosendcopy.rdstate() << '\n';
            DROPS::DiST::SendStreamCL tosendassg;
            tosendassg= tosend;
            cerr << " assg is: " << tosendassg.str() << " position: " << tosendassg.tellp() << " state: " << tosendassg.rdstate() << '\n';
        }

        DROPS::ProcCL::Barrier();
        cerr << debugOutputInsert.str();
        DROPS::ProcCL::Barrier();

//********************* S E N D I N G  P H A S E **********************
        if(myrank==0) cerr << "\n************************ S E N D I N G  P H A S E *************************\n";

        DROPS::ProcCL::RequestT req;
        req = tosend.Isend((myrank+1)%size);

//********************* R E C E I V I N G  P H A S E **********************
        if(myrank==0) cerr << "\n********************** R E C E I V I N G  P H A S E ***********************\n";

        torecv.Recv((myrank-1+size)%size);


//********************* E X T R A C T I O N  P H A S E **********************
        if(myrank==0) cerr << "\n********************* E X T R A C T I O N  P H A S E **********************\n";

        stringstream myReceivedData;
        myReceivedData << "\nI'm PROCESS " << myrank << " and I've just received from PROCESS " << (size+myrank-1)%size<< ":\n";

        torecv >> intValue;
        myReceivedData << "+ an integer (" << intValue << ")\n";

        torecv >> floatValue;
        myReceivedData << "+ a float (" << floatValue << ")\n";

        torecv >> doubleValue;
        myReceivedData << "+ a double (" << doubleValue << ")\n";

        torecv >> UintValue;
        myReceivedData << "+ an unsigned integer (" << UintValue << ")\n";

        torecv >> boolValue;
        myReceivedData << "+ a boolean (" << boolValue << ")\n";

        torecv >> UlintValue;
        myReceivedData << "+ an unsigned long integer (" << UlintValue << ")\n";

        torecv >> UsintValue;
        myReceivedData << "+ an unsigned short integer (" << UsintValue << ")\n";

        torecv >> byteValue;
        myReceivedData << "+ a signed char (" << byteValue << ")\n";

        torecv >> sepValue;
        myReceivedData << "+ the ascii item-terminator (" << sepValue << ")\n";

        torecv >> UbyteValue;
        myReceivedData << "+ an unsigned char (" << UbyteValue << ")\n";

        torecv >> SizetValue;
        myReceivedData << "+ a size_t value (" << SizetValue << ")\n";

        DROPS::ProcCL::Wait(req);
        cerr << myReceivedData.str();
        DROPS::ProcCL::Barrier();

//********************* V E R I F Y I N G  P H A S E **********************
        if(myrank==0) cerr << "\n*********************** V E R I F Y I N G  P H A S E ************************\n";
        DROPS::ProcCL::Barrier();

        stringstream myVerifiedData;
        int expectedInt = 24453 + (size+myrank-1)%size;
        float expectedFloat = 244.53+(size+myrank-1)%size;
        double expectedDouble = 0.53121+(size+myrank-1)%size;
        DROPS::Uint expectedUint = 244+(size+myrank-1)%size;
        bool expectedBool = true;
        DROPS::Ulint expectedUlint = 2442342+(size+myrank-1)%size;
        DROPS::Usint expectedUsint = 242+(size+myrank-1)%size;
        DROPS::byte  expectedByte = -24+(size+myrank-1)%size;
        char expectedsepValue= DROPS::DiST::SendRecvStreamAsciiTerminatorC;
        DROPS::Ubyte  expectedUbyte = 'x';
        size_t expectedSizet = 8+(size+myrank-1)%size;

        if (intValue!=expectedInt) {
            myVerifiedData << "\nI'm PROCESS " << myrank << ". The received integer and the expected integer don't match." << endl;
        }else if (floatValue!=expectedFloat) {
            myVerifiedData << "\nI'm PROCESS " << myrank << ". The received float and the expected float don't match." << endl;
        }else if (doubleValue!=expectedDouble) {
            myVerifiedData << "\nI'm PROCESS " << myrank << ". The received double and the expected double don't match.\t" << doubleValue << " != " << expectedDouble << ", diff = " << doubleValue-expectedDouble << endl;
        }else if (UintValue!=expectedUint) {
            myVerifiedData << "\nI'm PROCESS " << myrank << ". The received unsigned integer and the expected unsigned integer don't match." << endl;
        }else if (boolValue!=expectedBool) {
            myVerifiedData << "\nI'm PROCESS " << myrank << ". The received boolean and the expected boolean don't match." << endl;
        }else if (UlintValue!=expectedUlint) {
            myVerifiedData << "\nI'm PROCESS " << myrank << ". The received unsigned long integer and the expected unsigned long integer don't match." << endl;
        }else if (UsintValue!=expectedUsint) {
            myVerifiedData << "\nI'm PROCESS " << myrank << ". The received unsigned short integer and the expected unsigned short integer don't match." << endl;
        }else if (byteValue!=expectedByte) {
            myVerifiedData << "\nI'm PROCESS " << myrank << ". The received signed char and the expected signed char don't match." << endl;
        }else if (sepValue!=expectedsepValue) {
            myVerifiedData << "\nI'm PROCESS " << myrank << ". The received ascii item-terminator and the expected ascii item-terminator don't match." << endl;
        }else if (UbyteValue!=expectedUbyte) {
            myVerifiedData << "\nI'm PROCESS " << myrank << ". The received unsigned char and the expected unsigned char don't match." << endl;
        }else if (SizetValue!=expectedSizet) {
            myVerifiedData << "\nI'm PROCESS " << myrank << ". The received sizeof value and the expected sizeof value don't match." << endl;
        }else {
            myVerifiedData << "\nI'm PROCESS " << myrank << ". All the received and expected data match. Data transmission Succeed!" << endl;
        }

        cerr << myVerifiedData.str();
  }
    catch (DROPS::DROPSErrCL err) {err.handle();}
}
