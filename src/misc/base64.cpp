/// \file base64.cpp
/// \brief Base64 encoder and decoder.
/// \author LNM RWTH Aachen: Joerg Grande

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
 * Copyright 2014 LNM/SC RWTH Aachen, Germany
*/

#include "misc/base64.h"

namespace DROPS
{

namespace Base64Encoding
{

void encode (const unsigned char* x, const unsigned char* xend, std::ostream& os, bool wrap_lines)
{
    const unsigned char alphabet[]= "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                    "abcdefghijklmnopqrstuvwxyz"
                                    "0123456789+/",
                        sixbits=    0x3f;
//     if (sizeof(size_t) < 4)
//         std::cerr << "Base64Encoding::encode: sizeof(size_t): " << sizeof(size_t) << ".\n";

    const size_t rem= (xend - x)%3;
    size_t buf= 0;
    for ( ; x < xend - rem; x+= 3) {
        // Put x[0..2] into buf such that x[0] provides the high-order bits and x[2] provides the low order bits.
        buf=  size_t( x[0]) << 16;
        buf|= size_t( x[1]) <<  8;
        buf|= x[2];
        // Output sextets of buf from highest to lowest bits.
        os << alphabet[buf>>18 & sixbits]
           << alphabet[buf>>12 & sixbits]
           << alphabet[buf>> 6 & sixbits]
           << alphabet[buf     & sixbits];
        if (wrap_lines && (xend + 3 - rem - x)/3 % 19 == 0) // Output a newline after at most 19 octets written.
            os << '\n';
    }
    if (rem == 0)
        return;

    if (x > xend) // Guard against 0 iterations in the loop above.
        x-= 3;
    buf= size_t( x[0]) << 16;
    if (rem == 2)
         buf|= size_t( x[1]) << 8;
    os << alphabet[buf>>18 & sixbits]
       << alphabet[buf>>12 & sixbits];
    if (rem == 2)
        os << alphabet[buf>> 6 & sixbits];

    os << terminator + rem - 1; // Output padding (note rem > 0 at this point).
}

const unsigned char decode_alphabet[]= {
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255, 62,255,255,255, 63,
     52, 53, 54, 55, 56, 57, 58, 59, 60, 61,255,255,255,255,255,255,
    255,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
     15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,255,255,255,255,255,
    255, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
     41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255
};

} // end of namespace Base64Encoding

} // end of namespace DROPS
