/// \file base64.h
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

#ifndef DROPS_BASE64_H
#define DROPS_BASE64_H

#include <iostream>

namespace DROPS
{

namespace Base64Encoding
{

///\brief Write the base64 encoding of [x, xend) to os.
/// The output is padded with "=" as required by RFC 2045.
void encode (const unsigned char* x, const unsigned char* xend, std::ostream& os, bool wrap_lines= false);

///\brief Decode the base64-input from the character stream is to x.
/// Padding with "=" is not required; whitespace is ignored.
/// The OutputIter x is not checked for overflow (e.g, use a std::back_inserter).
template <typename OutputIter>
OutputIter decode (std::istream& is, OutputIter x);


/// =Private part of the header: definition of function templates.=

///\brief Map the ascii-value of the base64-alphabet to the original value.
extern const unsigned char decode_alphabet[];

///\brief Padding characters.
const unsigned char terminator[]= "==";

template <typename OutputIter>
OutputIter decode (std::istream& is, OutputIter x)
{
//     if (sizeof(size_t) < 4)
//         std::cerr << "Base64Encoding::decode: sizeof(size_t): " << sizeof(size_t) << ".\n";

    const size_t do_output_mask= size_t( 1) << 24,
                 empty_buf=      1;
    // The buffer is built by left shifting and or-ing in new data. So a buffer
    // containing the bytes 0 0 1 will have a leading 1 at bit 24. Hence, there is no ambiguity.
    size_t buf= empty_buf;
    for (unsigned char c= is.get(); !is.eof(); c= is.get()) {
        if (std::isspace(c)) // Ignore white space
            continue;
        else if (c == terminator[0]) // End terminator; all input data has been written to buf.
            break;
        else if (decode_alphabet[c] == 255) { // c is not in the alphabet of the encoder.
            std::cerr << "Skipping invalid input character '" << c << "'.\n";
            continue;
        }
        c= decode_alphabet[c];
        buf= buf<<6 | c;
        if ((buf & do_output_mask) == 0)
            continue;
        *x++= buf >> 16;
        *x++= buf >>  8;
        *x++= buf;
        buf= empty_buf;
    }
    // There could be 2, 1, or 0 chars remaining in the buffer, encoded in 18, 12, or 0 bits.
    if (buf & size_t( 1)<<18) {
            *x++= buf >> 10;
            *x++= buf >>  2;
    }
    else if (buf & size_t( 1)<<12)
        *x++= buf >> 4;
    return x;
}

} // end of namespace Base64NS

} // end of namespace DROPS

#endif