/// \file gauss.h
/// \brief Gauss solver with pivoting
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Helmut Jarausch, Volker Reichelt; SC RWTH Aachen:

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
 * Copyright 2012 LNM/SC RWTH Aachen, Germany
*/

#ifndef GAUSS_H_
#define GAUSS_H_

namespace DROPS{

template<typename Mat, typename Vec>
void
gauss_pivot(Mat& A, Vec& b)
{
    const size_t _Dim = b.size();
    double max;
    size_t ind_max;
    size_t* p = new size_t[_Dim];
    Vec b2(b);
    for (size_t i=0; i<_Dim; ++i) p[i]= i;

    for (size_t i=0; i<_Dim-1; ++i)
    {
        max= std::fabs(A(p[i], i));
        ind_max= i;
        for (size_t l=i+1; l<_Dim; ++l)
            if (std::fabs(A(p[l], i))>max)
            {
                max= std::fabs(A(p[l], i));
                ind_max= l;
            }
        if (max == 0.0) throw DROPSErrCL("gauss_pivot: Matrix is singular.");
        if (i!=ind_max) std::swap(p[i], p[ind_max]);
        const double piv= A(p[i], i);
        for (size_t j=i+1; j<_Dim; ++j)
        {
            const double fac= A(p[j], i)/piv;
            b2[p[j]]-= fac*b2[p[i]];
            for (size_t k=i+1; k<_Dim; ++k)
                A(p[j], k)-= fac* A(p[i], k);
        }
    }

    for (int i=_Dim-1; i>=0; --i)
    {
        for (int j=_Dim-1; j>i; --j)
            b2[p[i]]-= A(p[i], j)*b2[p[j]];
        b2[p[i]]/= A(p[i], i);
        b[i]= b2[p[i]];
    }
    delete[] p;
}

} // end of namespace DROPS


#endif /* GAUSS_H_ */
