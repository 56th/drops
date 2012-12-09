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

    void *handle;
    //double (*cosine)(double);
    char *error;
    handle = dlopen ( P.get<std::string>("LoadExternal","").c_str(), RTLD_LAZY);
    // load(P.get<std::string>("LoadExternal",""));
    if (!handle) {
        fprintf (stderr, "%s\n", dlerror());
        exit(1);
    }
    dlerror();    /* Clear any existing error */
    //cosine = dlsym(handle, "cos");
    if ((error = dlerror()) != NULL)  {
        fprintf (stderr, "%s\n", error);
        exit(1);
    }
    //printf ("%f\n", (*cosine)(2.0));