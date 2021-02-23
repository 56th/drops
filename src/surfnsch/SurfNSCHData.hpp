#ifndef SURF_NAVIER_STOKES_CAHN_HILLIARD_DATA_HPP
#define SURF_NAVIER_STOKES_CAHN_HILLIARD_DATA_HPP

#include "surfnavierstokes/SurfNavierStokesData.hpp"
#include "surfphasesep/SurfCahnHilliardData.hpp"

namespace DROPS {

    std::pair<SurfNavierStokesData, SurfCahnHilliardData> SurfNSCHDataFactory(std::string const & test, ParamCL const & param) {
        SurfNavierStokesData dataNS;
        SurfCahnHilliardData dataCH;
        // ... add KH, RT, and convergence test
        return data;
    }

}

#endif // SURF_NAVIER_STOKES_CAHN_HILLIARD_DATA_HPP