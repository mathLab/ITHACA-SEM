#include <MultiRegions/DisContField3D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        DisContField3D::DisContField3D()
        {
        }

        DisContField3D::DisContField3D( SpatialDomains::MeshGraph3D &graph3D,
                                        SpatialDomains::BoundaryConditions &bcs,
                                        const int bc_loc,
                                        bool SetUpJustDG)
        {
        }

        DisContField3D::DisContField3D( SpatialDomains::MeshGraph3D &graph3D,
                                        SpatialDomains::BoundaryConditions &bcs,
                                        const std::string variable,
                                        bool SetUpJustDG)
        {
        }

        DisContField3D::DisContField3D( SpatialDomains::MeshGraph3D &graph3D,
                                        const GlobalSysSolnType solnType,
                                        const bool constructMap)
        {
        }

        DisContField3D::DisContField3D(const DisContField3D &In)
        {
        }

        DisContField3D::~DisContField3D()
        {
        }

    } // end of namespace
} // end of namespace
