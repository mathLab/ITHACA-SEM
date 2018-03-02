#include <NekPyConfig.hpp>
#include <StdRegions/StdExpansion.h>

using namespace Nektar;
using namespace Nektar::StdRegions;

Array<OneD, NekDouble> StdExpansion_FwdTrans(
    StdExpansionSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetNcoeffs());
    exp->FwdTrans(in, out);
    return out;
}

Array<OneD, NekDouble> StdExpansion_BwdTrans(
    StdExpansionSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetTotPoints());
    exp->BwdTrans(in, out);
    return out;
}

Array<OneD, NekDouble> StdExpansion_IProductWRTBase(
    StdExpansionSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetNcoeffs());
    exp->BwdTrans(in, out);
    return out;
}

NekDouble StdExpansion_PhysEvaluate(
    StdExpansionSharedPtr exp,
    const Array<OneD, const NekDouble> &coords,
    const Array<OneD, const NekDouble> &physvals)
{
    return exp->PhysEvaluate(coords, physvals);
}

NekDouble StdExpansion_L2(
    StdExpansionSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    return exp->L2(in);
}

NekDouble StdExpansion_L2_Error(
    StdExpansionSharedPtr exp,
    const Array<OneD, const NekDouble> &in,
    const Array<OneD, const NekDouble> &err)
{
    return exp->L2(in, err);
}

py::tuple StdExpansion_GetCoords(StdExpansionSharedPtr exp)
{
    int nPhys = exp->GetTotPoints();
    int coordim = exp->GetCoordim();

    vector<Array<OneD, NekDouble> > coords(coordim);
    for (int i = 0; i < coordim; ++i)
    {
        coords[i] = Array<OneD, NekDouble>(nPhys);
    }

    switch (coordim)
    {
        case 1:
            exp->GetCoords(coords[0]);
            return py::make_tuple(coords[0]);
            break;
        case 2:
            exp->GetCoords(coords[0], coords[1]);
            return py::make_tuple(coords[0], coords[1]);
            break;
        case 3:
            exp->GetCoords(coords[0], coords[1], coords[2]);
            return py::make_tuple(coords[0], coords[1], coords[2]);
            break;
    }

    return py::tuple();
}

void export_StdExpansion()
{
    py::class_<StdExpansion,
               std::shared_ptr<StdExpansion>,
               boost::noncopyable>(
                   "StdExpansion", py::no_init)

        .def("GetNcoeffs", &StdExpansion::GetNcoeffs)
        .def("GetTotPoints", &StdExpansion::GetTotPoints)
        .def("GetBasisType", &StdExpansion::GetBasisType)
        .def("GetPointsType", &StdExpansion::GetPointsType)
        .def("GetNverts", &StdExpansion::GetNverts)
        .def("GetNedges", &StdExpansion::GetNedges)
        .def("GetNfaces", &StdExpansion::GetNfaces)
        .def("DetShapeType", &StdExpansion::DetShapeType)
        .def("GetShapeDimension", &StdExpansion::GetShapeDimension)
        .def("Integral", &StdExpansion::Integral)

        .def("GetBasis", &StdExpansion::GetBasis,
             py::return_value_policy<py::copy_const_reference>())

        .def("FwdTrans", &StdExpansion_FwdTrans)
        .def("BwdTrans", &StdExpansion_BwdTrans)
        .def("IProductWRTBase", &StdExpansion_IProductWRTBase)

        .def("PhysEvaluate", &StdExpansion_PhysEvaluate)
        .def("L2", &StdExpansion_L2)
        .def("L2", &StdExpansion_L2_Error)

        .def("GetCoords", &StdExpansion_GetCoords)
        ;
}
