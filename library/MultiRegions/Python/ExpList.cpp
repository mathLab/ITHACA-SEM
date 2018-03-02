#include <NekPyConfig.hpp>
#include <MultiRegions/ExpList.h>

#include <fstream>

using namespace Nektar;
using namespace Nektar::StdRegions;
using namespace Nektar::LocalRegions;
using namespace Nektar::SpatialDomains;
using namespace Nektar::MultiRegions;

int ExpList_GetNcoeffs(ExpListSharedPtr exp)
{
    return exp->GetNcoeffs();
}

ExpansionSharedPtr ExpList_GetExp(ExpListSharedPtr exp, int i)
{
    return exp->GetExp(i);
}

void ExpList_WriteVTK(ExpListSharedPtr exp, std::string filename)
{
    ofstream out(filename.c_str());
    exp->WriteVtkHeader(out);
    for (int i = 0; i < exp->GetExpSize(); ++i)
    {
        exp->WriteVtkPieceHeader(out, i);
        exp->WriteVtkPieceFooter(out, i);
    }
    exp->WriteVtkFooter(out);
}

Array<OneD, NekDouble> ExpList_FwdTrans(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetNcoeffs());
    exp->FwdTrans(in, out);
    return out;
}

Array<OneD, NekDouble> ExpList_BwdTrans(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetNpoints());
    exp->BwdTrans(in, out);
    return out;
}

Array<OneD, NekDouble> ExpList_IProductWRTBase(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetNcoeffs());
    exp->BwdTrans(in, out);
    return out;
}

NekDouble ExpList_L2(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    return exp->L2(in);
}

NekDouble ExpList_L2_Error(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in,
    const Array<OneD, const NekDouble> &err)
{
    return exp->L2(in, err);
}

py::tuple ExpList_GetCoords(ExpListSharedPtr exp)
{
    int nPhys = exp->GetNpoints();
    int coordim = exp->GetCoordim(0);

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


void export_ExpList()
{
    py::class_<ExpList,
               std::shared_ptr<ExpList>,
               boost::noncopyable>(
                   "ExpList", py::no_init)

        .def("GetNpoints", &ExpList::GetNpoints)
        .def("GetNcoeffs", &ExpList_GetNcoeffs)
        .def("GetExp", &ExpList_GetExp)
        .def("GetExpSize", &ExpList::GetExpSize)
        .def("WriteVTK", &ExpList_WriteVTK)
        .def("GetCoords", &ExpList_GetCoords)
        .def("FwdTrans", &ExpList_FwdTrans)
        .def("BwdTrans", &ExpList_BwdTrans)
        .def("IProductWRTBase", &ExpList_IProductWRTBase)
        .def("L2", &ExpList_L2)
        .def("L2", &ExpList_L2_Error)

        ;
}
