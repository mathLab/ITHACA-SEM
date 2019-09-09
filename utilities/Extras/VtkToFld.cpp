#include <LibUtilities/BasicUtils/VtkUtil.hpp>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCellDataToPointData.h>
#include <vtkContourFilter.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/HashUtils.hpp>
#include <SpatialDomains/MeshGraph.h>
#include <MultiRegions/ExpList2D.h>
#include <LocalRegions/Expansion2D.h>

using namespace std;
using namespace Nektar;

/**
 * @brief Represents a vertex in the mesh.
 *
 * Each vertex has a 3-component coordinate and a scalar value associated with
 * it. The factor provides a mechanism to specify precision of the coordinate
 * comparison. Although coordinates are provided as floating-point numbers, they
 * are stored as integers. The integer value is computed as
 *    floor (x * factor)
 * Therefore a factor of 100 would ensure a precision of 0.01 in the comparison.
 */
struct Vertex
{
    Vertex(double pX, double pY, double pZ, double pScalar, double factor)
        : x((int)floor(pX*factor)),
          y((int)floor(pY*factor)),
          z((int)floor(pZ*factor)),
          scalar(pScalar) {}

    int x;
    int y;
    int z;
    double scalar;

    bool operator==(const Vertex& v)
    {
        return (x == v.x && y == v.y && z == v.z);
    }
};
typedef std::shared_ptr<Vertex> VertexSharedPtr;

/// Define comparison operator for the vertex struct
bool operator==(const VertexSharedPtr& v1, const VertexSharedPtr& v2)
{
    return v1->operator==(*v2);
}


/**
 * @brief Hash function for the Vertex struct used for defining sets.
 */
struct VertexHash : std::unary_function<VertexSharedPtr, std::size_t>
{
    std::size_t operator()(VertexSharedPtr const& p) const
    {
        return hash_combine(p->x, p->y, p->z);
    }
};

/// Define a set of Vertices with associated hash function
typedef std::unordered_set<VertexSharedPtr, VertexHash> VertexSet;


/**
 * Main function.
 *
 * Usage: VtkToFld session.xml input.vtk output.fld [options]
 */
int main(int argc, char* argv[])
{
    // Set up available options
    po::options_description desc("Available options");
    desc.add_options()
        ("help,h",         "Produce this help message.")
        ("name,n", po::value<string>()->default_value("Intensity"),
                "Name of field in VTK file to use for intensity.")
        ("outname,m", po::value<string>()->default_value("intensity"),
                "Name of field in output FLD file.")
        ("precision,p",  po::value<double>()->default_value(1),
             "Precision of vertex matching.");

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("file",   po::value<vector<string> >(), "Input filename");

    po::options_description cmdline_options;
    cmdline_options.add(desc).add(hidden);

    po::positional_options_description p;
    p.add("file", -1);

    po::variables_map vm;

    // Parse command-line options
    try
    {
        po::store(po::command_line_parser(argc, argv).
                  options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (const std::exception& e)
    {
        cerr << e.what() << endl;
        cerr << desc;
        return 1;
    }

    if ( vm.count("help") || vm.count("file") == 0 ||
                             vm["file"].as<vector<string> >().size() != 3) {
        cerr << "Usage: VtkToFld session.xml intensity.vtk output.fld [options]"
             << endl;
        cerr << desc;
        return 1;
    }

    // Extract command-line argument values
    std::vector<std::string> vFiles = vm["file"].as<vector<string> >();
    const string infile  = vFiles[1];
    const string outfile = vFiles[2];
    const double factor  = vm["precision"].as<double>();
    const string name    = vm["name"].as<string>();
    const string outname = vm["outname"].as<string>();

    std::vector<std::string> vFilenames;
    LibUtilities::SessionReaderSharedPtr vSession;
    SpatialDomains::MeshGraphSharedPtr graph2D;
    MultiRegions::ExpList2DSharedPtr Exp;

    vFilenames.push_back(vFiles[0]);
    vSession = LibUtilities::SessionReader::CreateInstance(2, argv, vFilenames);

    try
    {
        //----------------------------------------------
        // Read in mesh from input file
        graph2D = SpatialDomains::MeshGraph::Read(vSession);
        //----------------------------------------------

        //----------------------------------------------
        // Define Expansion
        Exp = MemoryManager<MultiRegions::ExpList2D>::
                    AllocateSharedPtr(vSession,graph2D);
        //----------------------------------------------

        //----------------------------------------------
        // Set up coordinates of mesh
        int coordim = Exp->GetCoordim(0);
        int nq      = Exp->GetNpoints();

        Array<OneD, NekDouble> xc0(nq,0.0);
        Array<OneD, NekDouble> xc1(nq,0.0);
        Array<OneD, NekDouble> xc2(nq,0.0);

        switch(coordim)
        {
        case 2:
            Exp->GetCoords(xc0,xc1);
            break;
        case 3:
            Exp->GetCoords(xc0,xc1,xc2);
            break;
        default:
            ASSERTL0(false,"Coordim not valid");
            break;
        }
        //----------------------------------------------

        vtkPolyDataReader *vtkMeshReader = vtkPolyDataReader::New();
        vtkMeshReader->SetFileName(infile.c_str());
        vtkMeshReader->Update();

        vtkPolyData *vtkMesh = vtkMeshReader->GetOutput();
        vtkCellDataToPointData* c2p = vtkCellDataToPointData::New();
#if VTK_MAJOR_VERSION <= 5
        c2p->SetInput(vtkMesh);
#else
        c2p->SetInputData(vtkMesh);
#endif
        c2p->PassCellDataOn();
        c2p->Update();
        vtkPolyData *vtkDataAtPoints = c2p->GetPolyDataOutput();

        vtkPoints *vtkPoints = vtkMesh->GetPoints();
        ASSERTL0(vtkPoints, "ERROR: cannot get points from mesh.");

        vtkCellArray *vtkPolys = vtkMesh->GetPolys();
        ASSERTL0(vtkPolys,  "ERROR: cannot get polygons from mesh.");

        vtkPointData *vtkPData = vtkDataAtPoints->GetPointData();
        ASSERTL0(vtkPolys,  "ERROR: cannot get point data from file.");

        VertexSet points;
        VertexSet::iterator vIter;
        double p[3];
        double val;
        double x, y, z;
        int coeff_idx;
        int i,j,n;

        if (!vtkDataAtPoints->GetPointData()->HasArray(name.c_str())) {
            n = vtkDataAtPoints->GetPointData()->GetNumberOfArrays();
            cerr << "Input file '" << infile
                 << "' does not have a field named '"
                 << name << "'" << endl;
            cerr << "There are " << n << " arrays in this file." << endl;
            for (int i = 0; i < n; ++i)
            {
                cerr << "  "
                     << vtkDataAtPoints->GetPointData()->GetArray(i)->GetName()
                     << endl;
            }
            return 1;
        }

        // Build up an unordered set of vertices from the VTK file. For each
        // vertex a hashed value of the coordinates is generated to within a
        // given tolerance.
        n = vtkPoints->GetNumberOfPoints();
        for (i = 0; i < n; ++i)
        {
            vtkPoints->GetPoint(i,p);
            val = vtkPData->GetScalars(name.c_str())->GetTuple1(i);
            std::shared_ptr<Vertex> v(new Vertex(p[0],p[1],p[2],val,factor));
            points.insert(v);
        }

        // Now process each vertex of each element in the mesh
        SpatialDomains::PointGeomSharedPtr vert;
        for (i = 0; i < Exp->GetNumElmts(); ++i)
        {
            StdRegions::StdExpansionSharedPtr e = Exp->GetExp(i);
            for (j = 0; j < e->GetNverts(); ++j)
            {
                // Get the index of the coefficient corresponding to this vertex
                coeff_idx = Exp->GetCoeff_Offset(i) + e->GetVertexMap(j);

                // Get the coordinates of the vertex
                vert = e->as<LocalRegions::Expansion2D>()->GetGeom2D()
                                                         ->GetVertex(j);
                vert->GetCoords(x,y,z);

                // Look up the vertex in the VertexSet
                std::shared_ptr<Vertex> v(new Vertex(x,y,z,0.0,factor));
                vIter = points.find(v);

                // If not found, maybe the tolerance should be reduced?
                // If found, record the scalar value from the VTK file in the
                // corresponding coefficient.
                if (vIter == points.end())
                {
                    cerr << "Vertex " << i << " not found. Looking for ("
                            << x << ", " << y << ", " << z << ")" << endl;
                }
                else
                {
                    Exp->UpdateCoeffs()[coeff_idx] = (*vIter)->scalar;
                }
            }
        }
        Exp->SetPhysState(false);

        //-----------------------------------------------
        // Write solution to file
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                                                = Exp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        for(i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back(outname);
            Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        }

        LibUtilities::FieldIOSharedPtr vFld =
            LibUtilities::FieldIO::CreateDefault(vSession);
        vFld->Write(outfile, FieldDef, FieldData);
        //-----------------------------------------------
    }
    catch (...) {
        cout << "An error occurred." << endl;
    }
}
