#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <MultiRegions/ExpList2D.h>

#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCellDataToPointData.h>
#include <vtkContourFilter.h>
// Usage: VtkToFld session.xml input.vtk output.fld

#include <boost/unordered_set.hpp>


struct Vertex
{
    Vertex(double pX, double pY, double pZ, double pU, double pV, double pW, double factor)
        : x((int)floor(pX*factor)), y((int)floor(pY*factor)), z((int)floor(pZ*factor)), u(pU), v(pV), w(pW) {}
    int x;
    int y;
    int z;
    double u;
    double v;
    double w;

    bool operator=(const Vertex& v)
    {
        return (x == v.x && y == v.y && z == v.z);
    }
};
typedef boost::shared_ptr<Vertex> VertexSharedPtr;

struct VertexHash : std::unary_function<VertexSharedPtr, std::size_t>
{
    std::size_t operator()(VertexSharedPtr const& p) const
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, p -> x);
        boost::hash_combine(seed, p -> y);
        boost::hash_combine(seed, p -> z);
        return seed;
    }
};
typedef boost::unordered_set<VertexSharedPtr, VertexHash> VertexSet;

bool operator==(const VertexSharedPtr& v1, const VertexSharedPtr& v2)
{
    return v1->operator=(*v2);
}


int main(int argc, char* argv[])
{
    if (argc != 5)
    {
        cout << "Usage: FibreToNektar [session].xml [fibre].vtk [fibre-out].fld [scale]" << endl;
        exit(-1);
    }

    MultiRegions::ExpList2DSharedPtr Exp;

    std::vector<std::string> vFilenames;
    vFilenames.push_back(std::string(argv[1]));

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(2, argv, vFilenames);

    try
    {
        const double factor = atof(argv[4]);

        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraphSharedPtr graph2D = MemoryManager<SpatialDomains::MeshGraph2D>::AllocateSharedPtr(vSession);
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
        vtkMeshReader->SetFileName(argv[2]);
        vtkMeshReader->Update();
        vtkPolyData *vtkData = vtkMeshReader->GetOutput();

        vtkPoints *vtkPoints = vtkData->GetPoints();

        ASSERTL0(vtkPoints, "ERROR: cannot get points from mesh.");

        VertexSet points;
        VertexSet::iterator vIter;
        double p[3];
        double vel[3];
        double x, y, z;
        int coeff_idx;
        int i,j;
        Array<OneD, Array<OneD, NekDouble> > fibre(3);
        for (i = 0; i < 3; ++i)
        {
            fibre[i] = Array<OneD, NekDouble>(Exp->GetNcoeffs());
        }
        std::vector<std::string> fibrevar(3);
        fibrevar[0] = "fx";
        fibrevar[1] = "fy";
        fibrevar[2] = "fz";

        // Build up an unordered set of vertices from the VTK file. For each
        // vertex a hashed value of the coordinates is generated to within a
        // given tolerance.
        for (i = 0; i < vtkPoints->GetNumberOfPoints(); ++i)
        {
            vtkPoints->GetPoint(i,p);
            vtkData->GetPointData()->GetArray("FibreDirection")->GetTuple(i, vel);
            boost::shared_ptr<Vertex> v(new Vertex(p[0],p[1],p[2],vel[0],vel[1],vel[2],factor));
            points.insert(v);
        }

        // Now process each vertex of each element in the mesh
        for (i = 0; i < Exp->GetNumElmts(); ++i)
        {
            StdRegions::StdExpansionSharedPtr e = Exp->GetExp(i);
            for (j = 0; j < e->GetNverts(); ++j)
            {
                // Get the index of the coefficient corresponding to this vertex
                coeff_idx = Exp->GetCoeff_Offset(i) + e->GetVertexMap(j);

                // Get the coordinates of the vertex
                SpatialDomains::VertexComponentSharedPtr vert = e->GetGeom2D()->GetVertex(j);
                vert->GetCoords(x,y,z);

                // Look up the vertex in the VertexSet
                boost::shared_ptr<Vertex> v(new Vertex(x,y,z,0.0, 0.0, 0.0,factor));
                vIter = points.find(v);

                // If not found, maybe the tolerance should be reduced?
                // If found, record the scalar value from the VTK file in the
                // corresponding coefficient.
                if (vIter == points.end())
                {
                    cerr << "Vertex " << i << " not found."
                            << x << ", " << y << ", " << z << endl;
                    cerr << v->x << ", " << v->y << ", " << v->z << endl;
                }
                else
                {
                    double mag = sqrt(((*vIter)->u)*((*vIter)->u) + ((*vIter)->v)*((*vIter)->v) + ((*vIter)->w)*((*vIter)->w));
                    if (mag < 1e-06) mag = 1.0;
                    fibre[0][coeff_idx] = fabs((*vIter)->u) / mag;
                    fibre[1][coeff_idx] = fabs((*vIter)->v) / mag;
                    fibre[2][coeff_idx] = fabs((*vIter)->w) / mag;
                }
            }
        }

        //-----------------------------------------------
        // Write solution to file
        string   out(argv[3]);
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                                                    = Exp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        for (j = 0; j < 3; ++j)
        {
            for(i = 0; i < FieldDef.size(); ++i)
            {
                FieldDef[i]->m_fields.push_back(fibrevar[j]);
                Exp->AppendFieldData(FieldDef[i], FieldData[i], fibre[j]);
            }
        }
        LibUtilities::Write(out, FieldDef, FieldData);
        //-----------------------------------------------
    }
    catch (...) {
        cout << "An error occurred." << endl;
    }
}
