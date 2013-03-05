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
    Vertex(double pX, double pY, double pZ, double pScalar, double factor)
        : x((int)floor(pX*factor)), y((int)floor(pY*factor)), z((int)floor(pZ*factor)), scalar(pScalar) {}
    int x;
    int y;
    int z;
    double scalar;

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
        vtkPolyData *vtkMesh = vtkMeshReader->GetOutput();
        vtkCellDataToPointData* c2p = vtkCellDataToPointData::New();
        c2p->SetInput(vtkMesh);
        c2p->PassCellDataOn();
        c2p->Update();
        vtkPolyData *vtkDataAtPoints = c2p->GetPolyDataOutput();

        vtkPoints *vtkPoints = vtkMesh->GetPoints();
        vtkCellArray *vtkPolys = vtkMesh->GetPolys();

        ASSERTL0(vtkPoints, "ERROR: cannot get points from mesh.");
        ASSERTL0(vtkPolys,  "ERROR: cannot get polygons from mesh.");

        VertexSet points;
        VertexSet::iterator vIter;
        double p[3];
        double val;
        double x, y, z;
        int coeff_idx;
        int i,j;

        // Build up an unordered set of vertices from the VTK file. For each
        // vertex a hashed value of the coordinates is generated to within a
        // given tolerance.
        for (i = 0; i < vtkPoints->GetNumberOfPoints(); ++i)
        {
            vtkPoints->GetPoint(i,p);
            val = vtkDataAtPoints->GetPointData()->GetScalars("Image_Intensity")->GetTuple1(i);
            boost::shared_ptr<Vertex> v(new Vertex(p[0],p[1],p[2],val,factor));
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
                boost::shared_ptr<Vertex> v(new Vertex(x,y,z,0.0,factor));
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
                    Exp->UpdateCoeffs()[coeff_idx] = (*vIter)->scalar;
                }
            }
        }
        Exp->SetPhysState(false);

        //-----------------------------------------------
        // Write solution to file
        string   out(argv[3]);
        if (vSession->GetComm()->GetSize() > 1)
        {
            out += "." + boost::lexical_cast<string>(vSession->GetComm()->GetRank());
        }
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                                                    = Exp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        for(i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("intensity");
            Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        LibUtilities::Write(out, FieldDef, FieldData);
        //-----------------------------------------------
    }
    catch (...) {
        cout << "An error occurred." << endl;
    }
}
