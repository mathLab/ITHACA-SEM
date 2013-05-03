#include <string>
using namespace std;

#include "MeshElements.h"
#include "ProcessCyl.h"

#include <LocalRegions/SegExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/NodalTriExp.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessCyl::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "cyl"),
                ProcessCyl::create);
        /**
         * @brief Default constructor.
         */
        ProcessCyl::ProcessCyl(MeshSharedPtr m) : ProcessModule(m)
        {
            config["surf"] = ConfigOption(false, "-1",
                "Tag identifying surface to process.");
            config["r"] = ConfigOption(false, "0.0",
                "Radius of cylinder.");
            config["N"] = ConfigOption(false, "7",
                "Number of points along edge.");
        }
      
        /**
         * @brief Destructor.
         */
        ProcessCyl::~ProcessCyl()
        {
            
        }

        void ProcessCyl::Process()
        {
            int surfTag = config["surf"].as<int>();
            set<pair<int,int> > tmp;

            int prismedge[2][3] = {{0,5,4}, {2,6,7}};
            Node zax(0,0,0,1);
            
            for (int i = 0; i < m->element[m->expDim].size(); ++i)
            {
                ElementSharedPtr el = m->element[m->expDim][i];
                int nSurf = m->expDim == 3 ? el->GetFaceCount() : 
                    el->GetEdgeCount();
                        
                for (int j = 0; j < nSurf; ++j)
                {
                    int bl = el->GetBoundaryLink(j);
                    if (bl == -1)
                    {
                        continue;
                    }

                    ElementSharedPtr bEl  = m->element[m->expDim-1][bl];
                    vector<int>      tags = bEl->GetTagList();

                    if (find(tags.begin(), tags.end(), surfTag) ==
                        tags.end())
                    {
                        continue;
                    }

                    ASSERTL0(j == 1 || j == 3, "rofl");
                    
                    // Check all edge interior points.
                    for (int k = 0; k < 3; ++k)
                    {
                        EdgeSharedPtr edge = el->GetEdge(prismedge[(j-1)/2][k]);

                        // See if edge is parallel to z axis
                        /*
                        if (fabs((*(edge->n1) - *(edge->n2)).dot(zax)) < 1e-6)
                        {
                            continue;
                        }
                        */

                        GenerateEdgeNodes(edge);
                    }
                }
            }
        }

        void ProcessCyl::GenerateEdgeNodes(EdgeSharedPtr edge)
        {
            NodeSharedPtr n1 = edge->n1;
            NodeSharedPtr n2 = edge->n2;

            int nq = config["N"].as<int>();
            double r = config["r"].as<double>();
            double t1 = atan2(n1->y, n1->x);
            double t2 = atan2(n2->y, n2->x);
            double dt, dz;

            if (t1 < -M_PI/2.0 && t2 > 0.0)
            {
                t1 += 2*M_PI;
            }
            if (t2 < -M_PI/2.0 && t1 > 0.0)
            {
                t2 += 2*M_PI;
            }
            
            dt = (t2-t1) / (nq-1);
            dz = (n2->z - n1->z) / (nq-1);

            edge->edgeNodes.resize(nq-2);
            Node dx = (*n2-*n1) * (1.0/(nq-1.0));
            for (int i = 1; i < nq-1; ++i)
            {
                edge->edgeNodes[i-1] = NodeSharedPtr(
                    new Node(0, r*cos(t1 + i*dt), r*sin(t1 + i*dt), n1->z + i*dz));
            }
            edge->curveType = LibUtilities::ePolyEvenlySpaced;
        }
    }
}
