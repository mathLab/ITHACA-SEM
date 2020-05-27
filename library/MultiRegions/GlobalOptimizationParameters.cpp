///////////////////////////////////////////////////////////////////////////////
//
// File GlobalOptimizationParameters.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Source file of global optimisation parameters class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif
#include <tinyxml.h>
#include <MultiRegions/GlobalOptimizationParameters.h>

namespace Nektar
{
    namespace NekOptimize
    {
        /**
         * @class GlobalOptParam
         *
         * Global optimisation parameters determines how the various matrices
         * in the spectral/hp element formulation are evaluated. For details
         * see the page on \ref optimisation Optimisation.
         */

        /**
         * No global optimisation parameters present.
         */
        GlobalOptParam::GlobalOptParam(const int nel):
            m_doGlobalMatOp(SIZE_OptimizeOperationType,false),
            m_shapeList(1,LibUtilities::eNoShapeType),
            m_shapeNumElements(1,nel)
        {
            Array<OneD, bool> set_false(1,false);
            m_doBlockMatOp = Array<OneD, Array<OneD, bool > > 
                (SIZE_OptimizeOperationType,set_false);
        }

        /**
         * Read global optimisation parameters from a file and set up flags.
         * @param   fileName    File to read parameters from.
         */
        GlobalOptParam::GlobalOptParam(const LibUtilities::SessionReaderSharedPtr& pSession, const int dim,
                                         const Array<OneD, const int> &NumShapeElements):
            m_doGlobalMatOp(SIZE_OptimizeOperationType,false)
        {
            int i;
            int numShapes = 0;
            TiXmlDocument& doc = pSession->GetDocument();

            m_shapeNumElements = NumShapeElements;

            switch (dim)
            {
            case 1:
                numShapes = 1;
                ASSERTL0(false,"Needs setting up for dimension 1");
                break;
            case 2:
                numShapes = 2;
                m_shapeList = Array<OneD, LibUtilities::ShapeType>(numShapes);
                m_shapeList[0] = LibUtilities::eTriangle;
                m_shapeList[1] = LibUtilities::eQuadrilateral;
                break;
            case 3:
                numShapes = 4;
                m_shapeList = Array<OneD, LibUtilities::ShapeType>(numShapes);
                m_shapeList[0] = LibUtilities::eTetrahedron;
                m_shapeList[1] = LibUtilities::ePyramid;
                m_shapeList[2] = LibUtilities::ePrism;
                m_shapeList[3] = LibUtilities::eHexahedron;
                break;
            }

            m_doBlockMatOp = Array<OneD, Array<OneD,bool> > (SIZE_OptimizeOperationType);
            for(i = 0; i < SIZE_OptimizeOperationType; ++i)
            {
                m_doBlockMatOp[i] = Array<OneD, bool> (numShapes,false);
            }

            TiXmlHandle docHandle(&doc);
            TiXmlElement* master
                        = docHandle.FirstChildElement("NEKTAR").Element();
            ASSERTL0(master   , "Unable to find NEKTAR tag in file.");

            TiXmlElement* paramList
                        = docHandle.FirstChildElement("NEKTAR")
                            .FirstChildElement("GLOBALOPTIMIZATIONPARAMETERS")
                            .Element();

            // If no global optimisation parameters set, we're done
            if (!paramList)
            {
                return;
            }

            // Check if there is a reference an external file and if so, load it
            TiXmlElement* source
                        = paramList->FirstChildElement("SOURCE");
            if (source)
            {
                std::string sourceFile = source->Attribute("FILE");
                TiXmlDocument docSource;
                bool loadOkay = docSource.LoadFile(sourceFile);
                ASSERTL0(loadOkay, (std::string("Unable to load file: ") +
                                                sourceFile).c_str());
                TiXmlHandle docSourceHandle(&docSource);
                master = docHandle.FirstChildElement("NEKTAR").Element();
                ASSERTL0(master   , "Unable to find NEKTAR tag in file.");

                paramList = docHandle.FirstChildElement("NEKTAR")
                                .FirstChildElement("GLOBALOPTIMIZATIONPARAMETERS")
                                .Element();
                ASSERTL0(paramList, std::string("Specified source file '"
                        + sourceFile + "' is missing an "
                        "GLOBALOPTIMIZATIONPARAMETERS tag.").c_str());
            }

            int n;
            for(n = 0; n < SIZE_OptimizeOperationType; n++)
            {
                TiXmlElement* operationType = paramList->FirstChildElement(
                                      OptimizationOperationTypeMap[n]);

                if(operationType)
                {
                    TiXmlElement* arrayElement = operationType
                                        ->FirstChildElement("DO_GLOBAL_MAT_OP");
                    if(arrayElement)
                    {
                        int value;
                        int err;

                        err = arrayElement->QueryIntAttribute("VALUE", &value);
                        ASSERTL0(err == TIXML_SUCCESS,(
                           std::string("Unable to read DO_GLOBAL_MAT_OP "
                                       "attribute VALUE for ")
                         + std::string(OptimizationOperationTypeMap[n])
                         + std::string(".")
                        ));

                        m_doGlobalMatOp[n] = (bool) value;
                    }

                    arrayElement
                        = operationType->FirstChildElement("DO_BLOCK_MAT_OP");
                    if(arrayElement)
                    {
                        int value;
                        int err;

                        switch (dim)
                        {
                        case 1:
                            break;
                        case 2:
                            err = arrayElement->QueryIntAttribute("TRI", &value);
                            ASSERTL0(err == TIXML_SUCCESS, (
                               std::string("Unable to read DO_BLOCK_MAT_OP "
                                           "attribute TRI for ")
                             + std::string(OptimizationOperationTypeMap[n])
                               + std::string(".")));

                            m_doBlockMatOp[n][0] = (bool) value;

                            err = arrayElement->QueryIntAttribute("QUAD", &value);
                            ASSERTL0(err == TIXML_SUCCESS, (
                               std::string("Unable to read DO_BLOCK_MAT_OP "
                                           "attribute QUAD for ")
                             + std::string(OptimizationOperationTypeMap[n])
                             + std::string(".")));

                            m_doBlockMatOp[n][1] = (bool) value;
                            break;
                        case 3:
                            err = arrayElement->QueryIntAttribute("TET", &value);
                            ASSERTL0(err == TIXML_SUCCESS, (
                               std::string("Unable to read DO_BLOCK_MAT_OP "
                                           "attribute TET for ")
                             + std::string(OptimizationOperationTypeMap[n])
                             + std::string(".")));

                            m_doBlockMatOp[n][0] = (bool) value;

                            err = arrayElement->QueryIntAttribute("PYR", &value);
                            ASSERTL0(err == TIXML_SUCCESS, (
                               std::string("Unable to read DO_BLOCK_MAT_OP "
                                           "attribute PYR for ")
                             + std::string(OptimizationOperationTypeMap[n])
                             + std::string(".")));

                            m_doBlockMatOp[n][1] = (bool) value;

                            err = arrayElement->QueryIntAttribute("PRISM", &value);
                            ASSERTL0(err == TIXML_SUCCESS, (
                               std::string("Unable to read DO_BLOCK_MAT_OP "
                                           "attribute PRISM for ")
                             + std::string(OptimizationOperationTypeMap[n])
                             + std::string(".")));

                            m_doBlockMatOp[n][2] = (bool) value;

                            err = arrayElement->QueryIntAttribute("HEX", &value);
                            ASSERTL0(err == TIXML_SUCCESS, (
                               std::string("Unable to read DO_BLOCK_MAT_OP "
                                           "attribute HEX for ")
                             + std::string(OptimizationOperationTypeMap[n])
                             + std::string(".")));

                            m_doBlockMatOp[n][3] = (bool) value;
                            break;

                            break;

                        }
                    }
                }
            }
        }


    } // end of namespace
} // end of namespace
