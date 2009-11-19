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
// License for the specific language governing rights and limitations under
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
#include <tinyxml/tinyxml.h>
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
        GlobalOptParam::GlobalOptParam():
            m_doGlobalMatOp(SIZE_OptimizeOperationType,false),
            m_doBlockMatOp(SIZE_OptimizeOperationType,false)
        {
        }


        /**
         * Read global optimisation parameters from a file and set up flags.
         * @param   fileName    File to read parameters from.
         */
        GlobalOptParam::GlobalOptParam(const std::string& fileName):
            m_doGlobalMatOp(SIZE_OptimizeOperationType,false),
            m_doBlockMatOp(SIZE_OptimizeOperationType,false)
        {
            TiXmlDocument doc(fileName);
            bool loadOkay = doc.LoadFile();

            ASSERTL0(loadOkay, (std::string("Unable to load file: ") +
                                fileName).c_str());

            TiXmlHandle docHandle(&doc);
            TiXmlElement* master
                        = docHandle.FirstChildElement("NEKTAR").Element();
            TiXmlElement* paramList
                        = docHandle.FirstChildElement("NEKTAR")
                            .FirstChildElement("GLOBALOPTIMIZATIONPARAMETERS")
                            .Element();

            ASSERTL0(master   , "Unable to find NEKTAR tag in file.");
            ASSERTL0(paramList, "Unable to find ELEMENTALOPTIMIZATIONPARAMETERS"
                                " tag in file.");
            int n;
            for(n = 0; n < SIZE_OptimizeOperationType; n++)
            {
                TiXmlElement* operationType = paramList->FirstChildElement(
                                      ElementalOptimizationOperationTypeMap[n]);
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
                         + std::string(ElementalOptimizationOperationTypeMap[n])
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

                        err = arrayElement->QueryIntAttribute("VALUE", &value);
                        ASSERTL0(err == TIXML_SUCCESS, (
                           std::string("Unable to read DO_BLOCK_MAT_OP "
                                       "attribute VALUE for ")
                         + std::string(ElementalOptimizationOperationTypeMap[n])
                         + std::string(".")
                        ));

                        m_doBlockMatOp[n] = (bool) value;
                    }
                }
            }
        }

    } // end of namespace
} // end of namespace
