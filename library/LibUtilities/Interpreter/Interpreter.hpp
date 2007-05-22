///////////////////////////////////////////////////////////////////////////////
//
// File Interpreter.hpp
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
// Description: Wrapper around interp.y functions for function interpretation
//
///////////////////////////////////////////////////////////////////////////////

#ifndef INTERPRETER_H
#define INTERPRETER_H

#include <cstdio>

extern "C" 
{
    // -- Routines from initial.y:

    void   yy_initialize (void);
    double yy_interpret  (const char*);
    
    void   yy_vec_init   (const char*, const char*);
    void   yy_vec_interp (const int, ...);
    
    void   yy_help       (void);
    void   yy_show       (void);
}

namespace Nektar
{
    class Interpret 
    {
    public:
        
        Interpret()
        {
            yy_initialise();
        }

        ~Interpret()
        {

        }
        
        void SetVector (const char* v, const char* f)
        { 
            yy_vec_init (v, f); 
        }
        
        void EvaluateVector (const int n ... )
        { 
            yy_vec_interp (n ... ); 
        }
        
        double EvaluateDoubleValue (const char* s)
        { 
            return yy_interpret (s); 
        }

        int EvaluateIntValue (const char* s)
        { 
            return static_cast<int>(yy_interpret (s)); 
        }
        
    private:
        
    };

#endif
