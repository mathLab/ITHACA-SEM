#!/usr/bin/env python

"""Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

"""
This part of PyCml deals with converting CellML models into programming language code.
It is a thin executable wrapper around translators.py.
"""

import os #importing functions related to operating system
import sys #importing functions related to system (e.g. path directory etc)

# Make sure PyCml is on sys.path
pycml_path = os.path.dirname(os.path.realpath(__file__))
sys.path[0:0] = [pycml_path]

import translators #importing the main translator class
from translators import CellMLTranslator #again, just to make sure
import CellMLToNektarTranslator #importing the nektar sub-class

#This part actually runs the code

if __name__ == '__main__': #this part of the code is only run when CellMLToNektar.py is run directly, not imported
    CellMLTranslator.register(CellMLTranslator, 'C++') #this registers the main translator and it's name when called as an option in modified_translators
    CellMLTranslator.register(CellMLToNektarTranslator.CellMLToNektarTranslator, 'Nektar') #this registers the nektar subclass and it's name when called as an option in modified_translators
    if '--profile' in sys.argv: #sys.argv is the list of command line items passed to the script, this checks if a specific translator is requested
        import time, cProfile
        profile_name = '/tmp/pycml-profile-%f-%d' % (time.time(), os.getpid())
        cProfile.run('modified_translators.run()', profile_name) #the requested translator is run
    else:
        translators.run() #if no specific translator is requested, the default run function in modified_translators is run

#Below here seem to be variables that do some form of testing

    # For use in testing
    def euler(doc, t, nsteps=1000, dt=0.01):
        global tvar, state_vars, exprs #assigns tvar, state_vars and exprs to be global variables, meaning they can be used by any function
        tvar = t.free_vars[0] #takes the t-input and does "free_vars[0]" to it?
        state_vars = t.state_vars #defines state_vars as the state_vars of the t-input?
        for var in state_vars: #cycles through all the entries in "state_vars"
            var.set_value(float(var.initial_value)) #for all the entries in "state_vars" it sets the value to a float of the initial value
        tvar.set_value(0.0) #this sets the value of tvars to (0.0)?
        exprs = [e for e in doc.model.get_assignments() #I don't understand this part
                 if isinstance(e, modified_translators.mathml_apply)] #what is the mathml_apply function? Can't find it in modified_translators
        for _ in range(nsteps): #this is a for-loop through all the values up to the input of nsteps
            for expr in exprs: #this cycles through all the entries in exprs
                expr.evaluate()
            tvar.set_value(tvar.get_value() + dt)
            for var in state_vars:
                var.set_value(var.get_value() +
                              dt * var.get_value(ode=tvar))
        return

    def writefile(doc, outfn='test.cml'):
        # Write out CellML file
        st = modified_translators.open_output_stream(outfn)
        doc.xml(indent=1, stream=st)
        st.close()
        return

    def show_usage(doc):
        for comp in doc.model.component:
            for var in comp.variable:
                print var.fullname(), var._cml_usage_count


    def fix_divide_by_zero(doc):
        """
        Several models have equations of a form that may give rise to
        a divide by zero error on simulation, especially when lookup
        tables are used.  The general form is:
        
        (a * (V - v0)) / (exp(b * (V - v0)) - 1)
        
        When V = v0 this is undefined, however the limit of the
        function as V approaches v0 from either side is well-defined,
        and each limit is the same.  We approximate the limit by
        linear interpolation between values of the expression for
        (V-v0) = +/- 1e-10.
        """
        divides = [d.xml_parent
                   for d in doc.xml_xpath(u'//m:apply/m:divide')]
        for divide in divides:
            pass
        return
