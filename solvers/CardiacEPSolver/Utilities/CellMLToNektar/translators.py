# We want 1/2==0.5
from __future__ import division

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
This part of PyCml deals with converting CellML models into
programming language code, primarily C++ compatible with Chaste, but
supporting a few other languages also (and easily extensible).
"""

import optparse #adds options to the command line parsing stuff
import os #imports OS functions
import re #imports regular expression operations
import time #lets you do various time-related operations
import sys #imports system functions (e.g. which path we are on)
from cStringIO import StringIO #"This module implements a file-like class, StringIO, that reads and writes a string buffer (also known as memory files)."

# Common CellML processing stuff, all of these are in the /pycml folder
import pycml
from pycml import *  # Put contents in the local namespace as well
import optimize
import processors 
import validator


#this sets the version
__version__ = "$Revision: 25950 $"[11:-2]


#this function generates a version comment
def version_comment(note_time=True):
    """Generate a version comment, with optional time info."""
    if note_time:
        t = '\non ' + time.asctime()
    else:
        t = ''
    text = """Processed by pycml - CellML Tools in Python
    (translators: %s, pycml: %s, optimize: %s)%s""" % (
        __version__, pycml.__version__, optimize.__version__, t)
    return text

#some debugging function?
def debugexpr(e):
    "For debugging."
    v = None
    if isinstance(e, cellml_variable):
        v = e
    elif isinstance(e, mathml_apply):
        v = e.assigned_variable()
    if v:
        r = (v==e, v.name, v.get_usage_count())
    else:
        r = (False, '', -1)
    return r


class TranslationError(RuntimeError):
    """Error thrown if CellML translation fails."""
    pass

class ConfigurationError(ValueError):
    """Error thrown if configuration file is invalid."""
    pass



class CellMLTranslator(object):
    """
    Base class for translators from CellML to programming languages.

    Provides various methods & attributes that can be overridden to
    achieve the desired output language and style.

    Also contains a registration system for subclasses, so the
    command-line client can know what translators are available.  See
    the register method for more information.
    """
    
    translators = {} #creates an empty dictionary, into which the classes are registered (the registering is done in modified_translate)
    class NameAlreadyRegistered(ValueError):
        pass


#this is the function used by modified_translate to register each module into the translators dictionary
    @classmethod
    def register(cls, subclass, name):
        """Register a new translator subclass.

        Registers the subclass `subclass' with name `name' in the
        translators class attribute of CellMLTranslator.  If the name
        given is already in use, raises NameAlreadyRegistered.
        """
        if name in cls.translators:
            raise cls.NameAlreadyRegistered(name)
        cls.translators[name] = subclass
        return
    

    @staticmethod
    def generate_interface(doc, solver_info):
        """Generate an interface component connecting the model to whatever will use it.
        
        Stub method that subclasses can override to implement this functionality.
        """
        pass

    ###########################
    # Various language tokens #
    ###########################
    STMT_END = ';'        # End of statement
    EQ_ASSIGN = ' = '     # Assignment operator
    COMMENT_START = '// ' # Start of a 1 line comment
    DOXYGEN_COMMENT_START = '//! ' # Start of a 1 line Doxygen comment

    # Variable types
    TYPE_DOUBLE = 'NekDouble '
    TYPE_VOID = 'void '
    TYPE_CONST_DOUBLE = 'const NekDouble '
    TYPE_CONST_UNSIGNED = 'const unsigned '

    # Special constants
    TRUE = 'true'
    FALSE = 'false'
    PI = 'M_PI'
    E = 'M_E'
    NOT_A_NUMBER = 'NAN' # GNU extension, but fairly common
    
    # Whether the target language uses a subsidiary file, such as
    # a header file in C/C++
    USES_SUBSIDIARY_FILE = False
    # Mapping from primary file extension to subsidiary file extension
    FILE_EXTENSIONS = {'cpp': 'h',
                       'c': 'h',
                       'cxx': 'hxx'}

    def __init__(self, add_timestamp=True, options=None):
        """Create a translator."""
        self.options = options
        # Initially output should not be indented
        self.indent_level = 0
        # Character to indent with
        self.indent_char = ' '
        # No. of occurrences of indent_char per indent_level
        self.indent_factor = 4
        # Whether to use lookup tables where possible
        self.use_lookup_tables = True
        # Whether to add a timestamp comment to generated files
        self.add_timestamp = add_timestamp
        # Main output goes to the main file by default
        self._main_output_to_subsidiary = False

    def error(self, lines, xml=None):
        """Raise a translation error.

        lines is a list of strings describing what went wrong.
        A TranslationError with that message will be raised.

        If xml is given, it should be an element, which will be
        pretty-printed and included in the error.
        """
        if xml is not None:
            lines.extend(xml.xml(indent = u'yes',
                                 omitXmlDeclaration = u'yes').split('\n'))
        raise TranslationError('\n'.join(lines))

    @property
    def config(self):
        """Get the current document's configuration store."""
        return getattr(self.doc, '_cml_config', None)

    def translate(self, doc, model_filename, output_filename=None,
                  subsidiary_file_name=None,
                  class_name=None, v_variable=None,
                  continuation=None,
                  lookup_method_prefix='', row_lookup_method=False,
                  lt_index_uses_floor=True, constrain_table_indices=False):
        """Generate code for the given model.

        doc should be an instance of cellml_model representing a
        valid CellML model, such as might be produced from a call
        to
        >>> valid, doc = validator.CellMLValidator().validate(
        ...     model_filename, True)

        model_filename is the filename of the input model.
        The output program will by default be generated in the same
        folder, but with a different extension.  This can be
        overridden by supplying the output_filename keyword argument.

        By default the name of the class representing the model will
        be derived from the model name.  This can be overridden by
        passing an alternative as the class_name argument.

        The variable representing the transmembrane potential should
        be passed in using the v_variable argument.

        By default this method will perform some setup and then call
            self.output_top_boilerplate()
            self.output_mathematics()
            self.output_bottom_boilerplate()
        To alter this, pass a callable as the continuation parameter;
        this will then be called instead.

        lookup_method_prefix and row_lookup_method can be used to
        customise some aspects of lookup table usage.  The former is
        used by the Chaste translator to place lookup tables within a
        singleton class, while the latter can improve cache
        performance by looking up all tables in a single call, and
        returning an array of the results.
        
        lt_index_uses_floor specifies whether to use the floor()
        function to calculate the index into the lookup tables, or
        just cast to unsigned.

        constrain_table_indices specifies whether to throw an
        exception if lookup table index variables go outside the
        bounds specified (default), or just to cap them at the bounds.
        """
        self.doc = doc
        self.model = doc.model
        # Name of the class that will represent this model
        if class_name is None:
            self.class_name = u'CML_' + doc.model.name.replace('-', '_')
        else:
            self.class_name = class_name
        # Figure out the free & state vars in this model
        self.free_vars = doc.model.find_free_vars()
        self.state_vars = doc.model.find_state_vars()
        if len(self.free_vars) > 1:
            self.error(["Model has more than one free variable; exiting.",
                        "Free vars:" + str(self.free_vars)])
        if len(self.free_vars) == 0:
            if self.model.get_option('protocol'):
                # We may be dealing with an algebraic model; check for an Unknown variable
                for var in self.model.get_all_variables():
                    if var.get_type() == VarTypes.Unknown:
                        self.free_vars.append(var)
            if len(self.free_vars) != 1:
                self.error(["Model has no free variable; exiting."])
        # If only a single component, don't add it to variable names
        self.single_component = (len(getattr(self.model, u'component', [])) == 1)
        # Find the (index of the) transmembrane potential
        self.v_variable = v_variable
        if self.v_variable:
            self.v_variable_name = (v_variable.component.name, v_variable.name)
        else:
            self.v_variable = None
        for i, var in enumerate(self.state_vars):
            if var is v_variable:
                self.v_index = i
                break
        else:
            self.v_index = -1
        # Check to see if we're using lookup tables
        self.lookup_method_prefix = lookup_method_prefix
        self.row_lookup_method = row_lookup_method
        self.lt_index_uses_floor = lt_index_uses_floor
        self.constrain_table_indices = constrain_table_indices
        self.scan_for_lookup_tables()
        if not doc.lookup_tables:
            # No tables found
            self.use_lookup_tables = False

        # Extra configuration hook
        self.final_configuration_hook()

        # Open the output file(s)
        if output_filename is None:
            output_filename = self.output_file_name(model_filename)
        if self.USES_SUBSIDIARY_FILE:
            output_filename, self.subsidiary_filename = self.subsidiary_file_name(output_filename)
        self.output_filename = output_filename
        self.out = open_output_stream(output_filename)
        if self.USES_SUBSIDIARY_FILE:
            self.out2 = open_output_stream(self.subsidiary_filename)

        # Translate
        if continuation:
            continuation()
        else:
            self.output_top_boilerplate()
            self.output_mathematics()
            self.output_bottom_boilerplate()
        close_output_stream(self.out)
        if self.USES_SUBSIDIARY_FILE:
            close_output_stream(self.out2)
        return

    def final_configuration_hook(self):
        """A hook for subclasses to do some final configuration."""
        return

    def output_file_name(self, model_filename):
        """Generate a name for our output file, based on the input file."""
        return os.path.splitext(model_filename)[0] + '.cpp'
    
    def subsidiary_file_name(self, output_filename):
        """Generate a name for the subsidiary output file, based on the main one.
        
        Returns a pair (main_output_file_name, subsidiary_file_name).  This is in
        case the user specifies (e.g.) a .hpp file as the main output - we consider
        the main output to be the .cpp file.
        """
        base, ext = os.path.splitext(output_filename)
        ext = ext[1:] # Remove the '.'
        try:
            new_ext = self.FILE_EXTENSIONS[ext]
            swap = False
        except KeyError:
            swap = True
            for key, val in self.FILE_EXTENSIONS.iteritems():
                if val == ext:
                    new_ext = key
                    break
            else:
                # Turn off usage of subsidiary file
                self.USES_SUBSIDIARY_FILE = False
                return output_filename, None
        subsidiary_filename = base + '.' + new_ext
        if swap:
            output_filename, subsidiary_filename = subsidiary_filename, output_filename
        return output_filename, subsidiary_filename

    def send_main_output_to_subsidiary(self, to_subsidiary=True):
        """Set subsequent main-file writes to go to the subsidiary file instead.
        
        Supplying a False argument reverts this setting.
        
        Has no effect if not using a subsidiary file.
        """
        self._main_output_to_subsidiary = to_subsidiary

    def writeln(self, *args, **kwargs):
        """Write a line to our output file.

        Takes any number of strings as input, and writes them out to file.

        Unless the keyword argument indent is given as False, then the
        output will be indented to the level set by self.set_indent().
        Setting indent_level will override this value.
        Setting indent_offset will adjust the current level temporarily.

        If nl is set to False then a newline character will not be
        appended to the output.
        
        If subsidiary=True, then the line will be written to the subsidiary
        output file instead of the main one.  An error will be raised if
        there is no subsidiary output file.
        """
        if kwargs.get('subsidiary', False) or self._main_output_to_subsidiary:
            if not self.USES_SUBSIDIARY_FILE:
                self.error('Tried to write to non-existent subsidiary file')
            else:
                target = self.out2
        else:
            target = self.out
        indent = kwargs.get('indent', True)
        nl = kwargs.get('nl', True)
        if indent:
            level = kwargs.get('indent_level', self.indent_level)
            level += kwargs.get('indent_offset', 0)
            target.write(self.indent_char * self.indent_factor * level)
        target.write(''.join(map(str, args)))
        if nl:
            target.write('\n')

    def write(self, *args):
        """Write to our output file.

        This variant does not indent the output, or add a newline.
        """
        self.writeln(indent=False, nl=False, *args)
    
    def capture_output(self):
        """Make subsequent output operations write to a string buffer."""
        self._original_out = self.out
        self.out = StringIO()
    
    def get_captured_output(self):
        """Stop capturing output, and return what was captured as a string."""
        output = self.out.getvalue()
        self.out = self._original_out
        return output

    def output_comment(self, *args, **kwargs):
        """Output a (multi-line) string as a comment."""
        start = kwargs.get('start', self.COMMENT_START)
        if kwargs.get('pad', False):
            start = ' ' + start
        comment = ''.join(map(str, args))
        lines = comment.split('\n')
        for line in lines:
            self.writeln(start, line, **kwargs)

    def output_doxygen(self, *args, **kwargs):
        """Output a (multi-line) string as a Doxygen comment."""
        kwargs['start'] = self.DOXYGEN_COMMENT_START
        self.output_comment(*args, **kwargs)

    def set_indent(self, level=None, offset=None):
        """Set the indentation level for subsequent writes.

        If offset is given, adjust the level by that amount, otherwise
        set it to an absolute value.
        """
        if offset is not None:
            self.indent_level += offset
        else:
            self.indent_level = level

    def code_name(self, var, ode=False, prefix=None):
        """
        Return the full name of var in a form suitable for inclusion in a
        source file.
        
        If ode is True then return the name of the derivative of var
        instead.  We go directly to the source variable in this case,
        rather than including intermediate assignment statements as is
        done for connections.
        """
        if prefix is None:
            prefix = ['var_', 'd_dt_'][ode]
        if ode:
            var = var.get_source_variable(recurse=True)
        name = prefix + var.fullname(cellml=True)
        return name
    
    def varobj(self, varname):
        """Return the variable object that has code_name varname."""
        return cellml_variable.get_variable_object(self.model, varname)

    def var_display_name(self, var):
        """Return a display name for the given variable.
        
        If it has an oxmeta name, uses that.  Otherwise, looks first for a bqbiol:is annotation,
        or uses the cmeta:id if present, or the name attribute if not.  If there is an interface
        component, strip the name of it out of the display name.
        """
        if var.oxmeta_name:
            name = var.oxmeta_name
        else:
            for uri in var.get_rdf_annotations(('bqbiol:is', NSS['bqbiol'])):
                if '#' in uri:
                    name = uri[1 + uri.rfind('#'):]
                    break
            else:
                if hasattr(var, u'id') and var.id:
                    name = var.id
                else:
                    name = var.name
        iface = getattr(self.model, 'interface_component_name', '#N/A#')
        if name.startswith(iface):
            name = name[len(iface)+2:]
        return name

    @property
    def include_guard(self):
        """
        Get the include guard (for C/C++ output) for this cell model,
        based on the class name.
        """
        return self.class_name.upper() + '_HPP_'
    
    def output_top_boilerplate(self):
        """Output top boilerplate."""
        self.writeln('#ifndef _', self.include_guard, '_')
        self.writeln('#define _', self.include_guard, '_\n')
        self.output_comment('Model: ', self.model.name)
        self.output_comment(version_comment(self.add_timestamp))
        self.writeln()
        self.writeln('#include <cmath>')
        self.writeln('#include "AbstractOdeSystem.hpp"')
        self.writeln('#include "Exception.hpp"')
        self.writeln('#include "AbstractStimulusFunction.hpp"\n')
        self.writeln('class ', self.class_name, ' : public AbstractOdeSystem')
        self.writeln('{')
        self.writeln('private:')
        self.writeln('AbstractStimulusFunction *mpStimulus;\n',
                     indent_offset=1)
        self.writeln('public:')
        self.set_indent(1)
        self.writeln('const static unsigned _NB_OF_STATE_VARIABLES_ = ',
                   str(len(self.state_vars)), ';\n')
        self.writeln('//', ('-'*66))
        self.writeln('// Methods')
        self.writeln('//', ('-'*66), '\n')
        # Constructor
        self.writeln('', self.class_name,
                   '(AbstractStimulusFunction *stim)')
        self.writeln('    : AbstractOdeSystem(_NB_OF_STATE_VARIABLES_, ',
                     self.v_index, ')')
        self.open_block()
        self.writeln('mpStimulus = stim;\n')
        for var in self.state_vars:
            self.writeln('mVariableNames.push_back("', var.name, '");')
            self.writeln('mVariableUnits.push_back("', var.units, '");')
            init_val = getattr(var, u'initial_value', None)
            if init_val is None:
                init_comm = ' // Value not given in model'
                # Don't want compiler error, but shouldn't be a real number
                init_val = self.NOT_A_NUMBER
            else:
                init_comm = ''
            self.writeln('mInitialConditions.push_back(', init_val, ');',
                       init_comm, '\n')
        if self.use_lookup_tables: self.output_lut_generation()
        self.close_block()
        # Destructor
        self.writeln('~', self.class_name, '(void)')
        self.open_block()
        if self.use_lookup_tables: self.output_lut_deletion()
        self.close_block()
        # Lookup table declarations & methods
        if self.use_lookup_tables:
            self.output_lut_declarations()
            self.output_lut_methods()
        # Evaluation function
        self.writeln('void EvaluateYDerivatives (')
        self.writeln('        double ', self.code_name(self.free_vars[0]), ',')
        self.writeln('        const std::vector<double> &rY,')
        self.writeln('        std::vector<double> &rDY)')
        self.open_block()
        self.writeln('// Inputs:')
        self.writeln('// Time units: ', self.free_vars[0].units)
        for i, var in enumerate(self.state_vars):
            self.writeln('double ', self.code_name(var),
                         ' = rY[', str(i), '];')
            self.writeln('// Units: ', var.units, '; Initial value: ',
                         getattr(var, u'initial_value', 'Unknown'))
        self.writeln()
        if self.use_lookup_tables:
            self.output_table_index_generation()
        return

    def output_mathematics(self):
        """Output the mathematics in this model."""
        self.writeln(self.COMMENT_START, 'Mathematics')
        for expr in self.model.get_assignments():
            # Check this expression is actually used; don't output if not
            var = None
            if isinstance(expr, mathml_apply) and expr.is_assignment():
                var = expr.assigned_variable()
            elif isinstance(expr, cellml_variable):
                var = expr
            if not (var and var.get_usage_count() == 0):
                self.output_assignment(expr)
        return

    def output_bottom_boilerplate(self):
        """Output bottom boilerplate"""
        self.writeln('\n')
        for i, var in enumerate(self.state_vars):
            self.writeln('rDY[', str(i), '] = ', self.code_name(var, True),
                         ';')
        self.close_block()
        self.set_indent(offset=-1)
        self.writeln('};\n')
        self.writeln('#endif')
        return

    def output_assignment(self, expr):
        """Output an assignment expression."""
        if isinstance(expr, cellml_variable):
            # self.writeln('Test')
            # This may be the assignment of a mapped variable, or a constant
            t = expr.get_type()
            if t == VarTypes.Mapped:
                # self.writeln(self.code_name(expr))
                self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(expr),
                             self.EQ_ASSIGN,
                             self.code_name(expr.get_source_variable()),
                             self.STMT_END, nl=False)
                self.output_comment(expr.units, indent=False, pad=True)
            elif t == VarTypes.Constant:
                # self.writeln('Test')
                # self.writeln(self.code_name(expr))


                # setting various stim-parameters to zero, since Nektar uses its own stimulus
                if 'stim_amplitude' in str(self.code_name(expr)):
                    self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(expr),
                                 self.EQ_ASSIGN, nl=False)
                    self.output_number(0)
                    self.writeln(self.STMT_END, indent=False, nl=False)
                    self.output_comment(expr.units, indent=False, pad=True)
                else:
                    self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(expr),
                                 self.EQ_ASSIGN, nl=False)
                    self.output_number(expr.initial_value)
                    self.writeln(self.STMT_END, indent=False, nl=False)
                    self.output_comment(expr.units, indent=False, pad=True)
        else:
            # This is a mathematical expression
            self.writeln(self.TYPE_CONST_DOUBLE, nl=False)
            opers = expr.operands()
            self.output_lhs(opers.next())
            self.write(self.EQ_ASSIGN)
            self.output_expr(opers.next(), False)
            self.writeln(self.STMT_END, indent=False, nl=False)
            #1365: add a comment with the LHS units
            self.output_comment(expr._get_element_units(expr.eq.lhs, return_set=False).description(),
                                indent=False, pad=True)

    def output_lhs(self, expr):
        """Output the left hand side of an assignment expression."""
        if expr.localName == 'ci':
            self.output_variable(expr)
        elif expr.operator().localName == 'diff':
            self.write(self.code_name(expr.operator().dependent_variable, ode=True))

    def output_variable(self, ci_elt, ode=False):
        """Output a ci element, i.e. a variable lookup."""
        self.write(self.code_name(ci_elt.variable, ode=ode))

    def output_expr(self, expr, paren):
        """Output the expression expr.
        
        If paren is True then the context has requested parentheses around the
        output; if expr requires them then they will be added.
        """
        if self.use_lookup_tables and self.is_lookup_table(expr):
            self.output_table_lookup(expr, paren)
        elif isinstance(expr, mathml_apply):
            self.output_apply(expr, paren)
        elif isinstance(expr, mathml_piecewise):
            self.output_piecewise(expr, paren)
        elif isinstance(expr, mathml_ci):
            self.output_variable(expr)
        elif expr.localName == u'cn':
            self.output_number(expr)
        elif expr.localName == u'degree':
            # <degree> is just a wrapper around an expression
            self.output_expr(child_i(expr, 1), paren)
        elif expr.localName == u'logbase':
            # <logbase> is just a wrapper around an expression
            self.output_expr(child_i(expr, 1), paren)
        elif expr.localName == u'true':
            self.write(self.TRUE)
        elif expr.localName == u'false':
            self.write(self.FALSE)
        elif expr.localName == u'pi':
            self.write(self.PI)
        elif expr.localName == u'exponentiale':
            self.write(self.E)
        else:
            self.error(["Unsupported expression element " + expr.localName],
                       xml=expr)

    def output_number(self, expr):
        """Output the plain number expr.
        
        We make all constants parse as doubles to avoid problems with
        integer division or numbers too large for the int type.
        
        Negative numbers will be prefixed by a space to avoid unwanted
        decrement operations.
        """
        n = self.eval_number(expr)
        num = "%.17g" % n
        if num[0] == '-':
            num = ' ' + num
        if not '.' in num and not 'e' in num:
            num = num + '.0'
        self.write(num)

    def eval_number(self, expr):
        """Evaluate a number.

        If a (unicode) string, convert to float.
        If a cn element, call its evaluate method.
        """
        if isinstance(expr, mathml_cn):
            return expr.evaluate()
        else:
            return float(unicode(expr))

    # Map from operator element names to C++ function names,
    # where the translation is straightforward.
    function_map = {'power': 'pow', 'abs': 'fabs', 'ln': 'log', 'exp': 'exp',
                    'floor': 'floor', 'ceiling': 'ceil',
                    'factorial': 'factorial', # Needs external definition
                    'not': '!', 'rem': 'fmod',
                    'sin': 'sin', 'cos': 'cos', 'tan': 'tan',
                    'sec': '1/cos', 'csc': '1/sin', 'cot': '1/tan',
                    'sinh': 'sinh', 'cosh': 'cosh', 'tanh': 'tanh',
                    'sech': '1/cosh', 'csch': '1/sinh', 'coth': '1/tanh',
                    'arcsin': 'asin', 'arccos': 'acos', 'arctan': 'atan',
                    'arcsinh': 'asinh', 'arccosh': 'acosh', 'arctanh': 'atanh'}
    # Inverse reciprocal trig functions; these are represented by
    # key(x) = function_map[val](1/x)
    recip_trig = {'arcsec': 'arccos', 'arccsc': 'arcsin', 'arccot': 'arctan',
                  'arcsech': 'arccosh', 'arccsch': 'arcsinh', 'arccoth': 'arctanh'}
    # Operators
    nary_ops   = {'plus': '+', 'times': '*',
                  'and': '&&', 'or': '||'}
    binary_ops = {'divide': '/',
                  'xor': '!=', 'eq': '==', 'neq': '!=',
                  'geq': '>=','leq': '<=','gt': '>','lt': '<'}

    def output_apply(self, expr, paren):
        """Output an <apply> expression.
        
        paren is True if the context has requested parentheses.
        """
        op = expr.operator()
        if op.localName in self.function_map:
            self.output_function(self.function_map[op.localName],
                                 expr.operands(), paren)
        elif op.localName in self.recip_trig:
            self.output_function(self.function_map[self.recip_trig[op.localName]],
                                 expr.operands(), paren, reciprocal=True)
        elif op.localName == u'root':
            self.output_root(expr, paren)
        elif op.localName == u'log':
            self.output_log(expr, paren)
        elif op.localName in self.nary_ops:
            self.output_nary_operator(self.nary_ops[op.localName],
                                      expr.operands(), paren)
        elif op.localName in self.binary_ops:
            self.output_binary_operator(self.binary_ops[op.localName],
                                        expr.operands(), paren, expr)
        elif op.localName == u'minus':
            self.output_minus(expr, paren)
        elif op.localName == u'diff':
            # ODE occuring on the RHS
            self.write(self.code_name(op.dependent_variable, ode=True))
        else:
            # Unrecognised operator
            self.error(["Unsupported operator element " + str(op.localName)], xml=expr)

    def output_function(self, func_name, args, paren, reciprocal=False):
        """Output a function call with name func_name and arguments args.
        
        Parentheses are not required so paren is ignored.
        If reciprocal is True then pass the reciprocal of each arg to
        func_name.
        """
        self.write(func_name + '(')
        comma = False
        for arg in args:
            if comma: self.write(', ')
            else: comma = True
            if reciprocal:
                self.write('1/')
                self.output_expr(arg, True)
            else:
                self.output_expr(arg, False)
        self.write(')')

    def output_nary_operator(self, operator, operands, paren):
        """Output an n-ary operator (using infix notation).
        
        If paren is True, enclose the output in parentheses.
        """
        # TODO: Optimise - to use expm1(x) for computing exp(x)-1
        self.open_paren(paren)
        op = False
        for operand in operands:
            if op: self.write(' ' + operator + ' ')
            else: op = True
            self.output_expr(operand, True)
        self.close_paren(paren)

    def output_unary_operator(self, operator, operand, paren):
        """Output a unary operator (using prefix notation)."""
        self.open_paren(paren)
        self.write(operator)
        self.output_expr(operand, True)
        self.close_paren(paren)

    def output_binary_operator(self, operator, operands, paren, expr):
        """Output a binary operator.
        
        As output_nary_operator, but checks that len(list(operands)) == 2.
        """
        operands = list(operands)
        if len(operands) != 2:
            self.error(["Binary operator" + operator +
                        "does not have 2 operands."], xml=expr)
        self.output_nary_operator(operator, operands, paren)

    special_roots = {2: 'sqrt', 3: 'cbrt'}
    def output_root(self, expr, paren):
        """Output a root taken to some degree.

        If a degree qualifier element is not provided, uses default 2.
        """
        if hasattr(expr, u'degree'):
            # A degree is given.  Compute x^(1/b)
            # TODO: Optimise for when b==2 (sqrt) or b==3 (cbrt)
            # Try to evaluate expr.degree, and if the result is a key
            # of self.special_roots, use the value as the function to call.
            self.write('pow(')
            self.output_expr(expr.operands().next(), False)
            self.write(', 1/')
            self.output_expr(expr.degree, True)
            self.write(')')
        else:
            # Compute square root
            self.output_function('sqrt', expr.operands(), paren)

    def output_log(self, expr, paren):
        """Output a logarithm to the given base, which defaults to base 10."""
        if hasattr(expr, u'logbase'):
            # A base is provided.  Use the identity log_b(x) = log(x)/log(b)
            # TODO: Optimise for log2(x)
            self.open_paren(paren)
            self.output_function('log', expr.operands(), paren)
            self.write('/log(')
            self.output_expr(expr.logbase, False)
            self.write(')')
            self.close_paren(paren)
        else:
            # Use base 10
            self.output_function('log10', expr.operands(), paren)

    def output_minus(self, expr, paren):
        """Output either a unary or binary minus.

        Which is chosen depends on the number of operands.
        """
        operands = list(expr.operands())
        if len(operands) == 1:
            self.output_unary_operator('-', operands[0], paren)
        else:
            self.output_binary_operator('-', operands, paren, expr)

    def output_piecewise(self, expr, paren):
        """Output the piecewise expression expr.

        We use a cascading ternary if expression for simplicity.
        """
        self.open_paren(paren)
        for piece in getattr(expr, u'piece', []):
            self.output_expr(child_i(piece, 2), True) # Condition
            self.write(' ? ')
            self.output_expr(child_i(piece, 1), True) # Result
            self.write(' : ')
        if hasattr(expr, u'otherwise'):
            self.output_expr(child_i(expr.otherwise, 1), True) # Default case
        else:
            self.write(self.NOT_A_NUMBER)
        self.close_paren(paren)

    def open_paren(self, paren):
        if paren: self.write('(')
    def close_paren(self, paren):
        if paren: self.write(')')

    def open_block(self, **kwargs):
        """Open a new code block and increase indent."""
        self.writeln('{', **kwargs)
        self.set_indent(offset=1)
    def close_block(self, blank_line=True, **kwargs):
        """Close a code block and decrease indent."""
        self.set_indent(offset=-1)
        self.writeln('}', **kwargs)
        if blank_line:
            self.writeln(**kwargs)
        return

    ##############################
    # Dependency related methods #
    ##############################

    # These methods allow us to calculate which equations must be
    # output in order to compute a given set of quantities.
    def calculate_extended_dependencies(self, nodes, prune=[], prune_deps=[]):
        """Method moved to cellml_model."""
        
        # for var in nodes:
        #     self.writeln(var)
        # self.writeln()        

        # for var in self.model.calculate_extended_dependencies(nodes, prune, prune_deps):
        #     self.writeln(var)
        # self.writeln()


        return self.model.calculate_extended_dependencies(nodes, prune, prune_deps)

    def output_equations(self, nodeset):
        """Output the mathematics described by nodeset.

        nodeset represents a subset of the assignments in the model.
        Output assignments in the order given by a topological sort,
        but only include those in nodeset.

        Since a set of assignments is given, this method does not
        check whether variables are used - it is assumed that only
        assignments that are actually wanted are given in nodeset.
        """
        for expr in (e for e in self.model.get_assignments() if e in nodeset):
            self.output_assignment(expr)
        return

    def _vars_in(self, expr):
        """Return a set of variable objects used in the given expression.

        Will include state variables.  If the expression includes a derivative, the defining equation
        for that derivative will be included in the set.  Also if an expression is being
        replaced by a lookup table, this will only include the table key variable.
        """
        res = set()
        if self.use_lookup_tables and isinstance(expr, mathml) and self.is_lookup_table(expr):
            key_var = self.varobj(expr.getAttributeNS(NSS['lut'], u'var'))
            key_var = key_var.get_source_variable(recurse=True)
            res.add(key_var)
        elif isinstance(expr, mathml_ci):
            varobj = getattr(expr, '_cml_variable', None)
            if not varobj:
                varname = unicode(expr)
                varobj = self.varobj(varname.strip())
            if varobj:
                res.add(varobj)
        elif isinstance(expr, mathml_apply) and expr.operator().localName == u'diff':
            dep_varname = unicode(expr.ci)
            varobj = self.varobj(dep_varname.strip())
            res.add(varobj.get_ode_dependency(self.free_vars[0]))
        elif hasattr(expr, 'xml_children'):
            for child in expr.xml_children:
                res.update(self._vars_in(child))
        return res


    ########################
    # Lookup table methods #
    ########################

    # Lookup tables should be done in a cache- and memory-
    # efficient manner.
    #
    # Cache: Have one block of memory allocated for all tables with a
    # given index variable, such that entries are found at i*N+j,
    # where N is the no. of tables in the block, i is the index into a
    # table, and j is the table to read.  Change how lookups are done,
    # such that the lookup method is called once and returns a pointer
    # to the (i*N)'th entry.  Places where we now call the method then
    # index this pointer by j.
    #    The 'one block' part is done by default.
    #    The better lookup method is activated by --row-lookup-method.
    #
    # Memory: Extract the lookup tables into a separate class (in the
    # same .cpp file though).  This can then be made a singleton class
    # in a multi-cellular context.
    #    Chaste code generation has the option to do this, enabled by
    #    default.  Use --no-separate-lut-class to disable.

    def scan_for_lookup_tables(self):
        """Search for lookup tables used in this document.

        Store a list of suitable expressions on the document root.
        Generate a dictionary mapping tables to their index variables.
        """
        doc = self.doc
        # Get list of suitable expressions
        doc.lookup_tables = doc.xml_xpath(u"//*[@lut:possible='yes']")
        doc.lookup_tables.sort(cmp=element_path_cmp)
        # Map table keys (min, max, step, var) to an index variable
        doc.lookup_table_indexes = {}
        # Count the no. of lookup tables with each index variable
        doc.lookup_tables_num_per_index = {}
        if not doc.lookup_tables:
            # If no suitable expressions, we're done
            return
        # Search for table index variables already assigned
        table_indexes = [int(getattr(expr, u'table_index', -1))
                         for expr in doc.lookup_tables]
        tidx = max(table_indexes) + 1
        # Search for table names already assigned
        table_numbers = {}
        for expr in doc.lookup_tables:
            if hasattr(expr, u'table_name'):
                idx = expr.table_index
                table_numbers[idx] = max(table_numbers.get(idx, 0), 1 + int(expr.table_name))
        # Now assign new names, and construct mapping from tables to index variables
        for expr in doc.lookup_tables:
            # Get a suitable table index variable
            comp = expr.get_component()
            var = comp.get_variable_by_name(expr.var)
            var = var.get_source_variable(recurse=True)
            key = (expr.min, expr.max, expr.step, var)
            if not key in doc.lookup_table_indexes:
                var._cml_modifiable = True # Table index variables shouldn't be const, in case we constrain to table bounds
                if hasattr(expr, u'table_index'):
                    doc.lookup_table_indexes[key] = expr.table_index
                else:
                    doc.lookup_table_indexes[key] = unicode(tidx)
                    tidx += 1
                    expr.xml_set_attribute((u'lut:table_index', NSS['lut']),
                                           doc.lookup_table_indexes[key])
            # Get a table name, unique over all tables with this index variable
            if not hasattr(expr, u'table_name'):
                tnum = table_numbers.get(doc.lookup_table_indexes[key], 0)
                expr.xml_set_attribute((u'lut:table_name', NSS['lut']), unicode(tnum))
                table_numbers[doc.lookup_table_indexes[key]] = tnum + 1
        # Re-number table indices so they are contiguous starting from 0.
        table_index_map = {}
        table_name_map = {}
        tidx = 0
        for key in sorted(doc.lookup_table_indexes.keys()):
            idx = unicode(tidx)
            table_index_map[doc.lookup_table_indexes[key]] = idx
            table_name_map[idx] = {}
            doc.lookup_table_indexes[key] = idx
            doc.lookup_tables_num_per_index[idx] = 0
            tidx += 1
        # Make sure each lookup table is only listed once in doc.lookup_tables,
        # so we don't get 2 tables for the same expression!
        # Also re-number table names so they are contiguous starting from 0 for each table index.
        candidates = doc.lookup_tables[:]
        doc.lookup_tables = []
        listed = set()
        for expr in candidates:
            tid = (expr.table_index, expr.table_name)
            if not tid in listed:
                listed.add(tid)
                doc.lookup_tables.append(expr)
                # Renumber
                expr.table_index = table_index_map[expr.table_index]
                table_name_map[expr.table_index][expr.table_name] = unicode(doc.lookup_tables_num_per_index[expr.table_index])
                expr.table_name = table_name_map[expr.table_index][expr.table_name]
                # Increment count for this index variable
                doc.lookup_tables_num_per_index[expr.table_index] += 1
            else:
                # Just renumber to match the new id for this expression
                expr.table_index = table_index_map[expr.table_index]
                expr.table_name = table_name_map[expr.table_index][expr.table_name]
        return

    def lut_access_code(self, table_index, table_name, i):
        """Get the code for accessing the i'th element of the given table.
        """
        return '_lookup_table_%s[%s][%s]' % (table_index, i, table_name)
    
    def lut_parameters(self, key):
        """Get the bounds and step size for a particular table.
        
        key should be a key into self.lookup_table_indices.
        Returns (min, max, step, step_inverse) suitable for putting in generated code.
        """
        return key[0:3] + [unicode(1 / float(key[2]))]
    
    def lut_size_calculation(self, min, max, step):
        """Return the equivalent of '1 + (unsigned)((max-min)/step+0.5)'."""
        return '1 + (unsigned)((%s-%s)/%s+0.5)' % (max, min, step)

    def output_lut_generation(self, only_index=None):
        """Output code to generate lookup tables.

        There should be a list of suitable expressions available as self.doc.lookup_tables,
        to save having to search the whole model.
        
        If only_index is given, only generate tables using the given table index key.
        """
        # Don't use table lookups to generate the tables!
        self.use_lookup_tables = False
        # Allocate memory for tables
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            if only_index is None or only_index == idx:
                min, max, step, _ = self.lut_parameters(key)
                self.writeln(self.TYPE_CONST_UNSIGNED, '_table_size_', idx, self.EQ_ASSIGN,
                             self.lut_size_calculation(min, max, step), self.STMT_END)
                self.writeln('_lookup_table_', idx, self.EQ_ASSIGN, 'new double[_table_size_', idx,
                             '][', self.doc.lookup_tables_num_per_index[idx], ']', self.STMT_END)
        # Generate each table in a separate loop
        for expr in self.doc.lookup_tables:
            var = expr.component.get_variable_by_name(expr.var)
            key = (expr.min, expr.max, expr.step, var.get_source_variable(recurse=True))
            idx = self.doc.lookup_table_indexes[key]
            if only_index is not None and only_index != idx:
                continue
            min, max, step, _ = self.lut_parameters(key)
            j = expr.table_name
            self.writeln('for (unsigned i=0 ; i<_table_size_', idx, '; i++)')
            self.open_block()
            self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(var), self.EQ_ASSIGN, min,
                         ' + i*', step, self.STMT_END)
            self.writeln(self.lut_access_code(idx, j, 'i'), self.EQ_ASSIGN, nl=False)
            self.output_expr(expr, False)
            self.writeln(self.STMT_END, indent=False)
            self.close_block()
        self.use_lookup_tables = True

    def output_lut_deletion(self, only_index=None):
        """Output code to delete memory allocated for lookup tables."""
        for idx in self.doc.lookup_table_indexes.itervalues():
            if only_index is None or only_index == idx:
                self.writeln('if (_lookup_table_', idx, ')')
                self.open_block()
                self.writeln('delete[] _lookup_table_', idx, self.STMT_END)
                self.writeln('_lookup_table_', idx, self.EQ_ASSIGN, 'NULL', self.STMT_END)
                self.close_block(blank_line=False)

    def output_lut_declarations(self):
        """Output declarations for the lookup tables."""
        self.output_comment('Lookup tables')
        # Allocate memory, per index variable for cache efficiency
        for idx in self.doc.lookup_table_indexes.itervalues():
            num_tables = unicode(self.doc.lookup_tables_num_per_index[idx])
            self.writeln(self.TYPE_DOUBLE, '(*_lookup_table_', idx, ')[', num_tables, ']', self.STMT_END)
        self.writeln()

    def output_lut_index_declarations(self, idx):
        """Output declarations the variables used to index this table."""
        self.writeln('unsigned _table_index_', idx, self.STMT_END)
        factor = self.lut_factor(idx, include_type=True)
        if factor:
            self.writeln(factor, self.STMT_END)
        if self.row_lookup_method:
            self.writeln('double* _lt_', idx, '_row', self.STMT_END)

    def output_lut_indices(self):
        """Output declarations for the lookup table indices."""
        self.output_comment('Lookup table indices')
        for idx in self.doc.lookup_table_indexes.itervalues():
            self.output_lut_index_declarations(idx)
        self.writeln()
    
    def output_lut_methods(self):
        """Output the methods which look up values from lookup tables."""
        if self.row_lookup_method:
            self.output_lut_row_lookup_methods()
            return
        self.output_comment('Methods to look up values from lookup tables')
        self.output_comment('using ', self.config.options.lookup_type)
        for expr in self.doc.lookup_tables:
            j = expr.table_name
            idx = expr.table_index
            self.writeln('inline double _lookup_', j, '(unsigned i',
                         self.lut_factor('', include_type=True, include_comma=True), ')')
            self.open_block()
            self.output_single_lookup(idx, j, 'return ')
            self.close_block()
        self.writeln()
        return

    def output_single_lookup(self, tidx, tname, result):
        """Write the lookup calculation for a single entry.
        
        Used by output_lut_row_lookup_methods and output_lut_methods.
        """
        self.writeln(self.TYPE_CONST_DOUBLE, 'y1', self.EQ_ASSIGN,
                     self.lut_access_code(tidx, tname, 'i'), self.STMT_END)
        if self.config.options.lookup_type == 'linear-interpolation':
            self.writeln(self.TYPE_CONST_DOUBLE, 'y2', self.EQ_ASSIGN,
                         self.lut_access_code(tidx, tname, 'i+1'), self.STMT_END)
            self.writeln(result, 'y1 + (y2-y1)*', self.lut_factor(''), self.STMT_END)
        else:
            self.writeln(result, 'y1', self.STMT_END)

    def output_lut_row_lookup_methods(self):
        """Write methods that return a whole row of a lookup table.

        Note: assumes that table names are numbered sequentially from 0.
        """
        self.output_comment('Row lookup methods')
        self.output_comment('using ', self.config.options.lookup_type)
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            num_tables = unicode(self.doc.lookup_tables_num_per_index[idx])
            self.writeln('double* _lookup_', idx, '_row(unsigned i',
                         self.lut_factor('', include_type=True, include_comma=True), ')')
            self.open_block()
            self.writeln('for (unsigned j=0; j<', num_tables, '; j++)')
            self.open_block()
            self.output_single_lookup(idx, 'j', '_lookup_table_%s_row[j] = ' % idx)
            self.close_block(False)
            self.writeln('return _lookup_table_', idx, '_row;')
            self.close_block()
        self.writeln()
        return
    
    def output_lut_row_lookup_memory(self):
        """Output declarations for the memory used by the row lookup methods."""
        self.output_comment('Row lookup methods memory')
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            min, max, step, var = key
            num_tables = unicode(self.doc.lookup_tables_num_per_index[idx])
            self.writeln('double _lookup_table_', idx, '_row[', num_tables, '];')
        self.writeln()
        return

    def is_lookup_table(self, expr):
        """Return True iff expr can be replaced by a lookup table.

        Uses annotations from a previous analysis."""
        return expr.getAttributeNS(NSS['lut'], u'possible', '') == u'yes'
    
    def contained_table_indices(self, node):
        """Return all lookup tables used directly in computing this node.

        If this is an expression node, checks all its children for table
        lookups, and returns the set of table indices used.
        """
        result = set()
        if isinstance(node, amara.bindery.element_base):
            if self.is_lookup_table(node):
                result.add(node.table_index)
            else:
                for child in node.xml_children:
                    result.update(self.contained_table_indices(child))
        return result
    
    def lut_factor(self, idx, include_comma=False, include_type=False):
        """Return code for any extra factor needed to do a table lookup.
        
        Will return the empty string unless linear interpolation is being used.
        """
        if self.config.options.lookup_type == 'linear-interpolation':
            factor = '_factor_' + str(idx)
            if include_type: factor = self.TYPE_DOUBLE + factor
            if include_comma: factor = ', ' + factor
        else:
            factor = ''
        return factor

    def output_table_lookup(self, expr, paren):
        """Output code to look up expr in the appropriate table."""
        i = expr.table_index
        if self.row_lookup_method:
            self.write('_lt_', i, '_row[', expr.table_name, ']')
        else:
            self.write(self.lookup_method_prefix, '_lookup_', expr.table_name,
                       '(_table_index_', i, self.lut_factor(i, include_comma=True), ')')
        return

    def output_table_index_generation(self, time_name, nodeset=set()):
        """Output code to calculate indexes into any lookup tables.
        
        If time_name is given and table bounds are being checked, the time value will be included in the
        error message.  Note that we need to pass it in, since in some contexts the free variable is not
        defined.
        
        If nodeset is given, then filter the table indices calculated so
        that only those needed to compute the nodes in nodeset are defined.
        
        A nodeset is required if any table indices are computed variables rather than state variables.
        In this case, we use the equations within nodeset to calculate the values of the indices, and
        return a set containing just those nodes used, so we can avoid recalculating them later.
        """
        tables_to_index = set()
        nodes_used = set()
        for node in nodeset:
            tables_to_index.update(self.contained_table_indices(node))
        if tables_to_index or not nodeset:
            self.output_comment('Lookup table indexing')
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            if not nodeset or idx in tables_to_index:
                var = key[-1]
                if var.get_type() is VarTypes.Computed:
                    if not nodeset:
                        raise TranslationError('Unable to generate lookup table indexed on', var, 'as it is a computed variable')
                    var_nodes = self.calculate_extended_dependencies([var]) & nodeset
                    self.output_equations(var_nodes)
                    nodes_used.update(var_nodes)
                self.output_table_index_checking(key, idx)
                if self.config.options.check_lt_bounds:
                    self.writeln('#define COVERAGE_IGNORE', indent=False)
                    self.writeln('if (_oob_', idx, ')')
                    if time_name is None:
                        dump_state_args = 'rY'
                    else:
                        dump_state_args = 'rY, ' + time_name
                    self.writeln('EXCEPTION(DumpState("', self.var_display_name(key[-1]),
                                 ' outside lookup table range", ', dump_state_args,'));', indent_offset=1)
                    self.writeln('#undef COVERAGE_IGNORE', indent=False)
                self.output_table_index_generation_code(key, idx)
        self.writeln()
        return nodes_used
        
    def output_table_index_checking(self, key, idx):
        """Check whether a table index is out of bounds."""
        if self.config.options.check_lt_bounds:
            var = key[-1]
            min, max, _, _ = self.lut_parameters(key)
            varname = self.code_name(var)
            self.writeln('bool _oob_', idx, self.EQ_ASSIGN, 'false', self.STMT_END)
            self.writeln('if (', varname, '>', max, ' || ', varname, '<', min, ')')
            self.open_block()
            self.writeln('#define COVERAGE_IGNORE', indent=False)
            if self.constrain_table_indices:
                self.writeln('if (', varname, '>', max, ') ', varname, self.EQ_ASSIGN, max, self.STMT_END)
                self.writeln('else ', varname, self.EQ_ASSIGN, min, self.STMT_END)
            else:
                self.writeln('_oob_', idx, self.EQ_ASSIGN, 'true', self.STMT_END)
            self.writeln('#undef COVERAGE_IGNORE', indent=False)
            self.close_block(blank_line=False)
    
    def output_table_index_generation_code(self, key, idx):
        """Method called by output_table_index_generation to output the code for a single table."""
        index_type = 'const unsigned '
        factor_type = 'const double '
        row_type = 'const double* const '
        var = key[-1]
        min, max, _, step_inverse = self.lut_parameters(key)
        offset = '_offset_' + idx
        offset_over_step = offset + '_over_table_step'
        varname = self.code_name(var)
        self.writeln(self.TYPE_CONST_DOUBLE, offset, self.EQ_ASSIGN, varname, ' - ', min, self.STMT_END)
        self.writeln(self.TYPE_CONST_DOUBLE, offset_over_step, self.EQ_ASSIGN,
                     offset, ' * ', step_inverse, self.STMT_END)
        idx_var = '_table_index_' + str(idx)
        if self.config.options.lookup_type == 'nearest-neighbour':
            if self.lt_index_uses_floor:
                self.writeln(index_type, idx_var, ' = (unsigned) round(', offset_over_step, ');')
            else:
                self.writeln(index_type, idx_var, ' = (unsigned) (', offset_over_step, '+0.5);')
        else:
            if self.lt_index_uses_floor:
                self.writeln(index_type, idx_var, ' = (unsigned) floor(', offset_over_step, ');')
            else:
                self.writeln(index_type, idx_var, ' = (unsigned)(', offset_over_step, ');')
            factor = self.lut_factor(idx)
            if factor:
                self.writeln(factor_type, factor, ' = ', offset_over_step, ' - ', idx_var, self.STMT_END)
        if self.row_lookup_method:
            self.writeln(row_type, '_lt_', idx, '_row = ', self.lookup_method_prefix, '_lookup_', idx,
                         '_row(', idx_var, self.lut_factor(idx, include_comma=True), ');')




###############################################
# Register translation classes in this module #
###############################################


class SolverInfo(object):
    """Add information for specialised translator classes into a model."""
    def __init__(self, model, force=False):
        """Add information for the solvers as XML.

        The Jacobian and linearity analyses store their results in
        Python data structures as attributes of this object.
        Transcribe these into XML in a child <solver_info> element.

        If any of these elements exist in the model they will be left
        unaltered, unless force is set to True.
        
        This constructor just sets up the container element; call one
        of the add_* methods to actually add the information to it.
        """
        self._model = model
        if force and hasattr(model, u'solver_info'):
            model.xml_remove_child(model.solver_info)
        if hasattr(model, u'solver_info'):
            solver_info = model.solver_info
        else:
            solver_info = model.xml_create_element(u'solver_info', NSS[u'solver'])
            model.xml_append(solver_info)
        self._solver_info = solver_info
        self._component = None
        self._dt = None
    
    def add_all_info(self):
        """Actually add the info."""
        self.add_transmembrane_potential_name()
        self.add_membrane_ionic_current()
        self.add_linearised_odes()
        self.add_jacobian_matrix()
        self.add_dt_reference()
    
    def add_dt_reference(self):
        """Add a reference to the variable representing dt."""
        solver_info = self._solver_info
        model = self._model
        if not hasattr(solver_info, u'dt'):
            dt = self.get_dt()
            elt = model.xml_create_element(u'dt', NSS[u'solver'], content=dt.fullname(cellml=True))
            solver_info.xml_append(elt)
            self._model._add_sorted_assignment(dt)
    
    def add_transmembrane_potential_name(self):
        """The name of the transmembrane potential."""
        solver_info = self._solver_info
        model = self._model
        if not hasattr(solver_info, u'transmembrane_potential'):
            v_elt = model.xml_create_element(
                u'transmembrane_potential', NSS[u'solver'],
                content=model._cml_transmembrane_potential.fullname())
            solver_info.xml_append(v_elt)
    
    def add_linearised_odes(self):
        """Linearised ODEs - where du/dt = g + hu (and g, h are not functions of u).
        
        Structure looks like:
        <linear_odes>
            <math>
                <apply><eq/>
                    <apply><diff/>
                        <bvar><ci>t</ci></bvar>
                        <ci>u</ci>
                    </apply>
                    <apply><plus/>
                        g
                        <apply><times/>
                            h
                            <ci>u</ci>
                        </apply>
                    </apply>
                </apply>
                .
                .
                .
            </math>
        </linear_odes>
        """
        solver_info = self._solver_info
        model = self._model
        if not hasattr(solver_info, u'linear_odes') and model._cml_linear_update_exprs:
            odes_elt = model.xml_create_element(u'linear_odes', NSS[u'solver'])
            solver_info.xml_append(odes_elt)
            odes_math = model.xml_create_element(u'math', NSS[u'm'])
            odes_elt.xml_append(odes_math)
            linear_vars = model._cml_linear_update_exprs.keys()
            linear_vars.sort(key=lambda v: v.fullname())
            free_var = model._cml_free_var
            for var in linear_vars:
                g, h = model._cml_linear_update_exprs[var]
                hu = mathml_apply.create_new(model, u'times', [h, var.fullname()])
                rhs = mathml_apply.create_new(model, u'plus', [g, hu])
                odes_math.xml_append(mathml_diff.create_new(
                    model, free_var.fullname(), var.fullname(), rhs))
            # Ensure that the model has a special component
            self._get_special_component()
    
    def _fix_jac_var_name(self, vname):
        """
        If PE will be performed on a model with a single component, then we'll need full names in
        the variable attributes.
        """
        if vname[:4] == 'var_' and len(self._model.component) == 1 and not self._model.component.ignore_component_name:
            name = unicode('var_' + self._model.component.name + '__' + vname[4:])
        else:
            name = unicode(vname)
        return name
    
    def add_jacobian_matrix(self):
        """Jacobian matrix elements.
        
        Structure looks like:
        <jacobian>
            [<math> assignments of common sub-terms </math>]
            <entry var_i='varname' var_j='varname'>
                <math> apply|cn|ci ...</math>
            </entry>
        </jacobian>
        """
        solver_info = self._solver_info
        model = self._model
        if model._cml_jacobian and model._cml_jacobian_full:
            jac = model._cml_jacobian[1]
        else:
            # Old-style partial jacobian, or no jacobian
            jac = model._cml_jacobian
        if not hasattr(solver_info, u'jacobian') and jac:
            jac_elt = model.xml_create_element(u'jacobian', NSS[u'solver'])
            solver_info.xml_append(jac_elt)
            
            if model._cml_jacobian_full:
                # There may be temporaries
                temporaries = model._cml_jacobian[0]
                if temporaries:
                    jac_elt.xml_append(amara_parse_cellml(temporaries).math)

            jac_vars = jac.keys()
            jac_vars.sort() # Will sort by variable name
            for v_i, v_j in jac_vars:
                # Add (i,j)-th entry
                attrs = {u'var_i': self._fix_jac_var_name(v_i),
                         u'var_j': self._fix_jac_var_name(v_j)}
                entry = model.xml_create_element(u'entry', NSS[u'solver'], attributes=attrs)
                jac_elt.xml_append(entry)
                entry_doc = amara_parse_cellml(jac[(v_i, v_j)].xml())
                entry.xml_append(entry_doc.math)
            # Ensure that the model has a special component
            self._get_special_component()
        return
    
    def use_canonical_variable_names(self):
        """
        PE has just been performed, so we need to update variable names occurring outside
        the modifiable mathematics sections.
        """
        jac_elt = getattr(self._solver_info, u'jacobian', None)
        for entry in getattr(jac_elt, u'entry', []):
            for vlabel in ['var_i', 'var_j']:
                vname = getattr(entry, vlabel)
                var = self._get_variable(vname)
                new_name = var.get_source_variable(recurse=True).fullname()
                setattr(entry, vlabel, new_name)
        dt_elt = getattr(self._solver_info, u'dt', None)
        if dt_elt:
            var = self._get_variable(unicode(dt_elt))
            new_name = var.get_source_variable(recurse=True).fullname()
            dt_elt.xml_remove_child(unicode(dt_elt))
            dt_elt.xml_append(unicode(new_name))

    def add_membrane_ionic_current(self):
        """Add ionic current information as XML for solvers to use."""
        solver_info = self._solver_info
        model = self._model
        # The total ionic current.  This relies on having a configuration store.
        if hasattr(model.xml_parent, '_cml_config') and not hasattr(solver_info, u'ionic_current'):
            conf = model.xml_parent._cml_config
            if conf.i_ionic_vars:
                ionic_elt = model.xml_create_element(u'ionic_current', NSS[u'solver'])
                # Adds each ionic var to the xml doc from the config store
                for var in conf.i_ionic_vars:
                    varelt = model.xml_create_element(u'var', NSS[u'solver'],
                                                      content=var.fullname())
                    ionic_elt.xml_append(varelt)
                solver_info.xml_append(ionic_elt)
        return
    
    def add_linear_ode_update_equations(self):
        """Add the update equations for the linear ODEs.
        
        A linear ODE has the form du/dt = g+h.u where g & h are not functions of u.  The
        update expression then looks like u = (u + g.dt)/(1 - h.dt).
        
        This replaces the linear_odes block with the structure:
        <linear_odes>
            <math>
                <ci>u</ci>
                <ci>t</ci>
                <apply> <!-- (u + g.dt)/(1 - h.dt) --> </apply>
            </math>
            .
            .
            .
        </linear_odes>
        """
        block = getattr(self._solver_info, u'linear_odes', None)
        dt = self._model.get_config().dt_variable.fullname() # was dt = u'delta_t'
        # Add the new equations
        for u, t, gh in self.get_linearised_odes():
            g, h = gh
            g.safe_remove_child(g, g.xml_parent)
            g_dt = mathml_apply.create_new(block, u'times', [g, dt])
            numer = mathml_apply.create_new(block, u'plus', [u.fullname(), g_dt])
            h.safe_remove_child(h, h.xml_parent)
            h_dt = mathml_apply.create_new(block, u'times', [h, dt])
            denom = mathml_apply.create_new(block, u'minus', [(u'1', u'dimensionless'), h_dt])
            eqn = mathml_apply.create_new(block, u'divide', [numer, denom])
            math = block.xml_create_element(u'math', NSS[u'm'])
            math.xml_append(mathml_ci.create_new(block, u.fullname()))
            math.xml_append(mathml_ci.create_new(block, t.fullname()))
            math.xml_append(eqn)
            block.xml_append(math)
            self._add_variable_links(math)
        # Remove the old equations (first math element)
        block.xml_remove_child(block.math)
    
    def add_variable_links(self):
        """Link ci elements in the added XML to cellml_variable objects.
        
        This analyses the names in the ci elements to determine which variable in
        the model they refer to.
        """
        self._process_mathematics(self._add_variable_links)
        #1795 - classify temporary variables for the Jacobian matrix, and append
        # to the main list of assignments in the model
        solver_info = self._solver_info
        if hasattr(solver_info, u'jacobian') and hasattr(solver_info.jacobian, u'math'):
            for elt in solver_info.jacobian.math.apply:
                elt.classify_variables(root=True)
            for elt in solver_info.jacobian.math.apply:
                self._model.topological_sort(elt)
        #2418 - check if any state variables have been units-converted
        self._check_state_var_units_conversions()
    
    def _check_state_var_units_conversions(self):
        """Check if any Jacobian entries need to be altered because the units of state variables have changed.
        
        If any variable considered a state variable by the Jacobian is now of type Computed then it has been
        converted.  We figure out the conversion factor, update the Jacobian to reference the new state variable,
        and units-convert the derivative.
        """
        if not hasattr(self._solver_info, u'jacobian'):
            return
        # Helper methods
        def set_var_values(elt, vars=None):
            """Fake all variables appearing in the given expression being set to 1.0, and return them."""
            if vars is None:
                vars = []
            if isinstance(elt, mathml_ci):
                elt.variable.set_value(1.0)
                vars.append(elt.variable)
            else:
                for child in getattr(elt, 'xml_children', []):
                    set_var_values(child, vars)
            return vars
        # Find any converted state variables
        converted_state_vars = set()
        for entry in getattr(self._solver_info.jacobian, u'entry', []):
            var = self._get_variable(entry.var_i)
            if var.get_type() == VarTypes.Computed:
                converted_state_vars.add(var)
        if not converted_state_vars:
            return
        # Figure out the conversion factor in each case
        state_var_map = {}
        for var in converted_state_vars:
            defn = var.get_dependencies()[0]
            defn_vars = set_var_values(defn.eq.rhs)
            assert len(defn_vars) == 1, "Unexpected form of units conversion expression found"
            factor = defn.eq.rhs.evaluate()
            state_var_map[var] = (defn_vars[0], factor)
            defn_vars[0].unset_values()
        # Apply the conversion to relevant Jacobian entries
        for entry in getattr(self._solver_info.jacobian, u'entry', []):
            factor = 1
            var_i = self._get_variable(entry.var_i)
            if var_i in converted_state_vars:
                var_i, factor_i = state_var_map[var_i]
                var_i = var_i.get_source_variable(recurse=True)
                entry.var_i = unicode(var_i.fullname())
                factor /= factor_i
            var_j = self._get_variable(entry.var_j)
            if var_j in converted_state_vars:
                var_j, factor_j = state_var_map[var_j]
                var_j = var_j.get_source_variable(recurse=True)
                entry.var_j = unicode(var_j.fullname())
                factor *= factor_j
            if factor != 1:
                # Replace rhs with rhs * factor
                rhs = list(entry.math.xml_element_children())[0]
                entry.math.safe_remove_child(rhs)
                new_rhs = mathml_apply.create_new(entry, 'times', [(factor, 'dimensionless'), rhs])
                entry.math.xml_append(new_rhs)
    
    def do_binding_time_analysis(self):
        """Do a binding time analysis on the additional mathematics.
        
        This requires self.add_variable_links to have been called already.
        """
        self._process_mathematics(lambda elt: elt._get_binding_time())
        
    def _process_mathematics(self, func):
        """Apply func to each top-level mathematical construct in the solver info blocks.
        
        func must be able to accept mathml_piecewise, mathml_apply, mathml_ci and mathml_cn elements.
        """
        solver_info = self._solver_info
        # Jacobian
        if hasattr(solver_info, u'jacobian'):
            if hasattr(solver_info.jacobian, u'math'):
                for elt in solver_info.jacobian.math.apply:
                    func(elt)
            for entry in solver_info.jacobian.entry:
                for elt in entry.math.xml_element_children():
                    func(elt)
        # Linearised ODEs
        if hasattr(solver_info, u'linear_odes'):
            for math in solver_info.linear_odes.math:
                for elt in math.xml_element_children():
                    func(elt)
    
    def has_modifiable_mathematics(self):
        """Check if the solver info blocks contain any modifiable mathematics."""
        try:
            self.get_modifiable_mathematics().next()
            return True
        except StopIteration:
            return False
    
    def get_modifiable_mathematics(self):
        """Get an iterable over mathematical constructs in the solver info blocks that can be changed.
        
        Returned elements will be mathml_piecewise, mathml_apply, mathml_ci or mathml_cn instances.
        """
        solver_info = self._solver_info
        # Jacobian - entry definitions and temporaries can be changed
        if hasattr(solver_info, u'jacobian'):
            if hasattr(solver_info.jacobian, u'math'):
                for elt in solver_info.jacobian.math.apply:
                    yield elt
            for entry in solver_info.jacobian.entry:
                for elt in entry.math.xml_element_children():
                    yield elt
        # Linearised ODEs - only g & h can be changed
        if hasattr(solver_info, u'linear_odes'):
            for _, _, eqns in self.get_linearised_odes():
                for eqn in eqns:
                    yield eqn

    def get_linearised_odes(self):
        """Return an iterable over the linearised ODEs, i.e. ODEs of the form
        du/dt = g + hu (with g, h not functions of u).
        
        Yields tuples (u, t, eqns) where the form of eqns depends on whether
        add_linear_ode_update_equations has been called.  If so, it is a 1-tuple
        containing the update equation; if not, it is (g,h).
        """
        if hasattr(self._solver_info, u'linear_odes'):
            if hasattr(self._solver_info.linear_odes.math, u'ci'):
                for math in self._solver_info.linear_odes.math:
                    u, t, eqn = list(math.xml_element_children())
                    u = u.variable
                    t = t.variable
                    yield (u, t, (eqn,))
            else:
                for ode in self._solver_info.linear_odes.math.apply:
                    u = ode.apply.ci.variable
                    t = ode.apply.bvar.ci.variable
                    opers = ode.apply[1].operands()
                    g = opers.next()
                    h = opers.next().operands().next()
                    yield (u, t, (g,h))
    
    def _add_variable_links(self, elt):
        """Recursively link ci elements in the given XML tree to cellml_variable objects.
        
        Also sets component links: for ci elements, to the component containing the linked
        variable, and for cn elements, to the first component in the model.
        """
        if isinstance(elt, mathml_ci):
            var = self._get_variable(unicode(elt))
            elt._cml_variable = var
            elt._cml_component = var.component
        elif isinstance(elt, mathml_cn):
            # Fake a component, since it doesn't really have one
            elt._cml_component = elt.model.component
        elif hasattr(elt, 'xml_children'):
            for child in elt.xml_children:
                self._add_variable_links(child)

    _jac_temp_re = re.compile(r't[0-9]+')
    def _get_variable(self, varname):
        """Return the variable in the model with name varname."""
        try:
            if varname == 'delta_t':
                # Special case for the timestep in ComputeJacobian and elsewhere
                var = self.get_dt()
            elif self._jac_temp_re.match(varname):
                var = self._get_special_variable(varname, VarTypes.Unknown)
            else:
                var = cellml_variable.get_variable_object(self._model, varname)
        except KeyError:
            raise ValueError("Cannot find variable '%s' referenced in SolverInfo" % varname)
        return var
    
    def create_dt(self, modifier, comp, units):
        """Create the special 'dt' variable in the given component."""
        self._dt = modifier.add_variable(comp, modifier._uniquify_var_name(u'dt', comp), units)
        self._dt._set_type(VarTypes.Free)
        return self._dt
    
    def get_dt(self):
        """Get or create a special 'dt' variable."""
        if not self._dt:
            self._dt = self._get_special_variable(u'dt', VarTypes.Free)
        return self._dt

    def _get_special_variable(self, varname, ptype=VarTypes.Unknown):
        """Get or create a special variable object that doesn't really exist in the model."""
        comp = self._get_special_component()
        try:
            var = comp.get_variable_by_name(varname)
        except KeyError:
            var = cellml_variable.create_new(self._model, varname, u'dimensionless')
            comp._add_variable(var)
            var._set_type(ptype)
        return var

    def _get_special_component(self):
        """Get or create a special component for containing special variables."""
        if not self._component:
            self._component = cellml_component.create_new(self._model, u'')
            self._model._add_component(self._component, special=True)
        return self._component



class ConfigurationStore(object):
    """
    A container for configuration information, read in from XML
    configuration files.  The file structure is described in the
    read_configuration_file method.
    """
    def __init__(self, doc, options=None):
        """Create a new store.

        doc specifies a CellML document, the processing of which this configuration store will affect.

        If given, options should be an optparse.Values instance containing command-line options.
        """
        self.doc = doc
        doc._cml_config = self
        self.options = options
        self.unit_definitions = cellml_component.create_new(doc.model, '*lookup_table_units*')
        self.unit_definitions.xml_parent = doc.model # Needed for looking up standard units
        # Transmembrane potential
        self.V_definitions = [u'membrane,V']
        self.V_variable = None
        # Membrane capacitance
        self.Cm_definitions = []
        self.Cm_variable = None
        # Lookup table configuration
        self.lut_config = {}
        # Ionic currents configuration
        self.i_stim_definitions = [u'membrane,i_Stim']
        self.i_stim_var = None
        self.i_ionic_definitions = [u'membrane,i_.*']
        self.i_ionic_vars = []
        # Whether GetIIonic will need to negate the sum of i_ionic_vars
        self.i_ionic_negated = False
        # Whether the stimulus magnitude is positive, rather than negative
        self.i_stim_negated = False
        # Other variables that may be set by other code, for example an InterfaceGenerator
        self.dt_variable = None
        self.i_data_clamp_current = None
        self.i_data_clamp_conductance = None
        return

    def read_configuration_file(self, config_file):
        """Read configuration stored in config_file.

        The configuration file is expected to be XML, conforming to
        the following structure.  Currently little checking is done on
        the file format; incorrectly formatted files are unlikely to
        give particularly helpful error messages.

        The root element may contain a 'global' element, giving global
        configuration options.  These include:

         * 'lookup_tables'
           Contains one or more 'lookup_table' elements, one for each
           type of lookup table available.  These contain (a selection of)
           the elements:
           * 'var' - the variable to key on.  The component name given
             should be that from which the variable is exported.  Must be
             present.
           * 'min', 'max', 'step' - table bounds parameters.  Optional.
           Default values are used for parameters that are not present.
         * 'currents'
           Defines which variables hold the ionic and stimulus currents,
           if any.  It contains 2 elements:
           * 'stimulus' - the full name of the stimulus current variable
           * 'ionic_match' - a regular expression matching full names of
             ionic current variables.  It may also match the stimulus
             current, but the stimulus will never be considered an ionic
             current.  The value is split on ','; the first part is then
             matched against component names, and the second against
             variables in matching components.
             
             This is mostly redundant now, because the equation for dV/dt
             is used first to determine the ionic currents (see documentation
             for _find_transmembrane_currents_from_voltage_ode), and only
             if this fails to find suitable currents will the ionic_match
             definition be used.
         * 'transmembrane_potential'
           Defines which variable holds the transmembrane potential.
           Defaults to 'membrane,V' if not present.
         * 'membrane_capacitance'
           Defines which variable holds the cell membrane capacitance.
           
        The root element also contains 0 or more 'for_model' elements,
        which contain settings for individual models.  These must have
        at least one of an 'id' or 'name' attribute, which specify the
        model in question.  They can also contain anything allowable as
        global configuration options.  Options given here override those
        specified globally.

        Configuration which is identical for groups of models may be given
        using the 'for_models' element.  This has an 'ids' element as its
        first child, which contains 'id' elements holding either the name
        or id of a model.  The remaining contents of the 'for_models'
        element are as for 'for_model'.

        There are 3 ways of specifying variables:
        1. By name (var type='name')
           Variable names are given in full form, i.e. component_name,variable_name
        2. By standardised name (var type='oxmeta')
           Use the name from the oxmeta annotations
        3. By reference to a section of the config file (when defining lookup table keys)
           e.g. <var type='config-name'>transmembrane_potential</var>

        Within any element that specifies a variable, a list of <var> elements can be
        provided.  Each will be tried in turn to see if a match can be found in the model,
        and the first match wins.

        Some items are overridden if oxmeta annotations are present in the model, with
        the annotated variable taking precedence over the config file specification.
        """
        DEBUG('config', "Reading configuration from ", config_file)
        binder = amara.bindery.binder()
        binder.set_binding_class(None, "units", cellml_units)
        binder.set_binding_class(None, "unit", cellml_unit)
        rules = [bt.ws_strip_element_rule(u'*')]
        config_doc = amara_parse(config_file, rules=rules, binderobj=binder)
        # Store extra units definitions
        for defn in config_doc.xml_xpath(u'/*/units'):
            defn.xml_parent = self.unit_definitions # Needed for looking up the units this definition is derived from
            self.unit_definitions.add_units(defn.name, defn)
        # Overrides for command-line options
        if self.options and hasattr(config_doc.pycml_config, 'command_line_args'):
            args = map(str, config_doc.pycml_config.command_line_args.arg)
            args.append('dummy-file')
            get_options(args, self.options)
        # Sections to use in configuration; later sections take precedence
        sections = []
        # Use global configuration?
        glo = config_doc.xml_xpath(u'/*/global')
        if glo:
            sections.append(glo[0])
        # Get the config section(s) for our model.  Sections
        # specifically for this model come after sections covering
        # multiple models, so they take precedence.
        model_id = getattr(self.doc.model, u'id', self.doc.model.name)
        sections.extend(config_doc.xml_xpath(
            u'/*/for_models[ids/id="%s" or ids/id="%s"]'
            % (self.doc.model.name, model_id)))
        sections.extend(config_doc.xml_xpath(
            u'/*/for_model[@name="%s" or @id="%s"]'
            % (self.doc.model.name, model_id)))
        # Main items of configuration
        for section in sections:
            if hasattr(section, u'lookup_tables'):
                self._parse_lookup_tables(section.lookup_tables)
            if hasattr(section, u'currents'):
                self._parse_currents(section.currents)
            if hasattr(section, u'transmembrane_potential'):
                self._parse_Vm(section.transmembrane_potential)
            if hasattr(section, u'membrane_capacitance'):
                self._parse_Cm(section.membrane_capacitance)
    
    def finalize_config(self):
        """Having read all the configuration files, apply to the model."""
        # If no LT options given, add defaults
        if not self.lut_config:
            config_key = ('config-name', 'transmembrane_potential')
            self.lut_config[config_key] = {}
            self._set_lut_defaults(self.lut_config[config_key])
        # Identify the variables in the model
        self.find_transmembrane_potential()
        self.find_membrane_capacitance()
        if not self.options.protocol:
            self.find_current_vars()

    def _create_var_def(self, content, defn_type):
        """Create a variable definition object."""
        xml_fragment = '<var type="%s">%s</var>' % (defn_type, content)
        return amara.parse(str(xml_fragment)).var

    def _check_var_def(self, var_elt, var_desc):
        """Check a variable definition is syntactically valid.
        
        If type == 'name', it must have text content of the form "component_name,variable_name".
        If type == 'oxmeta', it must have text content drawn from METADATA_NAMES.
        If type == 'config-name', it must have text content either 'stimulus' or 'transmembrane_potential'.
        """
        defn_type = getattr(var_elt, u'type', u'name')
        if defn_type == u'name':
            name_parts = unicode(var_elt).strip().split(',')
            if len(name_parts) != 2:
                raise ConfigurationError('Invalid definition of ' + var_desc + ': '
                                         + unicode(var_elt))
        elif defn_type == u'oxmeta':
            if unicode(var_elt) not in cellml_metadata.METADATA_NAMES:
                raise ConfigurationError('"' + unicode(var_elt) + '" is not a valid oxmeta name')
        elif defn_type == u'config-name':
            if unicode(var_elt) not in [u'stimulus', u'transmembrane_potential', u'membrane_capacitance']:
                raise ConfigurationError('"' + unicode(var_elt) + '" is not a name known to the config file')
        else:
            raise ConfigurationError('"' + defn_type + '" is not a valid variable definition type')
        return

    def _parse_var(self, elt, name):
        """Parse definition of a special variable."""
        if hasattr(elt, 'var'):
            # List of possibilities
            defs = []
            for vardef in elt.var:
                self._check_var_def(vardef, name)
                defs.append(vardef)
        else:
            # Old style - single variable given by text content
            self._check_var_def(elt, name)
            defs = [elt]
        return defs

    def _parse_Vm(self, vm_elt):
        """Parse definition of variable holding the transmembrane potential."""
        self.V_definitions = self._parse_var(vm_elt, 'transmembrane_potential')
    
    def _parse_Cm(self, cm_elt):
        """Parse definition of variable holding the cell membrane capacitance."""
        self.Cm_definitions = self._parse_var(cm_elt, 'membrane_capacitance')

    def _parse_currents(self, currents):
        """Parse definitions of ionic and stimulus currents."""
        if hasattr(currents, u'stimulus'):
            self.i_stim_definitions = self._parse_var(currents.stimulus, 'stimulus current')
        if hasattr(currents, u'ionic_match'):
            self.i_ionic_definitions = self._parse_var(currents.ionic_match, 'ionic currents')
        return
    
    def _find_variable(self, defn, pe_done=False):
        """Find a variable matching the given definition.

        If pe_done is True, then partial evaluation has been performed
        on the model, so looking for variables by name needs to look for
        variables called compname__varname in the single component.
        """
        defn_type = getattr(defn, u'type', u'name')
        if defn_type == u'name':
            name_parts = unicode(defn).strip().split(',')
            if pe_done:
                try:
                    var = self.doc.model.component.get_variable_by_name(u'__'.join(name_parts))
                except KeyError:
                    var = None
            else:
                var = self.doc.model.xml_xpath(u'cml:component[@name="%s"]/cml:variable[@name="%s"]'
                                               % tuple(name_parts))
                if var:
                    var = var[0]
        elif defn_type == u'oxmeta':
            var = self.doc.model.get_variable_by_oxmeta_name(str(defn), throw=False)
        elif defn_type == u'config-name':
            if unicode(defn) == u'stimulus':
                var = self.i_stim_var
            elif unicode(defn) == u'transmembrane_potential':
                var = self.V_variable
            elif unicode(defn) == u'membrane_capacitance':
                var = self.Cm_variable
            else:
                raise ConfigurationError('"' + str(defn) + '" is not a valid configuration file variable name')
        else:
            raise ConfigurationError('"' + defn_type + '" is not a valid variable definition type')
        return var
    
    def _process_ci_elts(self, elt, func, **kwargs):
        """Recursively apply func to any ci elements in the tree rooted at elt."""
        if isinstance(elt, mathml_ci):
            func(elt, **kwargs)
        else:
            for child in getattr(elt, 'xml_children', []):
                self._process_ci_elts(child, func, **kwargs)
    
    def _find_transmembrane_currents_from_voltage_ode(self):
        """Analyse the expression for dV/dt to determine the transmembrane currents.
        
        Looks for an equation defining dV/d(something) and assumes the something is
        time; this will be checked during code generation for Chaste.  It then finds
        all variables on the RHS of this equation which have the same units as the
        stimulus current (self.i_stim_var) and identifies these as transmembrane
        currents.  Will automatically exclude the stimulus current.
        
        If self.V_variable is not set, returns the empty list.
        """
        if not self.V_variable:
            DEBUG('config', "Transmembrane potential not configured, so can't determine currents from its ODE")
            return []
        if self.i_stim_var:
            current_units = [self.i_stim_var.component.get_units_by_name(self.i_stim_var.units)]
        else:
            from CellMLToNektarTranslator import CellMLToNektarTranslator
            current_units = CellMLToNektarTranslator.get_current_units_options(self.doc.model)
        ionic_vars = []
        
        def find_units_match(test_units, units_list, remove_match=False, keep_only_match=False):
            """Look for a units definition dimensionally equivalent to test_units within units_list.
            
            If remove_match is True, remove the first match from the list.
            If keep_only_match is True, remove all entries except the first match from the list.
            Return the matching units, or None if there are no matches.
            """
            for units in units_list:
                if test_units.dimensionally_equivalent(units):
                    match = units
                    break
            else:
                match = None
            if match and remove_match:
                units_list.remove(match)
            if match and keep_only_match:
                units_list[:] = [match]
            return match

        def clear_values(expr, process_definitions=False):
            """Recursively clear saved values for variables in this expression.
            
            If process_definitions is True, recursively treat expressions defining variables
            used in this expression, too.
            """
            def process_var(var):
                var.unset_values()
                var._unset_binding_time(only_temporary=True)
                if process_definitions:
                    defn = var.get_dependencies()
                    if defn:
                        if isinstance(defn[0], mathml_apply):
                            clear_values(defn[0].eq.rhs, process_definitions=True)
                        elif isinstance(defn[0], cellml_variable):
                            process_var(defn[0])
            def process_ci(ci_elt):
                process_var(ci_elt.variable)
            self._process_ci_elts(expr, process_ci)
        
        def check_if_current(ci_elt, vars_found):
            """Check if this is a transmembrane current."""
            v = ci_elt.variable
            if v.get_source_variable(recurse=True) is not self.i_stim_var:
                vars_found.append(v)
                # Check units
                u = v.component.get_units_by_name(v.units)
                if find_units_match(u, current_units, keep_only_match=True):
                    ionic_vars.append(v.get_source_variable(recurse=True))
                    ionic_vars[-1]._cml_ref_in_dvdt = ci_elt # Hack for data clamp support (#2708)
            # Fake this variable being 1 so we can check the sign of GetIIonic
            if not v.is_statically_const(ignore_annotations=True):
                v.set_value(1.0)
        
        def bfs(func, vars, *args, **kwargs):
            """Do a breadth first search of the definitions of variables in vars.
            
            func is the recursive function to call.  It will be given the list of defining expressions
            as its first argument, and args and kwargs as remaining arguments.
            """
            def get_defn(var):
                defn = var.get_dependencies()
                if defn:
                    var._set_binding_time(BINDING_TIMES.static, temporary=True)
                    if isinstance(defn[0], cellml_variable):
                        defn = get_defn(defn[0])
                    else:
                        assert isinstance(defn[0], mathml_apply)
                        var.unset_values()
                        defn = defn[0].eq.rhs
                return defn
            defns = []
            for var in vars:
                defn = get_defn(var)
                if defn:
                    defns.append(defn)
            if defns:
                func(defns, *args, **kwargs)

        def find_currents(exprs, depth=0, maxdepth=2):
            """Find ionic currents by searching the given expressions.
            
            On the initial call, exprs should contain just the definition of dV/dt (i.e. the RHS).
            
            Uses breadth-first search of the equation dependency tree to find variables that
            have units dimensionally equivalent to one of the current formulations that Chaste
            can handle, or equivalent to the stimulus current's units if one is defined.
            
            Initially, A_per_F is removed from the list, since the RHS of dV/dt should always
            have equivalent dimensions.  If another option can't be found within maxdepth levels,
            we restart the search with A_per_F included.  The depth limit is intended to guard against
            unexpectedly finding something that isn't a current; it's somewhat dodgy, but won't
            break on any model I know, and I haven't thought of a better approach yet.
            
            When one variable with suitable units is found, further ionic currents must have units
            equivalent to its to be found.  Also once one ionic current is found, only the remaining
            expressions at its depth will be processed.
            """
            if depth == 0 and maxdepth > 0:
                dvdt_units = exprs[0].xml_parent.eq.lhs.get_units()
                A_per_F = find_units_match(dvdt_units, current_units, remove_match=True)
#                 # We could do this check, but actually it doesn't catch much and later checks will pick up on the problem
#                 if A_per_F is None and not self.i_stim_var:
#                     raise ConfigurationError('Units ' + dvdt_units.description() + ' of dV/dt are not equivalent to V/s - unable to continue.')
            # Process all expressions at this depth
            vars_found = []
            for expr in exprs:
                self._process_ci_elts(expr, check_if_current, vars_found=vars_found)
            if not ionic_vars and depth != maxdepth:
                # Process the definitions of expressions at this depth
                bfs(find_currents, vars_found, depth+1, maxdepth)
            # If we reached maxdepth unsuccessfully, try again with A_per_F included (if it was an option)
            if not ionic_vars and depth == 0 and maxdepth > 0 and A_per_F:
                current_units.append(A_per_F)
                find_currents(exprs, depth, maxdepth=-1)

        def assign_values_for_stimulus_check(exprs, found_stim=Sentinel()):
            """Assign temporary values to variables in order to check the stimulus sign.
            
            This will process defining expressions in a breadth first search until the stimulus
            current is found.  Each variable that doesn't have its definitions processed will
            be given a value as follows:
             - stimulus current = 1
             - other currents = 0
             - other variables = 1
            The stimulus current is then negated from the sign expected by Chaste if evaluating
            dV/dt gives a positive value.
            """
            assert len(current_units) == 1 # We are using the stimulus units
            vars = []
            def f(ci_elt):
                v = ci_elt.variable
                if v.get_source_variable(recurse=True) is self.i_stim_var:
                    v.set_value(1.0)
                    found_stim.set()
                else:
                    u = v.component.get_units_by_name(v.units)
                    if u.dimensionally_equivalent(current_units[0]):
                        v.set_value(0.0)
                    elif not v.is_statically_const(ignore_annotations=True):
                        v.set_value(1.0)
                    vars.append(v)
            for expr in exprs:
                self._process_ci_elts(expr, f)
            if not found_stim:
                bfs(assign_values_for_stimulus_check, vars, found_stim=found_stim)

        # Iterate over all expressions in the model, to find the one for dV/d(something)
        for expr in (e for e in self.doc.model.get_assignments() if isinstance(e, mathml_apply) and e.is_ode()):
            # Assume the independent variable is time; if it isn't, we'll catch this later
            (dep_var, time_var) = expr.assigned_variable()
            if dep_var.get_source_variable(recurse=True) is self.V_variable:
                # Recursively search for ionic currents
                find_currents([expr.eq.rhs])
                if not ionic_vars:
                    # The sign checks below will be nonsense in this case. An error will be raised later.
                    break
                # Check the sign of the RHS
                self.i_ionic_negated = expr.eq.rhs.evaluate() > 0.0
                clear_values(expr.eq.rhs, process_definitions=True)
                if self.i_stim_var:
                    # Check the sign of the stimulus current
                    assign_values_for_stimulus_check([expr.eq.rhs])
                    self.i_stim_negated = expr.eq.rhs.evaluate() > 0.0
                    clear_values(expr.eq.rhs, process_definitions=True)
                # Found dV/d(something); don't check any more expressions
                break
        DEBUG('config', "Found ionic currents from dV/dt: ", ionic_vars)
        call_if(self.i_ionic_negated, DEBUG, 'config', "Ionic current is negated")
        call_if(self.i_stim_negated, DEBUG, 'config', "Stimulus current is negated")
        return ionic_vars
    
    def _find_var(self, oxmeta_name, definitions):
        """Find the variable object in the model for a particular concept.
        
        Will look for a variable annotated with the given oxmeta_name first, then
        try the list of definitions from the configuration file in turn.
        """
        var = None
        # Prepend an oxmeta definition
        oxmeta_defn = self._create_var_def(oxmeta_name, 'oxmeta')
        for defn in [oxmeta_defn] + definitions:
            var = self._find_variable(defn)
            if var:
                break
        return var

    def find_current_vars(self):
        """Find the variables representing currents."""
        # Find the stimulus current, if it exists for this kind of model (some are self-excitatory)
        if not self.doc.model.is_self_excitatory():
            # self.i_stim_var = self._find_var('membrane_stimulus_current', self.i_stim_definitions)
            # DEBUG('config', 'Found stimulus', self.i_stim_var)
            if not self.i_stim_var:
                # No match :(
                msg = "No stimulus current found; you'll have trouble generating Nektar code"
                if self.options.fully_automatic:
                    raise ConfigurationError(msg)
                else:
                    print >>sys.stderr, msg
                    self.i_stim_var = None
        # For other ionic currents, try using the equation for dV/dt unless told otherwise
        if not self.options.use_i_ionic_regexp:
            self.i_ionic_vars = self._find_transmembrane_currents_from_voltage_ode()
        else:
            for defn in self.i_ionic_definitions:
                if getattr(defn, u'type', u'name') != u'name':
                    raise ConfigurationError('Ionic current definitions have to have type "name"')
                regexps = unicode(defn).strip().split(',')
                comp_re = re.compile(regexps[0] + '$')
                var_re = re.compile(regexps[1] + '$')
                for component in getattr(self.doc.model, u'component', []):
                    if comp_re.match(unicode(component.name).strip()):
                        for var in getattr(component, u'variable', []):
                            if (var is not self.i_stim_var and
                                var_re.match(unicode(var.name).strip())):
                                self.i_ionic_vars.append(var)
        if not self.i_ionic_vars:
            msg = "No ionic currents found; you'll have trouble generating Nektar code"
            if self.options.fully_automatic:
                raise ConfigurationError(msg)
            else:
                print >>sys.stderr, msg
        return

    def _parse_lookup_tables(self, lookup_tables):
        """Parse a lookup_tables element."""
        for lt in lookup_tables.lookup_table:
            var_type = getattr(lt.var, u'type', u'name')
            var_name = unicode(lt.var).strip()
            config_key = (var_type, var_name)
            if not config_key in self.lut_config:
                self.lut_config[config_key] = {}
                self._set_lut_defaults(self.lut_config[config_key])
            for elt in lt.xml_element_children():
                if elt.localName != u'var':
                    self.lut_config[config_key]['table_' + elt.localName] = unicode(elt).strip()
            if hasattr(lt, u'units'):
                try:
                    units = self.unit_definitions.get_units_by_name(lt.units)
                except KeyError:
                    raise ConfigurationError('The units "%s" referenced by the lookup table for "%s" do not exist'
                                             % (lt.units, var_name))
                self.lut_config[config_key]['table_units'] = units
        return

    def _set_lut_defaults(self, lut_dict):
        """Set default configuration for a lookup table."""
        def_dict = optimize.LookupTableAnalyser._LT_DEFAULTS
        for k, v in def_dict.iteritems():
            if k != 'table_var':
                lut_dict[k] = v
        lut_dict['table_units'] = None
        return

    def annotate_currents_for_pe(self):
        """Annotate ionic & stimulus current vars so PE doesn't remove them.
        Also annotate the membrane capacitance, if defined."""
        if self.i_stim_var:
            self.i_stim_var.set_pe_keep(True)
        for var in self.i_ionic_vars:
            var.set_pe_keep(True)
        if self.Cm_variable:
            self.Cm_variable.set_pe_keep(True)
        return
    
    def expose_variables(self):
        """Expose variables for access with GetAnyVariable if desired."""
        def annotate(var):
            t = var.get_type()
            if t == VarTypes.Constant:
                var.set_is_modifiable_parameter(True)
            elif t in [VarTypes.Computed, VarTypes.Free, VarTypes.Mapped]:
                var.set_is_derived_quantity(True)
        if self.options.expose_annotated_variables:
            for var in self.metadata_vars:
                if (not self.options.use_chaste_stimulus or
                    not var.oxmeta_name in cellml_metadata.STIMULUS_NAMES):
                    annotate(var)
            DEBUG('translate', "+++ Exposed annotated variables")
        if self.options.expose_all_variables:
            for var in self.doc.model.get_all_variables():
                annotate(var)
            DEBUG('translate', "+++ Exposed all variables")
    
    def annotate_metadata_for_pe(self):
        "Annotate all vars tagged with metadata so PE doesn't remove them."
        for var in self.metadata_vars:
            var.set_pe_keep(True)
        return

    def find_transmembrane_potential(self):
        """Find and store the variable object representing V.

        Tries metadata annotation first.  If that fails, uses the name given in
        the command line options, if present.  If that fails, uses the config file.
        """
        if not self.options:
            raise ValueError('No command line options given')
        # Check command line option before config file
        if self.options.transmembrane_potential:
            self.V_definitions[0:0] = [self.options.transmembrane_potential.strip().split(',')]
            if len(self.V_definitions[0]) != 2:
                raise ConfigurationError('The name of V must contain both component and variable name')
        self.V_variable = self._find_var('membrane_voltage', self.V_definitions)
        DEBUG('config', 'Found V', self.V_variable)
        if not self.V_variable and not self.options.protocol:
            raise ConfigurationError('No transmembrane potential found; check your configuration')
        return self.V_variable
    
    def find_membrane_capacitance(self):
        """Find and store the variable object representing the cell membrane capacitance.
        
        Uses first metadata, if present, then the configuration file."""
        self.Cm_variable = self._find_var('membrane_capacitance', self.Cm_definitions)
        DEBUG('config', 'Found capacitance', self.Cm_variable)

    def find_lookup_variables(self):
        """Find the variable objects used as lookup table keys.

        This method translates the variable names given in the configuration file into objects
        in the document, and then uses those objects as keys in our lut_config dictionary.
        The ultimate source variable for the variable specified is used, in order to avoid
        complications caused by intermediaries being removed (e.g. by PE).

        The table settings are also units-converted to match the units of the key variable.
        """
        new_config = {}
        for key in self.lut_config:
            defn_type, content = key
            defn = self._create_var_def(content, defn_type)
            var = self._find_variable(defn)
            if not var:
                # Variable doesn't exist, so we can't index on it
                LOG('lookup-tables', logging.WARNING, 'Variable', content, 'not found, so not using as table index.')
            else:
                var = var.get_source_variable(recurse=True)
                if not var in new_config:
                    new_config[var] = {}
                new_config[var].update(self.lut_config[key])
                # Apply units conversions to the table settings if required
                table_units = new_config[var]['table_units']
                if table_units:
                    var_units = var.get_units()
                    if not table_units.dimensionally_equivalent(var_units):
                        LOG('lookup-tables', logging.WARNING, 'Variable', content, 'is in units', var_units.description(),
                            'which are incompatible with', table_units.description(), 'so not using as table index.')
                    elif not table_units.equals(var_units):
                        # New setting[var_units] = m[var_units/table_units]*(setting-o1[table_units]) + o2[var_units]
                        # c.f. mathml_units_mixin._add_units_conversion
                        print 'LT conversion:', table_units.description(), 'to', var_units.description(), 'equal?', table_units.equals(var_units)
                        m = table_units.get_multiplicative_factor() / var_units.get_multiplicative_factor()
                        for setting in new_config[var]:
                            try:
                                old_value = float(new_config[var][setting])
                                new_value = m * (old_value - table_units.get_offset()) + var_units.get_offset()
                                new_config[var][setting] = unicode(new_value)
                                print 'LT conversion', setting, old_value, new_value
                            except (ValueError, TypeError):
                                pass
        self.lut_config = new_config
        DEBUG('config', 'Lookup tables configuration:', new_config)
        return

    # TODO - move into seperate metadata class?
    def validate_metadata(self, assume_valid=False):
        """Check that the metadata annotations are 'sensible'.
        
        Ensures that only names we know are used, and that the same name isn't used for multiple variables.
        """
        vars = cellml_metadata.find_variables(self.doc.model, ('bqbiol:is', NSS['bqbiol']))
        self.metadata_vars = filter(lambda v: v.oxmeta_name, vars)
        if assume_valid:
            return
        names_used = [var.oxmeta_name for var in self.metadata_vars]
        DEBUG('metadata', 'Names found: ', names_used)
        # Check all metadata is allowed
        unknown_names = frozenset(names_used) - cellml_metadata.METADATA_NAMES
        if unknown_names:
            msg = ['Unrecognised oxmeta variable names found (run with --assume-valid to ignore):']
            msg.extend(sorted(unknown_names))
            raise ConfigurationError('\n  '.join(msg))
        # Check for duplicates
        d = {}
        for name in names_used:
            if name in d:
                raise ConfigurationError(name + ' metadata attribute is duplicated in the cellml file.')
            else:
                d[name] = name


######################################################################
#                    For running as an executable                    #
######################################################################

def get_options(args, default_options=None):
    """get_options(args):
    Process our command-line options.

    args is a list of options & positional arguments.
    
    default_options, if given, is an instance of optparse.Values created by a
    previous call to this function.
    """
    usage = 'usage: %prog [options] <cellml file or URI>'
    parser = optparse.OptionParser(version="%%prog %s" % __version__,
                                   usage=usage)
    parser.add_option('-q', '--quiet', action='store_true', default=False,
                      help="don't show warning messages, only errors")
    # What type of translation is being performed
    parser.add_option('-T', '--translate',
                      dest='translate', action='store_true',
                      default=True,
                      help="output computer code [default]")
    parser.add_option('-C', '--output-cellml',
                      dest='translate', action='store_false',
                      help="output an annotated CellML file instead of translating, on stdout unless -o specified")
    translators = sorted(CellMLTranslator.translators)
    parser.add_option('-t', '--translate-type',
                      type='choice', choices=translators,
                      default='Nektar', metavar='TYPE',
                      help="the type of code to output [default: %default].  "
                      "Choices: " + str(translators))
    parser.add_option('-o', dest='outfilename', metavar='OUTFILE',
                      help="write program code to OUTFILE [default action is to use the input filename with a different extension]")
    # Global adjustment settings
    parser.add_option('--config-file',
                      action='append', default=[],
                      help="pathname of configuration file")
    parser.add_option('-A', '--fully-automatic',
                      action='store_true', default=False,
                      help="if human intervention is required, fail noisily")
    parser.add_option('--assume-valid',
                      action='store_true', default=False,
                      help="skip some of the model validation checks")
    parser.add_option('--warn-on-unit-conversions',
                      action='store_true', default=False,
                      help="generate a warning if unit conversions are required")
    parser.add_option('--Wu', '--warn-on-units-errors',
                      action='store_true', default=False,
                      dest='warn_on_units_errors',
                      help="give a warning instead of an error for dimensional inconsistencies")
    parser.add_option('-V', '--transmembrane-potential', default=None, metavar='POT_VAR',
                      help="POT_VAR is the full name of the variable representing the transmembrane potential."
                      "  If not specified here, the configuration file will be used, which is the prefered method."
                      "  Defaults to 'membrane,V'.")
    parser.add_option('-d', '--debug', action='store_true', default=False,
                      help="output debug info to stderr")
    parser.add_option('-D', '--debug-source', action='append',
                      help="only show debug info from the specified part of the code."
                      "  This option may appear more than once to select multiple sources.  Implies -d.")
    parser.add_option('--profile', action='store_true', default=False,
                      help="turn on profiling of PyCml")
    # To examine the profile do something like:
    #    import os,pstats
    #    os.chdir('/tmp')
    #    files = filter(lambda f: f.startswith('pycml'), os.listdir('.'))
    #    p = pstats.Stats(*files)
    #    p.strip_dirs().sort_stats('cum').print_stats(15)
    # What optimisations/transformations to do
    group = optparse.OptionGroup(parser, 'Transformations',
                                 "These options control which transformations (typically optimisations) are applied in the generated code")
    group.add_option('-l', '--lookup-tables',
                     dest='lut', action='store_true', default=False,
                     help="perform a lookup table analysis")
    group.add_option('-p', '--pe', '--partial-evaluation',
                     dest='pe', action='store_true', default=False,
                     help="partially evaluate the model")
    group.add_option('-u', '--units-conversions',
                     action='store_true', default=False,
                     help="add explicit units conversion mathematics")
    group.add_option('-j', '--maple-output',
                     metavar='FILENAME', default=None,
                     help="file containing output from a Maple script generated using -J.  The generated"
                     " code/CellML will then contain a symbolic Jacobian as computed by Maple.")
    group.add_option('-J', '--do-jacobian-analysis',
                     action='store_true', default=False,
                     help="generate code to perform Jacobian analysis for backward Euler & CVODE; implies -t Maple")
    group.add_option('--backward-euler',
                     action='store_true', default=False,
                     help="generate a specialised cell model that solves itself using a decoupled"
                     " backward Euler method.  Not compatible with --rush-larsen.  Implies -t Chaste."
                     "  Requires -j.")
    group.add_option('--rush-larsen',
                     action='store_true', default=False,
                     help="use the Rush-Larsen method to solve Hodgkin-Huxley style gating variable"
                     " equations.  Not compatible with --backward-euler.  Implies -t Chaste.")
    group.add_option('--grl1',
                     action='store_true', default=False,
                     help="use the GRL1 method to solve Hodgkin-Huxley style gating variable"
                     " equations.  Not compatible with the backward Euler transformation."
                     " Implies -t Chaste.")
    group.add_option('--grl2',
                     action='store_true', default=False,
                     help="use the GRL2 method to solve Hodgkin-Huxley style gating variable"
                     " equations.  Not compatible with the backward Euler transformation."
                     " Implies -t Chaste.")
    parser.add_option_group(group)
    # Settings tweaking the generated code
    group = optparse.OptionGroup(parser, 'Generated code options')
    group.add_option('-c', '--class-name', default=None,
                     help="explicitly set the name of the generated class")
    group.add_option('-a', '--augment-class-name',
                     dest='augment_class_name', action='store_true',
                     default=False,
                     help="alter the class name to show what transformations are used")
    group.add_option('--no-timestamp',
                     action='store_true', default=False,
                     help="don't add a timestamp comment to generated files")
    parser.add_option_group(group)
    # Options specific to Maple output
    group = optparse.OptionGroup(parser, 'Maple options', "Options specific to Maple code output")
    group.add_option('--dont-omit-constants',
                     dest='omit_constants', action='store_false', default=True,
                     help="when generating Maple code, include assignments of constants")
    group.add_option('--compute-partial-jacobian', dest='compute_full_jacobian',
                     action='store_false', default=True,
                     help="make generated Maple code compute a Jacobian specific to a Newton solve"
                     " of the nonlinear portion of the ODE system, rather than the full system Jacobian")
    parser.add_option_group(group)
    # Options specific to Python output
    group = optparse.OptionGroup(parser, 'Python options', "Options specific to Python code output")
    group.add_option('--no-numba', dest='numba', default=True, action='store_false',
                     help="turn off using Numba to optimise code on-the-fly")
    parser.add_option_group(group)
    # Options specific to Chaste output
    group = optparse.OptionGroup(parser, 'Chaste options', "Options specific to Chaste code output")
    group.add_option('-y', '--dll', '--dynamically-loadable',
                     dest='dynamically_loadable',
                     action='store_true', default=False,
                     help="add code to allow the model to be compiled to a shared library and dynamically loaded"
                     " (only works if -t Chaste is used)")
    group.add_option('--use-chaste-stimulus',
                     action='store_true', default=False,
                     help="when generating Chaste code, use Chaste's stimulus rather than that defined in the model")
    group.add_option('--no-use-chaste-stimulus', dest='use_chaste_stimulus',
                     action='store_false',
                     help="when generating Chaste code, use the model's stimulus, not Chaste's")
    group.add_option('-i', '--convert-interfaces',
                     action='store_true', default=False,
                     help="perform units conversions at interfaces to Chaste (only works if -t Chaste is used)")
    group.add_option('--use-i-ionic-regexp', dest='use_i_ionic_regexp',
                     action='store_true', default=False,
                     help="determine ionic currents from the regexp specified in the config file"
                     " rather than analysing the voltage derivative equation")
    group.add_option('--include-dt-in-tables',
                     action='store_true', default=False,
                     help="[experimental] allow timestep to be included in lookup tables.  By default"
                     " uses the timestep of the first cell created.  Requires support from external"
                     " code if timestep changes.  Only really useful for backward Euler cells.")
    group.add_option('-m', '--use-modifiers',
                     action='store_true', default=False,
                     help="[experimental] add modifier functions for certain"
                     " metadata-annotated variables for use in sensitivity analysis (only works if -t Chaste is used)")
    group.add_option('--use-data-clamp',
                     action='store_true', default=False,
                     help="[experimental] generate a data clamp subclass of CVODE cells"
                     " which contains data clamp currents for fitting experimental data (only works if -t CVODE is used)")
    group.add_option('--expose-annotated-variables',
                     action='store_true', default=False,
                     help="expose all oxmeta-annotated variables for access via the GetAnyVariable functionality")
    group.add_option('--expose-all-variables',
                     action='store_true', default=False,
                     help="expose all variables for access via the GetAnyVariable functionality")
    parser.add_option_group(group)
    # Options specific to Functional Curation
    group = optparse.OptionGroup(parser, 'Functional Curation options', "Options specific to use by Functional Curation")
    def protocol_callback(option, opt_str, value, parser):
        """
        Protocols don't always produce normal cardiac cell models.
        However, we want to allow a later option to override these changes.
        """
        parser.values.protocol = value
        parser.values.convert_interfaces = False
        parser.values.use_chaste_stimulus = False
    group.add_option('--protocol',
                     action='callback', callback=protocol_callback, type='string', nargs=1,
                     help="specify a simulation protocol to apply to the model prior to translation")
    group.add_option('--protocol-options', action='store', type='string',
                     help="extra options for the protocol")
    group.add_option('--expose-named-parameters',
                     action='store_true', default=False,
                     help="expose all constant variables with 'name' annotations for access as model parameters")
    parser.add_option_group(group)
    # Settings for lookup tables
    group = optparse.OptionGroup(parser, 'Lookup tables options', "Options specific to the lookup tables optimisation")
    lookup_type_choices = ['entry-below', 'nearest-neighbour', 'linear-interpolation']
    group.add_option('--lookup-type', choices=lookup_type_choices,
                     default='linear-interpolation',
                     help="the type of table lookup to perform [default: %default]."
                     " Choices: " + str(lookup_type_choices))
    group.add_option('--no-separate-lut-class', dest='separate_lut_class',
                     action='store_false', default=True,
                     help="don't put lookup tables in a separate class")
    group.add_option('--row-lookup-method',
                     action='store_true', default=True,
                     help="add and use a method to look up a whole row of a table")
    group.add_option('--no-row-lookup-method', dest='row_lookup_method',
                     action='store_false',
                     help="don't add and use a method to look up a whole row of a table")
    group.add_option('--combine-commutative-tables',
                     action='store_true', default=False,
                     help="optimise a special corner case to reduce the number of tables."
                     " See documentation for details.")
    group.add_option('--lt-index-uses-floor',
                     action='store_true', default=False,
                     help="use floor() to calculate LT indices, instead of just casting")
    group.add_option('--constrain-table-indices',
                     action='store_true', default=False,
                     help="constrain lookup table index variables to remain within the bounds specified,"
                     " rather than throwing an exception if they go outside the bounds")
    group.add_option('--no-check-lt-bounds', dest='check_lt_bounds',
                     action='store_false', default=True,
                     help="[unsafe] don't check for LT indexes going outside the table bounds")
    parser.add_option_group(group)
    # Settings for partial evaluation
    group = optparse.OptionGroup(parser, 'Partial evaluation options', "Options specific to the partial evaluation optimisation")
    group.add_option('--pe-convert-power',
                     action='store_true', default=False,
                     help="convert pow(x,3) to x*x*x; similarly for powers 2 & 4.")
    group.add_option('--no-partial-pe-commutative', dest='partial_pe_commutative',
                     action='store_false', default=True,
                     help="don't combine static operands of dynamic commutative associative applys")
    group.add_option('--no-pe-instantiate-tables', dest='pe_instantiate_tables',
                     action='store_false', default=True,
                     help="don't instantiate definitions that will be tables regardless of usage")
    parser.add_option_group(group)

    options, args = parser.parse_args(args, values=default_options)
    if len(args) != 1:
        parser.error("exactly one input CellML file must be specified")

    # Some options imply others
    if options.debug_source:
        options.debug = True
    if options.do_jacobian_analysis:
        options.translate_type = 'Maple'
        options.maple_output = False
        options.rush_larsen = False
        options.backward_euler = False
    if options.backward_euler:
        if not options.maple_output:
            parser.error("Backward Euler code generation requires maple output (-j)")
        options.rush_larsen = False
        options.grl1 = False
        options.grl2 = False
    if options.rush_larsen or options.backward_euler or options.grl1 or options.grl2:
        options.translate_type = 'Chaste'
    if options.use_data_clamp and not options.translate_type=='CVODE':
        parser.error("Data clamp option '--use-data-clamp' also requires CVODE ('-t CVODE'). If you are calling this via ConvertCellModel use '--cvode-data-clamp'.")
    # Numba may not be available
    if options.numba:
        try:
            import numba
        except:
            options.numba = False

    return options, args[0]


def load_model(model_file, options):
    """Load and validate a CellML model."""
    # Setup logging
    logging.thread = None # Hack: we're not multi-threaded, so be slightly quicker...
    if options.debug:
        formatter = logging.Formatter(fmt="%(name)s: %(message)s")
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(formatter)
        handler.addFilter(OnlyDebugFilter())
        if options.debug_source:
            handler.addFilter(OnlyTheseSourcesFilter(options.debug_source))
        logging.getLogger().addHandler(handler)
        logging.getLogger().setLevel(logging.DEBUG)

    # We can't translate if some warnings occur, as well as if the model is invalid
    notifier = NotifyHandler(level=logging.WARNING_TRANSLATE_ERROR)
    logging.getLogger('validator').addHandler(notifier)
    v = validator.CellMLValidator(create_relaxng_validator=not options.assume_valid)
    valid, doc = v.validate(model_file, return_doc=True, show_warnings=not options.quiet,
                            check_for_units_conversions=options.warn_on_unit_conversions,
                            warn_on_units_errors=options.warn_on_units_errors,
                            assume_valid=options.assume_valid)
    v.quit()
    del v

    if not valid or notifier.messages:
        print >>sys.stderr, model_file,
        if not valid:
            print >>sys.stderr, "is not a valid CellML file"
        else:
            print >>sys.stderr, "contains untranslatable constructs (see warnings above for details)"
        sys.exit(1)
    
    return doc

def run():
    """Translate the file given on the command line."""
    options, model_file = get_options(sys.argv[1:])
    doc = load_model(model_file, options)
    DEBUG('translate', "+++ Loaded model")

    config = ConfigurationStore(doc, options=options)
    for config_file in options.config_file:
        config.read_configuration_file(config_file)
    DEBUG('translate', "+++ Read config")
    
    # Apply protocol, if given
    if options.protocol:
        import protocol
        protocol.apply_protocol_file(doc, options.protocol)
        if options.debug:
            post_proto_cellml = options.outfilename or model_file
            post_proto_cellml = os.path.splitext(post_proto_cellml)[0] + '-proto.cellml.ppp'
            stream = open_output_stream(post_proto_cellml)
            doc.xml(indent=u'yes', stream=stream)
            close_output_stream(stream)
        DEBUG('translate', "+++ Applied protocol")

    config.finalize_config()
    DEBUG('translate', "+++ Processed config")

    solver_info = SolverInfo(doc.model)

    # Generate an interface component, if desired
    translator_klass = CellMLTranslator.translators[options.translate_type]
    if not options.protocol:
        translator_klass.generate_interface(doc, solver_info)
    config.validate_metadata(options.assume_valid)
    DEBUG('translate', "+++ Generated interface")
    
    if options.lut:
        config.find_lookup_variables()
        DEBUG('translate', "+++ Found LT keys")

    # These bits could do with improving, as they annotate more than is really needed!
    if options.pe:
        # We need to ensure PE doesn't remove ionic currents needed for GetIIonic
        config.annotate_currents_for_pe()
        # "Need" to ensure pe doesn't remove metadata-annotated variables (when using modifiers or default stimulus?)
        config.annotate_metadata_for_pe()
        DEBUG('translate', "+++ Annotated variables")
    # Deal with the 'expose' options
    config.expose_variables()

    class_name = options.class_name
    if not class_name:
        class_name = doc.model.name.replace('-', '_')
        if options.augment_class_name:
            class_name = u'CML_' + class_name
            if options.pe:
                class_name += '_pe'
            if options.lut:
                class_name += '_lut'
            if options.backward_euler:
                class_name += '_be'
            if options.use_modifiers:
                class_name += '_sens'
    if options.protocol:
        # Try to avoid OdeSystemInformation conflicts
        class_name += "_Proto_" + os.path.splitext(os.path.basename(options.protocol))[0]

    output_filename = getattr(options, 'outfilename', None)
    if not options.translate and not output_filename:
        output_filename = 'stdout'

    if options.units_conversions:
        doc.model.add_units_conversions()
        DEBUG('translate', "+++ Added units conversions")

    if options.do_jacobian_analysis:
        lin = optimize.LinearityAnalyser()
        lin.analyse_for_jacobian(doc, V=config.V_variable)
        DEBUG('translate', "+++ Analysed model for Jacobian")

    if options.maple_output:
        # Parse Jacobian matrix
        from maple_parser import MapleParser
        mp = MapleParser()
        jacobian_file = file(options.maple_output) # TODO: Error checking
        doc.model._cml_jacobian = mp.parse(jacobian_file)
        doc.model._cml_jacobian_full = mp.JacobianWasFullSize
        jacobian_file.close()
        if not options.backward_euler and doc.model._cml_jacobian_full:
            # Add full jacobian to XML
            solver_info.add_jacobian_matrix()
            solver_info.add_variable_links()

    if options.backward_euler:
        # Rearrange linear ODEs
        lin = optimize.LinearityAnalyser()
        lin.analyse_for_jacobian(doc, V=config.V_variable)
        lin.rearrange_linear_odes(doc)
        # Remove jacobian entries that don't correspond to nonlinear state variables
        jacobian = doc.model._cml_jacobian
        if isinstance(jacobian, tuple):
            assert doc.model._cml_jacobian_full
            jacobian = jacobian[1]
        nonlinear_vars = set([v.get_source_variable(recurse=True) for v in doc.model._cml_nonlinear_system_variables])
        def gv(vname):
            return cellml_variable.get_variable_object(doc.model, vname).get_source_variable(recurse=True)
        for var_i, var_j in jacobian.keys():
            if gv(var_i) not in nonlinear_vars or gv(var_j) not in nonlinear_vars:
                del jacobian[(var_i, var_j)]
        if doc.model._cml_jacobian_full:
            # Transform the Jacobian into the form needed by the Backward Euler code
            import maple_parser
            for key, expr in jacobian.iteritems():
                new_expr = None
                if key[0] == key[1]:
                    # 1 on the diagonal
                    new_expr = maple_parser.MNumber(['1'])
                if not (isinstance(expr, maple_parser.MNumber) and str(expr) == '0'):
                    # subtract delta_t * expr
                    args = []
                    if new_expr:
                        args.append(new_expr)
                    args.append(maple_parser.MOperator([maple_parser.MVariable(['delta_t']), expr], 'prod', 'times'))
                    new_expr = maple_parser.MOperator(args, '', 'minus')
                if new_expr:
                    jacobian[key] = new_expr
        # Add info as XML
        solver_info.add_all_info()
        # Analyse the XML, adding cellml_variable references, etc.
        solver_info.add_variable_links()
        solver_info.add_linear_ode_update_equations()
        DEBUG('translate', "+++ Parsed and incorporated Maple output")
    else:
        options.include_dt_in_tables = False

    if options.lut:
        # Create the analyser so PE knows which variables are table keys
        lut = optimize.LookupTableAnalyser()
    else:
        lut = None

    if options.pe:
        # Do partial evaluation
        pe = optimize.PartialEvaluator()
        pe.parteval(doc, solver_info, lut)
        DEBUG('translate', "+++ Done PE")

    if options.lut:
        # Do the lookup table analysis
        lut.analyse_model(doc, solver_info)
        DEBUG('translate', "+++ Done LT analysis")
    
    if options.rush_larsen:
        rl = optimize.RushLarsenAnalyser()
        rl.analyse_model(doc)
        DEBUG('translate', "+++ Done Rush-Larsen analysis")

    if options.translate:
        # Translate to code
        initargs = {'add_timestamp': not options.no_timestamp,
                    'options': options}
        transargs = {'v_variable': config.V_variable}
        transargs['row_lookup_method'] = options.row_lookup_method
        transargs['lt_index_uses_floor'] = options.lt_index_uses_floor
        transargs['constrain_table_indices'] = options.constrain_table_indices
        """if issubclass(translator_klass, CellMLToMapleTranslator):
            initargs['omit_constants'] = options.omit_constants
            initargs['compute_full_jacobian'] = options.compute_full_jacobian
        el"""
        from CellMLToNektarTranslator import CellMLToNektarTranslator
        if issubclass(translator_klass, CellMLToNektarTranslator):
            solver_info.add_membrane_ionic_current()
            transargs['use_chaste_stimulus'] = options.use_chaste_stimulus
            transargs['separate_lut_class'] = options.separate_lut_class
            transargs['convert_interfaces'] = options.convert_interfaces
            transargs['use_modifiers'] = options.use_modifiers
            transargs['use_data_clamp'] = options.use_data_clamp
            transargs['dynamically_loadable'] = options.dynamically_loadable
            transargs['use_protocol'] = bool(options.protocol)
        t = translator_klass(**initargs)
        t.translate(doc, model_file, output_filename, class_name=class_name, **transargs)
        cellml_metadata.remove_model(doc.model)
    else:
        # Add a comment element
        comment = pycml.comment_base(
            body=u'\n' + version_comment(not options.no_timestamp) + u'\n')
        doc.xml_insert_before(doc.model, comment)
        # Output annotated model
        stream = open_output_stream(output_filename)
        doc.xml(indent=u'yes', stream=stream)
        close_output_stream(stream)

        
    DEBUG('translate', "+++ Done translation")
