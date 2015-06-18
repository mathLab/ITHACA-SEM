////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldIOBenchmarker.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and laimitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Measure the performance of FieldIO
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/Communication/CommMpi.h>

using namespace Nektar;
using namespace LibUtilities;

namespace po = boost::program_options;

struct Experiment {
        bool write;
        bool hdf;
        bool verbose;
        int n;
        std::string dataSource;
        CommSharedPtr comm;
};

typedef std::vector<double> Results;

Results TestRead(Experiment& exp);
Results TestWrite(Experiment& exp);
void PrintResults(Experiment& exp, Results& results);

int main(int argc, char* argv[])
{
    Experiment exp;
    exp.write = false;
    exp.hdf = false;
    exp.verbose = false;
    exp.n = 3;
    exp.comm = GetCommFactory().CreateInstance("ParallelMPI", argc, argv);

    po::options_description desc("Available options");
    desc.add_options()
        ("help,h",
                "Produce this help message.")
        ("mode,m", po::value<char>(),
                "Choose r[ead] (default), x[ml write] or h[df5 write]")
        ("number,n", po::value<unsigned>(),
                "Number of iterations to perform, default 3")
        ("verbose,v",
                "Enable verbose mode.");

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-file",   po::value<std::string>(), "Input filename");

    po::options_description cmdline_options;
    cmdline_options.add(hidden).add(desc);

    po::options_description visible("Allowed options");
    visible.add(desc);

    po::positional_options_description p;
    p.add("input-file", -1);

    po::variables_map vm;

    try
    {
        po::store(po::command_line_parser(argc, argv).
                  options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << desc;
        return 1;
    }


    if (vm.count("help") || vm.count("input-file") != 1) {
        std::cerr << "Usage: FieldIOBenchmarker [options] inputfile"
             << endl;
        std::cout << desc;
        std::cout << endl;
        return 1;
    }

    ASSERTL0(vm.count("input-file"),
             "Must specify input.");

    exp.dataSource = vm["input-file"].as<std::string>();


    if (vm.count("verbose") && exp.comm->GetRank()==0)
    {
        exp.verbose = true;
    }

    if (vm.count("number"))
    {
        exp.n = vm["number"].as<unsigned>();
    }

    if (vm.count("mode"))
    {
        char mode = vm["mode"].as<char>();
        switch (mode) {
            case 'r':
                exp.write = false;
                break;
            case 'x':
                exp.write = true;
                exp.hdf = false;
                break;
            case 'h':
                exp.write = true;
                exp.hdf = true;
                break;
            default:
                std::cout << "Unrecognised mode: " << mode << std::endl;
                std::cout << desc << endl;
                return 1;
                break;
        }
    }

    Results res;
    if (exp.write)
    {
        res = TestWrite(exp);
    }
    else
    {
        res = TestRead(exp);
    }

    PrintResults(exp, res);
    exp.comm->Finalise();
}

Results TestRead(Experiment& exp)
{
    if (exp.verbose)
    {
        std::cout << "Beginning read experiment with " << exp.n << " loops." << std::endl;
        std::cout << "Determining file type... ";
    }

    const std::string ft = FieldIO::GetFileType(exp.dataSource, exp.comm);
    if (exp.verbose)
        std::cout << ft << endl;

    Results res(exp.n, 0.0);
    for (unsigned i=0; i<exp.n; ++i)
    {
        if (exp.verbose)
            std::cout << "Test " << i << " of " << exp.n;
	
	// Synchronise
	exp.comm->Block();
	
        double t0 = MPI_Wtime();

        FieldIOSharedPtr fio = GetFieldIOFactory().CreateInstance(ft, exp.comm);
        std::vector<FieldDefinitionsSharedPtr> fielddefs;
        std::vector<std::vector<NekDouble> > fielddata;

        fio->Import(exp.dataSource,
                 fielddefs,
                fielddata);

        double t1 = MPI_Wtime();
        t1 -= t0;

        if (exp.verbose)
            std::cout << ": t = " << t1 << " s" << std::endl;

        res[i] = t1;
    }
    return res;
}
Results TestWrite(Experiment& exp)
{
    if (exp.verbose)
        std::cout << "Reading in input: " << exp.dataSource << std::endl;

    const std::string intype = FieldIO::GetFileType(exp.dataSource, exp.comm);
    FieldIOSharedPtr fioRead = GetFieldIOFactory().CreateInstance(intype, exp.comm);
    std::vector<FieldDefinitionsSharedPtr> fielddefs;
    std::vector<std::vector<NekDouble> > fielddata;
    fioRead->Import(exp.dataSource,
            fielddefs,
            fielddata);

    std::string outfile = exp.dataSource + ".tmp";
    std::string outtype;
    if (exp.hdf)
        outtype = "Hdf5";
    else
        outtype = "Xml";

    if (exp.verbose)
    {
        std::cout << "Beginning write (" << outtype << ") experiment with " << exp.n << " loops." << std::endl;
        std::cout << "Writing to temp file: " << outfile << std::endl;
    }

    Results res(exp.n, 0);
    for (unsigned i=0; i<exp.n; ++i)
    {
        if (exp.verbose)
            std::cout << "Test " << i << " of " << exp.n << std::endl;
	
	// Synchronise
	exp.comm->Block();
        double t0 = MPI_Wtime();

        FieldIOSharedPtr fio = GetFieldIOFactory().CreateInstance(outtype, exp.comm);

        fio->Write(outfile,
                 fielddefs,
                 fielddata);

        double t1 = MPI_Wtime();
        t1 -= t0;
	
	
        if (exp.verbose)
            std::cout << ": t = " << t1 << " s" << std::endl;

        res[i] = t1;
    }
    return res;
}
void PrintResults(Experiment& exp, Results& results)
{
    double sum = 0.0;
    double sumSq = 0.0;

    for (Results::const_iterator it = results.begin(); it != results.end(); ++it)
    {
        double x = *it;
        sum += x;
        sumSq += x*x;
    }
    
    double mean = sum / exp.n;
    // double var = sumSq / exp.n - mean*mean;
    // double std = std::sqrt(var);
    
    if (exp.comm->GetSize() > 1) {
      // Use all version since reduce to specified rank isn't wrapped.
      exp.comm->AllReduce(mean, ReduceSum);
      mean  /= exp.comm->GetSize();
    }
    
    if (exp.comm->GetRank() == 0)
    {
      std::cout << "Mean: " << mean << std::endl;
      // std::cout << "Std: " << std << std::endl;
    }
}
