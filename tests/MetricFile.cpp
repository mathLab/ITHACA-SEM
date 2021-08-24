///////////////////////////////////////////////////////////////////////////////
//
// File: MetricFile.cpp
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
// Description: Implementation of the LInf metric.
//
///////////////////////////////////////////////////////////////////////////////

#include <iterator>
#include <fstream>
#include <vector>

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string.hpp>

#include <MetricFile.h>
#include <sha1.h>

namespace Nektar
{
    std::string MetricFile::type = GetMetricFactory().
        RegisterCreatorFunction("FILE", MetricFile::create);

    MetricFile::MetricFile(TiXmlElement *metric, bool generate) : 
        Metric(metric, generate)
    {
        TiXmlElement *file = metric->FirstChildElement("file");
        ASSERTL0(file, "Missing file tag for file metric!");
        
        while (file)
        {
            std::string filename, sha1hash;
            
            if (file->Attribute("filename"))
            {
                filename = file->Attribute("filename");
            }
            else
            {
                ASSERTL0(false, "Missing filename for file tag!");
            }
            
            if (!m_generate)
            {
                TiXmlElement *sha1 = file->FirstChildElement("sha1");
                ASSERTL0(sha1, "Missing SHA1 hash for file "+filename);
                sha1hash = sha1->GetText();
                ASSERTL0(sha1hash.size() == 40, 
                         "Incorrect length for SHA1 hash");
            }
            
            m_filehash[filename] = sha1hash;
            
            file = file->NextSiblingElement("file");
        }
    }
    
    std::string MetricFile::CalculateHash(std::string pfilename)
    {
        int    fdot   = pfilename.find_last_of('.');
        string ending = pfilename.substr(fdot);
        string filename; 
        
        if(ending == ".fld" || ending == ".chk" || ending == ".rst")
        {
            TiXmlDocument *xmlFldFile;
            fs::path pathfilename(pfilename);
            
            if(fs::is_directory(pathfilename))
            {
                std::vector<fs::path> dirfiles;
                
                // make list of all files in directoryj
                copy(fs::directory_iterator(pathfilename), 
                     fs::directory_iterator(),
                     back_inserter(dirfiles));
                
                xmlFldFile = new TiXmlDocument;
                
                // load all them into the 
                for(int i = 0; i < dirfiles.size(); ++i)
                {
                    std::string infile = PortablePath(dirfiles[i]);
                    std::ifstream file(infile.c_str());
                    ASSERTL0(file.good(), "Unable to open file: " + infile);
                    file >> (*xmlFldFile);
                }
            }
            else
            {
                xmlFldFile = new TiXmlDocument(pfilename);
                xmlFldFile->LoadFile(pfilename);
            }       

            // strip out meta data before check
            TiXmlElement* vNektar  = xmlFldFile->FirstChildElement("NEKTAR");
            while(vNektar)
            {
            
                TiXmlNode* vMetaData = vNektar->FirstChild("Metadata");
                
                // delete MetaData section
                if(vMetaData) 
                {
                    vNektar->RemoveChild(vMetaData);
                }
                vNektar = vNektar->NextSiblingElement("NEKTAR");
            }                

            filename = pfilename + ".tmp";
            xmlFldFile->SaveFile(filename);

        }
        else
        {
            filename = pfilename;
        }
            
        // Open file.
        std::ifstream testFile(filename.c_str(), std::ios::binary);
        ASSERTL0(testFile.is_open(), "Error opening file "+filename);

        // Read in file contents.
        std::vector<char> fileContents(
            (std::istreambuf_iterator<char>(testFile)),
            std::istreambuf_iterator<char>());

        // Calculate SHA1 hash.
        unsigned char hash[20];
        char strhash[41];
        
        sha1::calc((void *)&fileContents[0], fileContents.size(), hash);
        sha1::toHexString(hash, strhash);
        
        // Close file and return string containing SHA1 hash.
        testFile.close();
        
        return std::string(strhash);
    }
    
    bool MetricFile::v_Test(std::istream& pStdout, std::istream& pStderr)
    {
        boost::ignore_unused(pStdout, pStderr);

        std::map<std::string, std::string>::iterator it;
        bool success = true;
        
        for (it = m_filehash.begin(); it != m_filehash.end(); ++it)
        {
            std::string filehash = CalculateHash(it->first);
            if (!boost::iequals(filehash, it->second))
            {
                std::cerr << "Failed SHA1 hash test." << std::endl;
                std::cerr << "  Expected: " << it->second
                          << std::endl;
                std::cerr << "  Result:   " << filehash << std::endl;
                success = false;
            }
        }
        
        return success;
    }

    void MetricFile::v_Generate(std::istream& pStdout, std::istream& pStderr)
    {
        boost::ignore_unused(pStdout, pStderr);

        std::map<std::string, std::string>::iterator it;

        // Update SHA1 hashes.
        for (it = m_filehash.begin(); it != m_filehash.end(); ++it)
        {
            std::string filehash = CalculateHash(it->first);
            m_filehash[it->first] = filehash;
        }
        
        // Write new XML structure.
        TiXmlElement *file = m_metric->FirstChildElement("file");
        while (file)
        {
            std::string filename = file->Attribute("filename");
            file->Clear();

            TiXmlElement *sha1 = new TiXmlElement("sha1");

            ASSERTL0(m_filehash.count(filename) != 0,
                     "Couldn't find file " + filename + " in list of calculated"
                     "hashes");
            
            sha1->LinkEndChild(new TiXmlText(m_filehash[filename]));
            file->LinkEndChild(sha1);
            file = file->NextSiblingElement("file");
        }
    }
}
