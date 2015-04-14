#include <cstdlib>
#include <iostream>
#include <sstream>

#include <tinyxml.h>

#include <LibUtilities/BasicUtils/FileSystem.h>

int main(int argc, char *argv[]) 
{
    bool DeleteFiles = false;
    if (argc < 3)
    {
        std::cout << "Usage: CombinePartitions [DeleteFile] nproc outfile"
                  << std::endl;
        std::cout << "  [DeleteFiles]    = Delete partiion files (optional)" << std::endl;
        std::cout << "  nproc            = Number of partitions" << std::endl;
        std::cout << "  outfile          = Target output filename" << std::endl;
        exit(1);
    }

    if(argc == 4)
    {
        DeleteFiles = true;
    }

    TiXmlDocument docOutput;
    TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "utf-8", "");
    docOutput.LinkEndChild(decl);

    TiXmlElement *master = new TiXmlElement("NEKTAR");
    std::string basename = argv[argc-1];
    std::string extension = argv[argc-1];
    basename = basename.substr(0, basename.find_last_of("."));
    extension = extension.substr(extension.find_last_of(".") + 1);
    
    fs::path infile(argv[argc-1]);
    bool isdirectory = fs::is_directory(infile);
    
    for (int n = 0; n < atoi(argv[argc-2]); ++n)
    {
        std::stringstream filename;
        if(isdirectory)
        {
            filename << argv[argc-1] << "/P" << n << ".fld"; 
        }
        else
        {
            filename << basename << "_P" << n << "." << extension;
        }
        TiXmlDocument docInput;
        if (!docInput.LoadFile(filename.str()))
        {
            std::cerr << "Unable to open file '" << filename.str() << "'." << std::endl;
            exit(1);
        }
        TiXmlElement *nektar = docInput.FirstChildElement("NEKTAR");
        
        // load up root processor's meta data
        if(n == 0 && nektar->FirstChildElement("Metadata"))
        {
            TiXmlElement *metadata = nektar->FirstChildElement("Metadata");
            if(metadata)
            {
                master->LinkEndChild(new TiXmlElement(*metadata));
            }
        }

        // load the elements from seperate files. 
        TiXmlElement *elements = nektar->FirstChildElement("ELEMENTS");
        while (elements)
        {
            master->LinkEndChild(new TiXmlElement(*elements));
            elements = elements->NextSiblingElement();
        }
    }
    
    docOutput.LinkEndChild(master);
    
    if(isdirectory)
    {
        std::string outname = basename  + "_combined" + "." + extension; ;
        if (!docOutput.SaveFile(outname))
        {
            std::cerr << "Unable to write file '" << outname << "'." << std::endl;
        }
    }
    else if (!docOutput.SaveFile(argv[argc-1]))
    {
        std::cerr << "Unable to write file '" << argv[argc-1] << "'." << std::endl;
    }
    else
    {
        if(DeleteFiles)
        {

            for (int n = 0; n < atoi(argv[argc-2]); ++n)
            {
                std::string basename = argv[argc-1];
                std::string extension = argv[argc-1];
                basename = basename.substr(0, basename.find_last_of("."));
                extension = extension.substr(extension.find_last_of(".") + 1);
                std::stringstream filename;
                filename << basename << "_P" << n << "." << extension;
                std::string deloutput = "rm -rf ";
                deloutput = deloutput + filename.str();
                system(deloutput.c_str());
            }
        }
    }

    exit(0);
}
