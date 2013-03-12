#include <cstdlib>
#include <iostream>
#include <sstream>

#include <tinyxml/tinyxml.h>

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: CombinePartitions [nproc] [outfile]"
                  << std::endl;
        std::cout << "  [nproc]            = Number of partitions" << std::endl;
        std::cout << "  [outfile]          = Target output filename" << std::endl;
        exit(1);
    }

    TiXmlDocument docOutput;
    TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "utf-8", "");
    docOutput.LinkEndChild(decl);

    TiXmlElement *master = new TiXmlElement("NEKTAR");
    for (int n = 0; n < atoi(argv[1]); ++n)
    {
        std::string basename = argv[2];
        std::string extension = argv[2];
        basename = basename.substr(0, basename.find_last_of("."));
        extension = extension.substr(extension.find_last_of(".") + 1);
        std::stringstream filename;
        filename << basename << "_P" << n << "." << extension;
        TiXmlDocument docInput;
        if (!docInput.LoadFile(filename.str()))
        {
            std::cerr << "Unable to open file '" << filename.str() << "'." << std::endl;
            exit(1);
        }
        TiXmlElement *nektar = docInput.FirstChildElement("NEKTAR");
        
        // load up root processor's meta data
        if(n == 0)
        {
            TiXmlElement *metadata = nektar->FirstChildElement("FIELDMETADATA");
            master->LinkEndChild(new TiXmlElement(*metadata));
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
    if (!docOutput.SaveFile(argv[2]))
    {
        std::cerr << "Unable to write file '" << argv[1] << "'." << std::endl;
    }

    exit(0);
}
