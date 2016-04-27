README
Instructions for the CellML-to-Nektar converter

============================

The translator is run as a python script from the command line:
If the path working directory is already set as the Cell_Model_Translator:
python modified_translate.py [path to cell model]

Example: "fj412:Cell_Model_Translator/ $ python modified_translate.py ../../../../../CellML_models_nektar/LuoRudy91_my.cellml"

If input with a filename.cellml it outputs a filename.cpp and a filename.h to be used with Nektar.






