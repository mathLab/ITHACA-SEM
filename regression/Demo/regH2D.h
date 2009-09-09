#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "regDemo.h"

#ifndef REGH2D_H
#define REGH2D_H

using namespace std;
//header of class regDemo
class regH2D: public regDemo{

public:

	regH2D(std::string D, std::string bM);

private:

	std::string makeCommand(int);

};
#endif
