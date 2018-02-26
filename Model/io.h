#ifndef IO_H
#define IO_H
#include <fstream>
#include "statisticssampler.h"
class System;
using std::ofstream;

class IO
{
private:
    ofstream file;
public:
    IO(const char *filename);
    ~IO();

    void saveState(System &system, StatisticsSampler &statisticsSampler);
    void open(const char *filename);
    void close();

};


double * getPositions(System &system);


#endif
