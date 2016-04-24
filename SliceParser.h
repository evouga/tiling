#ifndef SLICEPARSER_H
#define SLICEPARSER_H

#include <vector>

class Slice;

void readSlicesFromFolder(const char *baseFilename, const char *objectname, std::vector<Slice *> &slices);

#endif
