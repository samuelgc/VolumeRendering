#include "Volume.h"

Volume::Volume() {}

Volume::~Volume() {}

double Volume::sample(double pos[3], int field) {
    return fields.at(field)->sample(pos);
}