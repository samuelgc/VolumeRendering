#include "Volume.h"

Volume::Volume() {}

Volume::Volume(double bot[], double top[], int count[]) {
    for(int i = 0; i < 3; i++) {
        min[i] = bot[i];
        max[i] = top[i];
        res[i] = count[i];
    }
    size = (max[0] - min[0]) / (double)res[0];
}

Volume::~Volume() {}

double Volume::sample(double pos[3], int field) {
    for(int i = 0; i < 3; i++) {
        pos[i] = (pos[i] - min[i]) / size;
    }
    return fields.at(field)->sample(pos);
}

string Volume::name(int field) {
    return field_names.at(field);
}

void Volume::addField(string name) {
    field_names.push_back(name);
    fields.push_back(new Field(res));
}