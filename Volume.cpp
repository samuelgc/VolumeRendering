#include "Volume.h"

Volume::Volume() {}

Volume::~Volume() {}

double Volume::sample(double pos[3], int field) {
    return fields.at(field)->sample(pos);
}

string Volume::name(int field) {
    return field_names.at(field);
}

void Volume::addField(string name) {
    field_names.push_back(name);
    fields.push_back(new Field(res));
}