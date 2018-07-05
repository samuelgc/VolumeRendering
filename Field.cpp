#include "Field.h"

Field::Field() {}

Field::Field(int res[]) {
    for (int i = 0; i < 3; i++)
        size[i] = res[i];
}

Field::~Field() {}

double Field::sample(double pos[3]) {
    int x = pos[0];
    int y = pos[1];
    int z = pos[2];

    if(x < 0 || y < 0 || z < 0)
        return -1;
    if(x > size[0] || y > size[1] || z > size[2])
        return -1;

    return values.at(x).at(y).at(z);
}