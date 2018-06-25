#include "Field.h"

Field::Field() {}

Field::Field(int res[]) {
    for (int i = 0; i < 3; i++)
        size[i] = res[i];
}

Field::~Field() {}

double Field::sample(double pos[3]) {
    return 0;
}