#include "Field.h"
#include <iostream>
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
    cout << " x = " << x << " y = " << y << " z = " << z << endl; 
    if(x < 0 || y < 0 || z < 0){
        // cout << "out of bounds below" << endl;
        return -1;
    }
    if(x > size[0] -1 || y > size[1] -1 || z > size[2] -1) // Brian added -1 b/c if size is say 16 then 0 - 15 not 0 - 16
    {
        // cout << "out of bounds above" << endl;
        return -1;
    }

    return values.at(x).at(y).at(z);
}
void Field::setValues(vector<vector<vector<double>>> v)
{
    values = v;
}