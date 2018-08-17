#include "Field.h"
#include <iostream>
Field::Field() {}

Field::Field(int res[]) {
    for (int i = 0; i < 3; i++)
        size[i] = res[i];
}

Field::~Field() {}
double Field::interpolate(double pos[3])
{
    // cout << "enter here " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    // cout << "enter here " << size[0] << " " << size[1] << " " << size[2] << "\n";
    //Bottom 4 corners
    double val_0_0_0 = values[int(pos[0])][int(pos[1])][int(pos[2])];
    double val_1_0_0 = values[int(pos[0] + 1)][int(pos[1])][int(pos[2])];
    double val_0_1_0 = values[int(pos[0])][int(pos[1] + 1)][int(pos[2])];
    double val_1_1_0 = values[int(pos[0] + 1)][int(pos[1] + 1)][int(pos[2])];
    //top 4 corners
    double val_0_0_1 = values[int(pos[0])][int(pos[1])][int(pos[2] + 1)];
    double val_1_0_1 = values[int(pos[0] + 1)][int(pos[1])][int(pos[2] + 1)];
    double val_0_1_1 = values[int(pos[0])][int(pos[1] + 1)][int(pos[2] + 1)];
    double val_1_1_1 = values[int(pos[0] + 1)][int(pos[1] + 1)][int(pos[2] + 1)];
    //distance for inter x-dir
    double dis0_x = pos[0] - int(pos[0]);
    double disx_1 = int(pos[0] + 1) - pos[0];
    //distance for inter y-dir
    double dis0_y = pos[1] - int(pos[1]);
    double disy_1 = int(pos[1] + 1) - pos[1];
    //distance for inter z-dir
    double dis0_z = pos[2] - int(pos[2]);
    double disz_1 = int(pos[2] + 1) - pos[2];
    // cout << dis0_x << " " << disx_1 << " " << dis0_y << " " << disy_1 << " " << dis0_z << " " << disz_1 << " " <<endl;
    //Time to interpolate x-dir
        // bottom x's
    double x0 = val_0_0_0 * disx_1 + val_1_0_0 * dis0_x;
    double x1 = val_0_1_0 * disx_1 + val_1_1_0 * dis0_x;
        //tops x's
    double x2 = val_0_0_1 * disx_1 + val_1_0_1 * dis0_x;
    double x3 = val_0_1_1 * disx_1 + val_1_1_1 * dis0_x;
    //Time to interpolate y-dir
    double by = x0 * disy_1 + x1 * dis0_y; //bottom y's
    double ty = x2 * disy_1 + x3 * dis0_y; //tops y's
    //time to interpolate z-dir
    double value_final = by * disz_1 + ty * dis0_z;
    // cout << "leave here\n";
    return value_final;
}
double Field::sample(double pos[3]) {
    int x = pos[0];
    int y = pos[1];
    int z = pos[2];
    // cout << " x = " << x << " y = " << y << " z = " << z << endl; 
    if(x < 0 || y < 0 || z < 0){
        // cout << "out of bounds below" << endl;
        return -1;
    }
    if(x > size[0] -1 || y > size[1] -1 || z > size[2] -1) // Brian added -1 b/c if size is say 16 then 0 - 15 not 0 - 16
    {
        // cout << "out of bounds above" << endl;
        return -1;
    }
    if(x > size[0] -2 || y > size[1] -2 || z > size[2] -2) // Brian added -1 b/c if size is say 16 then 0 - 15 not 0 - 16
    {
        return values.at(x).at(y).at(z);
    }
    
    return interpolate(pos);
    // return values.at(x).at(y).at(z);
}
void Field::setValues(vector<vector<vector<double>>> v)
{
    values = v;
}