#ifndef FIELD_H
#define FIELD_H

#include <vector>

using namespace std;

class Field {

public:

    Field();

    Field(int res[]);

    ~Field();
    double interpolate(double pos[3]);
    /**
     * Gives the value of this field at the specified index location
     *
     * @param pos - index of voxel to be sampled
     * @return
     */
    double sample(double pos[3]);
    /**
     * Used to set All the values in Values
     * @params v - all values to be set
     */
    void setValues(vector<vector<vector<double>>> v);

private:

    int size[3];
    vector<vector<vector<double>>> values;

};


#endif //FIELD_H
