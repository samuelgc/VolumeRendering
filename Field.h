#ifndef FIELD_H
#define FIELD_H

#include <vector>

using namespace std;

class Field {

public:

    Field();

    Field(int res[]);

    ~Field();

    /**
     * Gives the value of this field at the specified index location
     *
     * @param pos - index of voxel to be sampled
     * @return
     */
    double sample(double pos[3]);

private:

    int size[3];
    vector<vector<vector<double>>> values;

};


#endif //FIELD_H
