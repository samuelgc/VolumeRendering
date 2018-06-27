


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

    // Data structure that will hold the info
    // Likely a matrix

};