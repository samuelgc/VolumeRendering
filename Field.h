

using namespace std;

class Field {

public:

    Field();

    Field(int res[]);

    ~Field();

    double sample(double pos[3]);

private:

    int size[3];

    // Data structure that will hold the info
    // Likely a matrix

};