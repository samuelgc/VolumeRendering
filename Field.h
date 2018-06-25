

using namespace std;

class Field {

public:

    Field();

    ~Field();

    double sample(double pos[3]);

private:

    double max[3];
    double min[3];
    int res[3];

    // Data structure that will hold the info
    // Likely a matrix

};