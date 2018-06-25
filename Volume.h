#include <vector>
#include <string>

#include "Field.h"

using namespace std;

class Volume {

public:

    Volume();

    ~Volume();

    double sample(double pos[], int field);

    string name(int field);

    void addField(string name);

private:

    // Fields
    vector<string> field_names;
    vector<Field*> fields;

    // Bounds
    double max[3];
    double min[3];
    int res[3];

};