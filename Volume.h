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

private:

    // Fields
    vector<string> field_names;
    vector<Field*> fields;

};