#ifndef VOLUME_H
#define VOLUME_H


#include <vector>
#include <string>

#include "Field.h"

using namespace std;

class Volume {

public:

    Volume();

    Volume(double bot[], double top[], int count[]);

    ~Volume();

    /**
     * Returns the value for the field specified and the requested location
     *
     * @param pos - location in 3D world space
     *     ** pos is converted to an index before sampling the field **
     * @param field - name of the field being sampled
     * @return
     */
    double sample(double pos[], int field);

    /**
     * Gives the name of the field associated with the given index
     *
     * @param field - index of the field
     * @return
     */
    string name(int field);

    /**
     * Adds a new field with the given name
     *
     * NOTE: All fields have the same resolution and size
     *       This is ensured by the min, max, and res variables
     *
     * @param name - name of the new field
     */
    void addField(string name);

private:

    // Fields
    vector<string> field_names;
    vector<Field*> fields;

    // Bounds
    double max[3];
    double min[3];
    int res[3];
    double size;

};


#endif //VOLUME_H
