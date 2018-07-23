#include "Volume.h"
#include "Material.h"

Volume::Volume() {}

Volume::Volume(double bot[], double top[], int count[]) {
    for(int i = 0; i < 3; i++) {
        min[i] = bot[i];
        max[i] = top[i];
        res[i] = count[i];
    }
    size = (max[0] - min[0]) / (double)res[0];
    mat = new Material();
}

Volume::~Volume() {}

double Volume::sample(double pos[3], int field) {
    for(int i = 0; i < 3; i++) {
        pos[i] = (pos[i] - min[i]) / size;
    }
    return fields.at(field)->sample(pos);
}

string Volume::name(int field) {
    return field_names.at(field);
}

void Volume::addField(string name) {
    field_names.push_back(name);
    fields.push_back(new Field(res));
}

Material* Volume::getMat() {
    return mat;
}

double* Volume::getMin() {
    return min;
}

double* Volume::getMax() {
    return max;
}
void Volume::loadFireData(send_vol_data svd)
{
    field_names = svd.getVolNames();
    vector<volume_data> vd = svd.getVolData();
    for(size_t i = 0 ; i < vd.size();i++)
    {
        Field * f = new Field(vd[i].getDim());
        f->setValues(vd[i].getVolData());
        fields.push_back(f);
    }
    double * temp_bounds = svd.getBounds();
    cout <<"x" << temp_bounds[0] << "y" << temp_bounds[1] << "z" << temp_bounds[2] << "x" << temp_bounds[3] << "y" << temp_bounds[4] << "z" << temp_bounds[5] << endl;
    int * temp_res = svd.getReso();
    for(int i = 0 ; i < 3 ; i++)
    {
        min[i] = temp_bounds[i*2];
        max[i] = temp_bounds[i*2+1];
        res[i] = temp_res[i];
    }
    for(int i = 0 ; i < 3;i++)
        cout <<" min="<<min[i] << " max=" << max[i] << endl;
}

double Volume::getSize() {
    return size;
}

string Volume::toString() {
    string result = "";
    for(string name : field_names) {
        result += name;
        result += " ";
    }
    return result;
}
