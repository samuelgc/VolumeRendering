#include <iostream>

#include "Renderer.h"
#include "Volume.h"
#include "dataGetter/loadGeo2.0.h"

using namespace std;

// int main(int argc, char *argv[]) {
//     if (argc != 3) {
//         cout << "USAGE:\n"
//              << "   ./vrender [scene-descriptor.txt] [output-img.ppm]\n"
//              << "Error -- Incorrect number of arguments.\n";
//         return 1;
//     }

//     Renderer* vrender = new Renderer();
//     vrender->loadScene(argv[1]);
//     vrender->render();
//     vrender->write(argv[2]);

//     return 0;
// }
int main(int argc, char *argv[]) {
    // cout << "number of args " <<  argc << endl;
    // cout << "argv[0] = " << argv[0] << endl;
    string filename = argv[1];
    send_vol_data svd = getAllData(filename);
    Volume vol = Volume();
    vol.loadFireData(svd);
    double pos[] = {0,0,0};
    vol.sample(pos,0);
    // vector<string> field_names = svd.getVolNames();
    // vector<volume_data> vd = svd.getVolData();
    // // cout << "size " << field.size()  << " size " << vd.size() << endl;
    // // Renderer* vrender = new Renderer();
    // // vrender->loadScene(argv[1]);
    // // vrender->render();
    // // vrender->write(argv[2]);
    // Volume vol = Volume();
    // vol.setAllNames(field_names);

    return 0;
}