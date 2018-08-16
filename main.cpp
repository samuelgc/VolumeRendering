#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

#include "Renderer.h"

using namespace std;

int main(int argc, char *argv[]) {
    if (argc != 3) {
        cout << "USAGE:\n"
             << "   ./vrender [scene-descriptor.txt] [output-img.ppm] [samples-64]\n"
             << "Error -- Incorrect number of arguments.\n";
        return 1;
    }

    // srand(time(NULL));
    srand(1);
    Renderer* vrender = new Renderer();
    cout << "loading Scene\n";
    vrender->loadScene(argv[1]);
    cout << "Starting Render\n";
    vrender->render(64);//128.0);//argv[3]);

    vrender->write(argv[2]);

    return 0;
}

// int main(int argc, char *argv[]) {
    // cout << "number of args " <<  argc << endl;
    // cout << "argv[0] = " << argv[0] << endl;
    // string filename = argv[1];
    // send_vol_data svd = getAllData(filename);
    // double bot1[] = {0,0,0};
    // double top1[] = {34,29,34};
    // int res1[] = {35,30,35};
    // Volume vol = Volume(bot1,top1,res1);
    // vol.loadFireData(svd);
    // double pos[] = {16,0,0};
    // cout << "Field " << vol.name(0) << " value " << vol.sample(pos,0) << endl;
    // cout << "pos[" << pos[0] << "," << pos[1] << "," << pos[2] << "]" << endl;
    // cout << "Min = " << vol.getMin()[0] << "," << vol.getMin()[1] << "," << vol.getMin()[2] << endl;
    // vector<string> field_names = svd.getVolNames();
    // vector<volume_data> vd = svd.getVolData();
    // // cout << "size " << field.size()  << " size " << vd.size() << endl;
    // // Renderer* vrender = new Renderer();
    // // vrender->loadScene(argv[1]);
    // // vrender->render();
    // // vrender->write(argv[2]);
    // Volume vol = Volume();
    // vol.setAllNames(field_names);
//     runTest();
//     return 0;
// }

