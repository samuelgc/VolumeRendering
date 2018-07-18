#include <iostream>
#include <stdlib.h>
#include <time.h>

#include "Renderer.h"

using namespace std;

int main(int argc, char *argv[]) {
    if (argc != 3) {
        cout << "USAGE:\n"
             << "   ./vrender [scene-descriptor.txt] [output-img.ppm]\n"
             << "Error -- Incorrect number of arguments.\n";
        return 1;
    }

    srand(time(NULL));

    Renderer* vrender = new Renderer();
    vrender->loadScene(argv[1]);
    vrender->render();
    vrender->write(argv[2]);

    return 0;
}