#include <iostream>

#include "Renderer.h"

using namespace std;

int main(int argc, char *argv[]) {
    if (argc != 3) {
        cout << "USAGE:\n"
             << "   ./vrender [scene-descriptor.txt] [output-img.ppm]\n"
             << "Error -- Incorrect number of arguments.\n";
        return 1;
    }

    Renderer* vrender = new Renderer();
    vrender->newScene();
    vrender->render();
    vrender->write(argv[2]);

    return 0;
}