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
void runTest();
void test1();
void printFailed(string msg);
void printSuccess(string msg);
int main(int argc, char *argv[]) {
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
    runTest();
    return 0;
}

void runTest()
{
    test1();
}
void test1()
{
    cout << "Testing \033[1mhelp.json \033[0m\n";
    string testname_1 = "help.json";
    send_vol_data svd = getAllData(testname_1);
    double bot1[] = {0,0,0};
    double top1[] = {34,29,34};
    int res1[] = {35,30,35};
    Volume vol = Volume(bot1,top1,res1);
    vol.loadFireData(svd);
    double pos0[] = {0,0,0}; // first index
    double pos1[] = {1,0,0}; // second index
    double pos2[] = {16,0,0}; // first index 2nd tile
    double pos3[] = {31,15,15}; // last index 2nd tile
    double pos4[] = {34,29,34}; // last index all tiles
    double pos5[] = {0,0,-1}; // out of bounds below
    double pos6[] = {0,30,0}; // out of bounds above
    double pos7[] = {0,1,0}; // 16th index
    //Density Checking
    cout << "Testing Density \n";
    //Results
    double res0_0 = vol.sample(pos0,0);
    double res1_0 = vol.sample(pos1,0);
    double res2_0 = vol.sample(pos2,0);
    double res3_0 = vol.sample(pos3,0);
    double res4_0 = vol.sample(pos4,0);
    double res5_0 = vol.sample(pos5,0);
    double res6_0 = vol.sample(pos6,0);
    double res7_0 = vol.sample(pos7,0);
    //Expected
    double exp0_0 = 1.62603939e-18;
    double exp1_0 = 3.78290172e-21;
    double exp2_0 = 4.00170194e-07;
    double exp3_0 = 0.00473240018;
    double exp4_0 = 0;
    double exp5_0 = -1;
    double exp6_0 = -1;
    double exp7_0 = 6.06725192e-19;
    int density_good = 0;
    if(res0_0 == exp0_0) density_good++;
    else printFailed("Test 0.0 failed");
    if(res1_0 == exp1_0) density_good++;
    else printFailed("Test 0.1 failed");
    if(res2_0 == exp2_0) density_good++;
    else printFailed("Test 0.2 failed");
    if(res3_0 == exp3_0) density_good++;
    else printFailed("Test 0.3 failed");
    if(res4_0 == exp4_0) density_good++;
    else printFailed("Test 0.4 failed");
    if(res5_0 == exp5_0) density_good++;
    else printFailed("Test 0.5 failed");
    if(res6_0 == exp6_0) density_good++;
    else printFailed("Test 0.6 failed");
    if(res7_0 == exp7_0) density_good++;
    else printFailed("Test 0.7 failed");
    if(density_good == 8) printSuccess("8/8 test passed");
    else printFailed(to_string(density_good) + "/8 test passed");

    //testing Vel.x
    cout << "Testing Vel.x \n";
    //Results
    double res0_1 = vol.sample(pos0,1);
    double res1_1 = vol.sample(pos1,1);
    double res2_1 = vol.sample(pos2,1);
    double res3_1 = vol.sample(pos3,1);
    double res4_1 = vol.sample(pos4,1);
    double res5_1 = vol.sample(pos5,1);
    double res6_1 = vol.sample(pos6,1);
    double res7_1 = vol.sample(pos7,1);
    //Expected
    double exp0_1 = -0.193638742;
    double exp1_1 = -0.000301873311;
    double exp2_1 = 0.172829479;
    double exp3_1 = 0.309634209;
    double exp4_1 = -0.0597708039;
    double exp5_1 = -1;
    double exp6_1 = -1;
    double exp7_1 = -0.22961016;
    int vel_x_good = 0;
    if(res0_1 == exp0_1) vel_x_good++;
    else printFailed("Test 1.0 failed");
    if(res1_1 == exp1_1) vel_x_good++;
    else printFailed("Test 1.1 failed");
    if(res2_1 == exp2_1) vel_x_good++;
    else printFailed("Test 1.2 failed");
    if(res3_1 == exp3_1) vel_x_good++;
    else printFailed("Test 1.3 failed");
    if(res4_1 == exp4_1) vel_x_good++;
    else printFailed("Test 1.4 failed");
    if(res5_1 == exp5_1) vel_x_good++;
    else printFailed("Test 1.5 failed");
    if(res6_1 == exp6_1) vel_x_good++;
    else printFailed("Test 1.6 failed");
    if(res7_1 == exp7_1) vel_x_good++;
    else printFailed("Test 1.7 failed");
    if(vel_x_good == 8) printSuccess("8/8 test passed");
    else printFailed(to_string(vel_x_good) + "/8 test passed");


    //testing Fuel
    cout << "Testing Fuel \n";
    //Results
    double res0_6 = vol.sample(pos0,6);
    double res1_6 = vol.sample(pos1,6);
    double res2_6 = vol.sample(pos2,6);
    double res3_6 = vol.sample(pos3,6);
    double res4_6 = vol.sample(pos4,6);
    double res5_6 = vol.sample(pos5,6);
    double res6_6 = vol.sample(pos6,6);
    double res7_6 = vol.sample(pos7,6);
    //Expected // lol they are all zero XD
    double exp0_6 = 0;
    double exp1_6 = 0;
    double exp2_6 = 0;
    double exp3_6 = 0;
    double exp4_6 = 0;
    double exp5_6 = -1;
    double exp6_6 = -1;
    double exp7_6 = 0;
    int fuel_good = 0;
    if(res0_6 == exp0_6) fuel_good++;
    else printFailed("Test 6.0 failed");
    if(res1_6 == exp1_6) fuel_good++;
    else printFailed("Test 6.1 failed");
    if(res2_6 == exp2_6) fuel_good++;
    else printFailed("Test 6.2 failed");
    if(res3_6 == exp3_6) fuel_good++;
    else printFailed("Test 6.3 failed");
    if(res4_6 == exp4_6) fuel_good++;
    else printFailed("Test 6.4 failed");
    if(res5_6 == exp5_6) fuel_good++;
    else printFailed("Test 6.5 failed");
    if(res6_6 == exp6_6) fuel_good++;
    else printFailed("Test 6.6 failed");
    if(res7_6 == exp7_6) fuel_good++;
    else printFailed("Test 6.7 failed");
    if(fuel_good == 8) printSuccess("8/8 test passed");
    else printFailed(to_string(fuel_good) + "/8 test passed");
}

void printFailed(string msg)
{
    cout << "\033[31m" << msg << "\033[0m\n";
}
void printSuccess(string msg)
{
    cout << "\033[32m" << msg << "\033[0m\n";
}