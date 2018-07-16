#ifndef LOAD_GEO2_H
#define LOAD_GEO2_H

#include <string>
#include "dist/json/json.h"
#include "volume_data.h"
#include <vector>
#include <string>

class send_vol_data
{
    private:
        std::vector<volume_data> vol_dat;
        std::vector<std::string> vol_names;
        int res[3];
        double bounds[6];
    public:
        send_vol_data(std::vector<volume_data> vd , std::vector<std::string> vn, int resolution[3],double b[6])
        {
            vol_dat = vd;
            vol_names = vn;
            res[0] = resolution[0];
            res[1] = resolution[1];
            res[2] = resolution[2];
            bounds[0] = b[0];
            bounds[1] = b[1];
            bounds[2] = b[2];
            bounds[3] = b[3];
            bounds[4] = b[4];
            bounds[5] = b[5];
        }
        std::vector<volume_data> getVolData()
        {
            return vol_dat;
        }
        std::vector<std::string> getVolNames()
        {
            return vol_names;
        }
        int* getReso()
        {
            return res;
        }
        double* getBounds()
        {
            return bounds;
        }
};

//used for debugging
void print(std::string s);
//used for debugging
void print(int i);
//yay
send_vol_data getAllData(std::string filename);
void getBounds(Json::Value info,double bounds[6]);
/**
 * used to read and set the x,y,z dimensions
 * @param string info contains the line from the json file with key "volume_summary" 
 */
void loadDimensions(std::string info);
/**
 * used to get a vector with the name for each data type position in vector aligns with position vector<volume_data>
 * @param string info contains the line from the json file with key "volume_summary" 
 * @return vector<string> where string is the name of each data type e.g. density, fuel, etc.
 */
std::vector<std::string> getNameMapping(std::string info);
/**
 * @params Value allVoxelData the raw form of all stored data
 * @params vector<string> vDtype the name for each different data stored. 
 * @return vector<volume_data> where the vector contains the data for each different data type stored
 */
std::vector<volume_data> getAllVolumeData(Json::Value allVoxelData, std::vector<std::string> vDtypes);
/**  
 * @params Values vData -- contains all data for one volume type
 * return a volume_data with all needed informaiton
 */
volume_data getData4Volume(Json::Value vData);
/**
 * @parama Value tile -- a 16x16x16 (if full size) containing the data to be loaded
 * @param volume_data the volume in to which the data will be saved
 * @param xOff -- the off set of which tile in the x direction
 * @param yOff -- the off set of which tile in the y direction
 * @param zOff -- the off set of which tile in the z direction 
 * @param string compr -- If the tile is compressed or full etc...
 */
void setVolData(Json::Value tile, volume_data &v, int xOff, int yOff, int zOff, std::string compr);

#endif