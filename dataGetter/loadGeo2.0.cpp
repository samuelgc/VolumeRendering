#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "dist/json/json.h"
#include <iomanip>
#include "volume_data.h"
#include "loadGeo2.0.h"

using namespace std;

int xDim = 0; // THe dim in x direction
int yDim = 0; // THe dim in y direction
int zDim = 0; // THe dim in z direction



send_vol_data getAllData(string filename)
{
    // cout << "LOADING JSON..." << endl;
    ifstream jsonFile(filename,ifstream::binary);
    // cout << "PARSING JSON" << endl;
    Json::Value root;
    jsonFile >> root; 
    // cout << "Getting VOLUME INFO" << endl;
    string volume_info =  root[11]["volume_summary"].asString();
    double bounds[6];
    getBounds(root[11],bounds);
    loadDimensions(volume_info); // get the dimensions
    vector<string> vDtypes = getNameMapping(root[11]["volume_summary"].asString()); // gets the different kind of data stored as voxels.
    Json::Value allVoxelData = root[21];
    vector<volume_data> allVoxData = getAllVolumeData(allVoxelData,vDtypes);
    // cout << "DATA READ" << endl;
    // cout << "SAVING DATA" << endl;
    int res2send[] = {xDim,yDim,zDim};
    send_vol_data svd = send_vol_data(allVoxData,vDtypes,res2send,bounds);
    // allVoxData[0].writeSlice();
    return svd;
}

void getBounds(Json::Value info,double bounds[6])
{
    bounds[0] = info["bounds"][0].asDouble();
    bounds[1] = info["bounds"][1].asDouble();
    bounds[2] = info["bounds"][2].asDouble();
    bounds[3] = info["bounds"][3].asDouble();
    bounds[4] = info["bounds"][4].asDouble();
    bounds[5] = info["bounds"][5].asDouble();
}

void setVolData(Json::Value tile,volume_data &v,int xOff,int yOff,int zOff,string compr)
{
    int x_start = xOff * 16;
    int y_start = yOff * 16;
    int z_start = zOff * 16;

    int xEnd = min(xDim,(xOff+1)*16);
    int yEnd = min(yDim,(yOff+1)*16);
    int zEnd = min(zDim,(zOff+1)*16);

    int xLen = xEnd - x_start;
    int yLen = yEnd - y_start;
    int zLen = zEnd - z_start;
    if(xLen > 16)
        xLen = 16;
    if(yLen > 16)
        yLen = 16;
    if(zLen > 16)
        zLen = 16;
    bool full = true; // if no compression is used.
    if(compr == "2")
    {
        full = false;
    }
    int count = 0;
    // cout << "({" << xLen <<"},{" << yLen << "},{" << zLen << "})" << endl;
    for(int lz = 0 ; lz < zLen; lz++){
        for(int ly = 0 ; ly < yLen; ly++){
            for(int lx = 0 ; lx < xLen; lx++){
                if(full)
                {
                    v.setValue(tile[count].asDouble(),(x_start +lx),(y_start +ly),(z_start +lz));
                    count++;
                }
                else
                    v.setValue(tile.asDouble(),(x_start + lx),(y_start + ly),(z_start + lz));
            }
        }
    }
}

volume_data getData4Volume(Json::Value vData)
{
    volume_data vd1 = volume_data(xDim,yDim,zDim);
    Json::Value comprNames = vData[2][1][3];
    Json:: Value tiles = vData[2][1][5];
    int tileCount = 0;
    for(int tz = 0 ; tz < (zDim + 15) / 16 ;tz++){
        for(int ty = 0 ; ty < (yDim + 15) / 16 ;ty++){
            for(int tx = 0 ; tx < (xDim + 15) / 16 ;tx++){
                Json::Value tile = tiles[tileCount];
                setVolData(tile[3],vd1,tx,ty,tz,tile[1].asString());
                tileCount++;
            }
        }
    }
    return vd1;
}

vector<volume_data>  getAllVolumeData(Json::Value allVoxelData,vector<string> vDtypes)
{   
    vector<volume_data> volData;
    
    for(unsigned i = 0 ; i < vDtypes.size() ; i++)
    {
        volData.push_back(getData4Volume(allVoxelData[(i * 2) + 1]));
    }
    return volData;
}

void print(string s)
{
    cout << s << endl;
}

void print(int i)
{
    cout << i << endl;
}

vector<string> getNameMapping(string s)
{
    size_t found = s.find_first_of("(");
    size_t found2 = s.find_first_of(")");
    vector<string> names; 
    while (found!=string::npos)
    {  
        // cout << s.substr(found + 1,found2 - found - 1) << endl;
        names.push_back(s.substr(found + 1,found2 - found - 1));
        s = s.substr(found2 + 1,string::npos);
        found = s.find_first_of("(");
        found2 = s.find_first_of(")");
    }
    return names;    
}

//gets the X,Y,Z Dimensions
void loadDimensions(string info)
{
    int indexB1 = -1; // means index Bracket 1
    int indexB2 = -1; // means index Bracket 2
    for(unsigned int i = 0 ; i < info.length() ;i++)
    {
        if(info[i] == '[')
            indexB1 = i;
        if(info[i] == ']'){
            indexB2 = i;
            break;
        }
    }
    string num =  info.substr(indexB1 + 1,indexB2 - indexB1 - 1);
    num.erase(remove(num.begin(), num.end(),' '), num.end());
    int indexC1 = -1;
    int indexC2 = -1;
    //used to find commas
    for(unsigned int i = 0 ; i < num.length() ;i++){
        if(num[i] == ',')
        {
            if(indexC1 == -1)
                indexC1 = i;
            else{
                indexC2 = i;
                break;
            }
        }
    }
    string x = num.substr(0,indexC1);
    string y = num.substr(indexC1 + 1,indexC2 - indexC1 -1);
    string z = num.substr(indexC2 + 1,num.length() - indexC2 - 1);
    //convert to Integars
    try{
        xDim = stoi(x);
        yDim = stoi(y);
        zDim = stoi(z);
    }
    catch(exception& e)
    {
        cout << x << " " << y << " " << z <<  endl;
        cout << e.what() << endl;
    }
}