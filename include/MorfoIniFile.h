/*********************************************************************************
 * 
 *                  MorfoIniFile.h (Morfo70)
 * 
 *  Fucntions to read ini file and stores in OptionsParameters
 *  Implementation in: 
 *      - MorfoIniFile.cpp 
 *   Implementation with OptionsParametersConstants.h/cpp
 * 
 * *******************************************************************************/
#ifndef MORFOIO_H_
#define MORFOIO_H_

#include <fstream>
#include "OptionsParametersConstants.h"
#include "Morfo70.h"

using namespace std;

namespace MorfoIO{

/////////////////////////////////////////////////////////////////////////////////////
/// String functions to trim an compare
/////////////////////////////////////////////////////////////////////////////////////

/**
 *  Left trim
 *  @param str string to trim
 * */
inline std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    str.erase(0, str.find_first_not_of(chars));
    return str;
}
 
/**
 *  Right trim
 *  @param str string to trim
 * */
inline std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    str.erase(str.find_last_not_of(chars) + 1);
    return str;
}


/**
 *  Trim
 *  @param str string to trim
 * */
inline std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    return ltrim(rtrim(str, chars), chars);
}

/**
 *  Compare the first n =length(b) characters of strings  a and b, return the first
 *  vector position of a without b
 *  @param a larger string
 *  @param b small string
 * */
inline int compare(string a, string b)
{
    if (a.compare(0,b.length(),b)==0)
        return b.length();
    else return 0;
}
/////////////////////////////////////////////////////////////////////////////////////
/// END String functions to trim an compare
/////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////
/// Read Ini File (self understandble functions)
/////////////////////////////////////////////////////////////////////////////////////

/**
 *  Load Init file and fill Parameters and Options variables
 *  @param iniFile .ini file
 * */ 
 extern void loadIniFile(string IniFile);
 void readBottomParameter(string parameter);
 void readBottomOption(string option);
 void readOffshoreData(string line);
 
 void readWaveOption(string option);
 void readWaveParameter(string parameter);
 void readRollerOption(string option);
 void readRollerParameter(string parameter);
 void readHydroOption(string parameter);
 void readHydroParameter(string parameter);
 
 void readSedimentOption(string option);
 void readSedimentCWSParameter(string parameter);
 void readSedimentSwashParameter(string parameter);
 
 void readRunParameter(string parameter);
 void readRunOption(string options);
/////////////////////////////////////////////////////////////////////////////////////
/// END read Init file
/////////////////////////////////////////////////////////////////////////////////////

 
}
#endif