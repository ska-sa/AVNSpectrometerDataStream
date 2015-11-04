
//System includes

//Library includes

//Local includes
#include "SpectrometerDefinitions.h"

using namespace AVN::Spectrometer;
using namespace std;

string AVN::Spectrometer::digitiserTypeToString(uint16_t u16DigitiserType)
{
    std::string strDigitiserType;

    switch(u16DigitiserType)
    {
    case WB_SPECTROMETER_CFFT:
        strDigitiserType = std::string("WB spectrometer complex FFT");
        break;

    case WB_SPECTROMETER_LRQU:
        strDigitiserType = std::string("WB spectrometer left and right power, Stokes Q and U");
        break;

    case NB_SPECTROMETER_CFFT:
        strDigitiserType = std::string("NB spectrometer complex FFT");
        break;

    case NB_SPECTROMETER_LRQU:
        strDigitiserType = std::string("NB spectrometer left and right power, Stokes Q and U");
        break;

    default:
        strDigitiserType = std::string("Undefined");
        break;
    }

    return strDigitiserType;
}

string AVN::Spectrometer::getDigitiserChannelName(uint16_t u16DigitiserType, uint32_t u32ChanNo)
{
    std::string strDigitiserChannelName;

    switch(u16DigitiserType)
    {
    case WB_SPECTROMETER_LRQU:
    case NB_SPECTROMETER_LRQU:
        switch(u32ChanNo)
        {
        case 0: strDigitiserChannelName = std::string("Left Power"); break;
        case 1: strDigitiserChannelName = std::string("Right Power"); break;
        case 2: strDigitiserChannelName = std::string("Stoke Q"); break;
        case 3: strDigitiserChannelName = std::string("Stoke U"); break;
        default: strDigitiserChannelName = std::string("Undefined"); break;
        }
        break;

    case WB_SPECTROMETER_CFFT:
    case NB_SPECTROMETER_CFFT:
        switch(u32ChanNo)
        {
        case 0: strDigitiserChannelName = std::string("Channel 1 I"); break;
        case 1: strDigitiserChannelName = std::string("Channel 1 Q"); break;
        case 2: strDigitiserChannelName = std::string("Channel 2 I"); break;
        case 3: strDigitiserChannelName = std::string("Channel 2 Q"); break;
        default: strDigitiserChannelName = std::string("Undefined"); break;
        }
        break;

    default:
        strDigitiserChannelName = std::string("Undefined");
        break;
    }

    return strDigitiserChannelName;
}

uint32_t AVN::Spectrometer::getDigitiserFrameSize_nBins(uint16_t u16DigitiserType)
{
    uint32_t u32DigitiserFrameSize_nBins;

    switch(u16DigitiserType)
    {
    case WB_SPECTROMETER_CFFT:
        u32DigitiserFrameSize_nBins = 1024;
        break;

    case WB_SPECTROMETER_LRQU:
        u32DigitiserFrameSize_nBins = 1024;
        break;

    case NB_SPECTROMETER_CFFT:
        u32DigitiserFrameSize_nBins = 4096;
        break;

    case NB_SPECTROMETER_LRQU:
        u32DigitiserFrameSize_nBins = 4096;
        break;

    default:
        u32DigitiserFrameSize_nBins = 0;
        break;
    }

    return u32DigitiserFrameSize_nBins;
}
