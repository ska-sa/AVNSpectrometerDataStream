#ifndef SPECTROMETER_DEFINITIONS
#define SPECTROMETER_DEFINITIONS

//System includes
#include <string>
#include <inttypes.h>
#include <vector>

//Library includes

//Local includes

namespace AVN
{

namespace Spectrometer
{

static const uint32_t HEADER_SIZE_B = 16;
static const uint32_t SYNC_WORD = 0x1a2b3c4d;

enum digitiserType
{
    WB_SPECTROMETER_CFFT = 0,
    WB_SPECTROMETER_LRQU = 1,
    NB_SPECTROMETER_CFFT = 2,
    NB_SPECTROMETER_LRQU = 3,
    UNDEFINED = 0xffff
};

static std::string digitiserTypeToString(uint16_t u16DigitiserType)
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

static std::string getDigitiserChannelName(uint16_t u16DigitiserType, uint32_t u32ChanNo)
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

static uint32_t getDigitiserFrameSize_nBins(uint16_t u16DigitiserType)
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

} // namespace Spectrometer

} // namespace AVN

#endif //SPECTROMETER_DEFINITIONS
