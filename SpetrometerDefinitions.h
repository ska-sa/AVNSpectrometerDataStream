#ifndef SPECTROMETER_DEFINITIONS
#define SPECTROMETER_DEFINITIONS

//System includes
#include <string>
#include <inttypes.h>

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

}

}

#endif //SPECTROMETER_DEFINITIONS
