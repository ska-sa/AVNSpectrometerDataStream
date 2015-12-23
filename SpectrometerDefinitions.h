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
    WB_SPECTROMETER_LRPP = 0,
    WB_SPECTROMETER_LRQU = 1,
    NB_SPECTROMETER_LRPP = 2,
    NB_SPECTROMETER_LRQU = 3,
    UNDEFINED = 0xffff
};

std::string digitiserTypeToString(uint16_t u16DigitiserType);

std::string getDigitiserChannelName(uint16_t u16DigitiserType, uint32_t u32ChanNo);

uint32_t getDigitiserFrameSize_nBins(uint16_t u16DigitiserType);


} // namespace Spectrometer

} // namespace AVN

#endif //SPECTROMETER_DEFINITIONS
