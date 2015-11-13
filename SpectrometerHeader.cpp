
//System includes
#include <string>
#include <sstream>
#include <iostream>

//Library includes

//Local includes
#include "SpectrometerHeader.h"
#include "SpectrometerDefinitions.h"
#include "../../AVNUtilLibs/Timestamp/Timestamp.h"
#include "ByteSwap.h" //For old versions of GCC there are some byte swap functions missing

cSpectrometerHeader::cSpectrometerHeader(const char* pcData)
{
    deserialise(pcData);
}

cSpectrometerHeader::cSpectrometerHeader(const std::vector<char> &vcData)
{
    deserialise(vcData);
}

cSpectrometerHeader::cSpectrometerHeader():
    m_u32SyncWord(0),
    m_i64tTimestamp_us(0),
    m_u8SubframeNumber(0),
    m_u8NSubframes(0),
    m_u16DigitiserType(0xffff),
    m_bNoiseDiodeOn(false),
    m_bFlipEndianess(true) //Assume true by default Roach is big endian, most PCs are little endian.
{
}

cSpectrometerHeader::~cSpectrometerHeader()
{
}

std::string cSpectrometerHeader::toString()
{
    std::stringstream oSS;
    oSS << "Sync word:           " << std::hex << m_u32SyncWord << std::dec << std::endl;
    oSS << "Timestamp:           " << AVN::stringFromTimestamp_full(m_i64tTimestamp_us) << std::endl;
    oSS << "Subframe numer:      " << (uint32_t)m_u8SubframeNumber << std::endl;
    oSS << "Number of subframes: " << (uint32_t)m_u8NSubframes << std::endl;
    oSS << "Digitiser type:      " << m_u16DigitiserType << " (" << AVN::Spectrometer::digitiserTypeToString(m_u16DigitiserType) << ")" << std::endl;
    oSS << "Noise diode enabled: " << m_bNoiseDiodeOn << std::endl;

    return oSS.str();
}

bool cSpectrometerHeader::deserialise(const char* pcData)
{
    if(*(uint32_t*)pcData == AVN::Spectrometer::SYNC_WORD)
    {
        m_bFlipEndianess = false;
    }
#ifdef _WIN32
    else if(_byteswap_ulong( *(uint32_t*)pcData ) == AVN::Spectrometer::SYNC_WORD)
#else
    else if(__builtin_bswap32( *(uint32_t*)pcData ) == AVN::Spectrometer::SYNC_WORD)
#endif
    {
        m_bFlipEndianess = true;
    }
    else
    {
        //SYNC Error
        m_bValid = false;

        std::cout << "cSpectrometerHeader::deserialise(): Warning got wrong magic no: "
             << std::hex << *(uint32_t*)pcData << ". Expected " << AVN::Spectrometer::SYNC_WORD << std::dec << std::endl;

        return false;
    }

    if(m_bFlipEndianess)
    {
#ifdef _WIN32
        m_u32SyncWord = _byteswap_ulong( *(uint32_t*)(pcData) );
        m_i64tTimestamp_us = _byteswap_uint64( *(int64_t*)(pcData + 4) );
        m_u8SubframeNumber = *(uint8_t*)(pcData + 12);
        m_u8NSubframes = *(uint8_t*)(pcData + 13);
        m_u16DigitiserType = 0b0000000000001111 & _byteswap_ushort( *(uint16_t*)(pcData + 14) );
        m_bNoiseDiodeOn = ( 0b1000000000000000 & _byteswap_ushort( *(uint16_t*)(pcData + 14) ) ) >> 15;
#else
        m_u32SyncWord = (uint32_t)( __builtin_bswap32( *(uint32_t*)(pcData) ) );
        m_i64tTimestamp_us = (int64_t)( __builtin_bswap64( *(int64_t*)(pcData + 4) ) );
        m_u8SubframeNumber = *(uint8_t*)(pcData + 12);
        m_u8NSubframes = *(uint8_t*)(pcData + 13);
        m_u16DigitiserType = 0b0000000000001111 & (uint16_t)( __builtin_bswap16( *(uint16_t*)(pcData + 14) ) );
        m_bNoiseDiodeOn = ( 0b1000000000000000 & (uint16_t)(__builtin_bswap16( *(uint16_t*)(pcData + 14) ) ) ) >> 15;
#endif
    }
    else
    {
        m_u32SyncWord = *(uint32_t*)(pcData);
        m_i64tTimestamp_us = *(int64_t*)(pcData + 4);
        m_u8SubframeNumber = *(uint8_t*)(pcData + 12);
        m_u8NSubframes = *(uint8_t*)(pcData + 13);
        m_u16DigitiserType = 0b0000000000001111 & *(uint16_t*)(pcData + 14);
        m_bNoiseDiodeOn = ( 0b1000000000000000 & *(uint16_t*)(pcData + 14) ) >> 15;
    }

    m_bValid = true;
    return true;
}

bool cSpectrometerHeader::deserialise(const std::vector<char> &vcData)
{
    //Here we can check for size
    if(vcData.size() < AVN::Spectrometer::HEADER_SIZE_B)
    {
        std::cout << "cSpectrometerHeader::deserialise() Error provided vector is shorter that AVN::Spectrometer head size of " << AVN::Spectrometer::HEADER_SIZE_B << " bytes" << std::endl;
        std::cout << "Got length: " << vcData.size() << " bytes. Returning." << std::endl;
        return false;

        //In future this should also throw an exception
    }

    return deserialise(&vcData.front());
}

std::vector<char> cSpectrometerHeader::serialise()
{
    //TODO: Place holder. Implement on demand

    std::vector<char> vcSerialisedHeader;
    return vcSerialisedHeader;
}

//Accessors
uint32_t cSpectrometerHeader::getSyncWord() const
{
    return m_u32SyncWord;
}

int64_t cSpectrometerHeader::getTimestamp_us() const
{
    return m_i64tTimestamp_us;
}

uint8_t cSpectrometerHeader::getSubframeNumber() const
{
    return m_u8SubframeNumber;
}

uint8_t cSpectrometerHeader::getNSubframes() const
{
    return m_u8NSubframes;
}

uint16_t cSpectrometerHeader::getDigitiserType() const
{
    return m_u16DigitiserType;
}

bool cSpectrometerHeader::getNoiseDiodeOn() const
{
    return m_bNoiseDiodeOn;
}

//Mutators
void cSpectrometerHeader::setSyncWord(uint32_t u32SyncWord)
{
    m_u32SyncWord = u32SyncWord;
}

void cSpectrometerHeader::setTimestamp_us(int64_t i64Timestamp)
{
    m_i64tTimestamp_us = i64Timestamp;
}

void cSpectrometerHeader::setSubframeNumber(uint8_t u8SubframeNumber)
{
    m_u8SubframeNumber = u8SubframeNumber;
}

void cSpectrometerHeader::setNSubframes(uint8_t u8NSubframes)
{
    m_u8NSubframes = u8NSubframes;
}

void cSpectrometerHeader::setDigitiserType(uint16_t u16DigitiserType)
{
    m_u16DigitiserType = u16DigitiserType;
}

void cSpectrometerHeader::setNoiseDiodeOn(bool bNoiseDiodeOn)
{
    m_bNoiseDiodeOn = bNoiseDiodeOn;
}

bool cSpectrometerHeader::requiresEndianessFlip()
{
    return m_bFlipEndianess;
}

bool cSpectrometerHeader::isValid()
{
    return m_bValid;
}
