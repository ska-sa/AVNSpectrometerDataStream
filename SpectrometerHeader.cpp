
//System includes
#include <string>
#include <sstream>

//Library includes

//Local includes
#include "SpectrometerHeader.h"
#include "SpetrometerDefinitions.h"
#include "../../AVNUtilLibs/Timestamp/Timestamp.h"

cSpectrometerHeader::cSpectrometerHeader(const char*)
{

}

cSpectrometerHeader::cSpectrometerHeader():
    m_u32SyncWord(0),
    m_i64tTimestamp_us(0),
    m_u8SubframeNumber(0),
    m_u8NSubframes(0),
    m_u16DigitiserType(0xffff),
    m_bNoiseDiodeOn(false),
    m_bIsBigEndian(true)
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

void cSpectrometerHeader::deserialise(const char* pcData)
{
    if(*(uint32_t*)pcData == AVN::Spectrometer::SYNC_WORD)
    {
        m_bIsBigEndian = false;
    }
#ifdef _WIN32
    else if(_byteswap_ulong( *(uint32_t*)pcData ) == AVN::Spectrometer::SYNC_WORD)
#else
    else if(__builtin_bswap32( *(uint32_t*)pcData ) == AVN::Spectrometer::SYNC_WORD)
#endif
    {
        m_bIsBigEndian = true;
    }
    else
    {
        //SYNC Error
    }

    if(m_bIsBigEndian)
    {
        m_u32SyncWord = *(uint32_t*)(pcData);
        m_i64tTimestamp_us = *(int64_t*)(pcData + 4);
        m_u8SubframeNumber = *(uint8_t*)(pcData + 12);
        m_u8NSubframes = *(uint8_t*)(pcData + 13);
        m_u16DigitiserType = 0b0000000000001111 & *(uint16_t*)(pcData + 14);
        m_bNoiseDiodeOn = ( 0b1000000000000000 & *(uint16_t*)(pcData + 14) ) >> 15;
    }
    else
    {
#ifdef _WIN32
        m_u32SyncWord = _byteswap_ulong( *(uint32_t*)(pcData) );
        m_i64tTimestamp_us = _byteswap_uint64( *(int64_t*)(pcData + 4) );
        m_u8SubframeNumber = *(uint8_t*)(pcData + 12);
        m_u8NSubframes = *(uint8_t*)(pcData + 13);
        m_u16DigitiserType = 0b0000000000001111 & _byteswap_ushort( *(uint16_t*)(pcData + 14) );
        m_bNoiseDiodeOn = ( 0b1000000000000000 & _byteswap_ushort( *(uint16_t*)(pcData + 14) ) ) >> 15;
#else
        m_u32SyncWord = __builtin_bswap32( *(uint32_t*)(pcData) );
        m_i64tTimestamp_us = __builtin_bswap64( *(int64_t*)(pcData + 4) );
        m_u8SubframeNumber = *(uint8_t*)(pcData + 12);
        m_u8NSubframes = *(uint8_t*)(pcData + 13);
        m_u16DigitiserType = 0b0000000000001111 & __builtin_bswap16( *(uint16_t*)(pcData + 14) );
        m_bNoiseDiodeOn = ( 0b1000000000000000 & __builtin_bswap16( *(uint16_t*)(pcData + 14) ) ) >> 15;
#endif
    }
}

std::vector<char>   cSpectrometerHeader::serialise()
{

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

bool cSpectrometerHeader::isBigEndian()
{
    return m_bIsBigEndian;
}
