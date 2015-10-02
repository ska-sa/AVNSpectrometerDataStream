#ifndef SPECTROMETER_HEADER
#define SPECTORMETER_HEADER

//System includes
#include <string>
#include <vector>

//Library includes

//Local includes
#include "SpetrometerDefinitions.h"

class cSpectrometerHeader
{
public:
    cSpectrometerHeader(const char* pcData);
    cSpectrometerHeader();
    ~cSpectrometerHeader();

    std::string         toString();
    void                deserialise(const char* pcData);
    std::vector<char>   serialise();

    //Accessors
    uint32_t            getSyncWord() const;
    int64_t             getTimestamp_us() const;
    uint8_t             getSubframeNumber() const;
    uint8_t             getNSubframes() const;
    uint16_t            getDigitiserType() const;
    bool                getNoiseDiodeOn() const;

    //Mutators
    void                setSyncWord(uint32_t u32SyncWord = AVN::Spectrometer::SYNC_WORD);
    void                setTimestamp_us(int64_t i64Timestamp);
    void                setSubframeNumber(uint8_t u8SubframeNumber);
    void                setNSubframes(uint8_t u8NSubframes);
    void                setDigitiserType(uint16_t u16DigitiserType);
    void                setNoiseDiodeOn(bool bNoiseDiodeOn);

    bool                isBigEndian();

private:
    uint32_t            m_u32SyncWord;
    int64_t             m_i64tTimestamp_us;
    uint8_t             m_u8SubframeNumber;
    uint8_t             m_u8NSubframes;
    uint16_t            m_u16DigitiserType;
    bool                m_bNoiseDiodeOn;

    bool                m_bIsBigEndian;
};

#endif //SPECTROMETER_HEADER
