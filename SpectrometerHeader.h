#ifndef SPECTROMETER_HEADER
#define SPECTROMETER_HEADER

//System includes
#include <string>
#include <vector>

//Library includes

//Local includes
#include "SpectrometerDefinitions.h"

class cSpectrometerHeader
{
public:
    cSpectrometerHeader(const char* pcData);
    cSpectrometerHeader(const std::vector<char> &vcData);
    cSpectrometerHeader();
    ~cSpectrometerHeader();

    std::string         toString();
    bool                deserialise(const char* pcData);
    bool                deserialise(const std::vector<char> &vcData);
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

    bool                requiresEndianessFlip();

    bool                isValid();

private:
    uint32_t            m_u32SyncWord;
    int64_t             m_i64tTimestamp_us;
    uint8_t             m_u8SubframeNumber;
    uint8_t             m_u8NSubframes;
    uint16_t            m_u16DigitiserType;
    bool                m_bNoiseDiodeOn;

    bool                m_bFlipEndianess;

    bool                m_bValid;
};

#endif //SPECTROMETER_HEADER
