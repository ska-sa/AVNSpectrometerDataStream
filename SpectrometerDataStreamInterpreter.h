#ifndef SPECTROMETER_DATA_STREAM_INTERPRETER_H
#define SPECTROMETER_DATA_STREAM_INTERPRETER_H

//System includes
#include <vector>
#include <iostream>

#ifdef _WIN32
#include <stdint.h>

#ifndef int64_t
typedef __int64 int64_t;
#endif

#ifndef uint64_t
typedef unsigned __int64 uint64_t;
#endif

#else
#include <inttypes.h>
#endif

//Library includes
#ifndef Q_MOC_RUN //Qt's MOC and Boost have some issues don't let MOC process boost headers
#include <boost/thread/mutex.hpp>
#endif

//Local includes
#include "SpectrometerHeader.h"
#include "../../AVNAppLibs/SocketStreamers/SocketReceiverBase.h"

class cSpectrometerDataStreamInterpreter
{
public:
    cSpectrometerDataStreamInterpreter(boost::shared_ptr<cSocketReceiverBase>  pSocketReceiver);
    ~cSpectrometerDataStreamInterpreter();

    void                                    setUpdateRate(uint32_t u32UpdateRate_ms);

    bool                                    synchronise();
    bool                                    getNextFrame(float *pfChan0, float *pfChan1, float *pfChan2, float *pfChan3, uint32_t u32AllocateSize_nElements);
    bool                                    getNextFrame(int32_t *pi32Chan0, int32_t *pi32Chan1, int32_t *pi32Chan2, int32_t *pi32Chan3, uint32_t u32AllocateSize_nElements);

    const cSpectrometerHeader               &getLastHeader();

    uint32_t                                getNValuesPerChannelPerFrame();

    bool                                    isRunning();
    void                                    setIsRunning(bool bIsRunning);

private:
    boost::shared_ptr<cSocketReceiverBase>  m_pSocketReceiver;

    cSpectrometerHeader                     m_oCurrentHeader;
    cSpectrometerHeader                     m_oPreviousHeader;

    bool                                    m_bIsRunning;

    boost::shared_mutex                     m_oMutex;

    //Update frame rate. Output a frame every so many milliseconds.
    //Results in some data being discarded. 0 (default) implies no discarding
    uint32_t                                m_u32UpdateRate_ms;

    //Variables and vectors used interpretting:
    vector<char>                            m_vcPacket;

    uint32_t                                m_u32PacketSize_B;
    uint32_t                                m_u32NValuesPerFrame;
    uint32_t                                m_u32NValuesPerPacket;
    int64_t                                 m_i64LastUsedTimestamp_us;

    bool                                    m_bSkipFrame;

    //Inline functions

    inline bool                             headerConsistencyCheck(uint8_t u8ExpectedSubframeIndex)
    {
        //Some more consistency checks
        if(m_oCurrentHeader.getSubframeNumber() != u8ExpectedSubframeIndex)
        {
            std::cout << "Expected packet index " << (uint32_t)u8ExpectedSubframeIndex << ", got " << (uint32_t)m_oCurrentHeader.getSubframeNumber() << ". Resynchronising." << std::endl;
            return false;
        }

        if(m_oCurrentHeader.getNSubframes() != m_oPreviousHeader.getNSubframes())
        {
            std::cout << "Expected " << m_oPreviousHeader.getNSubframes() << " packets per frame, got " << (uint32_t)m_oCurrentHeader.getNSubframes() << ". Resynchronising." << std::endl;
            return false;
        }

        if(m_oCurrentHeader.getDigitiserType() != m_oPreviousHeader.getDigitiserType())
        {
            std::cout << "Expected plot type " << m_oPreviousHeader.getDigitiserType() << " , got " << m_oCurrentHeader.getDigitiserType() << ". Resynchronising." << std::endl;
            return false;
        }

        m_oPreviousHeader = m_oCurrentHeader;

        return true;
    }
};

#endif // SPECTROMETER_DATA_STREAM_INTERPRETER_H
