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

class cSpectrometerDataStreamInterpreter : public cSocketReceiverBase::cDataCallbackInterface
{
public:
    class cCallbackInterface
    {
    public:
        virtual void getNextFrame_callback(const std::vector<int> &vi32Chan0, const std::vector<int> &vi32Chan1, const std::vector<int> &vi32Chan2, std::vector<int> &vi32Chan3,
                                           const cSpectrometerHeader &oHeader) = 0;
    };

    cSpectrometerDataStreamInterpreter(boost::shared_ptr<cSocketReceiverBase>  pSocketReceiver);
    cSpectrometerDataStreamInterpreter();
    ~cSpectrometerDataStreamInterpreter();

    void                                                setUpdateRate(uint32_t u32UpdateRate_ms);

    bool                                                synchronise();
    bool                                                getNextFrame(float *pfChan0, float *pfChan1, float *pfChan2, float *pfChan3, uint32_t u32AllocateSize_nElements);
    bool                                                getNextFrame(int32_t *pi32Chan0, int32_t *pi32Chan1, int32_t *pi32Chan2, int32_t *pi32Chan3, uint32_t u32AllocateSize_nElements);

    const cSpectrometerHeader                           &getLastHeader();

    uint32_t                                            getNValuesPerChannelPerFrame();

    bool                                                isRunning();
    void                                                setIsRunning(bool bIsRunning);

    void                                                registerCallbackHandler(cCallbackInterface* pNewHandler);
    void                                                registerCallbackHandler(boost::shared_ptr<cCallbackInterface> pNewHandler);
    void                                                deregisterCallbackHandler(cCallbackInterface* pHandler);
    void                                                deregisterCallbackHandler(boost::shared_ptr<cCallbackInterface> pHandler);

    virtual void                                        offloadData_callback(char* pcData, uint32_t u32Size_B);

protected:
    boost::shared_ptr<cSocketReceiverBase>              m_pSocketReceiver;

    cSpectrometerHeader                                 m_oCurrentHeader;
    cSpectrometerHeader                                 m_oPreviousHeader;

    bool                                                m_bIsRunning;

    boost::shared_mutex                                 m_oMutex;

    //Update frame rate. Output a frame every so many milliseconds.
    //Results in some data being discarded. 0 (default) implies no discarding
    uint32_t                                            m_u32UpdateRate_ms;

    //Variables and vectors used interpretting:
    vector<char>                                        m_vcPacket;

    uint32_t                                            m_u32PacketSize_B;
    uint32_t                                            m_u32NValuesPerFrame;
    uint32_t                                            m_u32NValuesPerPacket;
    int64_t                                             m_i64LastUsedTimestamp_us;
    uint8_t                                             m_u8ExpectedSubframeIndex;

    bool                                                m_bSkipFrame;

    //Variables used for callback function mode
    std::vector< std::vector<int> >                     m_vviChannelData;

    int32_t                                             *m_pi32Chan0;
    int32_t                                             *m_pi32Chan1;
    int32_t                                             *m_pi32Chan2;
    int32_t                                             *m_pi32Chan3;

    bool                                                m_bSynchronised;

    std::vector<boost::shared_ptr<cCallbackInterface> > m_vpCallbackHandlers_shared;
    std::vector<cCallbackInterface*>                    m_vpCallbackHandlers; //For classes owning and instance of this class a regular pointer is simpler.
    boost::shared_mutex                                 m_oCallbackHandlersMutex;


    //Inline functions
    inline bool headerConsistencyCheck()
    {
        //Some more consistency checks
        if(m_oCurrentHeader.getSubframeNumber() != m_u8ExpectedSubframeIndex)
        {
            std::cout << "Expected packet index " << (uint32_t)m_u8ExpectedSubframeIndex << ", got " << (uint32_t)m_oCurrentHeader.getSubframeNumber() << ". Resynchronising." << std::endl;
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

    inline void deinterleaveInt32ToInt32(int32_t* pi32Data, int32_t *&pi32Chan0, int32_t *&pi32Chan1, int32_t *&pi32Chan2, int32_t *&pi32Chan3, uint32_t u32NValuesPerChan, bool bFlipEndianess)
    {
        if(bFlipEndianess)
        {
            for(uint32_t u32ValueNo = 0; u32ValueNo < u32NValuesPerChan; u32ValueNo++)
            {
#ifdef _WIN32
                *pi32Chan0++ = (int32_t)( _byteswap_ulong( *pi32Data++ ) );
                *pi32Chan1++ = (int32_t)( _byteswap_ulong( *pi32Data++ ) );
                *pi32Chan2++ = (int32_t)( _byteswap_ulong( *pi32Data++ ) );
                *pi32Chan3++ = (int32_t)( _byteswap_ulong( *pi32Data++ ) );
#else
                *pi32Chan0++ = (int32_t)( __builtin_bswap32( *pi32Data++ ) );
                *pi32Chan1++ = (int32_t)( __builtin_bswap32( *pi32Data++ ) );
                *pi32Chan2++ = (int32_t)( __builtin_bswap32( *pi32Data++ ) );
                *pi32Chan3++ = (int32_t)( __builtin_bswap32( *pi32Data++ ) );
#endif
            }
        }
        else
        {
            for(uint32_t u32ValueNo = 0; u32ValueNo < u32NValuesPerChan; u32ValueNo++)
            {
                *pi32Chan0++ =  (int32_t)( *pi32Data++ );
                *pi32Chan1++ =  (int32_t)( *pi32Data++ );
                *pi32Chan2++ =  (int32_t)( *pi32Data++ );
                *pi32Chan3++ =  (int32_t)( *pi32Data++ );
            }
        }
    }

    inline void deinterleaveInt32ToFloat(int32_t* pi32Data, float *&pfChan0, float *&pfChan1, float *&pfChan2, float *&pfChan3, uint32_t u32NValuesPerChan, bool bFlipEndianess)
    {
        if(bFlipEndianess)
        {
            for(uint32_t u32ValueNo = 0; u32ValueNo < u32NValuesPerChan; u32ValueNo++)
            {
#ifdef _WIN32
                *pfChan0++ = (int32_t)( _byteswap_ulong( *pi32Data++ ) );
                *pfChan1++ = (int32_t)( _byteswap_ulong( *pi32Data++ ) );
                *pfChan2++ = (int32_t)( _byteswap_ulong( *pi32Data++ ) );
                *pfChan3++ = (int32_t)( _byteswap_ulong( *pi32Data++ ) );
#else
                *pfChan0++ = (int32_t)( __builtin_bswap32( *pi32Data++ ) );
                *pfChan1++ = (int32_t)( __builtin_bswap32( *pi32Data++ ) );
                *pfChan2++ = (int32_t)( __builtin_bswap32( *pi32Data++ ) );
                *pfChan3++ = (int32_t)( __builtin_bswap32( *pi32Data++ ) );
#endif
            }
        }
        else
        {
            for(uint32_t u32ValueNo = 0; u32ValueNo < u32NValuesPerChan; u32ValueNo++)
            {
                *pfChan0++ = (int32_t)( *pi32Data++ );
                *pfChan1++ = (int32_t)( *pi32Data++ );
                *pfChan2++ = (int32_t)( *pi32Data++ );
                *pfChan3++ = (int32_t)( *pi32Data++ );
            }
        }
    }
};

#endif // SPECTROMETER_DATA_STREAM_INTERPRETER_H
