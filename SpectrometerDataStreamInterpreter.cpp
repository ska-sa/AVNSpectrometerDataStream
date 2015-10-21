//System includes

//Library includes

//Local includes
#include "SpectrometerDataStreamInterpreter.h"

using namespace std;

cSpectrometerDataStreamInterpreter::cSpectrometerDataStreamInterpreter(boost::shared_ptr<cSocketReceiverBase>  pSocketReceiver) :
    m_pSocketReceiver(pSocketReceiver),
    m_bIsRunning(true),
    m_u32UpdateRate_ms(0),
    m_u32PacketSize_B(0),
    m_u32NValuesPerFrame(0),
    m_u32NValuesPerPacket(0),
    m_i64LastUsedTimestamp_us(0),
    m_bSkipFrame(false)
{
}

cSpectrometerDataStreamInterpreter::~cSpectrometerDataStreamInterpreter()
{
    setIsRunning(false);
}

void cSpectrometerDataStreamInterpreter::setUpdateRate(uint32_t u32UpdateRate_ms)
{
    //Thread safe mutator
    boost::upgrade_lock<boost::shared_mutex>  oLock(m_oMutex);
    boost::upgrade_to_unique_lock<boost::shared_mutex>  oUniqueLock(oLock);

    m_u32UpdateRate_ms = u32UpdateRate_ms;
}

bool cSpectrometerDataStreamInterpreter::isRunning()
{
    //Thread safe accessor
    boost::upgrade_lock<boost::shared_mutex>  oLock(m_oMutex);

    return m_bIsRunning;
}

void cSpectrometerDataStreamInterpreter::setIsRunning(bool bIsRunning)
{
    //Thread safe flag mutator
    boost::upgrade_lock<boost::shared_mutex>  oLock(m_oMutex);
    boost::upgrade_to_unique_lock<boost::shared_mutex>  oUniqueLock(oLock);

    m_bIsRunning = bIsRunning;
}

bool cSpectrometerDataStreamInterpreter::synchronise()
{
    setIsRunning(true);

    int32_t i32NextPacketSize_B = 0;

    do
    {
        if(!isRunning())
            return false;

        i32NextPacketSize_B = m_pSocketReceiver->getNextPacketSize_B(500);
    }
    while(i32NextPacketSize_B == -1);

    //Resize array as required
    if(m_vcPacket.size() != (uint32_t)i32NextPacketSize_B)
    {
        m_vcPacket.resize(i32NextPacketSize_B);
        m_u32PacketSize_B = i32NextPacketSize_B;
    }

    //Synchronise: Find the last packet of the frame
    cout << "Resynchronising to frame border." << endl;
    do
    {
        if(!isRunning())
            return false;

        m_pSocketReceiver->getNextPacket(&m_vcPacket.front(), 500);

        //Check that we synced to the stream correctly
        if(!m_oCurrentHeader.deserialise(m_vcPacket))
        {
            cout << "cSpectrometerDataStreamInterpreter::getNextFrame(): Warning got wrong magic no: "
                 << std::hex << m_oCurrentHeader.getSyncWord() << ". Expected " << AVN::Spectrometer::SYNC_WORD << std::dec << endl;

            return false;
        }

        cout << "cSpectrometerDataStreamInterpreter::getNextFrame(): Synchronising: Got packet " << (uint32_t)m_oCurrentHeader.getSubframeNumber()
             << " of " << (uint32_t)m_oCurrentHeader.getNSubframes() << endl;
    }
    //Keep going until the next subframe will be subframe #0
    while(m_oCurrentHeader.getSubframeNumber() != m_oCurrentHeader.getNSubframes() - 1);

    cout << "cSpectrometerDataStreamInterpreter::getNextFrame(): Synchronisation successful. Now aquiring data for plotting..." << endl;

    m_u32NValuesPerPacket = (i32NextPacketSize_B - AVN::Spectrometer::HEADER_SIZE_B) / sizeof(int32_t);
    m_u32NValuesPerFrame = m_u32NValuesPerPacket * m_oCurrentHeader.getNSubframes();

    //Store current header
    m_oPreviousHeader = m_oCurrentHeader;

    return true;
}

bool cSpectrometerDataStreamInterpreter::getNextFrame(int32_t *pi32Chan0, int32_t *pi32Chan1, int32_t *pi32Chan2, int32_t *pi32Chan3, uint32_t u32AllocateSize_nElements)
{
    setIsRunning(true);

    int32_t *pi32Data = NULL;
    int32_t i32NextPacketSize_B = 0;

    //Check plot vectors sizes
    //4 channels of data (L,R,Q,U or I0, Q0, I1, Q1)
    //There number of samples per channel is total values per frame / 4
    if(u32AllocateSize_nElements != m_u32NValuesPerFrame / 4)
    {
        return false;
    }

    //Read the data

    for(uint8_t u8ExpectedSubframeIndex = 0; u8ExpectedSubframeIndex < m_oCurrentHeader.getNSubframes(); u8ExpectedSubframeIndex++)
    {

        do
        {
            if(!isRunning())
                return false;

            i32NextPacketSize_B = m_pSocketReceiver->getNextPacketSize_B(500);
        }
        while(i32NextPacketSize_B == -1);

        if((unsigned)i32NextPacketSize_B != m_u32PacketSize_B)
        {
            cout << "cPlotsWidget::getDataThreadFunction(): Got different packet size, resynchronising." << endl;
            return false;
        }

        //Read the next packet
        while(!m_pSocketReceiver->getNextPacket(&m_vcPacket.front(), 500))
        {
            if(!isRunning())
                return false;
        }

        //Check for data consistency
        if(!headerConsistencyCheck(u8ExpectedSubframeIndex))
            return false;

        if(!m_oCurrentHeader.deserialise(m_vcPacket))
        {
            cout << "cSpectrometerDataStreamInterpreter::getNextFrame(): Warning got wrong magic no: "
                 << std::hex << m_oCurrentHeader.getSyncWord() << ". Expected " << AVN::Spectrometer::SYNC_WORD << std::dec << endl;
            return false;
        }


        //Get timestamp on the first subframe
        if(!u8ExpectedSubframeIndex)
        {
            //Attempt to reach 30 frames per second.
            boost::upgrade_lock<boost::shared_mutex>  oLock(m_oMutex); //For for update rate variable

            m_bSkipFrame = (m_oCurrentHeader.getTimestamp_us() - m_i64LastUsedTimestamp_us) < m_u32UpdateRate_ms;
        }

        if(m_bSkipFrame)
        {
            continue;
        }

        pi32Data = (int32_t*)(&m_vcPacket.front() + AVN::Spectrometer::HEADER_SIZE_B); //Go to offset of first sample (header is 16 bytes).

        //Deinterleave data to output channels of type int32
        if(m_oCurrentHeader.requiresEndianessFlip())
        {
            for(uint32_t u32ValueNo = 0; u32ValueNo < m_u32NValuesPerPacket / 4; u32ValueNo++)
            {
#ifdef _WIN32
                *pi32Chan0++ = (int32_t)( _byteswap_long( *pi32Data++ ) );
                *pi32Chan1++ = (int32_t)( _byteswap_long( *pi32Data++ ) );
                *pi32Chan2++ = (int32_t)( _byteswap_long( *pi32Data++ ) );
                *pi32Chan3++ = (int32_t)( _byteswap_long( *pi32Data++ ) );
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
            for(uint32_t u32ValueNo = 0; u32ValueNo < m_u32NValuesPerPacket / 4; u32ValueNo++)
            {
                *pi32Chan0++ =  (int32_t)( *pi32Data++ );
                *pi32Chan1++ =  (int32_t)( *pi32Data++ );
                *pi32Chan2++ =  (int32_t)( *pi32Data++ );
                *pi32Chan3++ =  (int32_t)( *pi32Data++ );
            }
        }
    }

    m_i64LastUsedTimestamp_us = m_oCurrentHeader.getTimestamp_us();

    return true;
}

bool cSpectrometerDataStreamInterpreter::getNextFrame(float *pfChan0, float *pfChan1, float *pfChan2, float *pfChan3, uint32_t u32AllocateSize_nElements)
{
    setIsRunning(true);

    int32_t *pi32Data = NULL;
    int32_t i32NextPacketSize_B = 0;

    //Check plot vectors sizes
    //4 channels of data (L,R,Q,U or I0, Q0, I1, Q1)
    //There number of samples per channel is total values per frame / 4
    if(u32AllocateSize_nElements != m_u32NValuesPerFrame / 4)
    {
        return false;
    }

    //Read the data

    for(uint8_t u8ExpectedSubframeIndex = 0; u8ExpectedSubframeIndex < m_oCurrentHeader.getNSubframes(); u8ExpectedSubframeIndex++)
    {

        do
        {
            if(!isRunning())
                return false;

            i32NextPacketSize_B = m_pSocketReceiver->getNextPacketSize_B(500);
        }
        while(i32NextPacketSize_B == -1);

        if((unsigned)i32NextPacketSize_B != m_u32PacketSize_B)
        {
            cout << "cPlotsWidget::getDataThreadFunction(): Got different packet size, resynchronising." << endl;
            return false;
        }

        //Read the next packet
        while(!m_pSocketReceiver->getNextPacket(&m_vcPacket.front(), 500))
        {
            if(!isRunning())
                return false;
        }

        if(!m_oCurrentHeader.deserialise(m_vcPacket))
        {
            cout << "cSpectrometerDataStreamInterpreter::getNextFrame(): Warning got wrong magic no: "
                 << std::hex << m_oCurrentHeader.getSyncWord() << ". Expected " << AVN::Spectrometer::SYNC_WORD << std::dec << endl;
            return false;
        }

        //Check for data consistency
        if(!headerConsistencyCheck(u8ExpectedSubframeIndex))
            return false;

        //Get timestamp on the first subframe
        if(!u8ExpectedSubframeIndex)
        {
            //Attempt to reach 30 frames per second.
            boost::upgrade_lock<boost::shared_mutex>  oLock(m_oMutex); //For for update rate variable

            m_bSkipFrame = (m_oCurrentHeader.getTimestamp_us() - m_i64LastUsedTimestamp_us) < m_u32UpdateRate_ms;
        }

        if(m_bSkipFrame)
        {
            continue;
        }

        pi32Data = (int32_t*)(&m_vcPacket.front() + AVN::Spectrometer::HEADER_SIZE_B); //Go to offset of first sample (header is 16 bytes).

        //Deinterleave data to output channels of type float
        if(m_oCurrentHeader.requiresEndianessFlip())
        {
            for(uint32_t u32ValueNo = 0; u32ValueNo < m_u32NValuesPerPacket / 4; u32ValueNo++)
            {
#ifdef _WIN32
                *pfChan0++ = (int32_t)( _byteswap_long( *pi32Data++ ) );
                *pfChan1++ = (int32_t)( _byteswap_long( *pi32Data++ ) );
                *pfChan2++ = (int32_t)( _byteswap_long( *pi32Data++ ) );
                *pfChan3++ = (int32_t)( _byteswap_long( *pi32Data++ ) );
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
            for(uint32_t u32ValueNo = 0; u32ValueNo < m_u32NValuesPerPacket / 4; u32ValueNo++)
            {
                *pfChan0++ =  (int32_t)( *pi32Data++ );
                *pfChan1++ =  (int32_t)( *pi32Data++ );
                *pfChan2++ =  (int32_t)( *pi32Data++ );
                *pfChan3++ =  (int32_t)( *pi32Data++ );
            }
        }
    }

    m_i64LastUsedTimestamp_us = m_oCurrentHeader.getTimestamp_us();

    return true;
}

const cSpectrometerHeader &cSpectrometerDataStreamInterpreter::getLastHeader()
{
    return m_oCurrentHeader;
}

uint32_t cSpectrometerDataStreamInterpreter::getNValuesPerChannelPerFrame()
{
    return m_u32NValuesPerFrame / 4;
}
