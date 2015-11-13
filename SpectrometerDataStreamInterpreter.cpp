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
    m_bSkipFrame(false),
    m_vviChannelData(4),
    m_bSynchronised(false)
{
}

cSpectrometerDataStreamInterpreter::cSpectrometerDataStreamInterpreter() :
    m_bIsRunning(true),
    m_u32UpdateRate_ms(0),
    m_u32PacketSize_B(0),
    m_u32NValuesPerFrame(0),
    m_u32NValuesPerPacket(0),
    m_i64LastUsedTimestamp_us(0),
    m_bSkipFrame(false),
    m_vviChannelData(4),
    m_bSynchronised(false)
{
    m_pSocketReceiver.reset(); //Not used
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
             cout << "cSpectrometerDataStreamInterpreter::synchronise(): Deserialising header failed, resynchronising." << endl;

            return false;
        }

        cout << "cSpectrometerDataStreamInterpreter::synchronise(): Synchronising: Got packet " << (uint32_t)m_oCurrentHeader.getSubframeNumber()
             << " of " << (uint32_t)m_oCurrentHeader.getNSubframes() << endl;
    }
    //Keep going until the next subframe will be subframe #0
    while(m_oCurrentHeader.getSubframeNumber() != m_oCurrentHeader.getNSubframes() - 1);

    cout << "cSpectrometerDataStreamInterpreter::synchronise(): Synchronisation successful." << endl;

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

    for(m_u8ExpectedSubframeIndex = 0; m_u8ExpectedSubframeIndex < m_oCurrentHeader.getNSubframes(); m_u8ExpectedSubframeIndex++)
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
            cout << "cSpectrometerDataStreamInterpreter::getNextFrame(): Got different packet size, resynchronising." << endl;
            return false;
        }

        //Read the next packet
        while(!m_pSocketReceiver->getNextPacket(&m_vcPacket.front(), 500))
        {
            if(!isRunning())
                return false;
        }

        //Check for data consistency
        if(!headerConsistencyCheck())
            return false;

        if(!m_oCurrentHeader.deserialise(m_vcPacket))
        {
            cout << "cSpectrometerDataStreamInterpreter::getNextFrame(): Deserialising header failed, resynchronising." << endl;
            return false;
        }


        //Get timestamp on the first subframe
        if(!m_u8ExpectedSubframeIndex)
        {
            //Attempt to reach 30 frames per second.
            boost::upgrade_lock<boost::shared_mutex>  oLock(m_oMutex); //For for update rate variable

            m_bSkipFrame = (m_oCurrentHeader.getTimestamp_us() - m_i64LastUsedTimestamp_us) < m_u32UpdateRate_ms;
        }

        if(m_bSkipFrame)
        {
            continue;
        }

        pi32Data = (int32_t*)((char*)(&m_vcPacket.front()) + AVN::Spectrometer::HEADER_SIZE_B); //Go to offset of first sample (header is 16 bytes).

        //Deinterleave data to output channels of type float
        deinterleaveInt32ToInt32(pi32Data, pi32Chan0, pi32Chan1, pi32Chan2, pi32Chan3, m_u32NValuesPerPacket / 4, m_oCurrentHeader.requiresEndianessFlip());
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

    for(m_u8ExpectedSubframeIndex = 0; m_u8ExpectedSubframeIndex < m_oCurrentHeader.getNSubframes(); m_u8ExpectedSubframeIndex++)
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
            cout << "cSpectrometerDataStreamInterpreter::getNextFrame(): Got different packet size, resynchronising." << endl;
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
            cout << "cSpectrometerDataStreamInterpreter::getNextFrame(): Deserialising header failed, resynchronising." << endl;
            return false;
        }

        //Check for data consistency
        if(!headerConsistencyCheck())
            return false;

        //Get timestamp on the first subframe
        if(!m_u8ExpectedSubframeIndex)
        {
            //Attempt to reach 30 frames per second.
            boost::upgrade_lock<boost::shared_mutex>  oLock(m_oMutex); //For for update rate variable

            m_bSkipFrame = (m_oCurrentHeader.getTimestamp_us() - m_i64LastUsedTimestamp_us) < m_u32UpdateRate_ms;
        }

        if(m_bSkipFrame)
        {
            continue;
        }


        pi32Data = (int32_t*)((char*)(&m_vcPacket.front()) + AVN::Spectrometer::HEADER_SIZE_B); //Go to offset of first sample (header is 16 bytes). (NB: Note cast to char before adding bytes to pointer.)

        //Deinterleave data to output channels of type float
        deinterleaveInt32ToFloat(pi32Data, pfChan0, pfChan1, pfChan2, pfChan3, m_u32NValuesPerPacket / 4, m_oCurrentHeader.requiresEndianessFlip());
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

void cSpectrometerDataStreamInterpreter::offloadData_callback(char* pcData, uint32_t u32Size_B)
{    
    int32_t *pi32Data = NULL;

    if(!m_bSynchronised)
    {

        //Synchronise: Find the last packet of the frame
        cout << "Resynchronising to frame border." << endl;
        do
        {
            if(!isRunning())
                return;

            //Check that data is large enough for at least the header
            if(u32Size_B < AVN::Spectrometer::HEADER_SIZE_B)
            {
                cout << "cSpectrometerDataStreamInterpreter::offloadData_callback(): Warning got packet size of " << u32Size_B << "which is smaller than AVN spectrometer header size." << endl;
                return;
            }

            //Check that we synced to the stream correctly
            if(!m_oCurrentHeader.deserialise(pcData))
            {
                cout << "cSpectrometerDataStreamInterpreter::offloadData_callback(): Deserialising header failed, resynchronising." << endl;

                return;
            }

            cout << "cSpectrometerDataStreamInterpreter::offloadData_callback(): Synchronising: Got packet " << (uint32_t)m_oCurrentHeader.getSubframeNumber()
                 << " of " << (uint32_t)m_oCurrentHeader.getNSubframes() << endl;
        }
        //Keep going until the next subframe will be subframe #0
        while(m_oCurrentHeader.getSubframeNumber() != m_oCurrentHeader.getNSubframes() - 1);

        cout << "cSpectrometerDataStreamInterpreter::offloadData_callback(): Synchronisation successful." << endl;

        m_u32NValuesPerPacket = (u32Size_B - AVN::Spectrometer::HEADER_SIZE_B) / sizeof(int32_t);
        m_u32NValuesPerFrame = m_u32NValuesPerPacket * m_oCurrentHeader.getNSubframes();
        m_u32PacketSize_B = u32Size_B;

        //Store current header
        m_oPreviousHeader = m_oCurrentHeader;

        m_bSynchronised = true;

        //Resize channel vectors as necessary
        if((uint32_t)m_vviChannelData[0].size() != m_u32NValuesPerFrame / 4)
        {
            for(uint32_t ui = 0; ui <  m_vviChannelData.size(); ui++)
            {
                m_vviChannelData[ui].resize(m_u32NValuesPerFrame / 4);
            }
        }

        m_u8ExpectedSubframeIndex = 0;

        return;
    }

    //Check packet size
    if(u32Size_B != m_u32PacketSize_B)
    {
        cout << "cSpectrometerDataStreamInterpreter::offloadData_callback(): Got different packet size, resynchronising." << endl;
        m_bSynchronised = false;
        return;
    }

    //Unpack and check this frames header
    if(!m_oCurrentHeader.deserialise(pcData))
    {
         cout << "cSpectrometerDataStreamInterpreter::offloadData_callback(): Deserialising header failed, resynchronising." << endl;

        m_bSynchronised = false;
        return;
    }

    //Check for data consistency
    if(!headerConsistencyCheck())
    {
        m_bSynchronised = false;
        return;
    }

    //Get timestamp on the first subframe
    if(!m_u8ExpectedSubframeIndex)
    {
        //Attempt to reach 30 frames per second.
        boost::upgrade_lock<boost::shared_mutex>  oLock(m_oMutex); //For for update rate variable

        m_bSkipFrame = (m_oCurrentHeader.getTimestamp_us() - m_i64LastUsedTimestamp_us) < m_u32UpdateRate_ms;

        //Also the the pointers to the beginning of the deinterleaved channel data
        m_pi32Chan0 = &(m_vviChannelData[0].front());
        m_pi32Chan1 = &(m_vviChannelData[1].front());
        m_pi32Chan2 = &(m_vviChannelData[2].front());
        m_pi32Chan3 = &(m_vviChannelData[3].front());
    }

    m_u8ExpectedSubframeIndex++;

    if(m_bSkipFrame)
    {
        return;
    }

    pi32Data = (int32_t*)(pcData + AVN::Spectrometer::HEADER_SIZE_B); //Go to offset of first sample (header is 16 bytes).

    //Deinterleave data to output channels of type float
    deinterleaveInt32ToInt32(pi32Data, m_pi32Chan0, m_pi32Chan1, m_pi32Chan2, m_pi32Chan3, m_u32NValuesPerPacket / 4, m_oCurrentHeader.requiresEndianessFlip());

    m_i64LastUsedTimestamp_us = m_oCurrentHeader.getTimestamp_us();

    //If all data has been received for this data frame pass on to the callback handler(s)
    if(m_oCurrentHeader.getSubframeNumber() == m_oCurrentHeader.getNSubframes() - 1)
    {
        boost::shared_lock<boost::shared_mutex> oLock(m_oCallbackHandlersMutex);
        for(uint32_t ui = 0; ui < m_vpCallbackHandlers.size(); ui++)
        {
            m_vpCallbackHandlers[ui]->getNextFrame_callback(m_vviChannelData[0],  m_vviChannelData[1], m_vviChannelData[2], m_vviChannelData[3], m_oCurrentHeader);
        }

        //The next subframe will then be the first of the new complete data frame
        m_u8ExpectedSubframeIndex = 0;
    }
}

void cSpectrometerDataStreamInterpreter::registerCallbackHandler(boost::shared_ptr<cCallbackInterface> pNewHandler)
{
    boost::unique_lock<boost::shared_mutex> oLock(m_oCallbackHandlersMutex);

    m_vpCallbackHandlers.push_back(pNewHandler);

    cout << "cSpectrometerDataStreamInterpreter::registerCallbackHandler(): Successfully registered callback handler: " << pNewHandler.get() << endl;
}

void cSpectrometerDataStreamInterpreter::deregisterCallbackHandler(boost::shared_ptr<cCallbackInterface> pHandler)
{
    boost::unique_lock<boost::shared_mutex> oLock(m_oCallbackHandlersMutex);
    bool bSuccess = false;

    //Search for matching pointer values and erase
    for(uint32_t ui = 0; ui < m_vpCallbackHandlers.size();)
    {
        if(m_vpCallbackHandlers[ui].get() == pHandler.get())
        {
            m_vpCallbackHandlers.erase(m_vpCallbackHandlers.begin() + ui);

            cout << "cSpectrometerDataStreamInterpreter::deregisterCallbackHandler(): Deregistered callback handler: " << pHandler.get() << endl;
            bSuccess = true;
        }
        else
        {
            ui++;
        }
    }

    if(!bSuccess)
    {
        cout << "cSpectrometerDataStreamInterpreter::deregisterCallbackHandler(): Warning: Deregistering callback handler: " << pHandler.get() << " failed. Object instance not found." << endl;
    }
}
