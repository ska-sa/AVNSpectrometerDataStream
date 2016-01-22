#ifndef SPECTROMETER_HDF5_OUTPUT_FILE_H
#define SPECTROMETER_HDF5_OUTPUT_FILE_H

//System includes
#include <string>
#include <vector>

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

#include <vector>

//Library includes
extern "C" {
#include <hdf5.h>
}

#ifndef Q_MOC_RUN //Qt's MOC and Boost have some issues don't let MOC process boost headers
#include <boost/thread/shared_mutex.hpp>
#endif


//Local includes
#include "SpectrometerHeader.h"
#include "SpectrometerDefinitions.h"

class cSpectrometerHDF5OutputFile
{
    //Structures used for storing entries in HDF5 tables (NB: these must consist of POD! http://www.cplusplus.com/reference/type_traits/is_pod/)

    //Specific structures
    typedef struct cNoiseDiodeState
    {
        double              m_dTimestamp_s;
        char                m_caValue[1];
        char                m_caStatus[7];
    } cNoiseDiodeState;

    typedef struct cMotorTorques
    {
        double              m_dTimestamp_s;
        double              m_dAz0_Nm;
        double              m_dAz1_Nm;
        double              m_dEl0_Nm;
        double              m_dEl1_Nm;
        char                m_caStatus[7];
    } cMotorTorques;

    typedef struct cSourceSelection
    {
        double              m_dTimestamp_s;
        char                m_caSource[64];
        char                m_caStatus[7];
    } cSourceSelection;

    typedef struct cSizeOfFFTs
    {
        uint32_t            m_u32CoarseFFTSize;
        uint32_t            m_u32FineFFTSize;
        char                m_caStatus[7];
    } cSizeOfFFTs;

    typedef struct cAntennaStatus
    {
        double              m_dTimestamp_s;
        char                m_caAntennaStatus[8];
        char                m_caStatus[7];
    } cAntennaStatus;

    typedef struct cNoiseDiodeSource
    {
        double              m_dTimestamp_s;
        char                m_caSource[8];
        char                m_caStatus[7];
    } cNoiseDiodeSource;

    //General structures
    typedef struct cTimestampedInt
    {
        double              m_dTimestamp_s;
        int32_t             m_i32Value;
        char                m_caStatus[7];
    } cTimestampedInt;

    typedef struct cTimestampedUnsignedInt
    {
        double              m_dTimestamp_s;
        uint32_t            m_u32Value;
        char                m_caStatus[7];
    } cTimestampedUnsignedInt;

    typedef struct cTimestampedDouble
    {
        double              m_dTimestamp_s;
        double              m_dValue;
        char                m_caStatus[7];
    } cTimestampedDouble;

    typedef struct cTimestampedDualDouble
    {
        double              m_dTimestamp_s;
        double              m_dValue1;
        double              m_dValue2;
        char                m_caStatus[7];
    } cTimestampedDualDouble;




public:
    cSpectrometerHDF5OutputFile(const std::string &strFilename, AVN::Spectrometer::digitiserType eDigitiserType, uint32_t u32NFrequencyBins);
    ~cSpectrometerHDF5OutputFile();

    //Different functions can be called concurrently but any single function is not at present thread safe with regard to concurrent calls
    void                                    addFrame(const std::vector<int> &vi32Chan0, const std::vector<int> &vi32Chan1, const std::vector<int> &vi32Chan2, std::vector<int> &vi32Chan3,
                                                 const cSpectrometerHeader &oHeader);

    void                                    addRequestedAntennaAzEl(int64_t i64Timestamp_us, double dAzimuth_deg, double dElevation_deg);
    void                                    addActualAntennaAzEl(int64_t i64Timestamp_us, double dAzimuth_deg, double dElevation_deg);
    void                                    addActualSourceOffsetAzEl(int64_t i64Timestamp_us, double dAzimuthOffset_deg, double dElevationOffset_deg);
    void                                    addActualAntennaRADec(int64_t i64Timestamp_us, double dRighAscension_deg, double dDeclination_deg);

    void                                    addAntennaStatus(int64_t i64Timestamp_us, int32_t i32AntennaStatus, const std::string &strAntennaStatus);
    void                                    addMotorTorques(int64_t i64Timestamp_us, double dAz0_Nm, double dAz1_Nm, double dEl0_Nm, double dEl1_Nm);
    void                                    setAppliedPointingModel(const std::string &strModelName, const std::vector<double> &vdPointingModelParams);

    void                                    addNoiseDiodeSoftwareState(int64_t i64Timestamp_us, int32_t i32NoiseDiodeState);
    void                                    addNoiseDiodeSource(int64_t i64Timestamp_us, int32_t i32NoiseDiodeSource, const std::string &strNoiseSource);
    void                                    addNoiseDiodeCurrent(int64_t i64Timestamp_us, double dNoiseDiodeCurrent_A);

    void                                    addSourceSelection(int64_t i64Timestamp_us, const std::string &strSourceName, double dRighAscension_deg, double dDeclination_deg);

    void                                    addFrequencyRF(int64_t i64Timestamp_us, double dFreqencyRF_MHz);
    void                                    addFrequencyLOs(int64_t i64Timestamp_us, double dFrequencyLO1_MHz, double dFrequencyLO2_MHz);
    void                                    addBandwidthIF(int64_t i64Timestamp_us, double dBandwidthIF_MHz);

    void                                    addAccumulationLength(int64_t i64Timestamp_us, uint32_t u32NFrames);
    void                                    addNarrowBandChannelSelect(int64_t i64Timestamp_us, uint32_t u32ChannelNo);
    void                                    setFrequencyFs(double dFrequencyFs_MHz);
    void                                    setSizeOfFFTs(uint32_t u32CoarseSize_nSamp, uint32_t u32FineSize_nSamp);
    void                                    addCoarseFFTShiftMask(int64_t i64Timestamp_us, uint32_t u32ShiftMask);
    void                                    addAdcAttenuation(int64_t i64Timestamp_us, double dAttenuationChan0_dB, double dAttenuationChan1_dB);

private:
    std::string                             m_strFilename;
    AVN::Spectrometer::digitiserType        m_eDigitiserType;
    uint32_t                                m_u32NFrequencyBins;

    //HDF5:
    hid_t                                   m_iH5FileHandle;
    hid_t                                   m_iH5FileProperties;

    hid_t                                   m_iH5DataGroupHandle;
    hid_t                                   m_iH5MetaDataGroupHandle;
    hid_t                                   m_iH5SensorsGroupHandle;
    hid_t                                   m_iH5SensorsAntennasGroupHandle;
    hid_t                                   m_iH5SensorsRFEGroupHandle;
    hid_t                                   m_iH5SensorsDBEGroupHandle;
    hid_t                                   m_iH5SensorsAntennasAntenna1GroupHandle;
    hid_t                                   m_iH5ConfigurationGroupHandle;
    hid_t                                   m_iH5ConfigurationAntennasGroupHandle;
    hid_t                                   m_iH5ConfigurationAntennasAntenna1GroupHandle;

    hid_t                                   m_iH5DatasetVis;
    hid_t                                   m_iH5DatasetStokes;

    hsize_t                                 m_aChannelDatasetDims[3];
    hsize_t                                 m_aChannelDatasetExtensionDims[3];
    hsize_t                                 m_aChannelDataOffset[3];
    hsize_t                                 m_aMemspaceSize[3];

    //Values received / calculated from ROACH sample data stream
    std::vector<double>                     m_vdSampleDataTimestamps_s;
    std::vector<std::vector<float> >        m_vvfChannelAverages;
    std::vector<cNoiseDiodeState>           m_voROACHNoiseDiodeStateChanges;

    //Values received from station controller
    std::vector<cTimestampedDouble>         m_voRequestedAntennaAzs_deg;
    std::vector<cTimestampedDouble>         m_voRequestedAntennaEls_deg;
    std::vector<cTimestampedDouble>         m_voActualAntennaAzs_deg;
    std::vector<cTimestampedDouble>         m_voActualAntennaEls_deg;
    std::vector<cTimestampedDualDouble>     m_voActualSourceOffsetAzEls_deg;
    std::vector<cTimestampedDualDouble>     m_voActualAntennaRADecs_deg;

    std::vector<cAntennaStatus>             m_voAntennaStatuses;
    std::vector<cMotorTorques>              m_voMotorTorques_Nm;
    std::vector<double>                     m_vdPointingModelParams; //Only store most recent version
    std::string                             m_strPointModelName;

    std::vector<cTimestampedInt>            m_voNoiseDiodeSoftwareStates;
    std::vector<cNoiseDiodeSource>          m_voNoiseDiodeSources;
    std::vector<cTimestampedDouble>         m_voNoiseDiodeCurrents;

    std::vector<cSourceSelection>           m_voSelectedSources;

    std::vector<cTimestampedDouble>         m_voFrequenciesRF_MHz;
    std::vector<cTimestampedDualDouble>     m_voFrequenciesLOs_MHz;
    std::vector<cTimestampedDouble>         m_voBandwidthsIF_MHz;

    //Values received from ROACH TCPBorphServer
    std::vector<cTimestampedUnsignedInt>    m_voROACHAccumulationLengths_nFrames;
    std::vector<cTimestampedUnsignedInt>    m_voROACHNBChannelSelects;
    double                                  m_dROACHFrequencyFs_MHz; //Only store most recent version (shouldn't ever change from 800 MHz)
    cSizeOfFFTs                             m_oROACHSizesOfFFTs_nSamp; //Only store most recent version (shouldn't ever change for a given gateware)
    std::vector<cTimestampedUnsignedInt>    m_voROACHCoarseFFTShiftMasks;
    std::vector<cTimestampedDualDouble>     m_voROACHAdcAttenuations_dB;

    //Other
    cSpectrometerHeader                     m_oLastHeader;

    boost::shared_mutex                     m_oAppendDataMutex; //share-locked for augmenting all station controller and TCPBorph data. Unique locked when writing this data to file

    float                                   calculateFrameAverage(const std::vector<int32_t> &vi32ChannelData);

    void                                    addAttributeToDataSet(const std::string &strDescription, const std::string &strName, const std::string &strType, const std::string &strUnits, hid_t dataset);

    //Write logged data to file (after sample recording is completed)
    void                                    writeSampleDataTimestamps();
    void                                    writeChannelAverages();
    void                                    writeROACHNoiseDiodeStates();

    void                                    writeRequestedAntennaAzEls();
    void                                    writeActualAntennaAzEls();
    void                                    writeActualSourceOffsetAzEls();
    void                                    writeActualAntennaRADecs();

    void                                    writeAntennaStatuses();
    void                                    writeMotorTorques();
    void                                    writeAppliedPointingModel();

    void                                    writeNoiseDiodeSoftwareStates();
    void                                    writeNoiseDiodeSources();
    void                                    writeNoideDiodeCurrents();

    void                                    writeSelectedSources();

    void                                    writeRFFrequencies();
    void                                    writeLOFrequencies();
    void                                    writeIFBandwidths();

    void                                    writeROACHAccumulationLengths();
    void                                    writeROACHNBNarrowbandSelections();
    void                                    writeROACHSamplingFrequency();
    void                                    writeROACHSizeOfFFTs();
    void                                    writeROACHCoarseFFTShiftMask();
    void                                    writeROACHAdcAttentuations(); 

};

#endif // SPECTROMETER_HDF5_OUTPUT_FILE_H
