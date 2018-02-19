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
    typedef struct cMarkupLabels
    {
        double              m_dTimestamp_s;
        char                m_chaLabel[64];
    } cMarkupLabels;

    typedef struct cTimestampedChar
    {
        double              m_dTimestamp_s;
        char                m_chaValue[1];
        char                m_chaStatus[7];
    } cTimestampedChar;

    typedef struct cTimestampedBool
    {
        double              m_dTimestamp_s;
        int32_t             m_i32Value;
        char                m_chaStatus[7];
    } cTimestampedBool;

    typedef struct cSourceSelection
    {
        double              m_dTimestamp_s;
        char                m_chaSource[128];
        char                m_chaStatus[7];
    } cSourceSelection;

    typedef struct cAntennaStatus
    {
        double              m_dTimestamp_s;
        char                m_chaAntennaStatus[16];
        char                m_chaStatus[7];
    } cAntennaStatus;

    typedef struct cNoiseDiodeSource
    {
        double              m_dTimestamp_s;
        char                m_chaSource[8];
        char                m_chaStatus[7];
    } cNoiseDiodeSource;

    //General structures
    typedef struct cTimestampedInt
    {
        double              m_dTimestamp_s;
        int32_t             m_i32Value;
        char                m_chaStatus[7];
    } cTimestampedInt;

    typedef struct cTimestampedUnsignedInt
    {
        double              m_dTimestamp_s;
        uint32_t            m_u32Value;
        char                m_chaStatus[7];
    } cTimestampedUnsignedInt;

    typedef struct cTimestampedDouble
    {
        double              m_dTimestamp_s;
        double              m_dValue;
        char                m_chaStatus[7];
    } cTimestampedDouble;

    typedef struct cTimestampedDualDouble
    {
        double              m_dTimestamp_s;
        double              m_dValue1;
        double              m_dValue2;
        char                m_chaStatus[7];
    } cTimestampedDualDouble;

    typedef struct cAntennaConfiguration
    {
        char                m_chaAntennaName[64];
        char                m_chaAntennaDiameter_m[8];
        char                m_chaAntennaBeamwidth[8];
        char                m_chaAntennaLongitude_deg[16];
        char                m_chaAntennaLatitude_deg[16];
        char                m_chaAntennaAltitude_m[8];
    } cAntennaConfiguration;


public:
    cSpectrometerHDF5OutputFile(const std::string &strFilename, AVN::Spectrometer::digitiserType eDigitiserType, uint32_t u32NFrequencyBins);
    ~cSpectrometerHDF5OutputFile();

    //Different functions can be called concurrently but any single function is not at present thread safe with regard to concurrent calls
    void                                    addFrame(const std::vector<int> &vi32Chan0, const std::vector<int> &vi32Chan1, const std::vector<int> &vi32Chan2, std::vector<int> &vi32Chan3,
                                                     const cSpectrometerHeader &oHeader);

    void                                    addMarkupLabel(int64_t i64Timestamp_us, const std::string &strLabel);

    // Antenna-space values
    void                                    addAcsRequestedAz(int64_t i64Timestamp_us, double dAzimuth_deg, const std::string &strStatus);
    void                                    addAcsRequestedEl(int64_t i64Timestamp_us, double dElevation_deg, const std::string &strStatus);
    void                                    addAcsActualAz(int64_t i64Timestamp_us, double dAzimuth_deg, const std::string &strStatus);
    void                                    addAcsActualEl(int64_t i64Timestamp_us, double dElevation_deg, const std::string &strStatus);

    // Sky-space values
    void                                    addSkyRequestedAz(int64_t i64Timestamp_us, double dAzimuth_deg, const std::string &strStatus);
    void                                    addSkyRequestedEl(int64_t i64Timestamp_us, double dElevation_deg, const std::string &strStatus);
    void                                    addSkyActualAz(int64_t i64Timestamp_us, double dAzimuth_deg, const std::string &strStatus);
    void                                    addSkyActualEl(int64_t i64Timestamp_us, double dElevation_deg, const std::string &strStatus);

    void                                    addPointingModelParameter(uint8_t ui8ParameterNumber, double dParameterValue);

    void                                    addAntennaStatus(int64_t i64Timestamp_us, const std::string &strAntennaStatus, const std::string &strStatus);

    void                                    addNoiseDiodeInputSource(int64_t i64Timestamp_us, const std::string &strNoiseDiodeInputSource, const std::string &strStatus);
    void                                    addNoiseDiodeEnable(int64_t i64Timestamp_us, bool bNoiseDiodeEnable, const std::string &strStatus);
    void                                    addNoiseDiodeSelect(int64_t i64Timestamp_us, int32_t i32NoiseDiodeSelect, const std::string &strStatus);
    void                                    addNoiseDiodePWMMark(int64_t i64Timestamp_us, int32_t i32NoiseDiodePWMMark, const std::string &strStatus);
    void                                    addNoiseDiodePWMFrequency(int64_t i64Timestamp_us, double dNoiseDiodePWMFrequency, const std::string &strStatus);



    void                                    addSourceSelection(int64_t i64Timestamp_us, const std::string &strSourceName, const std::string &strStatus);

    void                                    addFrequencySelectLcp(int64_t i64Timestamp_us, bool bFrequencySelectLcp, const std::string &strStatus);
    void                                    addFrequencySelectRcp(int64_t i64Timestamp_us, bool bFrequencySelectRcp, const std::string &strStatus);
    void                                    addFrequencyLOIntermed5GHz(int64_t i64Timestamp_us, double dFrequencyLOIntermediate5GHz_Hz, const std::string &strStatus);
    void                                    addFrequencyLOIntermed6_7GHz(int64_t i64Timestamp_us, double dFrequencyLOIntermediate6_7GHz_MHz, const std::string &strStatus);
    void                                    addFrequencyLOFinal(int64_t i64Timestamp_us, double dFrequencyLOFinal_MHz, const std::string &strStatus);
    void                                    addReceiverBandwidthLcp(int64_t i64Timestamp_us, double dReceiverBandwidthLcp_MHz, const std::string &strStatus);
    void                                    addReceiverBandwidthRcp(int64_t i64Timestamp_us, double dReceiverBandwidthRcp_MHz, const std::string &strStatus);
    void                                    addReceiverLcpAttenuation(int64_t i64Timestamp_us, double dReceiverLcpAttenuation_dB, const std::string &strStatus);
    void                                    addReceiverRcpAttenuation(int64_t i64Timestamp_us, double dReceiverRcpAttenuation_dB, const std::string &strStatus);

    void                                    addWindSpeed(int64_t i64Timestamp_us, double dWindSpeed_mps, const std::string &strStatus);
    void                                    addWindDirection(int64_t i64Timestamp_us, double dWindDirection_degrees, const std::string &strStatus);
    void                                    addTemperature(int64_t i64Timestamp_us, double dTemperature_degreesC, const std::string &strStatus);
    void                                    addAbsolutePressure(int64_t i64Timestamp_us, double dPressure_mbar, const std::string &strStatus);
    void                                    addRelativeHumidity(int64_t i64Timestamp_us, double dHumidity_percent, const std::string &strStatus);

    void                                    addAccumulationLength(int64_t i64Timestamp_us, uint32_t u32NFrames);
    void                                    addCoarseChannelSelect(int64_t i64Timestamp_us, uint32_t u32ChannelNo);
    void                                    setFrequencyFs(double dFrequencyFs_Hz);
    void                                    setSizeOfCoarseFFT(uint32_t u32SizeOfCoarseFFT_nSamp);
    void                                    setSizeOfFineFFT(uint32_t u32SizeOfFineFFT_nSamp);
    void                                    addCoarseFFTShiftMask(int64_t i64Timestamp_us, uint32_t u32ShiftMask);
    void                                    addDspGain(int64_t i64Timestamp_us, double dDspGain);
    void                                    addAttenuationADCChan0(int64_t i64Timestamp_us, double dADCAttenuationChan0_dB);
    void                                    addAttenuationADCChan1(int64_t i64Timestamp_us, double dADCAttenuationChan1_dB);

    void                                    setAntennaName(const std::string &strAntennaName);
    void                                    setAntennaBeamwidth(const std::string &strAntennaBeamwidth_deg);
    void                                    setAntennaDiameter(const std::string &strAntennaDiameter_m);
    void                                    setAntennaLatitude(const std::string &strAntennaLatitude_deg);
    void                                    setAntennaAltitude(const std::string &strAntennaAltitude_m);
    void                                    setAntennaLongitude(const std::string &strAntennaLongitude_deg);

    void                                    setAntennaDelayModel(const std::vector<double> &vdDelayModelParams);

    std::string                             getFilename() const;

private:
    std::string                             m_strFilename;
    AVN::Spectrometer::digitiserType        m_eDigitiserType;
    uint32_t                                m_u32NFrequencyBins;

    //HDF5:
    hid_t                                   m_iH5FileHandle;
    hid_t                                   m_iH5FileProperties;

    hid_t                                   m_iH5DataGroupHandle;
    hid_t                                   m_iH5MetaDataGroupHandle;
    hid_t                                   m_iH5MarkupGroupHandle;
    hid_t                                   m_iH5SensorsGroupHandle;
    hid_t                                   m_iH5SensorsAntennasGroupHandle;
    hid_t                                   m_iH5SensorsRFEGroupHandle;
    hid_t                                   m_iH5SensorsEnvGroupHandle;
    hid_t                                   m_iH5SensorsDBEGroupHandle;
    hid_t                                   m_iH5SensorsAntennasAntenna1GroupHandle;
    hid_t                                   m_iH5ConfigurationGroupHandle;
    hid_t                                   m_iH5ConfigurationAntennasGroupHandle;
    hid_t                                   m_iH5ConfigurationAntennasAntenna1GroupHandle;
    hid_t                                   m_iH5ConfigurationObservationGroupHandle;
    hid_t                                   m_iH5ConfigurationDBEGroupHandle;

    hid_t                                   m_iH5DatasetVis;
    hid_t                                   m_iH5DatasetStokes;

    hsize_t                                 m_aChannelDatasetDims[3];
    hsize_t                                 m_aChannelDatasetExtensionDims[3];
    hsize_t                                 m_aChannelDataOffset[3];
    hsize_t                                 m_aMemspaceSize[3];

    //Values received / calculated from ROACH sample data stream
    std::vector<double>                     m_vdSampleDataTimestamps_s;
    std::vector<std::vector<float> >        m_vvfChannelAverages;
    std::vector<cTimestampedChar>           m_voROACHNoiseDiodeStateChanges;

    //Values received from station controller
    std::vector<cMarkupLabels>              m_voMarkupLabels;

    std::vector<cTimestampedDouble>         m_voAcsRequestedAntennaAzs_deg;
    std::vector<cTimestampedDouble>         m_voAcsRequestedAntennaEls_deg;
    std::vector<cTimestampedDouble>         m_voAcsActualAntennaAzs_deg;
    std::vector<cTimestampedDouble>         m_voAcsActualAntennaEls_deg;

    std::vector<cTimestampedDouble>         m_voSkyRequestedAntennaAzs_deg;
    std::vector<cTimestampedDouble>         m_voSkyRequestedAntennaEls_deg;
    std::vector<cTimestampedDouble>         m_voSkyActualAntennaAzs_deg;
    std::vector<cTimestampedDouble>         m_voSkyActualAntennaEls_deg;

    /* Marked for removal.
    std::vector<cTimestampedDouble>         m_voActualSourceOffsetAzs_deg;
    std::vector<cTimestampedDouble>         m_voActualSourceOffsetEls_deg;
    std::vector<cTimestampedDouble>         m_voActualAntennaRAs_deg;
    std::vector<cTimestampedDouble>         m_voActualAntennaDecs_deg;
    */

    std::vector<cAntennaStatus>             m_voAntennaStatuses;

    /* Marked for removal.
    std::vector<cTimestampedDouble>         m_voMotorTorquesAzMaster_mNm;
    std::vector<cTimestampedDouble>         m_voMotorTorquesAzSlave_mNm;
    std::vector<cTimestampedDouble>         m_voMotorTorquesElMaster_mNm;
    std::vector<cTimestampedDouble>         m_voMotorTorquesElSlave_mNm;
    */

    std::vector<cNoiseDiodeSource>          m_voNoiseDiodeInputSource;
    std::vector<cTimestampedBool>           m_voNoiseDiodeEnable;
    std::vector<cTimestampedInt>            m_voNoiseDiodeSelect;
    std::vector<cTimestampedInt>            m_voNoiseDiodePWMMark;
    std::vector<cTimestampedDouble>         m_voNoiseDiodePWMFrequency;

    std::vector<cSourceSelection>           m_voSelectedSources;

    std::vector<cTimestampedChar>           m_voFrequencySelectLcp;
    std::vector<cTimestampedChar>           m_voFrequencySelectRcp;
    std::vector<cTimestampedDouble>         m_voFrequenciesLOIntermediate5GHz_Hz;
    std::vector<cTimestampedDouble>         m_voFrequenciesLOIntermediate6_7GHz_Hz;
    std::vector<cTimestampedDouble>         m_voFrequenciesLOFinal_Hz;
    std::vector<cTimestampedDouble>         m_voReceiverBandwidthsLcp_Hz;
    std::vector<cTimestampedDouble>         m_voReceiverBandwidthsRcp_Hz;
    std::vector<cTimestampedDouble>         m_voReceiverAttenuationsLcp_dB;
    std::vector<cTimestampedDouble>         m_voReceiverAttenuationsRcp_dB;

    std::vector<cTimestampedDouble>         m_voWindSpeeds_mps;
    std::vector<cTimestampedDouble>         m_voWindDirections_degrees;
    std::vector<cTimestampedDouble>         m_voTemperatures_degreesC;
    std::vector<cTimestampedDouble>         m_voAbsolutePressures_mbar;
    std::vector<cTimestampedDouble>         m_voRelativeHumidities_percent;

    cAntennaConfiguration                   m_oAntennaConfiguration;
    std::vector<double>                     m_vdDelayModelParams;
    double                                  m_adPointingModelParams[30]; //Only store most recent version

    //Values received from ROACH TCPBorphServer
    std::vector<cTimestampedUnsignedInt>    m_voROACHAccumulationLengths_nFrames;
    std::vector<cTimestampedUnsignedInt>    m_voROACHNBChannelSelects;
    double                                  m_dROACHFrequencyFs_Hz; //Only store most recent version (shouldn't ever change from 800 MHz)
    uint32_t                                m_u32ROACHSizeOfCoarseFFT_nSamp; //Only store most recent version (shouldn't ever change for a given gateware)
    uint32_t                                m_u32ROACHSizeOfFineFFT_nSamp; //Only store most recent version (shouldn't ever change for a given gateware)
    std::vector<cTimestampedDouble>         m_voROACHDspGains;
    std::vector<cTimestampedUnsignedInt>    m_voROACHCoarseFFTShiftMasks;
    std::vector<cTimestampedDouble>         m_voROACHADCAttenuationsLcp_dB;
    std::vector<cTimestampedDouble>         m_voROACHADCAttenuationsRcp_dB;

    //Other
    cSpectrometerHeader                     m_oLastHeader;

    boost::shared_mutex                     m_oAppendDataMutex; //share-locked for augmenting all station controller and TCPBorph data. Unique locked when writing this data to file

    float                                   calculateFrameAverage(const std::vector<int32_t> &vi32ChannelData);

    void                                    addAttributesToFile(const std::string &strVersion, const std::string &strExperimentID, int64_t i64AugmentTimestamp_us, uint32_t u32AugmentErrors, hid_t fileHandle);
    void                                    addAttributeToDataSet(const std::string &strDescription, const std::string &strName, const std::string &strType, const std::string &strUnits, hid_t dataset);
    void                                    addAttributesToObservation(
                                                const std::string &strScriptName,
                                                const std::string &strScriptArguments,
                                                const std::string &strObserver,
                                                const std::string &strExperimentID,
                                                const std::string &strDescription,
                                                const std::string &strAntennas,
                                                const std::string &strStartTime,
                                                const std::string &strEndTime,
                                                const std::string &strNoiseDiodeParams,
                                                const std::string &strRFParams,
                                                const std::string &strStatus,
                                                hid_t observationGroup
                                            );
    void                                    addAttributesToDBE(const std::string &strVisOrdering, const std::string &strStokesOrdering, hid_t DBEGroup);

    //Write logged data to file (after sample recording is completed)
    void                                    writeSampleDataTimestamps();
    void                                    writeChannelAverages();
    void                                    writeROACHNoiseDiodeStates();

    void                                    writeMarkupLabels();

    void                                    writeAcsRequestedAntennaAzEls();
    void                                    writeAcsActualAntennaAzEls();
    void                                    writeSkyRequestedAntennaAzEls();
    void                                    writeSkyActualAntennaAzEls();

    /* Marked for removal.
    void                                    writeActualSourceOffsetAzEls();
    void                                    writeActualAntennaRADecs();
    */

    void                                    writeAntennaStatuses();
    void                                    writeMotorTorques();
    void                                    writeAntennaConfiguration();

    void                                    writeNoiseDiodeInformation();
    void                                    writeNoiseDiodeSources();
    void                                    writeNoideDiodeCurrents();

    void                                    writeSelectedSources();

    void                                    writeRFFrequencies();
    void                                    writeLOFrequencies();
    void                                    writeIFBandwidths();
    void                                    writeReceiverAttenuations();

    void                                    writeEnvironmentData();

    void                                    writeROACHNumberChannels();
    void                                    writeROACHAccumulationLengths();
    void                                    writeROACHNBNarrowbandSelections();
    void                                    writeROACHSamplingFreqBandwidth();
    void                                    writeROACHSizeOfFFTs();
    void                                    writeROACHCoarseFFTShiftMask();
    void                                    writeROACHDspGains();
    void                                    writeROACHAdcAttentuations();

};

#endif // SPECTROMETER_HDF5_OUTPUT_FILE_H
