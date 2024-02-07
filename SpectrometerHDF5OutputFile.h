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
        char                m_chaValue[2];
        char                m_chaStatus[8];
    } cTimestampedChar;

    typedef struct cTimestampedBool
    {
        double              m_dTimestamp_s;
        int32_t             m_i32Value;
        char                m_chaStatus[8];
    } cTimestampedBool;

    typedef struct cSourceSelection
    {
        double              m_dTimestamp_s;
        char                m_chaSource[128];
        char                m_chaStatus[8];
    } cSourceSelection;

    typedef struct cAntennaStatus
    {
        double              m_dTimestamp_s;
        char                m_chaAntennaStatus[16];
        char                m_chaStatus[8];
    } cAntennaStatus;

    typedef struct cNoiseDiodeSource
    {
        double              m_dTimestamp_s;
        char                m_chaSource[8];
        char                m_chaStatus[8];
    } cNoiseDiodeSource;

    //General structures
    typedef struct cTimestampedInt
    {
        double              m_dTimestamp_s;
        int32_t             m_i32Value;
        char                m_chaStatus[8];
    } cTimestampedInt;

    typedef struct cTimestampedUnsignedInt
    {
        double              m_dTimestamp_s;
        uint32_t            m_u32Value;
        char                m_chaStatus[8];
    } cTimestampedUnsignedInt;

    typedef struct cTimestampedDouble
    {
        double              m_dTimestamp_s;
        double              m_dValue;
        char                m_chaStatus[8];
    } cTimestampedDouble;

    typedef struct cTimestampedDualDouble
    {
        double              m_dTimestamp_s;
        double              m_dValue1;
        double              m_dValue2;
        char                m_chaStatus[8];
    } cTimestampedDualDouble;

    typedef struct cAntennaConfiguration
    {
        char                m_chaAntennaName[64];
        char                m_chaObserverName[64];
        double              m_dAntennaDiameter_m;
        double              m_dAntennaBeamwidth_deg;
        double              m_dAntennaLongitude_deg;
        double              m_dAntennaLatitude_deg;
        double              m_dAntennaAltitude_m;
    } cAntennaConfiguration;


public:
    cSpectrometerHDF5OutputFile(const std::string &strFilename, AVN::Spectrometer::digitiserType eDigitiserType, uint32_t u32NFrequencyBins);
    ~cSpectrometerHDF5OutputFile();

    //Different functions can be called concurrently but any single function is not at present thread safe with regard to concurrent calls
    void                                    addFrame(const std::vector<int> &vi32Chan0, const std::vector<int> &vi32Chan1, const std::vector<int> &vi32Chan2, std::vector<int> &vi32Chan3,
                                                     const cSpectrometerHeader &oHeader);

    void                                    addMarkupLabel(int64_t i64Timestamp_us, const std::string &strLabel);

    // Noise diode tables
    void                                    addNoiseDiodeData(const std::string &strPath);

    // Antenna-space values
    void                                    addAcsRequestedAz(int64_t i64Timestamp_us, double dAzimuth_deg, const std::string &strStatus);
    void                                    addAcsRequestedEl(int64_t i64Timestamp_us, double dElevation_deg, const std::string &strStatus);
    void                                    addAcsDesiredAz(int64_t i64Timestamp_us, double dAzimuth_deg, const std::string &strStatus);
    void                                    addAcsDesiredEl(int64_t i64Timestamp_us, double dElevation_deg, const std::string &strStatus);
    void                                    addAcsActualAz(int64_t i64Timestamp_us, double dAzimuth_deg, const std::string &strStatus);
    void                                    addAcsActualEl(int64_t i64Timestamp_us, double dElevation_deg, const std::string &strStatus);

    // Sky-space values
    void                                    addSkyRequestedAz(int64_t i64Timestamp_us, double dAzimuth_deg, const std::string &strStatus);
    void                                    addSkyRequestedEl(int64_t i64Timestamp_us, double dElevation_deg, const std::string &strStatus);
    void                                    addSkyDesiredAz(int64_t i64Timestamp_us, double dAzimuth_deg, const std::string &strStatus);
    void                                    addSkyDesiredEl(int64_t i64Timestamp_us, double dElevation_deg, const std::string &strStatus);
    void                                    addSkyActualAz(int64_t i64Timestamp_us, double dAzimuth_deg, const std::string &strStatus);
    void                                    addSkyActualEl(int64_t i64Timestamp_us, double dElevation_deg, const std::string &strStatus);

    void                                    addPointingModelParameter(uint8_t ui8ParameterNumber, double dParameterValue);

    void                                    addAntennaStatus(int64_t i64Timestamp_us, const std::string &strAntennaStatus, const std::string &strStatus);
    void                                    addObservationStatus(int64_t i64Timestamp_us, const std::string &strObservationStatus, const std::string &strStatus);

    void                                    addNoiseDiode5GHzInputSource(int64_t i64Timestamp_us, const std::string &strNoiseDiodeInputSource, const std::string &strStatus);
    void                                    addNoiseDiode5GHzLevel(int64_t i64Timestamp_us, int32_t i32NoiseDiodeLevel, const std::string &strStatus);
    void                                    addNoiseDiode5GHzPWMMark(int64_t i64Timestamp_us, int32_t i32NoiseDiodePWMMark, const std::string &strStatus);
    void                                    addNoiseDiode5GHzPWMFrequency(int64_t i64Timestamp_us, double dNoiseDiodePWMFrequency, const std::string &strStatus);
    void                                    addNoiseDiode6_7GHzInputSource(int64_t i64Timestamp_us, const std::string &strNoiseDiodeInputSource, const std::string &strStatus);
    void                                    addNoiseDiode6_7GHzLevel(int64_t i64Timestamp_us, int32_t i32NoiseDiodeLevel, const std::string &strStatus);
    void                                    addNoiseDiode6_7GHzPWMMark(int64_t i64Timestamp_us, int32_t i32NoiseDiodePWMMark, const std::string &strStatus);
    void                                    addNoiseDiode6_7GHzPWMFrequency(int64_t i64Timestamp_us, double dNoiseDiodePWMFrequency, const std::string &strStatus);


    void                                    addSourceSelection(int64_t i64Timestamp_us, const std::string &strSourceName, const std::string &strStatus);

    void                                    addBandSelectLcp(int64_t i64Timestamp_us, bool bBandSelectLcp, const std::string &strStatus);
    void                                    addBandSelectRcp(int64_t i64Timestamp_us, bool bBandSelectRcp, const std::string &strStatus);
    void                                    addFrequencySky5GHz(int64_t i64Timestamp_us, double dFrequencySky5GHz_Hz, const std::string &strStatus);
    void                                    addFrequencySky6_7GHz(int64_t i64Timestamp_us, double dFrequencySky6_7GHz_MHz, const std::string &strStatus);
    void                                    addReceiverGain5GHzLcp(int64_t i64Timestamp_us, double dGain_dB, const std::string &strStatus);
    void                                    addReceiverGain5GHzRcp(int64_t i64Timestamp_us, double dGain_dB, const std::string &strStatus);
    void                                    addReceiverGain6_7GHzLcp(int64_t i64Timestamp_us, double dGain_dB, const std::string &strStatus);
    void                                    addReceiverGain6_7GHzRcp(int64_t i64Timestamp_us, double dGain_dB, const std::string &strStatus);

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

    void                                    setObservationInfo(const std::string &strObservationInfo);
    void                                    setAntennaBeamwidth(const double &dAntennaBeamwidth_deg);

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
    hid_t                                   m_iH5NdGroupHandle;
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
    std::vector<cTimestampedDouble>         m_voAcsDesiredAntennaAzs_deg;
    std::vector<cTimestampedDouble>         m_voAcsDesiredAntennaEls_deg;
    std::vector<cTimestampedDouble>         m_voAcsActualAntennaAzs_deg;
    std::vector<cTimestampedDouble>         m_voAcsActualAntennaEls_deg;

    std::vector<cTimestampedDouble>         m_voSkyRequestedAntennaAzs_deg;
    std::vector<cTimestampedDouble>         m_voSkyRequestedAntennaEls_deg;
    std::vector<cTimestampedDouble>         m_voSkyDesiredAntennaAzs_deg;
    std::vector<cTimestampedDouble>         m_voSkyDesiredAntennaEls_deg;
    std::vector<cTimestampedDouble>         m_voSkyActualAntennaAzs_deg;
    std::vector<cTimestampedDouble>         m_voSkyActualAntennaEls_deg;

    std::vector<cAntennaStatus>             m_voAntennaStatuses;
    std::vector<cAntennaStatus>             m_voObservationStatuses; // Same datatype as above.

    std::vector<cNoiseDiodeSource>          m_voNoiseDiode5GHzInputSource;
    std::vector<cTimestampedInt>            m_voNoiseDiode5GHzLevel;
    std::vector<cTimestampedInt>            m_voNoiseDiode5GHzPWMMark;
    std::vector<cTimestampedDouble>         m_voNoiseDiode5GHzPWMFrequency;
    std::vector<cNoiseDiodeSource>          m_voNoiseDiode6_7GHzInputSource;
    std::vector<cTimestampedInt>            m_voNoiseDiode6_7GHzLevel;
    std::vector<cTimestampedInt>            m_voNoiseDiode6_7GHzPWMMark;
    std::vector<cTimestampedDouble>         m_voNoiseDiode6_7GHzPWMFrequency;

    std::vector<cSourceSelection>           m_voSelectedSources;

    std::vector<cTimestampedChar>           m_voBandSelectLcp;
    std::vector<cTimestampedChar>           m_voBandSelectRcp;
    std::vector<cTimestampedDouble>         m_voFrequenciesSky5GHz_Hz;
    std::vector<cTimestampedDouble>         m_voFrequenciesSky6_7GHz_Hz;
    std::vector<cTimestampedDouble>         m_voReceiverGain5GHzLcp_dB;
    std::vector<cTimestampedDouble>         m_voReceiverGain5GHzRcp_dB;
    std::vector<cTimestampedDouble>         m_voReceiverGain6_7GHzLcp_dB;
    std::vector<cTimestampedDouble>         m_voReceiverGain6_7GHzRcp_dB;

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
    void                                    writeAcsDesiredAntennaAzEls();
    void                                    writeAcsActualAntennaAzEls();
    void                                    writeSkyRequestedAntennaAzEls();
    void                                    writeSkyDesiredAntennaAzEls();
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

    void                                    writeRFBandSelects();
    void                                    writeSkyFrequencies();
    void                                    writeReceiverGains();

    void                                    writeEnvironmentData();

    void                                    writeROACHNumberChannels();
    void                                    writeROACHAccumulationLengths();
    void                                    writeROACHNBNarrowbandSelections();
    void                                    writeROACHSamplingFreqBandwidth();
    void                                    writeROACHSizeOfFFTs();
    void                                    writeROACHCoarseFFTShiftMask();
    void                                    writeROACHDspGains();
    void                                    writeROACHAdcAttentuations();

    // Noise diode csv files
    void                                    writeCsv2Hdf(const std::string &strPath, const std::string &strDataSetName);
    
// PJP
    double                                  DmsToDeg(std::string strDms);

};

#endif // SPECTROMETER_HDF5_OUTPUT_FILE_H
