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
//Local includes
#include "SpectrometerHeader.h"
#include "SpectrometerDefinitions.h"

class cSpectrometerHDF5OutputFile
{
    typedef struct cNoiseDiodeState
    {
        double             m_dTimeStamp_s;
        char               m_caState[1];
    } cNoiseDiodeState;

public:
    cSpectrometerHDF5OutputFile(const std::string &strFilename, AVN::Spectrometer::digitiserType eDigitiserType, uint32_t u32NFrequencyBins);
    ~cSpectrometerHDF5OutputFile();

    void                        addFrame(const std::vector<int> &vi32Chan0, const std::vector<int> &vi32Chan1, const std::vector<int> &vi32Chan2, std::vector<int> &vi32Chan3,
                                         const cSpectrometerHeader &oHeader);

private:
    std::string                         m_strFilename;
    AVN::Spectrometer::digitiserType    m_eDigitiserType;
    uint32_t                            m_u32NFrequencyBins;

    //HDF5:
    hid_t                               m_iH5FileHandle;
    hid_t                               m_iH5FileProperties;

    hid_t                               m_iH5DataGroupHandle;
    hid_t                               m_iH5MetaDataGroupHandle;
    hid_t                               m_iH5SensorsGroupHandle;
    hid_t                               m_iH5AntennasGroupHandle;
    hid_t                               m_iH5Antenna1GroupHandle;

    hid_t                               m_iH5DatasetVis;
    hid_t                               m_iH5DatasetStokes;

    hsize_t                             m_aChannelDatasetDims[3];
    hsize_t                             m_aChannelDatasetExtensionDims[3];
    hsize_t                             m_aChannelDataOffset[3];
    hsize_t                             m_aMemspaceSize[3];

    std::vector<double>                 m_vdTimestamps_s;
    std::vector<std::vector<float> >    m_vvfChannelAverages;
    std::vector<cNoiseDiodeState>       m_voNoiseDiodeStateChanges;

    cSpectrometerHeader                 m_oLastHeader;

    float                               calculateFrameAverage(const std::vector<int32_t> &vi32ChannelData, AVN::Spectrometer::digitiserType eDigitiserType);

    void                                writeTimestamps();
    void                                writeChannelAverages();
    void                                writeNoiseDiodeStates();

};

#endif // SPECTROMETER_HDF5_OUTPUT_FILE_H
