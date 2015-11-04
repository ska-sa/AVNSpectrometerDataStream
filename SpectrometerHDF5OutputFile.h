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

public:
    cSpectrometerHDF5OutputFile(const std::string &strFilename, AVN::Spectrometer::digitiserType eDigitiserType);
    ~cSpectrometerHDF5OutputFile();

    void                        addFrame(const std::vector<int> &vi32Chan0, const std::vector<int> &vi32Chan1, const std::vector<int> &vi32Chan2, std::vector<int> &vi32Chan3,
                                         const cSpectrometerHeader &oHeader);

    std::string                         m_strFilename;
    AVN::Spectrometer::digitiserType    m_eDigitiserType;

private:

    //HDF5:
    hid_t                               m_iH5FileHandle;
    hid_t                               m_iH5FileProperties;
    hid_t                               m_iH5DataGroupHandle;

    hid_t                               m_iH5Channel0Dataset;
    hid_t                               m_iH5Channel1Dataset;
    hid_t                               m_iH5Channel2Dataset;
    hid_t                               m_iH5Channel3Dataset;

    hsize_t                             m_aChannelDatasetDims[2];
    hsize_t                             m_aChannelDatasetExtensionDims[2];
    hsize_t                             m_aChannelDataOffset[2];

    std::vector<int64_t>                m_vi64Timestamps_us;
    std::vector<std::vector<float> >    m_vvfChannelAverages;

    cSpectrometerHeader                 m_oLastHeader;

    float                               calculateFrameAverage(const std::vector<int32_t> &vi32ChannelData, AVN::Spectrometer::digitiserType eDigitiserType);

    void                                writeTimestamps();
    void                                writeChannelAverages();

};

#endif // SPECTROMETER_HDF5_OUTPUT_FILE_H
