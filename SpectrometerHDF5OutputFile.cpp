
//System includes

//Library includes
extern "C" {
#include <hdf5_hl.h>
#include <H5LTpublic.h>
}

//Local includes
#include "SpectrometerHDF5OutputFile.h"

cSpectrometerHDF5OutputFile::cSpectrometerHDF5OutputFile(const std::string &strFilename, AVN::Spectrometer::digitiserType eDigitiserType) :
    m_strFilename(strFilename),
    m_eDigitiserType(eDigitiserType)
{
    //Create file (overwrite any existing one with the same name)
    m_iH5FileHandle = H5Fcreate(m_strFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t chunkDimensions[2] = {AVN::Spectrometer::getDigitiserFrameSize_nBins(eDigitiserType), 512}; //Each chuck is 512 frames


    m_iH5Channel0PTable = H5PTcreate_fl(m_iH5FileHandle, AVN::Spectrometer::getDigitiserChannelName(eDigitiserType, 0).c_str(), H5T_NATIVE_INT, chunkDimensions);
    m_iH5Channel1PTable = H5PTcreate_fl(m_iH5FileHandle, AVN::Spectrometer::getDigitiserChannelName(eDigitiserType, 1).c_str(), H5T_NATIVE_INT, chunkDimensions);
    m_iH5Channel2PTable = H5PTcreate_fl(m_iH5FileHandle, AVN::Spectrometer::getDigitiserChannelName(eDigitiserType, 2).c_str(), H5T_NATIVE_INT, chunkDimensions);
    m_iH5Channel3PTable = H5PTcreate_fl(m_iH5FileHandle, AVN::Spectrometer::getDigitiserChannelName(eDigitiserType, 3).c_str(), H5T_NATIVE_INT, chunkDimensions);

}

cSpectrometerHDF5OutputFile::~cSpectrometerHDF5OutputFile()
{
    //Close and free all HDF5 structures
    H5Fclose(m_iH5Channel0PTable);
    H5Fclose(m_iH5Channel1PTable);
    H5Fclose(m_iH5Channel2PTable);
    H5Fclose(m_iH5Channel3PTable);

    H5Fclose(m_iH5FileHandle);
}

void cSpectrometerHDF5OutputFile::addFrame(const std::vector<int> &vi32Chan0, const std::vector<int> &vi32Chan1, const std::vector<int> &vi32Chan2, std::vector<int> &vi32Chan3,
                                           const cSpectrometerHeader &oHeader)
{
    H5PTappend(m_iH5Channel0PTable, (hsize_t)1, &vi32Chan0.front());
    H5PTappend(m_iH5Channel1PTable, (hsize_t)1, &vi32Chan1.front());
    H5PTappend(m_iH5Channel2PTable, (hsize_t)1, &vi32Chan2.front());
    H5PTappend(m_iH5Channel3PTable, (hsize_t)1, &vi32Chan3.front());

    m_vi64Timestamps_us.push_back(oHeader.getTimestamp_us());
}

void cSpectrometerHDF5OutputFile::writeTimestamps()
{
    hsize_t dimension = m_vi64Timestamps_us.size();
    H5LTmake_dataset(m_iH5FileHandle, "Timestamps", 1, dimension, H5T_NATIVE_LLONG, &m_vi64Timestamps_us.front());
}

