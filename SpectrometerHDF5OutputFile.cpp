
//System includes
#include <iostream>
#include <sstream>
#include <complex>

//Library includes
extern "C" {
#include <hdf5.h> //Note must be included before other HDF5 libraries
#include <hdf5_hl.h>
}

//Local includes
#include "SpectrometerHDF5OutputFile.h"
#include "../../AVNUtilLibs/Timestamp/Timestamp.h"

using namespace std;

cSpectrometerHDF5OutputFile::cSpectrometerHDF5OutputFile(const std::string &strFilename, AVN::Spectrometer::digitiserType eDigitiserType) :
    m_strFilename(strFilename),
    m_eDigitiserType(eDigitiserType),
    m_vvfChannelAverages(4)
{
    //Create file (overwrite any existing one with the same name)
    m_iH5FileHandle = H5Fcreate(m_strFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if(m_iH5FileHandle == H5I_INVALID_HID)
    {
        cout << "Error opening HDF5 file." << endl;
    }

    m_iH5DataGroupHandle = H5Gcreate2(m_iH5FileHandle, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    m_aChannelDatasetDims[0] = 1;
    m_aChannelDatasetDims[1] = 1024;

    m_aChannelDatasetExtensionDims[0] = 1;
    m_aChannelDatasetExtensionDims[1] = 1024;

    m_aChannelDataOffset[0] = 0;
    m_aChannelDataOffset[1] = 0;

    hsize_t maxDims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};

    // Create a dataset creation property list and set it to use chunking
    hid_t datasetProperties = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(datasetProperties, 2, m_aChannelDatasetDims);

    /* Create the dataspace and the chunked dataset */
    hid_t dataSpace0 = H5Screate_simple(2, m_aChannelDatasetDims, maxDims);
    hid_t dataSpace1 = H5Screate_simple(2, m_aChannelDatasetDims, maxDims);
    hid_t dataSpace2 = H5Screate_simple(2, m_aChannelDatasetDims, maxDims);
    hid_t dataSpace3 = H5Screate_simple(2, m_aChannelDatasetDims, maxDims);

    m_iH5Channel0Dataset = H5Dcreate1(m_iH5DataGroupHandle, AVN::Spectrometer::getDigitiserChannelName(eDigitiserType, 0).c_str(), H5T_NATIVE_INT, dataSpace0, datasetProperties);
    m_iH5Channel1Dataset = H5Dcreate1(m_iH5DataGroupHandle, AVN::Spectrometer::getDigitiserChannelName(eDigitiserType, 1).c_str(), H5T_NATIVE_INT, dataSpace1, datasetProperties);
    m_iH5Channel2Dataset = H5Dcreate1(m_iH5DataGroupHandle, AVN::Spectrometer::getDigitiserChannelName(eDigitiserType, 2).c_str(), H5T_NATIVE_INT, dataSpace2, datasetProperties);
    m_iH5Channel3Dataset = H5Dcreate1(m_iH5DataGroupHandle, AVN::Spectrometer::getDigitiserChannelName(eDigitiserType, 3).c_str(), H5T_NATIVE_INT, dataSpace3, datasetProperties);


    if(m_iH5Channel0Dataset == H5I_INVALID_HID || m_iH5Channel1Dataset == H5I_INVALID_HID || m_iH5Channel2Dataset == H5I_INVALID_HID || m_iH5Channel3Dataset == H5I_INVALID_HID)
    {
        cout << "Error: Creating HDF5 datasets failed." << endl;
    }
    else
    {
        cout << "Successfully created HDF5 datasets." << endl;
    }

    H5Sclose(dataSpace0);
    H5Sclose(dataSpace1);
    H5Sclose(dataSpace2);
    H5Sclose(dataSpace3);
    H5Pclose(datasetProperties);
}

cSpectrometerHDF5OutputFile::~cSpectrometerHDF5OutputFile()
{
    cout << "cSpectrometerHDF5OutputFile::~cSpectrometerHDF5OutputFile(): Got close request, writing accumulated data to end of HDF5 file... " << endl;

    writeTimestamps();
    writeChannelAverages();

    cout << "cSpectrometerHDF5OutputFile::~cSpectrometerHDF5OutputFile(): Done writing accumulated data." << endl;

    //Close and free all HDF5 structures
    H5Dclose(m_iH5Channel0Dataset);
    H5Dclose(m_iH5Channel1Dataset);
    H5Dclose(m_iH5Channel2Dataset);
    H5Dclose(m_iH5Channel3Dataset);

    H5Gclose(m_iH5DataGroupHandle);

    H5Fclose(m_iH5FileHandle);
}

void cSpectrometerHDF5OutputFile::addFrame(const std::vector<int> &vi32Chan0, const std::vector<int> &vi32Chan1, const std::vector<int> &vi32Chan2, std::vector<int> &vi32Chan3,
                                           const cSpectrometerHeader &oHeader)
{
    m_aChannelDatasetDims[0] += m_aChannelDatasetExtensionDims[0]; //Extend in the row dimension
    m_aChannelDatasetDims[1] = m_aChannelDatasetDims[1];

    H5Dset_extent(m_iH5Channel0Dataset, m_aChannelDatasetDims);
    hid_t filespace0 = H5Dget_space(m_iH5Channel0Dataset);
    H5Sselect_hyperslab (filespace0, H5S_SELECT_SET, m_aChannelDataOffset, NULL, m_aChannelDatasetExtensionDims, NULL);
    hid_t memspace0 = H5Screate_simple (2, m_aChannelDatasetExtensionDims, NULL);

    herr_t err0 = H5Dwrite(m_iH5Channel0Dataset, H5T_NATIVE_INT, memspace0, filespace0, H5P_DEFAULT, &vi32Chan0.front());
    if(err0 < 0)
    {
        cout << "HDF5 chunk extend error on 1st channel: " << err0 << endl;
    }

    H5Dset_extent(m_iH5Channel1Dataset, m_aChannelDatasetDims);
    hid_t filespace1 = H5Dget_space(m_iH5Channel1Dataset);
    H5Sselect_hyperslab (filespace1, H5S_SELECT_SET, m_aChannelDataOffset, NULL, m_aChannelDatasetExtensionDims, NULL);
    hid_t memspace1 = H5Screate_simple (2, m_aChannelDatasetExtensionDims, NULL);

    herr_t err1 = H5Dwrite(m_iH5Channel1Dataset, H5T_NATIVE_INT, memspace1, filespace1, H5P_DEFAULT, &vi32Chan1.front());
    if(err1 < 0)
    {
        cout << "HDF5 chunk extend error on 2nd channel: " << err1 << endl;
    }

    H5Dset_extent(m_iH5Channel2Dataset, m_aChannelDatasetDims);
    hid_t filespace2 = H5Dget_space(m_iH5Channel2Dataset);
    H5Sselect_hyperslab (filespace2, H5S_SELECT_SET, m_aChannelDataOffset, NULL, m_aChannelDatasetExtensionDims, NULL);
    hid_t memspace2 = H5Screate_simple (2, m_aChannelDatasetExtensionDims, NULL);

    herr_t err2 = H5Dwrite(m_iH5Channel2Dataset, H5T_NATIVE_INT, memspace2, filespace2, H5P_DEFAULT, &vi32Chan2.front());
    if(err2 < 0)
    {
        cout << "HDF5 chunk extend error on 3rd channel: " << err2 << endl;
    }

    H5Dset_extent(m_iH5Channel3Dataset, m_aChannelDatasetDims);
    hid_t filespace3 = H5Dget_space(m_iH5Channel3Dataset);
    H5Sselect_hyperslab (filespace3, H5S_SELECT_SET, m_aChannelDataOffset, NULL, m_aChannelDatasetExtensionDims, NULL);
    hid_t memspace3 = H5Screate_simple (2, m_aChannelDatasetExtensionDims, NULL);

    herr_t err3 = H5Dwrite(m_iH5Channel3Dataset, H5T_NATIVE_INT, memspace3, filespace3, H5P_DEFAULT, &vi32Chan3.front());
    if(err3 < 0)
    {
        cout << "HDF5 chunk extend error on 4th channel: " << err3 << endl;
    }

    H5Sclose (memspace0);
    H5Sclose (memspace1);
    H5Sclose (memspace2);
    H5Sclose (memspace3);

    H5Sclose (filespace0);
    H5Sclose (filespace1);
    H5Sclose (filespace2);
    H5Sclose (filespace3);

    //Store values to be written after channel data in memory
    m_vi64Timestamps_us.push_back(oHeader.getTimestamp_us());
    m_vvfChannelAverages[0].push_back( calculateFrameAverage(vi32Chan0, (AVN::Spectrometer::digitiserType)oHeader.getDigitiserType()) );
    m_vvfChannelAverages[1].push_back( calculateFrameAverage(vi32Chan1, (AVN::Spectrometer::digitiserType)oHeader.getDigitiserType()) );
    m_vvfChannelAverages[2].push_back( calculateFrameAverage(vi32Chan2, (AVN::Spectrometer::digitiserType)oHeader.getDigitiserType()) );
    m_vvfChannelAverages[3].push_back( calculateFrameAverage(vi32Chan3, (AVN::Spectrometer::digitiserType)oHeader.getDigitiserType()) );

    m_aChannelDataOffset[0] += m_aChannelDatasetExtensionDims[0];

    m_oLastHeader = oHeader;

}

void cSpectrometerHDF5OutputFile::writeTimestamps()
{
    hsize_t dimension = m_vi64Timestamps_us.size();
    herr_t err = H5LTmake_dataset(m_iH5DataGroupHandle, "Timestamps_us", 1, &dimension, H5T_NATIVE_LLONG, &m_vi64Timestamps_us.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeTimestamps(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeTimestamps(): Wrote " << m_vi64Timestamps_us.size() << " timestamps to dataset." << endl;
    }
}

void cSpectrometerHDF5OutputFile::writeChannelAverages()
{
    hsize_t dimension;

    for(uint64_t ui = 0; ui <  m_vvfChannelAverages.size(); ui++)
    {
        dimension = m_vvfChannelAverages[ui].size();

        stringstream oSS;
        oSS << AVN::Spectrometer::getDigitiserChannelName(m_oLastHeader.getDigitiserType(), ui);
        oSS << string(" time average");

        herr_t err = H5LTmake_dataset(m_iH5DataGroupHandle, oSS.str().c_str(), 1, &dimension, H5T_NATIVE_FLOAT, &m_vvfChannelAverages[ui].front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeChannelAverages(): HDF5 make dataset error for channel " << ui << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeChannelAverages(): Wrote " << m_vvfChannelAverages[ui].size() << " time averages to dataset for channel " << ui << "." << endl;
        }
    }
}

float cSpectrometerHDF5OutputFile::calculateFrameAverage(const vector<int32_t> &vi32ChannelData, AVN::Spectrometer::digitiserType eDigitiserType)
{
    double dAverage = 0.0;

    switch(eDigitiserType)
    {
    case AVN::Spectrometer::NB_SPECTROMETER_LRQU:
    case AVN::Spectrometer::WB_SPECTROMETER_LRQU:

        for(uint32_t ui = 0; ui < vi32ChannelData.size(); ui++)
        {
            dAverage += vi32ChannelData[ui];
        }

        dAverage /= (double)vi32ChannelData.size();

        break;

    case AVN::Spectrometer::NB_SPECTROMETER_CFFT:
    case AVN::Spectrometer::WB_SPECTROMETER_CFFT:

        for(uint32_t ui = 0; ui < vi32ChannelData.size() / 2; ui += 2)
        {
            dAverage += abs( *( (complex<int32_t>*)&vi32ChannelData[ui] ) );
        }

        dAverage /= (double)(vi32ChannelData.size() / 2);

        break;

    default:
        return 0.0f;
    }

    return (float)dAverage;
}
