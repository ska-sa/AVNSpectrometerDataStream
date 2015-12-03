
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

cSpectrometerHDF5OutputFile::cSpectrometerHDF5OutputFile(const std::string &strFilename, AVN::Spectrometer::digitiserType eDigitiserType, uint32_t u32NFrequencyBins) :
    m_strFilename(strFilename),
    m_eDigitiserType(eDigitiserType),
    m_u32NFrequencyBins(u32NFrequencyBins),
    m_vvfChannelAverages(4)
{
    //Create file (overwrite any existing one with the same name)
    m_iH5FileHandle = H5Fcreate(m_strFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if(m_iH5FileHandle == H5I_INVALID_HID)
    {
        cout << "Error opening HDF5 file." << endl;
    }

    //Setup data groups:
    //Level 1:
    m_iH5DataGroupHandle = H5Gcreate2(m_iH5FileHandle, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5MetaDataGroupHandle = H5Gcreate2(m_iH5FileHandle, "/MetaData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //Level 2:
    m_iH5SensorsGroupHandle = H5Gcreate2(m_iH5MetaDataGroupHandle, "/MetaData/Sensors", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //Level 3:
    m_iH5AntennasGroupHandle  = H5Gcreate2(m_iH5SensorsGroupHandle, "/MetaData/Sensors/Antennas", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //Level 4:
    m_iH5Antenna1GroupHandle  = H5Gcreate2(m_iH5AntennasGroupHandle, "/MetaData/Sensors/Antennas/ant1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //Note data dimensions are [ time x frequency bin x data channel ]
    m_aChannelDatasetDims[0] = 1;
    m_aChannelDatasetDims[1] = m_u32NFrequencyBins;
    m_aChannelDatasetDims[2] = 2;

    m_aChannelDatasetExtensionDims[0] = 1;
    m_aChannelDatasetExtensionDims[1] = m_aChannelDatasetDims[1];
    m_aChannelDatasetExtensionDims[2] = m_aChannelDatasetDims[2];

    m_aMemspaceSize[0] = 1;
    m_aMemspaceSize[1] = m_aChannelDatasetExtensionDims[1];
    m_aMemspaceSize[2] = 1;

    m_aChannelDataOffset[0] = 0;
    m_aChannelDataOffset[1] = 0;
    m_aChannelDataOffset[2] = 0;

    hsize_t maxDims[3] = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED};

    // Create a dataset creation property list and set it to use chunking
    hid_t datasetPropertiesVis      = H5Pcreate(H5P_DATASET_CREATE);
    hid_t datasetPropertiesStokes   = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(datasetPropertiesVis, 3, m_aChannelDatasetDims);
    H5Pset_chunk(datasetPropertiesStokes, 3, m_aChannelDatasetDims);

    /* Create the dataspace and the chunked dataset */
    hid_t dataSpaceVis = H5Screate_simple(3, m_aChannelDatasetDims, maxDims);
    hid_t dataSpaceStokes = H5Screate_simple(3, m_aChannelDatasetDims, maxDims);


    m_iH5DatasetVis     = H5Dcreate1(m_iH5DataGroupHandle, "VisData", H5T_NATIVE_INT, dataSpaceVis, datasetPropertiesVis);
    m_iH5DatasetStokes  = H5Dcreate1(m_iH5DataGroupHandle, "StokesData", H5T_NATIVE_INT, dataSpaceStokes, datasetPropertiesStokes);


    if(m_iH5DatasetVis == H5I_INVALID_HID || m_iH5DatasetStokes == H5I_INVALID_HID)
    {
        cout << "Error: Creating HDF5 datasets failed." << endl;
    }
    else
    {
        cout << "Successfully created HDF5 datasets." << endl;
    }

    H5Sclose(dataSpaceVis);
    H5Sclose(dataSpaceStokes);
    H5Pclose(datasetPropertiesVis);
    H5Pclose(datasetPropertiesStokes);
}

cSpectrometerHDF5OutputFile::~cSpectrometerHDF5OutputFile()
{
    cout << "cSpectrometerHDF5OutputFile::~cSpectrometerHDF5OutputFile(): Got close request, writing accumulated data to end of HDF5 file... " << endl;

    writeTimestamps();
    writeChannelAverages();
    writeNoiseDiodeStates();

    cout << "cSpectrometerHDF5OutputFile::~cSpectrometerHDF5OutputFile(): Done writing accumulated data." << endl;

    //Close and free all HDF5 structures
    H5Dclose(m_iH5DatasetVis);
    H5Dclose(m_iH5DatasetStokes);

    //HDF5 groups:
    //Level 4
    H5Gclose(m_iH5Antenna1GroupHandle);
    //Level 3
    H5Gclose(m_iH5AntennasGroupHandle);
    //Level 2
    H5Gclose(m_iH5SensorsGroupHandle);
    //Level 1
    H5Gclose(m_iH5MetaDataGroupHandle);
    H5Gclose(m_iH5DataGroupHandle);

    //Level 0 (file)
    H5Fclose(m_iH5FileHandle);
}

void cSpectrometerHDF5OutputFile::addFrame(const std::vector<int> &vi32Chan0, const std::vector<int> &vi32Chan1, const std::vector<int> &vi32Chan2, std::vector<int> &vi32Chan3,
                                           const cSpectrometerHeader &oHeader)
{
    m_aChannelDatasetDims[0] += m_aChannelDatasetExtensionDims[0]; //Extend in the time dimension

    //Visibilty dataset (LL + RR)
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    H5Dset_extent(m_iH5DatasetVis, m_aChannelDatasetDims);
    hid_t filespaceVis = H5Dget_space(m_iH5DatasetVis);

    m_aChannelDataOffset[2] = 0;
    H5Sselect_hyperslab (filespaceVis, H5S_SELECT_SET, m_aChannelDataOffset, NULL, m_aMemspaceSize, NULL);
    hid_t memspace0 = H5Screate_simple (3, m_aMemspaceSize, NULL);

    herr_t err0 = H5Dwrite(m_iH5DatasetVis, H5T_NATIVE_INT, memspace0, filespaceVis, H5P_DEFAULT, &vi32Chan0.front());
    if(err0 < 0)
    {
        cout << "HDF5 chunk extend error on 1st channel: " << err0 << endl;
    }

    m_aChannelDataOffset[2] = 1;
    H5Sselect_hyperslab (filespaceVis, H5S_SELECT_SET, m_aChannelDataOffset, NULL, m_aMemspaceSize, NULL);
    hid_t memspace1 = H5Screate_simple (3, m_aMemspaceSize, NULL);

    herr_t err1 = H5Dwrite(m_iH5DatasetVis, H5T_NATIVE_INT, memspace1, filespaceVis, H5P_DEFAULT, &vi32Chan1.front());
    if(err1 < 0)
    {
        cout << "HDF5 chunk extend error on 2nd channel: " << err1 << endl;
    }

    //Stokes dataset (Q + U)
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    H5Dset_extent(m_iH5DatasetStokes, m_aChannelDatasetDims);
    hid_t filespaceStokes = H5Dget_space(m_iH5DatasetStokes);

    m_aChannelDataOffset[2] = 0;
    H5Sselect_hyperslab (filespaceStokes, H5S_SELECT_SET, m_aChannelDataOffset, NULL, m_aMemspaceSize, NULL);
    hid_t memspace2 = H5Screate_simple (3, m_aMemspaceSize, NULL);

    herr_t err2 = H5Dwrite(m_iH5DatasetStokes, H5T_NATIVE_INT, memspace2, filespaceStokes, H5P_DEFAULT, &vi32Chan2.front());
    if(err2 < 0)
    {
        cout << "HDF5 chunk extend error on 3rd channel: " << err2 << endl;
    }

    m_aChannelDataOffset[2] = 1;
    H5Sselect_hyperslab (filespaceStokes, H5S_SELECT_SET, m_aChannelDataOffset, NULL, m_aMemspaceSize, NULL);
    hid_t memspace3 = H5Screate_simple (3, m_aMemspaceSize, NULL);

    herr_t err3 = H5Dwrite(m_iH5DatasetStokes, H5T_NATIVE_INT, memspace3, filespaceStokes, H5P_DEFAULT, &vi32Chan3.front());
    if(err3 < 0)
    {
        cout << "HDF5 chunk extend error on 4th channel: " << err3 << endl;
    }

    H5Sclose (memspace0);
    H5Sclose (memspace1);
    H5Sclose (memspace2);
    H5Sclose (memspace3);

    H5Sclose (filespaceVis);
    H5Sclose (filespaceStokes);

    //Store values to be written after channel data in memory:

    //Data entry timestamps
    double dTimestamp_s = (double)oHeader.getTimestamp_us() / 1e6;
    m_vdTimestamps_s.push_back(dTimestamp_s); //As per KAT7 data: Seconds since Unix Epoch stored as double.

    //Data averages
    m_vvfChannelAverages[0].push_back( calculateFrameAverage(vi32Chan0, (AVN::Spectrometer::digitiserType)oHeader.getDigitiserType()) );
    m_vvfChannelAverages[1].push_back( calculateFrameAverage(vi32Chan1, (AVN::Spectrometer::digitiserType)oHeader.getDigitiserType()) );
    m_vvfChannelAverages[2].push_back( calculateFrameAverage(vi32Chan2, (AVN::Spectrometer::digitiserType)oHeader.getDigitiserType()) );
    m_vvfChannelAverages[3].push_back( calculateFrameAverage(vi32Chan3, (AVN::Spectrometer::digitiserType)oHeader.getDigitiserType()) );

    //Noise diode state changes
    if(oHeader.getNoiseDiodeOn() != m_oLastHeader.getNoiseDiodeOn() || !m_voNoiseDiodeStateChanges.size())
    {
        cNoiseDiodeState oState;

        oState.m_dTimeStamp_s = dTimestamp_s;

        if(oHeader.getNoiseDiodeOn())
        {
            oState.m_caValue[0] = '1';
        }
        else
        {
            oState.m_caValue[0] = '0';
        }

        sprintf(oState.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal.
        //At present AVN doesn't make hardware provision for any sort of feedback with the noise diode

        m_voNoiseDiodeStateChanges.push_back(oState);

        cout << "cSpectrometerHDF5OutputFile::addFrame(): Logged noise state change to " << string(oState.m_caValue, 1) << " at " << AVN::stringFromTimestamp_full(oHeader.getTimestamp_us()) << endl;
    }

    m_aChannelDataOffset[0] += m_aChannelDatasetExtensionDims[0];

    m_oLastHeader = oHeader;

}

void cSpectrometerHDF5OutputFile::writeTimestamps()
{
    hsize_t dimension = m_vdTimestamps_s.size();
    herr_t err = H5LTmake_dataset(m_iH5DataGroupHandle, "Timestamps", 1, &dimension, H5T_NATIVE_DOUBLE, &m_vdTimestamps_s.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeTimestamps(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeTimestamps(): Wrote " << m_vdTimestamps_s.size() << " timestamps to dataset." << endl;
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

void cSpectrometerHDF5OutputFile::writeNoiseDiodeStates()
{
    string strDatasetName("roach.noise.diode.on");

    //Create the data space
    hsize_t dimension[] = { m_voNoiseDiodeStateChanges.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cNoiseDiodeState));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cNoiseDiodeState, m_dTimeStamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the noise diode state (string of 1 character "0" or "1")
    hid_t stringTypeValue = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeValue, sizeof(cNoiseDiodeState::m_caValue));
    H5Tinsert(compoundDataType, "value", HOFFSET(cNoiseDiodeState, m_caValue), stringTypeValue);

    //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cNoiseDiodeState::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cNoiseDiodeState, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5Antenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voNoiseDiodeStateChanges.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeStates(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeStates(): Wrote " << m_voNoiseDiodeStateChanges.size() << " noise diode states to dataset." << endl;
    }

    addAttributeToDataSet(string("AVN frontend noise diode"), strDatasetName, string("boolean"), string(""), dataset);

    H5Tclose(stringTypeValue);
    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
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

void cSpectrometerHDF5OutputFile::addAttributeToDataSet(const std::string &strDescription, const std::string &strName, const std::string &strType, const std::string &strUnits, hid_t dataset)
{
    //Adds the 4 attributes as per KAT7 standard

    //Note: dataset must be open before calling this function and will remove open on function return.

    //Add attributes for noise diode state dataset
    hid_t variableLengthStringType;
    variableLengthStringType = H5Tcopy (H5T_C_S1);
    H5Tset_size (variableLengthStringType, H5T_VARIABLE);

    hid_t attrDataspace = H5Screate(H5S_SCALAR);

    hid_t attrDescription = H5Acreate(dataset, "description", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcDescription = &strDescription[0];
    H5Awrite(attrDescription, variableLengthStringType, &pcDescription);

    hid_t attrName = H5Acreate(dataset, "name", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcName = &strName[0];
    H5Awrite(attrName, variableLengthStringType, &pcName);

    hid_t attrType = H5Acreate(dataset, "type", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcType = &strType[0];
    H5Awrite(attrType, variableLengthStringType, &pcType);

    hid_t attrUnits = H5Acreate(dataset, "units", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcUnits = &strUnits[0];
    H5Awrite(attrUnits, variableLengthStringType, &pcUnits);

    H5Aclose (attrDescription);
    H5Aclose (attrName);
    H5Aclose (attrType);
    H5Aclose (attrUnits);
    H5Tclose(variableLengthStringType);
    H5Sclose(attrDataspace);
}
