
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
#include "../../AVNUtilLibs/CoordinatePosition/CoordinatePosition.h"

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

    //Setup data groups (hierachy like a directory structure within the HDF5 container):
    //Level 1:
    m_iH5DataGroupHandle                            = H5Gcreate2(m_iH5FileHandle, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5MetaDataGroupHandle                        = H5Gcreate2(m_iH5FileHandle, "/MetaData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //Level 2:
    m_iH5SensorsGroupHandle                         = H5Gcreate2(m_iH5MetaDataGroupHandle, "/MetaData/Sensors", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5ConfigurationGroupHandle                   = H5Gcreate2(m_iH5MetaDataGroupHandle, "/MetaData/Configuration", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //Level 3:
    m_iH5SensorsAntennasGroupHandle                 = H5Gcreate2(m_iH5SensorsGroupHandle, "/MetaData/Sensors/Antennas", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5SensorsRFEGroupHandle                      = H5Gcreate2(m_iH5SensorsGroupHandle, "/MetaData/Sensors/RFE", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5SensorsDBEGroupHandle                      = H5Gcreate2(m_iH5SensorsGroupHandle, "/MetaData/Sensors/DBE", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5ConfigurationAntennasGroupHandle           = H5Gcreate2(m_iH5ConfigurationGroupHandle, "/MetaData/Configuration/Antennas", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //Level 4:
    m_iH5SensorsAntennasAntenna1GroupHandle         = H5Gcreate2(m_iH5SensorsAntennasGroupHandle, "/MetaData/Sensors/Antennas/ant1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5ConfigurationAntennasAntenna1GroupHandle   = H5Gcreate2(m_iH5ConfigurationAntennasGroupHandle, "/MetaData/Configuration/Antennas/ant1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

    writeSampleDataTimestamps();
    writeChannelAverages();
    writeROACHNoiseDiodeStates();

    writeRequestedAntennaAzEls();
    writeActualAntennaAzEls();
    writeActualSourceOffsetAzEls();
    writeActualAntennaRADecs();

    writeAntennaStatuses();
    writeMotorTorques();
    writeAppliedPointingModel();

    writeNoiseDiodeSoftwareStates();
    writeNoiseDiodeSources();
    writeNoideDiodeCurrents();

    writeSelectedSources();

    writeRFFrequencies();
    writeLOFrequencies();
    writeIFBandwidths();

    writeROACHAccumulationLengths();
    writeROACHNBNarrowbandSelections();
    writeROACHSamplingFrequency();
    writeROACHSizeOfFFTs();
    writeROACHCoarseFFTShiftMask();
    writeROACHAdcAttentuations();

    cout << "cSpectrometerHDF5OutputFile::~cSpectrometerHDF5OutputFile(): Done writing accumulated data." << endl;

    //Close and free all HDF5 structures
    H5Dclose(m_iH5DatasetVis);
    H5Dclose(m_iH5DatasetStokes);

    //HDF5 groups close deepest first:
    //Level 4
    H5Gclose(m_iH5SensorsAntennasAntenna1GroupHandle);
    H5Gclose(m_iH5ConfigurationAntennasAntenna1GroupHandle);
    //Level 3
    H5Gclose(m_iH5SensorsAntennasGroupHandle);
    H5Gclose(m_iH5SensorsRFEGroupHandle);
    H5Gclose(m_iH5SensorsDBEGroupHandle);
    H5Gclose(m_iH5ConfigurationAntennasGroupHandle);
    //Level 2
    H5Gclose(m_iH5SensorsGroupHandle);
    H5Gclose(m_iH5ConfigurationGroupHandle);
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
    m_vdSampleDataTimestamps_s.push_back(dTimestamp_s); //As per KAT7 data: Seconds since Unix Epoch stored as double.

    //Data averages
    m_vvfChannelAverages[0].push_back( calculateFrameAverage(vi32Chan0) );
    m_vvfChannelAverages[1].push_back( calculateFrameAverage(vi32Chan1) );
    m_vvfChannelAverages[2].push_back( calculateFrameAverage(vi32Chan2) );
    m_vvfChannelAverages[3].push_back( calculateFrameAverage(vi32Chan3) );

    //Noise diode state changes
    if(oHeader.getNoiseDiodeOn() != m_oLastHeader.getNoiseDiodeOn() || !m_voROACHNoiseDiodeStateChanges.size())
    {
        cNoiseDiodeState oState;

        oState.m_dTimestamp_s = dTimestamp_s;

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

        m_voROACHNoiseDiodeStateChanges.push_back(oState);

        cout << "cSpectrometerHDF5OutputFile::addFrame(): Logged noise state change to " << string(oState.m_caValue, 1) << " at " << AVN::stringFromTimestamp_full(oHeader.getTimestamp_us()) << endl;
    }

    m_aChannelDataOffset[0] += m_aChannelDatasetExtensionDims[0];

    m_oLastHeader = oHeader;

}

void cSpectrometerHDF5OutputFile::writeSampleDataTimestamps()
{
    hsize_t dimension = m_vdSampleDataTimestamps_s.size();
    herr_t err = H5LTmake_dataset(m_iH5DataGroupHandle, "Timestamps", 1, &dimension, H5T_NATIVE_DOUBLE, &m_vdSampleDataTimestamps_s.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeTimestamps(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeTimestamps(): Wrote " << m_vdSampleDataTimestamps_s.size() << " timestamps to dataset." << endl;
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

void cSpectrometerHDF5OutputFile::writeROACHNoiseDiodeStates()
{
    string strDatasetName("noise-diode.roach.on");

    //Create the data space
    hsize_t dimension[] = { m_voROACHNoiseDiodeStateChanges.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cNoiseDiodeState));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cNoiseDiodeState, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the noise diode state (string of 1 character "0" or "1")
    hid_t stringTypeValue = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeValue, sizeof(cNoiseDiodeState::m_caValue));
    H5Tinsert(compoundDataType, "value", HOFFSET(cNoiseDiodeState, m_caValue), stringTypeValue);

    //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cNoiseDiodeState::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cNoiseDiodeState, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voROACHNoiseDiodeStateChanges.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeStates(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeStates(): Wrote " << m_voROACHNoiseDiodeStateChanges.size() << " noise diode states to dataset." << endl;
    }

    addAttributeToDataSet(string("AVN frontend noise diode"), strDatasetName, string("boolean"), string(""), dataset);

    H5Tclose(stringTypeValue);
    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls()
{
    //Requested values to motor drive after pointing model
    //As per KAT7 Azimuth and elevation are in seperate datasets

    //Azimuth:
    ////////////////////////////////////////////////////////////////////////////////////////////
    {
        string strDatasetName("pos.actual-pointm-azim");

        //Create the data space
        hsize_t dimension[] = { m_voRequestedAntennaAzs_deg.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the azimuth value (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the elevation sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_caStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_caStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voRequestedAntennaAzs_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls(): Wrote " << m_voRequestedAntennaAzs_deg.size() << " requested antenna azimuths to dataset." << endl;
        }

        addAttributeToDataSet(string("Requested azimuth after pointing model"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    //Elevation:
    ////////////////////////////////////////////////////////////////////////////////////////////
    {
        string strDatasetName("pos.actual-pointm-elev");

        //Create the data space
        hsize_t dimension[] = { m_voRequestedAntennaEls_deg.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the elevation value (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the elevation sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_caStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_caStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voRequestedAntennaEls_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls(): Wrote " << m_voRequestedAntennaEls_deg.size() << " requested antenna elevations to dataset." << endl;
        }

        addAttributeToDataSet(string("Requested elevation after pointing model"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeActualAntennaAzEls()
{
    //Values return from drive about current Az/El
    //Apply inverse pointing model to get to ideal Az/El

    //Requested values to motor drive after pointing model
    //As per KAT7 Azimuth and elevation are in seperate datasets

    //Azimuth:
    ////////////////////////////////////////////////////////////////////////////////////////////
    {
        string strDatasetName("pos.request-pointm-azim");

        //Create the data space
        hsize_t dimension[] = { m_voActualAntennaAzs_deg.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the azimuth value (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the elevation sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_caStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_caStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voActualAntennaAzs_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaAzEls(): Wrote " << m_voActualAntennaAzs_deg.size() << " actual antenna azimuths to dataset." << endl;
        }

        addAttributeToDataSet(string("Actual azimuth after pointing model"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    //Elevation:
    ////////////////////////////////////////////////////////////////////////////////////////////
    {
        string strDatasetName("pos.request-pointm-elev");

        //Create the data space
        hsize_t dimension[] = { m_voActualAntennaEls_deg.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the elevation value (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the elevation sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_caStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_caStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voActualAntennaEls_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaAzEls(): Wrote " << m_voActualAntennaEls_deg.size() << " actual antenna elevations to dataset." << endl;
        }

        addAttributeToDataSet(string("Actual elevation after pointing model"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeActualSourceOffsetAzEls()
{
    string strDatasetName("pos.source-offset-azim-elev");

    //Create the data space
    hsize_t dimension[] = { m_voActualSourceOffsetAzEls_deg.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDualDouble));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDualDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the azimuth value (double)
    H5Tinsert(compoundDataType, "azimuth offset", HOFFSET(cTimestampedDualDouble, m_dValue1), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the azimuth value (double)
    H5Tinsert(compoundDataType, "elevation offset", HOFFSET(cTimestampedDualDouble, m_dValue2), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedDualDouble::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDualDouble, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voActualSourceOffsetAzEls_deg.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeActualSourceOffsetAzEls(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeActualSourceOffsetAzEls(): Wrote " << m_voActualSourceOffsetAzEls_deg.size() << " actual source offsets to dataset." << endl;
    }

    addAttributeToDataSet(string("Azimuth/elevation offset from the source after pointing model"), strDatasetName, string("double"), string("deg"), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeActualAntennaRADecs()
{
    string strDatasetName("pos.actual-ra-dec");

    //Create the data space
    hsize_t dimension[] = { m_voActualAntennaRADecs_deg.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDualDouble));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDualDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the RA value (double)
    H5Tinsert(compoundDataType, "right ascension value", HOFFSET(cTimestampedDualDouble, m_dValue1), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the Dec value (double)
    H5Tinsert(compoundDataType, "declination value", HOFFSET(cTimestampedDualDouble, m_dValue2), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedDualDouble::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDualDouble, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voActualAntennaRADecs_deg.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeActualAntennaRADecs(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeActualAntennaRADecs(): Wrote " << m_voActualAntennaRADecs_deg.size() << " actual RA and Decs to dataset." << endl;
    }

    addAttributeToDataSet(string("Actual right ascention/declination derived from actual azimuth/elevation after pointing model"), strDatasetName, string("double"), string("deg"), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeAntennaStatuses()
{
    string strDatasetName("activity");

    //Create the data space
    hsize_t dimension[] = { m_voAntennaStatuses.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cAntennaStatus));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cAntennaStatus, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the antenna state (string)
    hid_t stringTypeValue = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeValue, sizeof(cAntennaStatus::m_caAntennaStatus));
    H5Tinsert(compoundDataType, "value", HOFFSET(cAntennaStatus, m_caAntennaStatus), stringTypeValue);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cAntennaStatus::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cAntennaStatus, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voAntennaStatuses.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeAntennaStatuses(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeAntennaStatuses(): Wrote " << m_voAntennaStatuses.size() << " antenna statuses to dataset." << endl;
    }

    addAttributeToDataSet(string("Combined state of the antenna proxy"), strDatasetName, string("discrete"), string(""), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeMotorTorques()
{
    string strDatasetName("motor-torques");

    //Create the data space
    hsize_t dimension[] = { m_voMotorTorques_Nm.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cMotorTorques));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cMotorTorques, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data types: Torque for azimuth motors 0, 1 and elevation motors 0, 1 (all doubles)
    H5Tinsert(compoundDataType, "azim motor 0", HOFFSET(cMotorTorques, m_dAz0_Nm), H5T_NATIVE_DOUBLE);
    H5Tinsert(compoundDataType, "azim motor 1", HOFFSET(cMotorTorques, m_dAz1_Nm), H5T_NATIVE_DOUBLE);
    H5Tinsert(compoundDataType, "elev motor 0", HOFFSET(cMotorTorques, m_dEl0_Nm), H5T_NATIVE_DOUBLE);
    H5Tinsert(compoundDataType, "elev motor 1", HOFFSET(cMotorTorques, m_dEl1_Nm), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cMotorTorques::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cMotorTorques, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voMotorTorques_Nm.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeMotorTorques(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeMotorTorques(): Wrote " << m_voMotorTorques_Nm.size() << " antenna statuses to dataset." << endl;
    }

    addAttributeToDataSet(string("Torques of motors steering the antenna"), strDatasetName, string("double"), string("Nm"), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeAppliedPointingModel()
{
    hsize_t dimension = m_vdPointingModelParams.size();
    string strDatasetName("pointing-model-params");

    herr_t err = H5LTmake_dataset(m_iH5ConfigurationAntennasAntenna1GroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_DOUBLE, &m_vdPointingModelParams.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeAppliedPointingModel(): HDF5 make dataset error." << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeAppliedPointingModel(): Wrote " << m_vdPointingModelParams.size() << " pointing model parameters to file." << endl;
    }

    //Add pointing model name as data attribute
    stringstream oSS;
    oSS << "Pointing model name = ";
    oSS << m_strPointModelName;

    //Need to open the dataset again here as the H5LTmake_dataset used above does leave an open handle.
    hid_t dataset_id = H5Dopen2(m_iH5ConfigurationAntennasAntenna1GroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
    addAttributeToDataSet(oSS.str(), strDatasetName, string("string"), string(""), dataset_id);

    H5Dclose(dataset_id);
}

void cSpectrometerHDF5OutputFile::writeNoiseDiodeSoftwareStates()
{
    string strDatasetName("noise-diode.software.on");

    //Create the data space
    hsize_t dimension[] = { m_voNoiseDiodeSoftwareStates.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cNoiseDiodeState));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cNoiseDiodeState, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the noise diode state (string of 1 character "0" or "1")
    hid_t stringTypeValue = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeValue, sizeof(cNoiseDiodeState::m_caValue));
    H5Tinsert(compoundDataType, "value", HOFFSET(cNoiseDiodeState, m_caValue), stringTypeValue);

    //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cNoiseDiodeState::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cNoiseDiodeState, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voNoiseDiodeSoftwareStates.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeSoftwareStates(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeSoftwareStates(): Wrote " << m_voNoiseDiodeSoftwareStates.size() << " noise diode states to dataset." << endl;
    }

    addAttributeToDataSet(string("software controlled noise diode state"), strDatasetName, string("boolean"), string(""), dataset);

    H5Tclose(stringTypeValue);
    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeNoiseDiodeSources()
{
    string strDatasetName("noise-diode.source");

    //Create the data space
    hsize_t dimension[] = { m_voNoiseDiodeSources.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cNoiseDiodeSource));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cNoiseDiodeSource, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the noise diode source (c string)
    hid_t stringTypeValue = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeValue, sizeof(cNoiseDiodeSource::m_caSource));
    H5Tinsert(compoundDataType, "value", HOFFSET(cNoiseDiodeSource, m_caSource), stringTypeValue);

    //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cNoiseDiodeSource::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cNoiseDiodeSource, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voNoiseDiodeSources.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeSources(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeSources(): Wrote " << m_voNoiseDiodeSources.size() << " noise diode sources to dataset." << endl;
    }

    addAttributeToDataSet(string("source of noise diode control"), strDatasetName, string("discrete"), string(""), dataset);

    H5Tclose(stringTypeValue);
    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeNoideDiodeCurrents()
{
    string strDatasetName("noise-diode.current");

    //Create the data space
    hsize_t dimension[] = { m_voNoiseDiodeCurrents.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the noise diode current (double)
    H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cNoiseDiodeState::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voNoiseDiodeCurrents.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeNoideDiodeCurrents(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeNoideDiodeCurrents(): Wrote " << m_voNoiseDiodeCurrents.size() << " noise diode currents to dataset." << endl;
    }

    addAttributeToDataSet(string("current draw of noise diode"), strDatasetName, string("double"), string("A"), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeSelectedSources()
{
    string strDatasetName("target");

    //Create the data space
    hsize_t dimension[] = { m_voSelectedSources.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cSourceSelection));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cSourceSelection, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the source name and ra/dec (c string)
    hid_t stringTypeValue = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeValue, sizeof(cSourceSelection::m_caSource));
    H5Tinsert(compoundDataType, "value", HOFFSET(cSourceSelection, m_caSource), stringTypeValue);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cSourceSelection::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cSourceSelection, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voSelectedSources.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeSelectedSources(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeSelectedSources(): Wrote " << m_voSelectedSources.size() << " astronomical sources to dataset." << endl;
    }

    addAttributeToDataSet(string("target"), strDatasetName, string("string"), string(""), dataset);

    H5Tclose(stringTypeValue);
    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeRFFrequencies()
{
    string strDatasetName("rfe.rf.frequency");

    //Create the data space
    hsize_t dimension[] = { m_voFrequenciesRF_MHz.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the RF frequency (double)
    H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voFrequenciesRF_MHz.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeRFFrequencies(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeRFFrequencies(): Wrote " << m_voFrequenciesLOs_MHz.size() << " RF frequencies." << endl;
    }

    addAttributeToDataSet(string("RF centre frequency"), strDatasetName, string("double"), string("MHz"), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeLOFrequencies()
{
    string strDatasetName("rfe.lo.frequencies");

    //Create the data space
    hsize_t dimension[] = { m_voFrequenciesLOs_MHz.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDualDouble));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDualDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the LO frequencies (double)
    H5Tinsert(compoundDataType, "LO1", HOFFSET(cTimestampedDualDouble, m_dValue1), H5T_NATIVE_DOUBLE);
    H5Tinsert(compoundDataType, "LO2", HOFFSET(cTimestampedDualDouble, m_dValue2), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedDualDouble::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDualDouble, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voFrequenciesLOs_MHz.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): Wrote " << m_voFrequenciesLOs_MHz.size() << " LO frequency pairs." << endl;
    }

    addAttributeToDataSet(string("frequencies of the LOs"), strDatasetName, string("double"), string("MHz"), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeIFBandwidths()
{
    string strDatasetName("rfe.if.bandwidth");

    //Create the data space
    hsize_t dimension[] = { m_voBandwidthsIF_MHz.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the If bandwidth (double)
    H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voBandwidthsIF_MHz.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeIFBandwidths(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeIFBandwidths(): Wrote " << m_voBandwidthsIF_MHz.size() << " IF bandwidths." << endl;
    }

    addAttributeToDataSet(string("IF bandwidth"), strDatasetName, string("double"), string("MHz"), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeROACHAccumulationLengths()
{
    string strDatasetName("dbe.if.bandwidth");

    //Create the data space
    hsize_t dimension[] = { m_voROACHAccumulationLengths_nFrames.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedUnsignedInt));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedUnsignedInt, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the accumulation length (unsigned int)
    H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedUnsignedInt, m_u32Value), H5T_NATIVE_UINT32);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedUnsignedInt::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedUnsignedInt, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voROACHAccumulationLengths_nFrames.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHAccumulationLengths(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHAccumulationLengths(): Wrote " << m_voROACHAccumulationLengths_nFrames.size() << " accumulation lengths." << endl;
    }

    addAttributeToDataSet(string("accumulation length"), strDatasetName, string("unsigned int"), string("no. of FFT frames"), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeROACHNBNarrowbandSelections()
{
    string strDatasetName("dbe.nb-chan");

    //Create the data space
    hsize_t dimension[] = { m_voROACHNBChannelSelects.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedUnsignedInt));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedUnsignedInt, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the narrowband channel no (unsigned int)
    H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedUnsignedInt, m_u32Value), H5T_NATIVE_UINT32);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedUnsignedInt::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedUnsignedInt, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voROACHNBChannelSelects.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHNBNarrowbandSelections(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHNBNarrowbandSelections(): Wrote " << m_voROACHNBChannelSelects.size() << " narrow band channel selections." << endl;
    }

    addAttributeToDataSet(string("narrowband channel selection"), strDatasetName, string("unsigned int"), string("cannonical bin no."), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency()
{
    hsize_t dimension = 1;
    string strDatasetName("dbe.fs");

    herr_t err = H5LTmake_dataset(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_DOUBLE, &m_dROACHFrequencyFs_MHz);

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): HDF5 make dataset error." << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): Wrote ROACH Fs." << endl;
    }

    //Need to open the dataset again here for attribute as the H5LTmake_dataset used above does leave an open handle.
    hid_t dataset_id = H5Dopen2(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
    addAttributeToDataSet(string("sampling frequency of ROACH"), strDatasetName, string("double"), string("MHz"), dataset_id);

    H5Dclose(dataset_id);
}

void cSpectrometerHDF5OutputFile::writeROACHSizeOfFFTs()
{
    string strDatasetName("rfe.fft.sizes");

    //Create the data space
    hsize_t dimension[] = { 1 };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cSizeOfFFTs));

    //Add to compound data type: the FFT sizes (unsigned int)
    H5Tinsert(compoundDataType, "Coarse", HOFFSET(cSizeOfFFTs, m_u32CoarseFFTSize), H5T_NATIVE_UINT32);
    H5Tinsert(compoundDataType, "Fine", HOFFSET(cSizeOfFFTs, m_u32FineFFTSize), H5T_NATIVE_UINT32);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cSizeOfFFTs::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cSizeOfFFTs, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_oROACHSizesOfFFTs_nSamp);

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSizeOfFFTs(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSizeOfFFTs(): Wrote size of FFTs." << endl;
    }

    addAttributeToDataSet(string("size of FFTs"), strDatasetName, string("unsigned int"), string("no. of input samples"), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeROACHCoarseFFTShiftMask()
{
    string strDatasetName("dbe.coarse-fft-shift-mask");

    //Create the data space
    hsize_t dimension[] = { m_voROACHCoarseFFTShiftMasks.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedUnsignedInt));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedUnsignedInt, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the shift mask applied to the course FFT no (unsigned int)
    H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedUnsignedInt, m_u32Value), H5T_NATIVE_UINT32);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedUnsignedInt::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedUnsignedInt, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voROACHCoarseFFTShiftMasks.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHCoarseFFTShiftMask(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHCoarseFFTShiftMask(): Wrote " << m_voROACHCoarseFFTShiftMasks.size() << " coarse FFT shift masks." << endl;
    }

    addAttributeToDataSet(string("the shift mask applied to the course FFT"), strDatasetName, string("unsigned int"), string(""), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeROACHAdcAttentuations()
{
    string strDatasetName("dbe.adc.attenuations");

    //Create the data space
    hsize_t dimension[] = { m_voROACHAdcAttenuations_dB.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDualDouble));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDualDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the ADC attenuations (double)
    H5Tinsert(compoundDataType, "Chan0", HOFFSET(cTimestampedDualDouble, m_dValue1), H5T_NATIVE_DOUBLE);
    H5Tinsert(compoundDataType, "Chan1", HOFFSET(cTimestampedDualDouble, m_dValue2), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedDualDouble::m_caStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDualDouble, m_caStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voROACHAdcAttenuations_dB.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): Wrote " << m_voROACHAdcAttenuations_dB.size() << " ADC attenuation pairs." << endl;
    }

    addAttributeToDataSet(string("attenuation of input to ADC"), strDatasetName, string("double"), string("dB"), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

float cSpectrometerHDF5OutputFile::calculateFrameAverage(const vector<int32_t> &vi32ChannelData)
{
    double dAverage = 0.0;

    for(uint32_t ui = 0; ui < vi32ChannelData.size(); ui++)
    {
        dAverage += vi32ChannelData[ui];
    }

    dAverage /= (double)vi32ChannelData.size();

    return (float)dAverage;
}

void cSpectrometerHDF5OutputFile::addAttributeToDataSet(const std::string &strDescription, const std::string &strName, const std::string &strType, const std::string &strUnits, hid_t dataset)
{
    //Adds the 4 attributes as per KAT7 standard

    //Note: dataset must be open before calling this function and will remain open on function return.

    //Add attributes for dataset
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

//Functions for adding logged data from outside this class
//Each function aquires a shared lock of a common shared mutex for adding data to its respective vector.
//When all of this data is written to file a unique lock is obtained from that mutex so that none of the data be altered at that point
void cSpectrometerHDF5OutputFile::addRequestedAntennaAzEl(int64_t i64Timestamp_us, double dAzimuth_deg, double dElevation_deg)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewRequestedAntennaAz;
    oNewRequestedAntennaAz.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewRequestedAntennaAz.m_dValue = dAzimuth_deg;
    sprintf(oNewRequestedAntennaAz.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    cTimestampedDouble oNewRequestedAntennaEl;
    oNewRequestedAntennaEl.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewRequestedAntennaEl.m_dValue = dElevation_deg;
    sprintf(oNewRequestedAntennaEl.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voRequestedAntennaAzs_deg.push_back(oNewRequestedAntennaAz);
    m_voRequestedAntennaEls_deg.push_back(oNewRequestedAntennaEl);
}

void cSpectrometerHDF5OutputFile::addActualAntennaAzEl(int64_t i64Timestamp_us, double dAzimuth_deg, double dElevation_deg)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewActualAntennaAz;
    oNewActualAntennaAz.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualAntennaAz.m_dValue = dAzimuth_deg;
    sprintf(oNewActualAntennaAz.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.


    cTimestampedDouble oNewActualAntennaEl;
    oNewActualAntennaEl.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualAntennaEl.m_dValue = dElevation_deg;
    sprintf(oNewActualAntennaAz.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.


    m_voActualAntennaAzs_deg.push_back(oNewActualAntennaAz);
    m_voActualAntennaEls_deg.push_back(oNewActualAntennaEl);
}

void cSpectrometerHDF5OutputFile::addActualSourceOffsetAzEl(int64_t i64Timestamp_us, double dAzimuthOffset_deg, double dElevationOffset_deg)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDualDouble oNewActualSourceOffsetAzEl;
    oNewActualSourceOffsetAzEl.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualSourceOffsetAzEl.m_dValue1 = dAzimuthOffset_deg;
    oNewActualSourceOffsetAzEl.m_dValue2 = dElevationOffset_deg;
    sprintf(oNewActualSourceOffsetAzEl.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voActualSourceOffsetAzEls_deg.push_back(oNewActualSourceOffsetAzEl);
}

void cSpectrometerHDF5OutputFile::addActualAntennaRADec(int64_t i64Timestamp_us, double dRighAscension_deg, double dDeclination_deg)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDualDouble oNewActualAntennaRADec;
    oNewActualAntennaRADec.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualAntennaRADec.m_dValue1 = dRighAscension_deg;
    oNewActualAntennaRADec.m_dValue2 = dDeclination_deg;
    sprintf(oNewActualAntennaRADec.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voActualAntennaRADecs_deg.push_back(oNewActualAntennaRADec);
}

void cSpectrometerHDF5OutputFile::addAntennaStatus(int64_t i64Timestamp_us, int32_t i32AntennaStatus, const string &strAntennaStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cAntennaStatus oNewAntennaStatus;
    oNewAntennaStatus.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    sprintf( oNewAntennaStatus.m_caAntennaStatus, strAntennaStatus.substr(0, sizeof(oNewAntennaStatus.m_caAntennaStatus)).c_str() ); //Limit to size of the char array
    sprintf( oNewAntennaStatus.m_caStatus, "nominal" ); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    //The numerical value is not used here

    m_voAntennaStatuses.push_back(oNewAntennaStatus);
}

void cSpectrometerHDF5OutputFile::addMotorTorques(int64_t i64Timestamp_us, double dAz0_Nm, double dAz1_Nm, double dEl0_Nm, double dEl1_Nm)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cMotorTorques oNewMotorTorques;
    oNewMotorTorques.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewMotorTorques.m_dAz0_Nm = dAz0_Nm;
    oNewMotorTorques.m_dAz0_Nm = dAz1_Nm;
    oNewMotorTorques.m_dAz0_Nm = dEl0_Nm;
    oNewMotorTorques.m_dAz0_Nm = dEl1_Nm;
    sprintf(oNewMotorTorques.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voMotorTorques_Nm.push_back(oNewMotorTorques);
}

void cSpectrometerHDF5OutputFile::setAppliedPointingModel(const string &strModelName, const vector<double> &vdPointingModelParams)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    //This only a single set of values, not a history of values.
    //Update the current record whenever a new value is received

    m_vdPointingModelParams = vdPointingModelParams;
    m_strPointModelName = strModelName;
}

void cSpectrometerHDF5OutputFile::addNoiseDiodeSoftwareState(int64_t i64Timestamp_us, int32_t i32NoiseDiodeState)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedInt oNewNoiseDiodeState;
    oNewNoiseDiodeState.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewNoiseDiodeState.m_i32Value = i32NoiseDiodeState;
    sprintf(oNewNoiseDiodeState.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voNoiseDiodeSoftwareStates.push_back(oNewNoiseDiodeState);
}

void cSpectrometerHDF5OutputFile::addNoiseDiodeSource(int64_t i64Timestamp_us, int32_t i32NoiseDiodeSource, const string &strNoiseSource)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cNoiseDiodeSource oNewNoiseDiodeSource;
    oNewNoiseDiodeSource.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    sprintf( oNewNoiseDiodeSource.m_caSource, strNoiseSource.substr(0, sizeof(oNewNoiseDiodeSource.m_caSource)).c_str() ); //Limit to size of the char array
    sprintf( oNewNoiseDiodeSource.m_caStatus, "nominal" ); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voNoiseDiodeSources.push_back(oNewNoiseDiodeSource);
}

void cSpectrometerHDF5OutputFile::addNoiseDiodeCurrent(int64_t i64Timestamp_us, double dNoiseDiodeCurrent_A)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewNoiseDiodeCurrent;
    oNewNoiseDiodeCurrent.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewNoiseDiodeCurrent.m_dValue = dNoiseDiodeCurrent_A;
    sprintf(oNewNoiseDiodeCurrent.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.


    m_voNoiseDiodeCurrents.push_back(oNewNoiseDiodeCurrent);
}

void cSpectrometerHDF5OutputFile::addSourceSelection(int64_t i64Timestamp_us, const string &strSourceName, double dRighAscension_deg, double dDeclination_deg)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    //Construct as per KAT7 standard
    //Source name if available, RA and Dec in DMS

    int32_t i32RA_deg;
    int32_t i32RA_min;
    double dRA_s;
    AVN::cCoordinatePosition::decimalDegreesToDMS(dRighAscension_deg, i32RA_deg, i32RA_min, dRA_s);

    int32_t i32Dec_deg;
    int32_t i32Dec_min;
    double dDec_s;
    AVN::cCoordinatePosition::decimalDegreesToDMS(dDeclination_deg, i32Dec_deg, i32Dec_min, dDec_s);

    stringstream oSS;
    if(strSourceName.length())
    {
        oSS << strSourceName;
        oSS << ", ";
    }

    oSS << "radec, ";
    oSS << i32RA_deg << ":" << i32RA_min << ":" << dRA_s << ", ";
    oSS << i32Dec_deg << ":" << i32Dec_min << ":" << dDec_s;

    cSourceSelection oNewSourceSelection;

    oNewSourceSelection.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    sprintf(oNewSourceSelection.m_caSource, oSS.str().substr(0, sizeof(oNewSourceSelection.m_caSource)).c_str() ); //Limit to size of the char array
    sprintf(oNewSourceSelection.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voSelectedSources.push_back(oNewSourceSelection);
}


void cSpectrometerHDF5OutputFile::addFrequencyRF(int64_t i64Timestamp_us, double dFreqencyRF_MHz)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewRFFrequency;
    oNewRFFrequency.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewRFFrequency.m_dValue = dFreqencyRF_MHz;
    sprintf(oNewRFFrequency.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voFrequenciesRF_MHz.push_back(oNewRFFrequency);
}

void cSpectrometerHDF5OutputFile::addFrequencyLOs(int64_t i64Timestamp_us, double dFrequencyLO1_MHz, double dFrequencyLO2_MHz)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDualDouble oNewLOFrequencies;
    oNewLOFrequencies.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewLOFrequencies.m_dValue1 = dFrequencyLO1_MHz;
    oNewLOFrequencies.m_dValue2 = dFrequencyLO2_MHz;
    sprintf(oNewLOFrequencies.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voFrequenciesLOs_MHz.push_back(oNewLOFrequencies);
}

void cSpectrometerHDF5OutputFile::addBandwidthIF(int64_t i64Timestamp_us, double dBandwidthIF_MHz)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewBandwidthIF;
    oNewBandwidthIF.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewBandwidthIF.m_dValue = dBandwidthIF_MHz;
    sprintf(oNewBandwidthIF.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voBandwidthsIF_MHz.push_back(oNewBandwidthIF);
}

void cSpectrometerHDF5OutputFile::addAccumulationLength(int64_t i64Timestamp_us, uint32_t u32NFrames)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedUnsignedInt oNewAccumulationLength;
    oNewAccumulationLength.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewAccumulationLength.m_u32Value = u32NFrames;
    sprintf(oNewAccumulationLength.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voROACHAccumulationLengths_nFrames.push_back(oNewAccumulationLength);
}

void cSpectrometerHDF5OutputFile::addNarrowBandChannelSelect(int64_t i64Timestamp_us, uint32_t u32ChannelNo)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedUnsignedInt oNewNarrowBandChannelSelection;
    oNewNarrowBandChannelSelection.m_dTimestamp_s = i64Timestamp_us;
    oNewNarrowBandChannelSelection.m_u32Value = u32ChannelNo;
    sprintf(oNewNarrowBandChannelSelection.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voROACHNBChannelSelects.push_back(oNewNarrowBandChannelSelection);
}

void cSpectrometerHDF5OutputFile::setFrequencyFs(double dFrequencyFs_MHz)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    //This only a single value, not a history of values.
    //Update the current record whenever a new value is received

    m_dROACHFrequencyFs_MHz = dFrequencyFs_MHz;
}

void cSpectrometerHDF5OutputFile::setSizeOfFFTs(uint32_t u32CoarseSize_nSamp, uint32_t u32FineSize_nSamp)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_oROACHSizesOfFFTs_nSamp.m_u32CoarseFFTSize = u32CoarseSize_nSamp;
    m_oROACHSizesOfFFTs_nSamp.m_u32FineFFTSize = u32FineSize_nSamp;
}

void cSpectrometerHDF5OutputFile::addCoarseFFTShiftMask(int64_t i64Timestamp_us, uint32_t u32ShiftMask)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedUnsignedInt oNewCoarseFFTShift;
    oNewCoarseFFTShift.m_dTimestamp_s = i64Timestamp_us;
    oNewCoarseFFTShift.m_u32Value = u32ShiftMask;
    sprintf(oNewCoarseFFTShift.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voROACHCoarseFFTShiftMasks.push_back(oNewCoarseFFTShift);
}

void cSpectrometerHDF5OutputFile::addAdcAttenuation(int64_t i64Timestamp_us, double dAttenuationChan0_dB, double dAttenuationChan1_dB)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDualDouble oNewAdcAttenuantion;
    oNewAdcAttenuantion.m_dTimestamp_s = i64Timestamp_us;
    oNewAdcAttenuantion.m_dValue1 = dAttenuationChan0_dB;
    oNewAdcAttenuantion.m_dValue2 = dAttenuationChan1_dB;
    sprintf(oNewAdcAttenuantion.m_caStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voROACHAdcAttenuations_dB.push_back(oNewAdcAttenuantion);
}
