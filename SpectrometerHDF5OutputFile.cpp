
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
    m_iH5MarkupGroupHandle                          = H5Gcreate2(m_iH5FileHandle, "/Markup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //Level 2:
    m_iH5SensorsGroupHandle                         = H5Gcreate2(m_iH5MetaDataGroupHandle, "/MetaData/Sensors", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5ConfigurationGroupHandle                   = H5Gcreate2(m_iH5MetaDataGroupHandle, "/MetaData/Configuration", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //Level 3:
    m_iH5SensorsAntennasGroupHandle                 = H5Gcreate2(m_iH5SensorsGroupHandle, "/MetaData/Sensors/Antennas", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5SensorsRFEGroupHandle                      = H5Gcreate2(m_iH5SensorsGroupHandle, "/MetaData/Sensors/RFE", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5SensorsDBEGroupHandle                      = H5Gcreate2(m_iH5SensorsGroupHandle, "/MetaData/Sensors/DBE", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5ConfigurationAntennasGroupHandle           = H5Gcreate2(m_iH5ConfigurationGroupHandle, "/MetaData/Configuration/Antennas", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5ConfigurationObservationGroupHandle        = H5Gcreate2(m_iH5MetaDataGroupHandle, "/MetaData/Configuration/Observation", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_iH5ConfigurationDBEGroupHandle                = H5Gcreate2(m_iH5MetaDataGroupHandle, "/MetaData/Configuration/DBE", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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

    // TODO: This data needs to be retrieved from the KatCP server. Somewhere.
    string strExperimentID("test_experiment");
    addAttributesToFile("2.5", strExperimentID, AVN::getTimeNow_us(), 0, m_iH5FileHandle); // Zero errors, because our software is perfect!
    addAttributesToObservation("script name", "script arguments", "Observer name", strExperimentID, "description", "ant1", "start time", "end time", "noise diode params", "rf params", "status", m_iH5ConfigurationObservationGroupHandle);
    addAttributesToDBE("ll,rr", "q,u", m_iH5ConfigurationDBEGroupHandle);

    writeMarkupLabels();

    writeSampleDataTimestamps();
    writeChannelAverages();
    writeROACHNoiseDiodeStates();

    writeRequestedAntennaAzEls();
    writeActualAntennaAzEls();
    writeActualSourceOffsetAzEls();
    writeActualAntennaRADecs();

    writeAntennaStatuses();
    writeMotorTorques();
    writeAntennaConfiguration();

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
    H5Gclose(m_iH5ConfigurationDBEGroupHandle);
    //Level 3
    H5Gclose(m_iH5SensorsAntennasGroupHandle);
    H5Gclose(m_iH5SensorsRFEGroupHandle);
    H5Gclose(m_iH5SensorsDBEGroupHandle);
    H5Gclose(m_iH5ConfigurationAntennasGroupHandle);
    H5Gclose(m_iH5ConfigurationObservationGroupHandle);
    //Level 2
    H5Gclose(m_iH5SensorsGroupHandle);
    H5Gclose(m_iH5ConfigurationGroupHandle);
    //Level 1
    H5Gclose(m_iH5MarkupGroupHandle);
    H5Gclose(m_iH5MetaDataGroupHandle);
    H5Gclose(m_iH5DataGroupHandle);

    //Level 0 (file)
    H5Fclose(m_iH5FileHandle);
}

string cSpectrometerHDF5OutputFile::getFilename() const
{
    return m_strFilename;
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
    //(Not written directly to disk because they will be interleaved on the disk and cause long read times.)

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
            oState.m_chaValue[0] = '1';
        }
        else
        {
            oState.m_chaValue[0] = '0';
        }

        sprintf(oState.m_chaStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal.
        //At present AVN doesn't make hardware provision for any sort of feedback with the noise diode

        m_voROACHNoiseDiodeStateChanges.push_back(oState);
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
    //This vector should always have at least one vector in it, and the noise diode
    //should by default be off (handled elsewhere).
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
    H5Tset_size(stringTypeValue, sizeof(cNoiseDiodeState::m_chaValue));
    H5Tinsert(compoundDataType, "value", HOFFSET(cNoiseDiodeState, m_chaValue), stringTypeValue);

    //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cNoiseDiodeState::m_chaStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cNoiseDiodeState, m_chaStatus), stringTypeStatus);

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

void cSpectrometerHDF5OutputFile::writeMarkupLabels()
{
    //If the vector is empty, no point in adding it to the HDF5 file.
    if (m_voMarkupLabels.size())
    {
        string strDatasetName("labels");

        //Create the data space TODO
        hsize_t dimension[] = { m_voMarkupLabels.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cMarkupLabels));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cMarkupLabels, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the label (string)
        hid_t stringTypeLabel = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeLabel, sizeof(cMarkupLabels::m_chaLabel));
        H5Tinsert(compoundDataType, "label", HOFFSET(cMarkupLabels, m_chaLabel), stringTypeLabel);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5MarkupGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voMarkupLabels.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeMarkupLabels(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeMarkupLabels(): Wrote " << m_voMarkupLabels.size() << " markup labels to dataset." << endl;
        }

        H5Tclose(stringTypeLabel);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls()
{
    //Requested values to motor drive after pointing model
    //As per KAT7 Azimuth and elevation are in seperate datasets

    //Azimuth:
    ////////////////////////////////////////////////////////////////////////////////////////////
    if (m_voRequestedAntennaAzs_deg.size())
    {
        string strDatasetName("pos.request-scan-azim");

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
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

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

        addAttributeToDataSet(string("Requested azimuth after scan offset"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    //Elevation:
    ////////////////////////////////////////////////////////////////////////////////////////////
    if (m_voRequestedAntennaEls_deg.size())
    {
        string strDatasetName("pos.request-scan-elev");

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
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

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

        addAttributeToDataSet(string("Requested elevation after scan offset"), strDatasetName, string("double"), string("deg"), dataset);

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
    if (m_voActualAntennaAzs_deg.size())
    {
        string strDatasetName("pos.actual-scan-azim");

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
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

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

        addAttributeToDataSet(string("Actual azimuth after scan offset"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    //Elevation:
    ////////////////////////////////////////////////////////////////////////////////////////////
    if (m_voActualAntennaEls_deg.size())
    {
        string strDatasetName("pos.actual-scan-elev");

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
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

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

        addAttributeToDataSet(string("Actual elevation after scan offset"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeActualSourceOffsetAzEls()
{
    if (m_voActualSourceOffsetAzs_deg.size())
    {
        string strDatasetName("pos.source-offset-azim");

        //Create the data space
        hsize_t dimension[] = { m_voActualSourceOffsetAzs_deg.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the azimuth value (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voActualSourceOffsetAzs_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualSourceOffsetAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualSourceOffsetAzEls(): Wrote " << m_voActualSourceOffsetAzs_deg.size() << " actual source azimuth offsets to dataset." << endl;
        }

        addAttributeToDataSet(string("Azimuth offset from the source after pointing model"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voActualSourceOffsetEls_deg.size())
    {
        string strDatasetName("pos.source-offset-elev");

        //Create the data space
        hsize_t dimension[] = { m_voActualSourceOffsetEls_deg.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the elevation value (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voActualSourceOffsetEls_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualSourceOffsetAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualSourceOffsetAzEls(): Wrote " << m_voActualSourceOffsetEls_deg.size() << " actual source elevation offsets to dataset." << endl;
        }

        addAttributeToDataSet(string("Elevation offset from the source after pointing model"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeActualAntennaRADecs()
{
    if (m_voActualAntennaRAs_deg.size())
    {
        string strDatasetName("pos.actual-ra");

        //Create the data space
        hsize_t dimension[] = { m_voActualAntennaRAs_deg.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the RA value (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voActualAntennaRAs_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaRADecs(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaRADecs(): Wrote " << m_voActualAntennaRAs_deg.size() << " actual RAs to dataset." << endl;
        }

        addAttributeToDataSet(string("Actual right ascention derived from actual azimuth/elevation after pointing model"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voActualAntennaRAs_deg.size())
    {
        string strDatasetName("pos.actual-dec");

        //Create the data space
        hsize_t dimension[] = { m_voActualAntennaRAs_deg.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the Dec value (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voActualAntennaRAs_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaRADecs(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaRADecs(): Wrote " << m_voActualAntennaRAs_deg.size() << " actual Decs to dataset." << endl;
        }

        addAttributeToDataSet(string("Actual declination derived from actual azimuth/elevation after pointing model"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeAntennaStatuses()
{
    if (m_voAntennaStatuses.size())
    {
        string strDatasetName("activity");

        //Create the data space
        hsize_t dimension[] = { m_voAntennaStatuses.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cAntennaStatus));
        //hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, H5T_VARIABLE);

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cAntennaStatus, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the antenna state (string)
        hid_t stringTypeValue = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeValue, sizeof(cAntennaStatus::m_chaAntennaStatus));
        H5Tinsert(compoundDataType, "value", HOFFSET(cAntennaStatus, m_chaAntennaStatus), stringTypeValue);
        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cAntennaStatus::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cAntennaStatus, m_chaStatus), stringTypeStatus);

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
}

void cSpectrometerHDF5OutputFile::writeMotorTorques()
{
    if (m_voMotorTorquesAzMaster_mNm.size())
    {
        string strDatasetName("motor-torque.az-master");

        //Create the data space
        hsize_t dimension[] = { m_voMotorTorquesAzMaster_mNm.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data types: Torque for azimuth-master motor
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voMotorTorquesAzMaster_mNm.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeMotorTorques(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeMotorTorques(): Wrote " << m_voMotorTorquesAzMaster_mNm.size() << " azimuth-master motor torques to dataset." << endl;
        }

        addAttributeToDataSet(string("Torques of azimuth-master motor"), strDatasetName, string("double"), string("Nm"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voMotorTorquesAzSlave_mNm.size())
    {
        string strDatasetName("motor-torque.az-slave");

        //Create the data space
        hsize_t dimension[] = { m_voMotorTorquesAzSlave_mNm.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data types: Torque for azimuth-slave motor
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voMotorTorquesAzSlave_mNm.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeMotorTorques(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeMotorTorques(): Wrote " << m_voMotorTorquesAzSlave_mNm.size() << " azimuth-slave motor torques to dataset." << endl;
        }

        addAttributeToDataSet(string("Torques of azimuth-slave motor"), strDatasetName, string("double"), string("Nm"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voMotorTorquesElMaster_mNm.size())
    {
        string strDatasetName("motor-torque.el-master");

        //Create the data space
        hsize_t dimension[] = { m_voMotorTorquesElMaster_mNm.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data types: Torque for elevation-master motor
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voMotorTorquesElMaster_mNm.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeMotorTorques(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeMotorTorques(): Wrote " << m_voMotorTorquesElMaster_mNm.size() << " elevation-master motor torques to dataset." << endl;
        }

        addAttributeToDataSet(string("Torques of elevation-master motor"), strDatasetName, string("double"), string("Nm"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voMotorTorquesElSlave_mNm.size())
    {
        string strDatasetName("motor-torque.el-slave");

        //Create the data space
        hsize_t dimension[] = { m_voMotorTorquesElSlave_mNm.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data types: Torque for elevation-slave motor
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voMotorTorquesElSlave_mNm.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeMotorTorques(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeMotorTorques(): Wrote " << m_voMotorTorquesElSlave_mNm.size() << " elevation-slave motor torques to dataset." << endl;
        }

        addAttributeToDataSet(string("Torques of elevation-slave motor"), strDatasetName, string("double"), string("Nm"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeAntennaConfiguration()
{
    //TODO: this might be a bit naive. I suspect that there's more complexity in this than the other datasets.
    if (m_vdPointingModelParams.size())
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
        oSS << m_oAntennaConfiguration.m_chaPointModelName;

        //Need to open the dataset again here as the H5LTmake_dataset used above does leave an open handle.
        hid_t dataset_id = H5Dopen2(m_iH5ConfigurationAntennasAntenna1GroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
        addAttributeToDataSet(oSS.str(), strDatasetName, string("string"), string(""), dataset_id);

        H5Dclose(dataset_id);
    }

    if (m_vdDelayModelParams.size())
    {
        hsize_t dimension = m_vdDelayModelParams.size();
        string strDatasetName("delay-model-params");

        herr_t err = H5LTmake_dataset(m_iH5ConfigurationAntennasAntenna1GroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_DOUBLE, &m_vdDelayModelParams.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeAntennaConfiguration(): HDF5 make dataset error." << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeAntennaConfiguration(): Wrote " << m_vdDelayModelParams.size() << " delay model parameters to file." << endl;
        }
    }

    //TODO: This stuff is coming out funky. Is it actually being assigned somewhere?
    {
        //Add the rest of the data as attributes to the group.
        hid_t variableLengthStringType;
        variableLengthStringType = H5Tcopy (H5T_C_S1);

        H5Tset_size (variableLengthStringType, H5T_VARIABLE);

        hid_t attrDataspace = H5Screate(H5S_SCALAR);
        hid_t groupHandle = m_iH5ConfigurationAntennasAntenna1GroupHandle; // Just to keep the lines a bit shorter.

        hid_t attrAntennaName = H5Acreate(groupHandle, "name", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
        const char* pcAntennaName = m_oAntennaConfiguration.m_chaAntennaName;
        H5Awrite(attrAntennaName, variableLengthStringType, &pcAntennaName);

        hid_t attrAntennaDiameter = H5Acreate(groupHandle, "diameter", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
        const char* pcAntennaDiameter_m = m_oAntennaConfiguration.m_chaAntennaDiameter_m;
        H5Awrite(attrAntennaDiameter, variableLengthStringType, &pcAntennaDiameter_m);

        hid_t attrAntennaBeamwidth = H5Acreate(groupHandle, "beamwidth", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
        const char* pcAntennaBeamwidth_deg = m_oAntennaConfiguration.m_chaAntennaBeamwidth_deg;
        H5Awrite(attrAntennaBeamwidth, variableLengthStringType, &pcAntennaBeamwidth_deg);

        hid_t attrAntennaLatitude = H5Acreate(groupHandle, "latitude", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
        const char* pcAntennaLatitude_deg = m_oAntennaConfiguration.m_chaAntennaLatitude_deg;
        H5Awrite(attrAntennaLatitude, variableLengthStringType, &pcAntennaLatitude_deg);

        hid_t attrAntennaLongitude = H5Acreate(groupHandle, "longitude", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
        const char* pcAntennaLongitude_deg = m_oAntennaConfiguration.m_chaAntennaLongitude_deg;
        H5Awrite(attrAntennaLongitude, variableLengthStringType, &pcAntennaLongitude_deg);

        H5Aclose(attrAntennaBeamwidth);
        H5Aclose(attrAntennaDiameter);
        H5Aclose(attrAntennaLatitude);
        H5Aclose(attrAntennaLongitude);
        H5Aclose(attrAntennaName);
        H5Tclose(variableLengthStringType);
        H5Sclose(attrDataspace);
    }
}

void cSpectrometerHDF5OutputFile::writeNoiseDiodeSoftwareStates()
{
    //TODO: Figure out what to do with this.
    if (m_voNoiseDiodeSoftwareStates.size())
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
        H5Tset_size(stringTypeValue, sizeof(cNoiseDiodeState::m_chaValue));
        H5Tinsert(compoundDataType, "value", HOFFSET(cNoiseDiodeState, m_chaValue), stringTypeValue);

        //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cNoiseDiodeState::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cNoiseDiodeState, m_chaStatus), stringTypeStatus);

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
}

void cSpectrometerHDF5OutputFile::writeNoiseDiodeSources()
{
    if (m_voNoiseDiodeSources.size())
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
        H5Tset_size(stringTypeValue, sizeof(cNoiseDiodeSource::m_chaSource));
        H5Tinsert(compoundDataType, "value", HOFFSET(cNoiseDiodeSource, m_chaSource), stringTypeValue);

        //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cNoiseDiodeSource::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cNoiseDiodeSource, m_chaStatus), stringTypeStatus);

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
}

void cSpectrometerHDF5OutputFile::writeNoideDiodeCurrents()
{
    if (m_voNoiseDiodeCurrents.size())
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
        H5Tset_size(stringTypeStatus, sizeof(cNoiseDiodeState::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

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
}

void cSpectrometerHDF5OutputFile::writeSelectedSources()
{
    if (m_voSelectedSources.size())
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
        H5Tset_size(stringTypeValue, sizeof(cSourceSelection::m_chaSource));
        H5Tinsert(compoundDataType, "value", HOFFSET(cSourceSelection, m_chaSource), stringTypeValue);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cSourceSelection::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cSourceSelection, m_chaStatus), stringTypeStatus);

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
}

void cSpectrometerHDF5OutputFile::writeRFFrequencies()
{
    if (m_voFrequenciesRFChan0_MHz.size())
    {
        string strDatasetName("rfe.rf.chan0.frequency");

        //Create the data space
        hsize_t dimension[] = { m_voFrequenciesRFChan0_MHz.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the RF frequency for channel 0 (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voFrequenciesRFChan0_MHz.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeRFFrequencies(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeRFFrequencies(): Wrote " << m_voFrequenciesRFChan0_MHz.size() << " RF frequencies for channel 0." << endl;
        }

        addAttributeToDataSet(string("RF centre frequency for final IF channel 0"), strDatasetName, string("double"), string("MHz"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voFrequenciesRFChan1_MHz.size())
    {
        string strDatasetName("rfe.rf.chan1.frequency");

        //Create the data space
        hsize_t dimension[] = { m_voFrequenciesRFChan1_MHz.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the RF frequency for channel 1 (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voFrequenciesRFChan1_MHz.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeRFFrequencies(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeRFFrequencies(): Wrote " << m_voFrequenciesRFChan1_MHz.size() << " RF frequencies for channel 1." << endl;
        }

        addAttributeToDataSet(string("RF centre frequency for final IF channel 1"), strDatasetName, string("double"), string("MHz"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeLOFrequencies()
{
    if (m_voFrequenciesLO0Chan0_MHz.size())
    {
        string strDatasetName("rfe.lo0.chan0.frequency");

        //Create the data space
        hsize_t dimension[] = { m_voFrequenciesLO0Chan0_MHz.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the chan 0 LO0 frequency (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voFrequenciesLO0Chan0_MHz.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): Wrote " << m_voFrequenciesLO0Chan0_MHz.size() << " LO0 chan 0 frequencies." << endl;
        }

        addAttributeToDataSet(string("frequency of channel 0, LO 0"), strDatasetName, string("double"), string("MHz"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voFrequenciesLO0Chan1_MHz.size())
    {
        string strDatasetName("rfe.lo0.chan1.frequency");

        //Create the data space
        hsize_t dimension[] = { m_voFrequenciesLO0Chan1_MHz.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the chan 1 LO0 frequency (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voFrequenciesLO0Chan1_MHz.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): Wrote " << m_voFrequenciesLO0Chan1_MHz.size() << " LO0 chan 1 frequencies." << endl;
        }

        addAttributeToDataSet(string("frequency of channel 0, LO 0"), strDatasetName, string("double"), string("MHz"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voFrequenciesLO1_MHz.size())
    {
        string strDatasetName("rfe.lo1.frequency");

        //Create the data space
        hsize_t dimension[] = { m_voFrequenciesLO1_MHz.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the LO1 frequency (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voFrequenciesLO1_MHz.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): Wrote " << m_voFrequenciesLO1_MHz.size() << " LO1 frequencies." << endl;
        }

        addAttributeToDataSet(string("frequency of LO 1"), strDatasetName, string("double"), string("MHz"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeIFBandwidths()
{
    if (m_voReceiverBandwidthsChan0_MHz.size())
    {
        string strDatasetName("rfe.if.chan0.receiver_bandwidth");

        //Create the data space
        hsize_t dimension[] = { m_voReceiverBandwidthsChan0_MHz.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the receiver bandwidth (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voReceiverBandwidthsChan0_MHz.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeIFBandwidths(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeIFBandwidths(): Wrote " << m_voReceiverBandwidthsChan0_MHz.size() << " chan0 receiver bandwidths." << endl;
        }

        addAttributeToDataSet(string("Analogue 3 dB bandwidth available to the ADC chan0"), strDatasetName, string("double"), string("MHz"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voReceiverBandwidthsChan1_MHz.size())
    {
        string strDatasetName("rfe.if.chan1.receiver_bandwidth");

        //Create the data space
        hsize_t dimension[] = { m_voReceiverBandwidthsChan1_MHz.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the receiver bandwidth (double)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voReceiverBandwidthsChan1_MHz.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeIFBandwidths(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeIFBandwidths(): Wrote " << m_voReceiverBandwidthsChan1_MHz.size() << " chan1 receiver bandwidths." << endl;
        }

        addAttributeToDataSet(string("Analogue 3 dB bandwidth available to the ADC chan1"), strDatasetName, string("double"), string("MHz"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeROACHAccumulationLengths()
{
    //Ideally this dataset should only be one element long. So we'll trim.
    //TODO: This might not be the case. Should we be able to adjust dynamically? Check with Ludwig.
    vector<cTimestampedUnsignedInt> voUniqueAccumulationLengths_nFrames;
    uint32_t ui32CurrentAccumulationLength_nFrames = 0;
    for (uint32_t ui = 0; ui < m_voROACHAccumulationLengths_nFrames.size(); ui++)
    {
        if (m_voROACHAccumulationLengths_nFrames[ui].m_u32Value != ui32CurrentAccumulationLength_nFrames)
        {
            ui32CurrentAccumulationLength_nFrames = m_voROACHAccumulationLengths_nFrames[ui].m_u32Value;
            voUniqueAccumulationLengths_nFrames.push_back(m_voROACHAccumulationLengths_nFrames[ui]);
        }
    }

    hsize_t dimension = 1;
    string strDatasetName("accum_length");

    herr_t err = H5LTmake_dataset(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_UINT32, &ui32CurrentAccumulationLength_nFrames);

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHAccumulationLengths(): HDF5 make dataset error." << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHAccumulationLengths(): Wrote accumulation length." << endl;
    }

    //Need to open the dataset again here for attribute as the H5LTmake_dataset used above does leave an open handle.
    hid_t dataset_id = H5Dopen2(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
    addAttributeToDataSet(string("number of PFB/FFT frames accumulated by the ROACH."), strDatasetName, string("double"), string("frames"), dataset_id);

    H5Dclose(dataset_id);

    if (0) // (voUniqueAccumulationLengths_nFrames.size() != 1)
    {
        cout << "Warning! Unexpected change in accumulation length during recording. Saving changes in additional sensor dataset." << endl;

        string strDatasetName("accum_length_history");

        //Create the data space
        hsize_t dimension[] = { voUniqueAccumulationLengths_nFrames.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedUnsignedInt));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedUnsignedInt, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the accumulation length (unsigned int)
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedUnsignedInt, m_u32Value), H5T_NATIVE_UINT32);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedUnsignedInt::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedUnsignedInt, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &voUniqueAccumulationLengths_nFrames.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeROACHAccumulationLengths(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeROACHAccumulationLengths(): Wrote " << voUniqueAccumulationLengths_nFrames.size() << " accumulation lengths." << endl;
        }

        addAttributeToDataSet(string("accumulation length"), strDatasetName, string("unsigned int"), string("no. of FFT frames"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}


void cSpectrometerHDF5OutputFile::writeROACHNBNarrowbandSelections()
{
    if (m_voROACHNBChannelSelects.size())
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
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedUnsignedInt::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedUnsignedInt, m_chaStatus), stringTypeStatus);

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
}

void cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency()
{
    hsize_t dimension = 1;
    string strDatasetName("adc_clk");

    herr_t err = H5LTmake_dataset(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_DOUBLE, &m_dROACHFrequencyFs_MHz);

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): HDF5 make dataset error." << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): Wrote ADC Sampling frequency." << endl;
    }

    //Need to open the dataset again here for attribute as the H5LTmake_dataset used above does leave an open handle.
    hid_t dataset_id = H5Dopen2(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
    addAttributeToDataSet(string("sampling frequency of ROACH / KatADC"), strDatasetName, string("double"), string("MHz"), dataset_id);

    H5Dclose(dataset_id);
}

void cSpectrometerHDF5OutputFile::writeROACHSizeOfFFTs()
{
    {
    hsize_t dimension = 1;
    string strDatasetName("dbe.fft.coarse.size");

    herr_t err = H5LTmake_dataset(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_UINT32, &m_dROACHSizeOfCoarseFFT_nSamp);

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): HDF5 make dataset error." << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): Wrote coarse FFT size: " << m_dROACHSizeOfCoarseFFT_nSamp << endl;
    }

    //Need to open the dataset again here for attribute as the H5LTmake_dataset used above does leave an open handle.
    hid_t dataset_id = H5Dopen2(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
    addAttributeToDataSet(string("size of coarse FFT"), strDatasetName, string("unsigned int"), string("no. of time domain input samples"), dataset_id);

    H5Dclose(dataset_id);
    }

    {
    hsize_t dimension = 1;
    string strDatasetName("dbe.fft.fine.size");

    herr_t err = H5LTmake_dataset(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_UINT32, &m_dROACHSizeOfFineFFT_nSamp);

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): HDF5 make dataset error." << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): Wrote fine FFT size: " << m_dROACHSizeOfFineFFT_nSamp << endl;
    }

    //Need to open the dataset again here for attribute as the H5LTmake_dataset used above does leave an open handle.
    hid_t dataset_id = H5Dopen2(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
    addAttributeToDataSet(string("size of fine FFT"), strDatasetName, string("unsigned int"), string("no. of time domain input samples"), dataset_id);

    H5Dclose(dataset_id);
    }
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

    //Add to compound data type: the shift mask applied to the coarse FFT no (unsigned int)
    H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedUnsignedInt, m_u32Value), H5T_NATIVE_UINT32);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedUnsignedInt::m_chaStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedUnsignedInt, m_chaStatus), stringTypeStatus);

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

    addAttributeToDataSet(string("the shift mask applied to the coarse FFT"), strDatasetName, string("unsigned int"), string(""), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

void cSpectrometerHDF5OutputFile::writeROACHAdcAttentuations()
{
    {
    string strDatasetName("dbe.adc.chan0.attenuation");

    //Create the data space
    hsize_t dimension[] = { m_voROACHADCAttenuationsChan0_dB.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the ADC attenuation (double)
    H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voROACHADCAttenuationsChan0_dB.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): Wrote " << m_voROACHADCAttenuationsChan0_dB.size() << " ADC attenuations for chan0." << endl;
    }

    addAttributeToDataSet(string("attenuation at input to ADC chan0"), strDatasetName, string("double"), string("dB"), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    }

    {
    string strDatasetName("dbe.adc.chan1.attenuation");

    //Create the data space
    hsize_t dimension[] = { m_voROACHADCAttenuationsChan1_dB.size() };
    hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

    //Create a compound data type consisting of different native types per entry:
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the ADC attenuation (double)
    H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the status of the sensor (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsDBEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voROACHADCAttenuationsChan1_dB.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): Wrote " << m_voROACHADCAttenuationsChan1_dB.size() << " ADC attenuations for chan1." << endl;
    }

    addAttributeToDataSet(string("attenuation at input to ADC chan1"), strDatasetName, string("double"), string("dB"), dataset);

    H5Tclose(stringTypeStatus);
    H5Tclose(compoundDataType);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    }
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


void cSpectrometerHDF5OutputFile::addAttributesToFile(const std::string &strVersion, const std::string &strExperimentID, int64_t i64AugmentTimestamp_us, uint32_t u32AugmentErrors, hid_t fileHandle)
{
    //Adds the attributes needed for the file as a whole.

    hid_t variableLengthStringType;
    variableLengthStringType = H5Tcopy (H5T_C_S1);

    H5Tset_size (variableLengthStringType, H5T_VARIABLE);

    hid_t attrDataspace = H5Screate(H5S_SCALAR);

    hid_t attrVersion = H5Acreate(fileHandle, "version", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcVersion = &strVersion[0];
    H5Awrite(attrVersion, variableLengthStringType, &pcVersion);

    hid_t attrExperimentID = H5Acreate(fileHandle, "experiment_id", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcExperimentID = &strExperimentID[0];
    H5Awrite(attrExperimentID, variableLengthStringType, &pcExperimentID);

    double dTimestamp_s = (double)i64AugmentTimestamp_us / 1e6;
    hid_t attrAugmentTimestamp = H5Acreate(fileHandle, "augment_ts", H5T_NATIVE_DOUBLE, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrAugmentTimestamp, H5T_NATIVE_DOUBLE, &dTimestamp_s);

    hid_t attrAugmentErrors = H5Acreate(fileHandle, "augment_errors", H5T_NATIVE_UINT32, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrAugmentErrors, H5T_NATIVE_UINT32, &u32AugmentErrors);

    H5Aclose (attrVersion);
    H5Aclose (attrExperimentID);
    H5Aclose (attrAugmentTimestamp);
    H5Aclose (attrAugmentErrors);
    H5Tclose(variableLengthStringType);
    H5Sclose(attrDataspace);
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

void cSpectrometerHDF5OutputFile::addAttributesToObservation(const std::string &strScriptName, const std::string &strScriptArguments, const std::string &strObserver, const std::string &strExperimentID, const std::string &strDescription, const std::string &strAntennas, const std::string &strStartTime, const std::string &strEndTime, const std::string &strNoiseDiodeParams, const std::string &strRFParams, const std::string &strStatus, hid_t observationGroup)
{
    //Adds attributes to the Configuration/Observation group

    hid_t variableLengthStringType;
    variableLengthStringType = H5Tcopy (H5T_C_S1);

    H5Tset_size (variableLengthStringType, H5T_VARIABLE);

    hid_t attrDataspace = H5Screate(H5S_SCALAR);

    hid_t attrScriptName = H5Acreate(observationGroup, "script_name", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcScriptName = &strScriptName[0];
    H5Awrite(attrScriptName, variableLengthStringType, &pcScriptName);

    hid_t attrScriptArguments = H5Acreate(observationGroup, "script_arguments", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcScriptArguments = &strScriptArguments[0];
    H5Awrite(attrScriptArguments, variableLengthStringType, &pcScriptArguments);

    hid_t attrObserver = H5Acreate(observationGroup, "observer", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcObserver = &strObserver[0];
    H5Awrite(attrObserver, variableLengthStringType, &pcObserver);

    hid_t attrExperimentID = H5Acreate(observationGroup, "experiment_id", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcExperimentID = &strExperimentID[0];
    H5Awrite(attrExperimentID, variableLengthStringType, &pcExperimentID);

    hid_t attrDescription = H5Acreate(observationGroup, "description", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcDescription = &strDescription[0];
    H5Awrite(attrDescription, variableLengthStringType, &pcDescription);

    hid_t attrAntennas = H5Acreate(observationGroup, "ants", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcAntennas = &strAntennas[0];
    H5Awrite(attrAntennas, variableLengthStringType, &pcAntennas);

    hid_t attrStartTime = H5Acreate(observationGroup, "starttime", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcStartTime = &strStartTime[0];
    H5Awrite(attrStartTime, variableLengthStringType, &pcStartTime);

    hid_t attrEndTime = H5Acreate(observationGroup, "endtime", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcEndTime = &strEndTime[0];
    H5Awrite(attrEndTime, variableLengthStringType, &pcEndTime);

    hid_t attrNoiseDiodeParams = H5Acreate(observationGroup, "nd_params", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcNoiseDiodeParams = &strNoiseDiodeParams[0];
    H5Awrite(attrNoiseDiodeParams, variableLengthStringType, &pcNoiseDiodeParams);

    hid_t attrRFParams = H5Acreate(observationGroup, "rf_params", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcRFParams = &strRFParams[0];
    H5Awrite(attrRFParams, variableLengthStringType, &pcRFParams);

    hid_t attrStatus = H5Acreate(observationGroup, "status", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcStatus = &strStatus[0];
    H5Awrite(attrStatus, variableLengthStringType, &pcStatus);

    H5Aclose(attrScriptName);
    H5Aclose(attrScriptArguments);
    H5Aclose(attrObserver);
    H5Aclose(attrExperimentID);
    H5Aclose(attrDescription);
    H5Aclose(attrAntennas);
    H5Aclose(attrStartTime);
    H5Aclose(attrEndTime);
    H5Aclose(attrNoiseDiodeParams);
    H5Aclose(attrRFParams);
    H5Aclose(attrStatus);
    H5Tclose(variableLengthStringType);
    H5Sclose(attrDataspace);
}

void cSpectrometerHDF5OutputFile::addAttributesToDBE(const std::string &strVisOrdering, const std::string &strStokesOrdering, hid_t DBEGroup)
{
    hid_t variableLengthStringType;
    variableLengthStringType = H5Tcopy (H5T_C_S1);

    H5Tset_size (variableLengthStringType, H5T_VARIABLE);

    hid_t attrDataspace = H5Screate(H5S_SCALAR);

    hid_t attrVisOrdering = H5Acreate(DBEGroup, "vis_ordering", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcVisOrdering = &strVisOrdering[0];
    H5Awrite(attrVisOrdering, variableLengthStringType, &pcVisOrdering);

    hid_t attrStokesOrdering = H5Acreate(DBEGroup, "stokes_ordering", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
    const char* pcStokesOrdering = &strStokesOrdering[0];
    H5Awrite(attrStokesOrdering, variableLengthStringType, &pcStokesOrdering);

    H5Aclose(attrVisOrdering);
    H5Aclose(attrStokesOrdering);
    H5Tclose(variableLengthStringType);
    H5Sclose(attrDataspace);
}

//Functions for adding logged data from outside this class
//Each function aquires a shared lock of a common shared mutex for adding data to its respective vector.
//When all of this data is written to file a unique lock is obtained from that mutex so that none of the data be altered at that point
void cSpectrometerHDF5OutputFile::addMarkupLabel(int64_t i64Timestamp_us, const string &strLabel)
{
    cMarkupLabels oNewMarkupLabel;
    oNewMarkupLabel.m_dTimestamp_s = double(i64Timestamp_us) / 1e6;
    sprintf(oNewMarkupLabel.m_chaLabel, strLabel.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voMarkupLabels.push_back(oNewMarkupLabel);
}

void cSpectrometerHDF5OutputFile::addRequestedAntennaAz(int64_t i64Timestamp_us, double dAzimuth_deg, const string &strStatus)
{
    cTimestampedDouble oNewRequestedAntennaAz;
    oNewRequestedAntennaAz.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewRequestedAntennaAz.m_dValue = dAzimuth_deg;
    sprintf(oNewRequestedAntennaAz.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voRequestedAntennaAzs_deg.push_back(oNewRequestedAntennaAz);
}

void cSpectrometerHDF5OutputFile::addRequestedAntennaEl(int64_t i64Timestamp_us, double dElevation_deg, const string &strStatus)
{
    cTimestampedDouble oNewRequestedAntennaEl;
    oNewRequestedAntennaEl.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewRequestedAntennaEl.m_dValue = dElevation_deg;
    sprintf(oNewRequestedAntennaEl.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voRequestedAntennaEls_deg.push_back(oNewRequestedAntennaEl);
}

void cSpectrometerHDF5OutputFile::addActualAntennaAz(int64_t i64Timestamp_us, double dAzimuth_deg, const string &strStatus)
{
    cTimestampedDouble oNewActualAntennaAz;
    oNewActualAntennaAz.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualAntennaAz.m_dValue = dAzimuth_deg;
    sprintf(oNewActualAntennaAz.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voActualAntennaAzs_deg.push_back(oNewActualAntennaAz);
}

void cSpectrometerHDF5OutputFile::addActualAntennaEl(int64_t i64Timestamp_us, double dElevation_deg, const string &strStatus)
{
    cTimestampedDouble oNewActualAntennaEl;
    oNewActualAntennaEl.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualAntennaEl.m_dValue = dElevation_deg;
    sprintf(oNewActualAntennaEl.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voActualAntennaEls_deg.push_back(oNewActualAntennaEl);
}

void cSpectrometerHDF5OutputFile::addActualSourceOffsetAz(int64_t i64Timestamp_us, double dAzimuthOffset_deg, const string &strStatus)
{
    cTimestampedDouble oNewActualSourceOffsetAz;
    oNewActualSourceOffsetAz.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualSourceOffsetAz.m_dValue = dAzimuthOffset_deg;
    sprintf(oNewActualSourceOffsetAz.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voActualSourceOffsetAzs_deg.push_back(oNewActualSourceOffsetAz);
}

void cSpectrometerHDF5OutputFile::addActualSourceOffsetEl(int64_t i64Timestamp_us, double dElevationOffset_deg, const string &strStatus)
{
    cTimestampedDouble oNewActualSourceOffsetEl;
    oNewActualSourceOffsetEl.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualSourceOffsetEl.m_dValue = dElevationOffset_deg;
    sprintf(oNewActualSourceOffsetEl.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voActualSourceOffsetEls_deg.push_back(oNewActualSourceOffsetEl);
}

void cSpectrometerHDF5OutputFile::addActualAntennaRA(int64_t i64Timestamp_us, double dRighAscension_deg, const string &strStatus)
{
    cTimestampedDouble oNewActualAntennaRA;
    oNewActualAntennaRA.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualAntennaRA.m_dValue = dRighAscension_deg;
    sprintf(oNewActualAntennaRA.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voActualAntennaRAs_deg.push_back(oNewActualAntennaRA);
}

void cSpectrometerHDF5OutputFile::addActualAntennaDec(int64_t i64Timestamp_us, double dDeclination_deg, const string &strStatus)
{
    cTimestampedDouble oNewActualAntennaDec;
    oNewActualAntennaDec.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualAntennaDec.m_dValue = dDeclination_deg;
    sprintf(oNewActualAntennaDec.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voActualAntennaDecs_deg.push_back(oNewActualAntennaDec);
}

void cSpectrometerHDF5OutputFile::addAntennaStatus(int64_t i64Timestamp_us, const string &strAntennaStatus, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cAntennaStatus oNewAntennaStatus;
    oNewAntennaStatus.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    sprintf( oNewAntennaStatus.m_chaAntennaStatus, strAntennaStatus.substr(0, sizeof(oNewAntennaStatus.m_chaAntennaStatus)).c_str() ); //Limit to size of the char array
    sprintf( oNewAntennaStatus.m_chaStatus, strStatus.c_str());

    //The numerical value is not used here

    m_voAntennaStatuses.push_back(oNewAntennaStatus);
}

void cSpectrometerHDF5OutputFile::motorTorqueAzMaster(int64_t i64Timestamp_us, double dAzMaster_mNm, const string &strStatus)
{
    cTimestampedDouble oNewMotorTorque;
    oNewMotorTorque.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewMotorTorque.m_dValue = dAzMaster_mNm;
    sprintf(oNewMotorTorque.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voMotorTorquesAzMaster_mNm.push_back(oNewMotorTorque);
}

void cSpectrometerHDF5OutputFile::motorTorqueAzSlave(int64_t i64Timestamp_us, double dAzSlave_mNm, const string &strStatus)
{
    cTimestampedDouble oNewMotorTorque;
    oNewMotorTorque.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewMotorTorque.m_dValue = dAzSlave_mNm;
    sprintf(oNewMotorTorque.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voMotorTorquesAzSlave_mNm.push_back(oNewMotorTorque);
}

void cSpectrometerHDF5OutputFile::motorTorqueElMaster(int64_t i64Timestamp_us, double dElMaster_mNm, const string &strStatus)
{
    cTimestampedDouble oNewMotorTorque;
    oNewMotorTorque.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewMotorTorque.m_dValue = dElMaster_mNm;
    sprintf(oNewMotorTorque.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voMotorTorquesElMaster_mNm.push_back(oNewMotorTorque);
}

void cSpectrometerHDF5OutputFile::motorTorqueElSlave(int64_t i64Timestamp_us, double dElSlave_mNm, const string &strStatus)
{
    cTimestampedDouble oNewMotorTorque;
    oNewMotorTorque.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewMotorTorque.m_dValue = dElSlave_mNm;
    sprintf(oNewMotorTorque.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voMotorTorquesElSlave_mNm.push_back(oNewMotorTorque);
}

void cSpectrometerHDF5OutputFile::setAppliedPointingModel(const string &strModelName, const vector<double> &vdPointingModelParams)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    //This only a single set of values, not a history of values.
    //Update the current record whenever a new value is received

    m_vdPointingModelParams = vdPointingModelParams;
    sprintf(m_oAntennaConfiguration.m_chaPointModelName, strModelName.c_str());
}

void cSpectrometerHDF5OutputFile::setAntennaName(const string &strAntennaName)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //Only a single set of values
    sprintf(m_oAntennaConfiguration.m_chaAntennaName, strAntennaName.c_str());
}

void cSpectrometerHDF5OutputFile::setAntennaDiameter(const string &strAntennaDiameter_m)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //Only a single set of values
    sprintf(m_oAntennaConfiguration.m_chaAntennaDiameter_m, strAntennaDiameter_m.c_str());
}

void cSpectrometerHDF5OutputFile::setAntennaBeamwidth(const string &strAntennaBeamwidth_deg)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //Only a single set of values
    sprintf(m_oAntennaConfiguration.m_chaAntennaBeamwidth_deg, strAntennaBeamwidth_deg.c_str());
}

void cSpectrometerHDF5OutputFile::setAntennaLongitude(const string &strAntennaLongitude_deg)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //Only a single set of values
    sprintf(m_oAntennaConfiguration.m_chaAntennaLongitude_deg, strAntennaLongitude_deg.c_str());
}

void cSpectrometerHDF5OutputFile::setAntennaLatitude(const string &strAntennaLatitude_deg)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //Only a single set of values
    sprintf(m_oAntennaConfiguration.m_chaAntennaLatitude_deg, strAntennaLatitude_deg.c_str());
}

/*void cSpectrometerHDF5OutputFile::addAntennaDelayModel(const string &strAntennaDelayModel)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //Only a single set of values
    sprintf(m_oAntennaConfiguration.m_chaAntennaDelayModel, strAntennaDelayModel.c_str());
}*/

void cSpectrometerHDF5OutputFile::addNoiseDiodeSoftwareState(int64_t i64Timestamp_us, int32_t i32NoiseDiodeState, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedInt oNewNoiseDiodeState;
    oNewNoiseDiodeState.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewNoiseDiodeState.m_i32Value = i32NoiseDiodeState;
    sprintf(oNewNoiseDiodeState.m_chaStatus, strStatus.c_str());

    m_voNoiseDiodeSoftwareStates.push_back(oNewNoiseDiodeState);
}

void cSpectrometerHDF5OutputFile::addNoiseDiodeSource(int64_t i64Timestamp_us, const string &strNoiseSource, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cNoiseDiodeSource oNewNoiseDiodeSource;
    oNewNoiseDiodeSource.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    sprintf( oNewNoiseDiodeSource.m_chaSource, strNoiseSource.substr(0, sizeof(oNewNoiseDiodeSource.m_chaSource)).c_str() ); //Limit to size of the char array
    sprintf( oNewNoiseDiodeSource.m_chaStatus, strStatus.c_str());

    m_voNoiseDiodeSources.push_back(oNewNoiseDiodeSource);
}

void cSpectrometerHDF5OutputFile::addNoiseDiodeCurrent(int64_t i64Timestamp_us, double dNoiseDiodeCurrent_A, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewNoiseDiodeCurrent;
    oNewNoiseDiodeCurrent.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewNoiseDiodeCurrent.m_dValue = dNoiseDiodeCurrent_A;
    sprintf(oNewNoiseDiodeCurrent.m_chaStatus, strStatus.c_str());


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
    sprintf(oNewSourceSelection.m_chaSource, oSS.str().substr(0, sizeof(oNewSourceSelection.m_chaSource)).c_str() ); //Limit to size of the char array
    sprintf(oNewSourceSelection.m_chaStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voSelectedSources.push_back(oNewSourceSelection);
}


void cSpectrometerHDF5OutputFile::addFrequencyRFChan0(int64_t i64Timestamp_us, double dFreqencyRFChan0_MHz, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewRFFrequency;
    oNewRFFrequency.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewRFFrequency.m_dValue = dFreqencyRFChan0_MHz;
    sprintf(oNewRFFrequency.m_chaStatus, strStatus.c_str());

    m_voFrequenciesRFChan0_MHz.push_back(oNewRFFrequency);
}

void cSpectrometerHDF5OutputFile::addFrequencyRFChan1(int64_t i64Timestamp_us, double dFreqencyRFChan1_MHz, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewRFFrequency;
    oNewRFFrequency.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewRFFrequency.m_dValue = dFreqencyRFChan1_MHz;
    sprintf(oNewRFFrequency.m_chaStatus, strStatus.c_str());

    m_voFrequenciesRFChan1_MHz.push_back(oNewRFFrequency);
}

void cSpectrometerHDF5OutputFile::addFrequencyLO0Chan0(int64_t i64Timestamp_us, double dFrequencyLO0Chan0_MHz, const string &strStatus)
{
    cTimestampedDouble oNewLOFrequency;
    oNewLOFrequency.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewLOFrequency.m_dValue = dFrequencyLO0Chan0_MHz;
    sprintf(oNewLOFrequency.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voFrequenciesLO0Chan0_MHz.push_back(oNewLOFrequency);
}

void cSpectrometerHDF5OutputFile::addFrequencyLO0Chan1(int64_t i64Timestamp_us, double dFrequencyLO0Chan1_MHz, const string &strStatus)
{
    cTimestampedDouble oNewLOFrequency;
    oNewLOFrequency.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewLOFrequency.m_dValue = dFrequencyLO0Chan1_MHz;
    sprintf(oNewLOFrequency.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voFrequenciesLO0Chan1_MHz.push_back(oNewLOFrequency);
}

void cSpectrometerHDF5OutputFile::addFrequencyLO1(int64_t i64Timestamp_us, double dFrequencyLO1_MHz, const string &strStatus)
{
    cTimestampedDouble oNewLOFrequency;
    oNewLOFrequency.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewLOFrequency.m_dValue = dFrequencyLO1_MHz;
    sprintf(oNewLOFrequency.m_chaStatus, strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voFrequenciesLO1_MHz.push_back(oNewLOFrequency);
}

void cSpectrometerHDF5OutputFile::addReceiverBandwidthChan0(int64_t i64Timestamp_us, double dReceiverBandwidthChan0_MHz, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewBandwidthIF;
    oNewBandwidthIF.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewBandwidthIF.m_dValue = dReceiverBandwidthChan0_MHz;
    sprintf(oNewBandwidthIF.m_chaStatus, strStatus.c_str());

    m_voReceiverBandwidthsChan0_MHz.push_back(oNewBandwidthIF);
}

void cSpectrometerHDF5OutputFile::addReceiverBandwidthChan1(int64_t i64Timestamp_us, double dReceiverBandwidthChan1_MHz, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewBandwidthIF;
    oNewBandwidthIF.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewBandwidthIF.m_dValue = dReceiverBandwidthChan1_MHz;
    sprintf(oNewBandwidthIF.m_chaStatus, strStatus.c_str());

    m_voReceiverBandwidthsChan1_MHz.push_back(oNewBandwidthIF);
}

void cSpectrometerHDF5OutputFile::addAccumulationLength(int64_t i64Timestamp_us, uint32_t u32NFrames)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedUnsignedInt oNewAccumulationLength;
    oNewAccumulationLength.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewAccumulationLength.m_u32Value = u32NFrames;
    sprintf(oNewAccumulationLength.m_chaStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voROACHAccumulationLengths_nFrames.push_back(oNewAccumulationLength);
}

void cSpectrometerHDF5OutputFile::addCoarseChannelSelect(int64_t i64Timestamp_us, uint32_t u32ChannelNo)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    if (!m_voROACHNBChannelSelects.size() || u32ChannelNo != m_voROACHNBChannelSelects[m_voROACHNBChannelSelects.size() - 1].m_u32Value)
    {
        cTimestampedUnsignedInt oNewCoarseChannelSelection;
        oNewCoarseChannelSelection.m_dTimestamp_s = i64Timestamp_us;
        oNewCoarseChannelSelection.m_u32Value = u32ChannelNo;
        sprintf(oNewCoarseChannelSelection.m_chaStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

        m_voROACHNBChannelSelects.push_back(oNewCoarseChannelSelection);
    }
}

void cSpectrometerHDF5OutputFile::setFrequencyFs(double dFrequencyFs_MHz)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    //This only a single value, not a history of values.
    //Update the current record whenever a new value is received

    m_dROACHFrequencyFs_MHz = dFrequencyFs_MHz;
}

void cSpectrometerHDF5OutputFile::setSizeOfCoarseFFT(uint32_t u32SizeOfCoarseFFT_nSamp)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_dROACHSizeOfCoarseFFT_nSamp = u32SizeOfCoarseFFT_nSamp;
    //cout << "cSpectrometerHDF5OutputFile::setSizeOfCoarseFFT(): wrote size of coarse FFT " << u32SizeOfCoarseFFT_nSamp << endl;
}

void cSpectrometerHDF5OutputFile::setSizeOfFineFFT(uint32_t u32SizeOfFineFFT_nSamp)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_dROACHSizeOfFineFFT_nSamp = u32SizeOfFineFFT_nSamp;
    //cout << "cSpectrometerHDF5OutputFile::setSizeOfFineFFT(): wrote size of fine FFT " << u32SizeOfFineFFT_nSamp << endl;
}

void cSpectrometerHDF5OutputFile::addCoarseFFTShiftMask(int64_t i64Timestamp_us, uint32_t u32ShiftMask)
{
    cTimestampedUnsignedInt oNewCoarseFFTShift;
    oNewCoarseFFTShift.m_dTimestamp_s = i64Timestamp_us;
    oNewCoarseFFTShift.m_u32Value = u32ShiftMask;
    sprintf(oNewCoarseFFTShift.m_chaStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voROACHCoarseFFTShiftMasks.push_back(oNewCoarseFFTShift);
}

void cSpectrometerHDF5OutputFile::addAttenuationADCChan0(int64_t i64Timestamp_us, double dADCAttenuationChan0_dB)
{
    cTimestampedDouble oNewAdcAttenuantion;
    oNewAdcAttenuantion.m_dTimestamp_s = i64Timestamp_us;
    oNewAdcAttenuantion.m_dValue = dADCAttenuationChan0_dB;
    sprintf(oNewAdcAttenuantion.m_chaStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voROACHADCAttenuationsChan0_dB.push_back(oNewAdcAttenuantion);
}

void cSpectrometerHDF5OutputFile::addAttenuationADCChan1(int64_t i64Timestamp_us, double dADCAttenuationChan1_dB)
{
    cTimestampedDouble oNewAdcAttenuantion;
    oNewAdcAttenuantion.m_dTimestamp_s = i64Timestamp_us;
    oNewAdcAttenuantion.m_dValue = dADCAttenuationChan1_dB;
    sprintf(oNewAdcAttenuantion.m_chaStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voROACHADCAttenuationsChan1_dB.push_back(oNewAdcAttenuantion);
}

