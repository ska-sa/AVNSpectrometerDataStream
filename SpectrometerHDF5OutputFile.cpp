
//System includes
#include <iostream>
#include <sstream>
#include <complex>

//Library includes
extern "C" {
#include <hdf5.h> //Note must be included before other HDF5 libraries
#include <hdf5_hl.h>
}

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

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
    m_iH5SensorsEnvGroupHandle                      = H5Gcreate2(m_iH5SensorsGroupHandle, "/MetaData/Sensors/Enviro", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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

    //Read antenna config information from ini file
    pt::ptree ObsInfo;
    pt::ini_parser::read_ini("ObservatoryInfo.ini", ObsInfo); // TODO: Think about not hardcoding this.
    setAntennaName(ObsInfo.get<string>("Antenna.name"));
    setAntennaLatitude(ObsInfo.get<string>("Antenna.latitude"));
    setAntennaLongitude(ObsInfo.get<string>("Antenna.longitude"));
    setAntennaAltitude(ObsInfo.get<string>("Antenna.altitude"));
    setAntennaDiameter(ObsInfo.get<string>("Antenna.diameter"));
    setAntennaBeamwidth(ObsInfo.get<string>("Antenna.beamwidth"));
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

    writeAcsRequestedAntennaAzEls();
    writeAcsActualAntennaAzEls();
    writeSkyRequestedAntennaAzEls();
    writeSkyActualAntennaAzEls();

    /* Marked for removal.
    writeActualSourceOffsetAzEls();
    writeActualAntennaRADecs();*/

    writeAntennaStatuses();
    /* Marked for removal.
    writeMotorTorques();
    */
    writeAntennaConfiguration();

    writeNoiseDiodeInformation();

    writeSelectedSources();

    writeRFFrequencies();
    writeLOFrequencies();
    writeIFBandwidths();
    writeReceiverAttenuations();

    writeEnvironmentData();

    writeROACHAccumulationLengths();
    writeROACHNBNarrowbandSelections();
    writeROACHSamplingFreqBandwidth();
    writeROACHNumberChannels();
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
    H5Gclose(m_iH5SensorsEnvGroupHandle);
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
        cTimestampedChar oState;

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
    hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedChar));

    //Add to compound data type: a timestamp (double)
    H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedChar, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

    //Add to compound data type: the noise diode state (string of 1 character "0" or "1")
    hid_t stringTypeValue = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeValue, sizeof(cTimestampedChar::m_chaValue));
    H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedChar, m_chaValue), stringTypeValue);

    //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
    hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
    H5Tset_size(stringTypeStatus, sizeof(cTimestampedChar::m_chaStatus));
    H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedChar, m_chaStatus), stringTypeStatus);

    //Create the data set of of the new compound datatype
    hid_t dataset = H5Dcreate1(m_iH5SensorsAntennasAntenna1GroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voROACHNoiseDiodeStateChanges.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHNoiseDiodeStates(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHNoiseDiodeStates(): Wrote " << m_voROACHNoiseDiodeStateChanges.size() << " noise diode states to dataset." << endl;
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

void cSpectrometerHDF5OutputFile::writeAcsRequestedAntennaAzEls()
{
    //Requested values to motor drive after pointing model
    //As per KAT7 Azimuth and elevation are in seperate datasets

    //Azimuth:
    ////////////////////////////////////////////////////////////////////////////////////////////
    if (m_voAcsRequestedAntennaAzs_deg.size())
    {
        string strDatasetName("pos.request-pointm-azim");

        //Create the data space
        hsize_t dimension[] = { m_voAcsRequestedAntennaAzs_deg.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voAcsRequestedAntennaAzs_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls(): Wrote " << m_voAcsRequestedAntennaAzs_deg.size() << " requested antenna azimuths to dataset." << endl;
        }

        addAttributeToDataSet(string("Requested antenna-space azimuth"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    //Elevation:
    ////////////////////////////////////////////////////////////////////////////////////////////
    if (m_voAcsRequestedAntennaEls_deg.size())
    {
        string strDatasetName("pos.request-pointm-elev");

        //Create the data space
        hsize_t dimension[] = { m_voAcsRequestedAntennaEls_deg.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voAcsRequestedAntennaEls_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls(): Wrote " << m_voAcsRequestedAntennaEls_deg.size() << " requested antenna elevations to dataset." << endl;
        }

        addAttributeToDataSet(string("Requested antenna-space elevation"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeAcsActualAntennaAzEls()
{
    //Values return from drive about current Az/El
    //Apply inverse pointing model to get to ideal Az/El

    //Requested values to motor drive after pointing model
    //As per KAT7 Azimuth and elevation are in seperate datasets

    //Azimuth:
    ////////////////////////////////////////////////////////////////////////////////////////////
    if (m_voAcsActualAntennaAzs_deg.size())
    {
        string strDatasetName("pos.actual-pointm-azim");

        //Create the data space
        hsize_t dimension[] = { m_voAcsActualAntennaAzs_deg.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voAcsActualAntennaAzs_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaAzEls(): Wrote " << m_voAcsActualAntennaAzs_deg.size() << " actual antenna azimuths to dataset." << endl;
        }

        addAttributeToDataSet(string("Actual antenna-space azimuth"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    //Elevation:
    ////////////////////////////////////////////////////////////////////////////////////////////
    if (m_voAcsActualAntennaEls_deg.size())
    {
        string strDatasetName("pos.actual-pointm-elev");

        //Create the data space
        hsize_t dimension[] = { m_voAcsActualAntennaEls_deg.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voAcsActualAntennaEls_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaAzEls(): Wrote " << m_voAcsActualAntennaEls_deg.size() << " actual antenna elevations to dataset." << endl;
        }

        addAttributeToDataSet(string("Actual antenna-space elevation"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeSkyRequestedAntennaAzEls()
{
    //Requested values to motor drive after pointing model
    //As per KAT7 Azimuth and elevation are in seperate datasets

    //Azimuth:
    ////////////////////////////////////////////////////////////////////////////////////////////
    if (m_voSkyRequestedAntennaAzs_deg.size())
    {
        string strDatasetName("pos.request-scan-azim");

        //Create the data space
        hsize_t dimension[] = { m_voSkyRequestedAntennaAzs_deg.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voSkyRequestedAntennaAzs_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls(): Wrote " << m_voSkyRequestedAntennaAzs_deg.size() << " requested sky azimuths to dataset." << endl;
        }

        addAttributeToDataSet(string("Requested sky-space azimuth"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    //Elevation:
    ////////////////////////////////////////////////////////////////////////////////////////////
    if (m_voSkyRequestedAntennaEls_deg.size())
    {
        string strDatasetName("pos.request-scan-elev");

        //Create the data space
        hsize_t dimension[] = { m_voSkyRequestedAntennaEls_deg.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voSkyRequestedAntennaEls_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeRequestedAntennaAzEls(): Wrote " << m_voSkyRequestedAntennaEls_deg.size() << " requested sky elevations to dataset." << endl;
        }

        addAttributeToDataSet(string("Requested sky-space elevation"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeSkyActualAntennaAzEls()
{
    //Values return from drive about current Az/El
    //Apply inverse pointing model to get to ideal Az/El

    //Requested values to motor drive after pointing model
    //As per KAT7 Azimuth and elevation are in seperate datasets

    //Azimuth:
    ////////////////////////////////////////////////////////////////////////////////////////////
    if (m_voSkyActualAntennaAzs_deg.size())
    {
        string strDatasetName("pos.actual-scan-azim");

        //Create the data space
        hsize_t dimension[] = { m_voSkyActualAntennaAzs_deg.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voSkyActualAntennaAzs_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaAzEls(): Wrote " << m_voSkyActualAntennaAzs_deg.size() << " actual sky azimuths to dataset." << endl;
        }

        addAttributeToDataSet(string("Actual sky-space azimuth"), strDatasetName, string("double"), string("deg"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    //Elevation:
    ////////////////////////////////////////////////////////////////////////////////////////////
    if (m_voSkyActualAntennaEls_deg.size())
    {
        string strDatasetName("pos.actual-scan-elev");

        //Create the data space
        hsize_t dimension[] = { m_voSkyActualAntennaEls_deg.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voSkyActualAntennaEls_deg.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaAzEls(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeActualAntennaAzEls(): Wrote " << m_voSkyActualAntennaEls_deg.size() << " actual sky elevations to dataset." << endl;
        }

        addAttributeToDataSet(string("Actual sky-space elevation"), strDatasetName, string("double"), string("deg"), dataset);

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


void cSpectrometerHDF5OutputFile::writeAntennaConfiguration()
{
    //TODO: this might be a bit naive. I suspect that there's more complexity in this than the other datasets.
    {
        hsize_t dimension = 30;
        string strDatasetName("pointing-model-params");

        herr_t err = H5LTmake_dataset(m_iH5ConfigurationAntennasAntenna1GroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_DOUBLE, m_adPointingModelParams);

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeAppliedPointingModel(): HDF5 make dataset error." << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeAppliedPointingModel(): Wrote 30 pointing model parameters to file." << endl;
        }

        /*//Add pointing model name as data attribute
        stringstream oSS;
        oSS << "Pointing model name = ";
        oSS << m_oAntennaConfiguration.m_chaPointModelName;

        //Need to open the dataset again here as the H5LTmake_dataset used above does leave an open handle.
        hid_t dataset_id = H5Dopen2(m_iH5ConfigurationAntennasAntenna1GroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
        addAttributeToDataSet(oSS.str(), strDatasetName, string("string"), string(""), dataset_id);

        H5Dclose(dataset_id);*/
    }

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
        const char* pcAntennaBeamwidth_deg = m_oAntennaConfiguration.m_chaAntennaBeamwidth;
        H5Awrite(attrAntennaBeamwidth, variableLengthStringType, &pcAntennaBeamwidth_deg);

        hid_t attrAntennaLatitude = H5Acreate(groupHandle, "latitude", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
        const char* pcAntennaLatitude_deg = m_oAntennaConfiguration.m_chaAntennaLatitude_deg;
        H5Awrite(attrAntennaLatitude, variableLengthStringType, &pcAntennaLatitude_deg);

        hid_t attrAntennaLongitude = H5Acreate(groupHandle, "longitude", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
        const char* pcAntennaLongitude_deg = m_oAntennaConfiguration.m_chaAntennaLongitude_deg;
        H5Awrite(attrAntennaLongitude, variableLengthStringType, &pcAntennaLongitude_deg);

        hid_t attrAntennaAltitude = H5Acreate(groupHandle, "altitude", variableLengthStringType, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
        const char* pcAntennaAltitude_m = m_oAntennaConfiguration.m_chaAntennaAltitude_m;
        H5Awrite(attrAntennaAltitude, variableLengthStringType, &pcAntennaAltitude_m);

        H5Aclose(attrAntennaAltitude);
        H5Aclose(attrAntennaBeamwidth);
        H5Aclose(attrAntennaDiameter);
        H5Aclose(attrAntennaLatitude);
        H5Aclose(attrAntennaLongitude);
        H5Aclose(attrAntennaName);
        H5Tclose(variableLengthStringType);
        H5Sclose(attrDataspace);
    }
}

void cSpectrometerHDF5OutputFile::writeNoiseDiodeInformation()
{
    if (m_voNoiseDiodeInputSource.size())
    {
        string strDatasetName("noise-diode.control-source");

        //Create the data space
        hsize_t dimension[] = { m_voNoiseDiodeInputSource.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cNoiseDiodeSource));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cNoiseDiodeSource, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the noise diode source (an 8-character string).
        hid_t stringTypeValue = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeValue, sizeof(cNoiseDiodeSource::m_chaSource));
        H5Tinsert(compoundDataType, "value", HOFFSET(cNoiseDiodeSource, m_chaSource), stringTypeValue);

        //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cNoiseDiodeSource::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cNoiseDiodeSource, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voNoiseDiodeInputSource.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeInformation(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeInformation(): Wrote " << m_voNoiseDiodeInputSource.size() << " noise diode input sources to dataset." << endl;
        }

        addAttributeToDataSet(string("noise diode input source"), strDatasetName, string("string"), string(""), dataset);

        H5Tclose(stringTypeValue);
        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voNoiseDiodeEnable.size())
    {
        string strDatasetName("noise-diode.enabled");

        //Create the data space
        hsize_t dimension[] = { m_voNoiseDiodeEnable.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedBool));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedBool, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        // Add to compound data type: noise diode enabled boolean state.
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedBool, m_i32Value), H5T_NATIVE_INT32);

        //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedBool::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedBool, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voNoiseDiodeEnable.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeInformation(): HDF5 make dataset error: " << strDatasetName  << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeInformation(): Wrote " << m_voNoiseDiodeEnable.size() << " noise diode enableds to dataset." << endl;
        }

        addAttributeToDataSet(string("noise diode enabled"), strDatasetName, string("boolean"), string(""), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voNoiseDiodeSelect.size())
    {
        string strDatasetName("noise-diode.select");

        //Create the data space
        hsize_t dimension[] = { m_voNoiseDiodeSelect.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedInt));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedInt, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the noise diode selection (integer).
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedInt, m_i32Value), H5T_NATIVE_INT32);

        //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedInt::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedInt, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voNoiseDiodeSelect.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeInformation(): HDF5 make dataset error: " << strDatasetName << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeInformation(): Wrote " << m_voNoiseDiodeSelect.size() << " noise diode selections to dataset." << endl;
        }

        addAttributeToDataSet(string("noise diode selection"), strDatasetName, string("int"), string(""), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voNoiseDiodePWMMark.size())
    {
        string strDatasetName("noise-diode.pwm-mark");

        //Create the data space
        hsize_t dimension[] = { m_voNoiseDiodePWMMark.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedInt));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedInt, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the noise diode pwm mark (integer).
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedInt, m_i32Value), H5T_NATIVE_INT32);

        //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedInt::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedInt, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voNoiseDiodePWMMark.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeInformation(): HDF5 make dataset error: " << strDatasetName << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeInformation(): Wrote " << m_voNoiseDiodePWMMark.size() << " noise diode PWM marks to dataset." << endl;
        }

        addAttributeToDataSet(string("noise diode pwm mark"), strDatasetName, string("int"), string("duty cycle, fraction out of 64"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voNoiseDiodePWMFrequency.size())
    {
        string strDatasetName("noise-diode.pwm-frequency");

        //Create the data space
        hsize_t dimension[] = { m_voNoiseDiodePWMFrequency.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedDouble));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedDouble, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the noise diode pwm frequency (double).
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedDouble, m_dValue), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the status of the noise diode equiptment (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedDouble::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedDouble, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voNoiseDiodePWMFrequency.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeInformation(): HDF5 make dataset error: " << strDatasetName << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeNoiseDiodeInformation(): Wrote " << m_voNoiseDiodePWMFrequency.size() << " noise diode frequencies to dataset." << endl;
        }

        addAttributeToDataSet(string("noise diode pwm frequency"), strDatasetName, string("double"), string("Hz"), dataset);

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
    if (m_voFrequencySelectLcp.size())
    {
        string strDatasetName("rfe.band.select.LCP");

        //Create the data space
        hsize_t dimension[] = { m_voFrequencySelectLcp.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedChar));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedChar, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the FrequencySelectLCP state (string of 1 character "0" or "1")
        // Zero means 5 GHz, 1 means 6.7 GHz.
        hid_t stringTypeValue = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeValue, sizeof(cTimestampedChar::m_chaValue));
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedChar, m_chaValue), stringTypeValue);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedChar::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedChar, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voFrequencySelectLcp.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeRFFrequencies(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeRFFrequencies(): Wrote " << m_voFrequencySelectLcp.size() << " RF band selections for LCP." << endl;
        }

        addAttributeToDataSet(string("RF band selection for LCP (0 - 5GHz, 1 - 6.7 GHz)"), strDatasetName, string("boolean"), string(""), dataset);

        H5Tclose(stringTypeValue);
        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voFrequencySelectRcp.size())
    {
        string strDatasetName("rfe.band.select.RCP");

        //Create the data space
        hsize_t dimension[] = { m_voFrequencySelectRcp.size() };
        hid_t dataspace = H5Screate_simple(1, dimension, NULL); // 1 = 1 dimensional

        //Create a compound data type consisting of different native types per entry:
        hid_t compoundDataType = H5Tcreate (H5T_COMPOUND, sizeof(cTimestampedChar));

        //Add to compound data type: a timestamp (double)
        H5Tinsert(compoundDataType, "timestamp", HOFFSET(cTimestampedChar, m_dTimestamp_s), H5T_NATIVE_DOUBLE);

        //Add to compound data type: the FrequencySelectLCP state (string of 1 character "0" or "1")
        // Zero means 5 GHz, 1 means 6.7 GHz.
        hid_t stringTypeValue = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeValue, sizeof(cTimestampedChar::m_chaValue));
        H5Tinsert(compoundDataType, "value", HOFFSET(cTimestampedChar, m_chaValue), stringTypeValue);

        //Add to compound data type: the status of the sensor (string typically containing "nominal")
        hid_t stringTypeStatus = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringTypeStatus, sizeof(cTimestampedChar::m_chaStatus));
        H5Tinsert(compoundDataType, "status", HOFFSET(cTimestampedChar, m_chaStatus), stringTypeStatus);

        //Create the data set of of the new compound datatype
        hid_t dataset = H5Dcreate1(m_iH5SensorsRFEGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voFrequencySelectRcp.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeRFFrequencies(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeRFFrequencies(): Wrote " << m_voFrequencySelectLcp.size() << " RF band selections for RCP." << endl;
        }

        addAttributeToDataSet(string("RF band selection for LCP (0 - 5GHz, 1 - 6.7 GHz)"), strDatasetName, string("boolean"), string(""), dataset);

        H5Tclose(stringTypeValue);
        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeLOFrequencies()
{
    if (m_voFrequenciesLOIntermediate5GHz_Hz.size())
    {
        string strDatasetName("rfe.lo-intermediate.5GHz.frequency"); // chan 0 is the 5 GHz receiver

        //Create the data space
        hsize_t dimension[] = { m_voFrequenciesLOIntermediate5GHz_Hz.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voFrequenciesLOIntermediate5GHz_Hz.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): Wrote " << m_voFrequenciesLOIntermediate5GHz_Hz.size() << " LO0 chan 0 frequencies." << endl;
        }

        addAttributeToDataSet(string("frequency of synth a, valon 0 (5 GHz LO)"), strDatasetName, string("double"), string("Hz"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voFrequenciesLOIntermediate6_7GHz_Hz.size())
    {
        string strDatasetName("rfe.lo-intermediate.6_7GHz.frequency"); // chan1 is the 6.7 GHz receiver

        //Create the data space
        hsize_t dimension[] = { m_voFrequenciesLOIntermediate6_7GHz_Hz.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voFrequenciesLOIntermediate6_7GHz_Hz.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): Wrote " << m_voFrequenciesLOIntermediate6_7GHz_Hz.size() << " LO0 chan 1 frequencies." << endl;
        }

        addAttributeToDataSet(string("frequency of synth b, valon 0 (6.7 GHz LO)"), strDatasetName, string("double"), string("Hz"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voFrequenciesLOFinal_Hz.size())
    {
        string strDatasetName("rfe.lo-final.frequency");

        //Create the data space
        hsize_t dimension[] = { m_voFrequenciesLOFinal_Hz.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voFrequenciesLOFinal_Hz.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): Wrote " << m_voFrequenciesLOFinal_Hz.size() << " LO1 frequencies." << endl;
        }

        addAttributeToDataSet(string("frequency of valon 1 synth a (Intermediate LO)"), strDatasetName, string("double"), string("Hz"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}


void cSpectrometerHDF5OutputFile::writeIFBandwidths()
{
    if (m_voReceiverBandwidthsLcp_Hz.size())
    {
        string strDatasetName("rfe.if.chan0.receiver_bandwidth");

        //Create the data space
        hsize_t dimension[] = { m_voReceiverBandwidthsLcp_Hz.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voReceiverBandwidthsLcp_Hz.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeIFBandwidths(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeIFBandwidths(): Wrote " << m_voReceiverBandwidthsLcp_Hz.size() << " chan0 receiver bandwidths." << endl;
        }

        addAttributeToDataSet(string("Analogue 3 dB bandwidth available to the ADC chan0"), strDatasetName, string("double"), string("Hz"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voReceiverBandwidthsRcp_Hz.size())
    {
        string strDatasetName("rfe.if.chan1.receiver_bandwidth");

        //Create the data space
        hsize_t dimension[] = { m_voReceiverBandwidthsRcp_Hz.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voReceiverBandwidthsRcp_Hz.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeIFBandwidths(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeIFBandwidths(): Wrote " << m_voReceiverBandwidthsRcp_Hz.size() << " chan1 receiver bandwidths." << endl;
        }

        addAttributeToDataSet(string("Analogue 3 dB bandwidth available to the ADC chan1"), strDatasetName, string("double"), string("Hz"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}


void cSpectrometerHDF5OutputFile::writeReceiverAttenuations()
{
    if (m_voReceiverAttenuationsLcp_dB.size())
    {
        string strDatasetName("rfe.LCP.attenuation");

        //Create the data space
        hsize_t dimension[] = { m_voReceiverAttenuationsLcp_dB.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voReceiverAttenuationsLcp_dB.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeReceiverAttenuations(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeReceiverAttenuations(): Wrote " << m_voReceiverAttenuationsLcp_dB.size() << " receiver LCP attenuation values." << endl;
        }

        addAttributeToDataSet(string("Receiver chain LCP variable attenuator setting."), strDatasetName, string("double"), string("dB"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }

    if (m_voReceiverAttenuationsRcp_dB.size())
    {
        string strDatasetName("rfe.RCP.attenuation");

        //Create the data space
        hsize_t dimension[] = { m_voReceiverAttenuationsRcp_dB.size() };
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

        herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voReceiverAttenuationsRcp_dB.front());

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeReceiverAttenuations(): HDF5 make dataset error" << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeReceiverAttenuations(): Wrote " << m_voReceiverAttenuationsRcp_dB.size() << " receiver RCP attenuation values." << endl;
        }

        addAttributeToDataSet(string("Receiver chain RCP variable attenuator setting."), strDatasetName, string("double"), string("dB"), dataset);

        H5Tclose(stringTypeStatus);
        H5Tclose(compoundDataType);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
}

void cSpectrometerHDF5OutputFile::writeEnvironmentData()
{
    if (m_voWindSpeeds_mps.size())
        {
            string strDatasetName("wind.speed");

            //Create the data space
            hsize_t dimension[] = { m_voWindSpeeds_mps.size() };
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
            hid_t dataset = H5Dcreate1(m_iH5SensorsEnvGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

            herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voWindSpeeds_mps.front());

            if(err < 0)
            {
                cout << "cSpectrometerHDF5OutputFile::writeEnvironmentData(): HDF5 make dataset error" << endl;
            }
            else
            {
                cout << "cSpectrometerHDF5OutputFile::writeEnvironmentData(): Wrote " << m_voWindSpeeds_mps.size() << " wind speed values." << endl;
            }

            addAttributeToDataSet(string("Wind speed."), strDatasetName, string("double"), string("m/s"), dataset);

            H5Tclose(stringTypeStatus);
            H5Tclose(compoundDataType);
            H5Sclose(dataspace);
            H5Dclose(dataset);
        }

    if (m_voWindDirections_degrees.size())
        {
            string strDatasetName("wind.direction");

            //Create the data space
            hsize_t dimension[] = { m_voWindDirections_degrees.size() };
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
            hid_t dataset = H5Dcreate1(m_iH5SensorsEnvGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

            herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voWindDirections_degrees.front());

            if(err < 0)
            {
                cout << "cSpectrometerHDF5OutputFile::writeEnvironmentData(): HDF5 make dataset error" << endl;
            }
            else
            {
                cout << "cSpectrometerHDF5OutputFile::writeEnvironmentData(): Wrote " << m_voWindDirections_degrees.size() << " wind direction values." << endl;
            }

            addAttributeToDataSet(string("Wind direction."), strDatasetName, string("double"), string("deg"), dataset);

            H5Tclose(stringTypeStatus);
            H5Tclose(compoundDataType);
            H5Sclose(dataspace);
            H5Dclose(dataset);
        }

    if (m_voTemperatures_degreesC.size())
        {
            string strDatasetName("air.temperature");

            //Create the data space
            hsize_t dimension[] = { m_voTemperatures_degreesC.size() };
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
            hid_t dataset = H5Dcreate1(m_iH5SensorsEnvGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

            herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voTemperatures_degreesC.front());

            if(err < 0)
            {
                cout << "cSpectrometerHDF5OutputFile::writeEnvironmentData(): HDF5 make dataset error" << endl;
            }
            else
            {
                cout << "cSpectrometerHDF5OutputFile::writeEnvironmentData(): Wrote " << m_voTemperatures_degreesC.size() << " air temperature values." << endl;
            }

            addAttributeToDataSet(string("Air temperature."), strDatasetName, string("double"), string("degC"), dataset);

            H5Tclose(stringTypeStatus);
            H5Tclose(compoundDataType);
            H5Sclose(dataspace);
            H5Dclose(dataset);
        }

    if (m_voAbsolutePressures_mbar.size())
        {
            string strDatasetName("air.pressure");

            //Create the data space
            hsize_t dimension[] = { m_voAbsolutePressures_mbar.size() };
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
            hid_t dataset = H5Dcreate1(m_iH5SensorsEnvGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

            herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voAbsolutePressures_mbar.front());

            if(err < 0)
            {
                cout << "cSpectrometerHDF5OutputFile::writeEnvironmentData(): HDF5 make dataset error" << endl;
            }
            else
            {
                cout << "cSpectrometerHDF5OutputFile::writeEnvironmentData(): Wrote " << m_voAbsolutePressures_mbar.size() << " air pressure values." << endl;
            }

            addAttributeToDataSet(string("Air pressure."), strDatasetName, string("double"), string("mbar"), dataset);

            H5Tclose(stringTypeStatus);
            H5Tclose(compoundDataType);
            H5Sclose(dataspace);
            H5Dclose(dataset);
        }

    if (m_voRelativeHumidities_percent.size())
        {
            string strDatasetName("relative.humidity");

            //Create the data space
            hsize_t dimension[] = { m_voRelativeHumidities_percent.size() };
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
            hid_t dataset = H5Dcreate1(m_iH5SensorsEnvGroupHandle, strDatasetName.c_str(), compoundDataType, dataspace, H5P_DEFAULT);

            herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voRelativeHumidities_percent.front());

            if(err < 0)
            {
                cout << "cSpectrometerHDF5OutputFile::writeEnvironmentData(): HDF5 make dataset error" << endl;
            }
            else
            {
                cout << "cSpectrometerHDF5OutputFile::writeEnvironmentData(): Wrote " << m_voRelativeHumidities_percent.size() << " relative humidity values." << endl;
            }

            addAttributeToDataSet(string("Relative humidity."), strDatasetName, string("double"), string("percent"), dataset);

            H5Tclose(stringTypeStatus);
            H5Tclose(compoundDataType);
            H5Sclose(dataspace);
            H5Dclose(dataset);
        }
}

void cSpectrometerHDF5OutputFile::writeROACHNumberChannels()
{
    hsize_t dimension = 1;
    string strDatasetName("n_chans");

    herr_t err = H5LTmake_dataset(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_UINT32, &m_u32NFrequencyBins);

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): HDF5 make dataset error." << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): Wrote number of channels." << endl;
    }

    //Need to open the dataset again here for attribute as the H5LTmake_dataset used above does leave an open handle.
    hid_t dataset_id = H5Dopen2(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
    addAttributeToDataSet(string("Number of frequency channels"), strDatasetName, string("uint32"), string(""), dataset_id);

    H5Dclose(dataset_id);
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

void cSpectrometerHDF5OutputFile::writeROACHSamplingFreqBandwidth()
{
    {
        hsize_t dimension = 1;
        string strDatasetName("adc_clk");

        herr_t err = H5LTmake_dataset(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_DOUBLE, &m_dROACHFrequencyFs_Hz);

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFreqBandwidth(): HDF5 make dataset error." << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFreqBandwidth(): Wrote ADC Sampling frequency." << endl;
        }

        //Need to open the dataset again here for attribute as the H5LTmake_dataset used above does leave an open handle.
        hid_t dataset_id = H5Dopen2(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
        addAttributeToDataSet(string("sampling frequency of ROACH / KatADC"), strDatasetName, string("double"), string("Hz"), dataset_id);

        H5Dclose(dataset_id);
    }

    {
        //Calculate the 'digital' bandwidth of the file. Wideband / radiometer mode assumed to be default.
        double iBandwidth_Hz = m_dROACHFrequencyFs_Hz / 2.0;
        if (m_u32ROACHSizeOfFineFFT_nSamp)
        {
            //If the fineFFT size is nonzero, then it's a narrowband mode.
            //Bandwidth is reduced by a factor of the number of coarse channels,
            //which is in turn the number of points divided by 2.
            int iNumCoarseChannels_nChan = m_u32ROACHSizeOfCoarseFFT_nSamp / 2;
            iBandwidth_Hz = iBandwidth_Hz / iNumCoarseChannels_nChan;
        }

        hsize_t dimension = 1;
        string strDatasetName("bandwidth");

        herr_t err = H5LTmake_dataset(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_DOUBLE, &iBandwidth_Hz);

        if(err < 0)
        {
            cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFreqBandwidth(): HDF5 make dataset error." << endl;
        }
        else
        {
            cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFreqBandwidth(): Wrote observation bandwidth." << endl;
        }

        //Need to open the dataset again here for attribute as the H5LTmake_dataset used above does leave an open handle.
        hid_t dataset_id = H5Dopen2(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
        addAttributeToDataSet(string("Bandwidth of observation"), strDatasetName, string("double"), string("Hz"), dataset_id);

        H5Dclose(dataset_id);
    }
}

void cSpectrometerHDF5OutputFile::writeROACHSizeOfFFTs()
{
    {
    hsize_t dimension = 1;
    string strDatasetName("dbe.fft.coarse.size");

    herr_t err = H5LTmake_dataset(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_UINT32, &m_u32ROACHSizeOfCoarseFFT_nSamp);

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): HDF5 make dataset error." << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): Wrote coarse FFT size: " << m_u32ROACHSizeOfCoarseFFT_nSamp << endl;
    }

    //Need to open the dataset again here for attribute as the H5LTmake_dataset used above does leave an open handle.
    hid_t dataset_id = H5Dopen2(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
    addAttributeToDataSet(string("size of coarse FFT"), strDatasetName, string("unsigned int"), string("no. of time domain input samples"), dataset_id);

    H5Dclose(dataset_id);
    }

    {
    hsize_t dimension = 1;
    string strDatasetName("dbe.fft.fine.size");

    herr_t err = H5LTmake_dataset(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), 1, &dimension, H5T_NATIVE_UINT32, &m_u32ROACHSizeOfFineFFT_nSamp);

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): HDF5 make dataset error." << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeROACHSamplingFrequency(): Wrote fine FFT size: " << m_u32ROACHSizeOfFineFFT_nSamp << endl;
    }

    //Need to open the dataset again here for attribute as the H5LTmake_dataset used above does leave an open handle.
    hid_t dataset_id = H5Dopen2(m_iH5ConfigurationDBEGroupHandle, strDatasetName.c_str(), H5P_DEFAULT);
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
    hsize_t dimension[] = { m_voROACHADCAttenuationsLcp_dB.size() };
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

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voROACHADCAttenuationsLcp_dB.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): Wrote " << m_voROACHADCAttenuationsLcp_dB.size() << " ADC attenuations for chan0." << endl;
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
    hsize_t dimension[] = { m_voROACHADCAttenuationsRcp_dB.size() };
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

    herr_t err = H5Dwrite(dataset, compoundDataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_voROACHADCAttenuationsRcp_dB.front());

    if(err < 0)
    {
        cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): HDF5 make dataset error" << endl;
    }
    else
    {
        cout << "cSpectrometerHDF5OutputFile::writeLOFrequencies(): Wrote " << m_voROACHADCAttenuationsRcp_dB.size() << " ADC attenuations for chan1." << endl;
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
    sprintf(oNewMarkupLabel.m_chaLabel, "%s", strLabel.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voMarkupLabels.push_back(oNewMarkupLabel);
}

void cSpectrometerHDF5OutputFile::addAcsRequestedAz(int64_t i64Timestamp_us, double dAzimuth_deg, const string &strStatus)
{
    cTimestampedDouble oNewRequestedAntennaAz;
    oNewRequestedAntennaAz.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewRequestedAntennaAz.m_dValue = dAzimuth_deg;
    sprintf(oNewRequestedAntennaAz.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voAcsRequestedAntennaAzs_deg.push_back(oNewRequestedAntennaAz);
}

void cSpectrometerHDF5OutputFile::addAcsRequestedEl(int64_t i64Timestamp_us, double dElevation_deg, const string &strStatus)
{
    cTimestampedDouble oNewRequestedAntennaEl;
    oNewRequestedAntennaEl.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewRequestedAntennaEl.m_dValue = dElevation_deg;
    sprintf(oNewRequestedAntennaEl.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voAcsRequestedAntennaEls_deg.push_back(oNewRequestedAntennaEl);
}

void cSpectrometerHDF5OutputFile::addAcsActualAz(int64_t i64Timestamp_us, double dAzimuth_deg, const string &strStatus)
{
    cTimestampedDouble oNewActualAntennaAz;
    oNewActualAntennaAz.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualAntennaAz.m_dValue = dAzimuth_deg;
    sprintf(oNewActualAntennaAz.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voAcsActualAntennaAzs_deg.push_back(oNewActualAntennaAz);
}

void cSpectrometerHDF5OutputFile::addAcsActualEl(int64_t i64Timestamp_us, double dElevation_deg, const string &strStatus)
{
    cTimestampedDouble oNewActualAntennaEl;
    oNewActualAntennaEl.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualAntennaEl.m_dValue = dElevation_deg;
    sprintf(oNewActualAntennaEl.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voAcsActualAntennaEls_deg.push_back(oNewActualAntennaEl);
}

void cSpectrometerHDF5OutputFile::addSkyRequestedAz(int64_t i64Timestamp_us, double dAzimuth_deg, const string &strStatus)
{
    cTimestampedDouble oNewRequestedAntennaAz;
    oNewRequestedAntennaAz.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewRequestedAntennaAz.m_dValue = dAzimuth_deg;
    sprintf(oNewRequestedAntennaAz.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voSkyRequestedAntennaAzs_deg.push_back(oNewRequestedAntennaAz);
}

void cSpectrometerHDF5OutputFile::addSkyRequestedEl(int64_t i64Timestamp_us, double dElevation_deg, const string &strStatus)
{
    cTimestampedDouble oNewRequestedAntennaEl;
    oNewRequestedAntennaEl.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewRequestedAntennaEl.m_dValue = dElevation_deg;
    sprintf(oNewRequestedAntennaEl.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voSkyRequestedAntennaEls_deg.push_back(oNewRequestedAntennaEl);
}

void cSpectrometerHDF5OutputFile::addSkyActualAz(int64_t i64Timestamp_us, double dAzimuth_deg, const string &strStatus)
{
    cTimestampedDouble oNewActualAntennaAz;
    oNewActualAntennaAz.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualAntennaAz.m_dValue = dAzimuth_deg;
    sprintf(oNewActualAntennaAz.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voSkyActualAntennaAzs_deg.push_back(oNewActualAntennaAz);
}

void cSpectrometerHDF5OutputFile::addSkyActualEl(int64_t i64Timestamp_us, double dElevation_deg, const string &strStatus)
{
    cTimestampedDouble oNewActualAntennaEl;
    oNewActualAntennaEl.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewActualAntennaEl.m_dValue = dElevation_deg;
    sprintf(oNewActualAntennaEl.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voSkyActualAntennaEls_deg.push_back(oNewActualAntennaEl);
}


void cSpectrometerHDF5OutputFile::addPointingModelParameter(uint8_t ui8ParameterNumber, double dParameterValue)
{
    boost::unique_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_adPointingModelParams[ui8ParameterNumber] = dParameterValue;
}


void cSpectrometerHDF5OutputFile::addAntennaStatus(int64_t i64Timestamp_us, const string &strAntennaStatus, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cAntennaStatus oNewAntennaStatus;
    oNewAntennaStatus.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    sprintf( oNewAntennaStatus.m_chaAntennaStatus, "%s", strAntennaStatus.substr(0, sizeof(oNewAntennaStatus.m_chaAntennaStatus)).c_str() ); //Limit to size of the char array
    sprintf( oNewAntennaStatus.m_chaStatus, "%s", strStatus.c_str());

    //The numerical value is not used here

    m_voAntennaStatuses.push_back(oNewAntennaStatus);
}

/* Marked for removal.
void cSpectrometerHDF5OutputFile::motorTorqueAzMaster(int64_t i64Timestamp_us, double dAzMaster_mNm, const string &strStatus)
{
    cTimestampedDouble oNewMotorTorque;
    oNewMotorTorque.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewMotorTorque.m_dValue = dAzMaster_mNm;
    sprintf(oNewMotorTorque.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voMotorTorquesAzMaster_mNm.push_back(oNewMotorTorque);
}

void cSpectrometerHDF5OutputFile::motorTorqueAzSlave(int64_t i64Timestamp_us, double dAzSlave_mNm, const string &strStatus)
{
    cTimestampedDouble oNewMotorTorque;
    oNewMotorTorque.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewMotorTorque.m_dValue = dAzSlave_mNm;
    sprintf(oNewMotorTorque.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voMotorTorquesAzSlave_mNm.push_back(oNewMotorTorque);
}

void cSpectrometerHDF5OutputFile::motorTorqueElMaster(int64_t i64Timestamp_us, double dElMaster_mNm, const string &strStatus)
{
    cTimestampedDouble oNewMotorTorque;
    oNewMotorTorque.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewMotorTorque.m_dValue = dElMaster_mNm;
    sprintf(oNewMotorTorque.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voMotorTorquesElMaster_mNm.push_back(oNewMotorTorque);
}

void cSpectrometerHDF5OutputFile::motorTorqueElSlave(int64_t i64Timestamp_us, double dElSlave_mNm, const string &strStatus)
{
    cTimestampedDouble oNewMotorTorque;
    oNewMotorTorque.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewMotorTorque.m_dValue = dElSlave_mNm;
    sprintf(oNewMotorTorque.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voMotorTorquesElSlave_mNm.push_back(oNewMotorTorque);
}
*/

void cSpectrometerHDF5OutputFile::setAntennaName(const string &strAntennaName)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //Only a single set of values
    sprintf(m_oAntennaConfiguration.m_chaAntennaName, "%s", strAntennaName.c_str());
}

void cSpectrometerHDF5OutputFile::setAntennaDiameter(const string &strAntennaDiameter_m)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //Only a single set of values
    sprintf(m_oAntennaConfiguration.m_chaAntennaDiameter_m, "%s", strAntennaDiameter_m.c_str());
}

void cSpectrometerHDF5OutputFile::setAntennaBeamwidth(const string &strAntennaBeamwidth_deg)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //Only a single set of values
    sprintf(m_oAntennaConfiguration.m_chaAntennaBeamwidth, "%s", strAntennaBeamwidth_deg.c_str());
}

void cSpectrometerHDF5OutputFile::setAntennaLongitude(const string &strAntennaLongitude_deg)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //Only a single set of values
    sprintf(m_oAntennaConfiguration.m_chaAntennaLongitude_deg, "%s", strAntennaLongitude_deg.c_str());
}

void cSpectrometerHDF5OutputFile::setAntennaLatitude(const string &strAntennaLatitude_deg)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //Only a single set of values
    sprintf(m_oAntennaConfiguration.m_chaAntennaLatitude_deg, "%s", strAntennaLatitude_deg.c_str());
}

void cSpectrometerHDF5OutputFile::setAntennaAltitude(const string &strAntennAltitude_m)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //Only a single set of values
    sprintf(m_oAntennaConfiguration.m_chaAntennaAltitude_m, "%s", strAntennAltitude_m.c_str());
}

void cSpectrometerHDF5OutputFile::setAntennaDelayModel(const vector<double> &vdDelayModelParams)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);
    //This only a single set of values, not a history of values.
    //Update the current record whenever a new value is received

    m_vdDelayModelParams = vdDelayModelParams;
}

void cSpectrometerHDF5OutputFile::addNoiseDiodeInputSource(int64_t i64Timestamp_us, const string &strNoiseDiodeInputSource, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cNoiseDiodeSource oNewNoiseDiodeInputSource;
    oNewNoiseDiodeInputSource.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    sprintf(oNewNoiseDiodeInputSource.m_chaSource, "%s", strNoiseDiodeInputSource.c_str());
    sprintf(oNewNoiseDiodeInputSource.m_chaStatus, "%s", strStatus.c_str());

    m_voNoiseDiodeInputSource.push_back(oNewNoiseDiodeInputSource);
}

void cSpectrometerHDF5OutputFile::addNoiseDiodeEnable(int64_t i64Timestamp_us, bool bNoiseDiodeEnable, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedBool oNewNoiseDiodeEnable;
    oNewNoiseDiodeEnable.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewNoiseDiodeEnable.m_i32Value = bNoiseDiodeEnable;
    sprintf(oNewNoiseDiodeEnable.m_chaStatus, "%s", strStatus.c_str());

    m_voNoiseDiodeEnable.push_back(oNewNoiseDiodeEnable);
}

void cSpectrometerHDF5OutputFile::addNoiseDiodeSelect(int64_t i64Timestamp_us, int32_t i32NoiseDiodeSelect, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedInt oNewNoiseDiodeSelect;
    oNewNoiseDiodeSelect.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewNoiseDiodeSelect.m_i32Value = i32NoiseDiodeSelect & 0x0000000f;
    sprintf(oNewNoiseDiodeSelect.m_chaStatus, "%s", strStatus.c_str());

    m_voNoiseDiodeSelect.push_back(oNewNoiseDiodeSelect);
}

void cSpectrometerHDF5OutputFile::addNoiseDiodePWMMark(int64_t i64Timestamp_us, int32_t i32NoiseDiodePWMMark, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedInt oNewNoiseDiodePWMMark;
    oNewNoiseDiodePWMMark.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewNoiseDiodePWMMark.m_i32Value = i32NoiseDiodePWMMark & 0x003f;
    sprintf(oNewNoiseDiodePWMMark.m_chaStatus, "%s", strStatus.c_str());

    m_voNoiseDiodePWMMark.push_back(oNewNoiseDiodePWMMark);
}

void cSpectrometerHDF5OutputFile::addNoiseDiodePWMFrequency(int64_t i64Timestamp_us, double dNoiseDiodePWMFrequency, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewNoiseDiodePWMFrequency;
    oNewNoiseDiodePWMFrequency.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewNoiseDiodePWMFrequency.m_dValue = dNoiseDiodePWMFrequency;
    sprintf(oNewNoiseDiodePWMFrequency.m_chaStatus, "%s", strStatus.c_str());

    m_voNoiseDiodePWMFrequency.push_back(oNewNoiseDiodePWMFrequency);
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
    sprintf(oNewSourceSelection.m_chaSource, "%s", oSS.str().substr(0, sizeof(oNewSourceSelection.m_chaSource)).c_str() ); //Limit to size of the char array
    sprintf(oNewSourceSelection.m_chaStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

    m_voSelectedSources.push_back(oNewSourceSelection);
}

// TODO: Think about whether to keep this as a char or as a bool. I have the noise diode enabled sensor as a bool.
void cSpectrometerHDF5OutputFile::addFrequencySelectLcp(int64_t i64Timestamp_us, bool bFrequencySelectLcp, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedChar oNewFrequencyLCP;
    oNewFrequencyLCP.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    //oNewFrequencyLCP.m_chaValue = bFreqencySelectLCP;
    sprintf(oNewFrequencyLCP.m_chaValue, "%s", (bFrequencySelectLcp)?"1":"0");
    sprintf(oNewFrequencyLCP.m_chaStatus, "%s", strStatus.c_str());

    m_voFrequencySelectLcp.push_back(oNewFrequencyLCP);
}

void cSpectrometerHDF5OutputFile::addFrequencySelectRcp(int64_t i64Timestamp_us, bool bFrequencySelectRcp, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedChar oNewRFFrequency;
    oNewRFFrequency.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    //oNewRFFrequency.m_chaValue = bFrequencySelectRCP;
    sprintf(oNewRFFrequency.m_chaValue, "%s", (bFrequencySelectRcp)?"1":"0");
    sprintf(oNewRFFrequency.m_chaStatus, "%s", strStatus.c_str());

    m_voFrequencySelectRcp.push_back(oNewRFFrequency);
}

void cSpectrometerHDF5OutputFile::addFrequencyLOIntermed5GHz(int64_t i64Timestamp_us, double dFrequencyLOIntermediate5GHz_Hz, const string &strStatus)
{
    cTimestampedDouble oNewLOFrequency;
    oNewLOFrequency.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewLOFrequency.m_dValue = dFrequencyLOIntermediate5GHz_Hz;
    sprintf(oNewLOFrequency.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voFrequenciesLOIntermediate5GHz_Hz.push_back(oNewLOFrequency);
}

void cSpectrometerHDF5OutputFile::addFrequencyLOIntermed6_7GHz(int64_t i64Timestamp_us, double dFrequencyLOIntermediate6_7GHz_MHz, const string &strStatus)
{
    cTimestampedDouble oNewLOFrequency;
    oNewLOFrequency.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewLOFrequency.m_dValue = dFrequencyLOIntermediate6_7GHz_MHz;
    sprintf(oNewLOFrequency.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voFrequenciesLOIntermediate6_7GHz_Hz.push_back(oNewLOFrequency);
}

void cSpectrometerHDF5OutputFile::addFrequencyLOFinal(int64_t i64Timestamp_us, double dFrequencyLOFinal_MHz, const string &strStatus)
{
    cTimestampedDouble oNewLOFrequency;
    oNewLOFrequency.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewLOFrequency.m_dValue = dFrequencyLOFinal_MHz;
    sprintf(oNewLOFrequency.m_chaStatus, "%s", strStatus.c_str());

    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_voFrequenciesLOFinal_Hz.push_back(oNewLOFrequency);
}

void cSpectrometerHDF5OutputFile::addReceiverBandwidthLcp(int64_t i64Timestamp_us, double dReceiverBandwidthLcp_MHz, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewBandwidthIF;
    oNewBandwidthIF.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewBandwidthIF.m_dValue = dReceiverBandwidthLcp_MHz;
    sprintf(oNewBandwidthIF.m_chaStatus, "%s", strStatus.c_str());

    m_voReceiverBandwidthsLcp_Hz.push_back(oNewBandwidthIF);
}

void cSpectrometerHDF5OutputFile::addReceiverBandwidthRcp(int64_t i64Timestamp_us, double dReceiverBandwidthRcp_MHz, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewBandwidthIF;
    oNewBandwidthIF.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewBandwidthIF.m_dValue = dReceiverBandwidthRcp_MHz;
    sprintf(oNewBandwidthIF.m_chaStatus, "%s", strStatus.c_str());

    m_voReceiverBandwidthsRcp_Hz.push_back(oNewBandwidthIF);
}

void cSpectrometerHDF5OutputFile::addReceiverLcpAttenuation(int64_t i64Timestamp_us, double dReceiverLcpAttenuation_dB, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewAttenuation;
    oNewAttenuation.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewAttenuation.m_dValue = dReceiverLcpAttenuation_dB;
    sprintf(oNewAttenuation.m_chaStatus, "%s", strStatus.c_str());

    m_voReceiverAttenuationsLcp_dB.push_back(oNewAttenuation);
}

void cSpectrometerHDF5OutputFile::addReceiverRcpAttenuation(int64_t i64Timestamp_us, double dReceiverRcpAttenuation_dB, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewAttenuation;
    oNewAttenuation.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewAttenuation.m_dValue = dReceiverRcpAttenuation_dB;
    sprintf(oNewAttenuation.m_chaStatus, "%s", strStatus.c_str());

    m_voReceiverAttenuationsRcp_dB.push_back(oNewAttenuation);
}

// Env values

void cSpectrometerHDF5OutputFile::addWindSpeed(int64_t i64Timestamp_us, double dWindSpeed_mps, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewWindSpeed;
    oNewWindSpeed.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewWindSpeed.m_dValue = dWindSpeed_mps;
    sprintf(oNewWindSpeed.m_chaStatus, "%s", strStatus.c_str());

    m_voWindSpeeds_mps.push_back(oNewWindSpeed);
}

void cSpectrometerHDF5OutputFile::addWindDirection(int64_t i64Timestamp_us, double dWindDirection_degrees, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewWindDirection;
    oNewWindDirection.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewWindDirection.m_dValue = dWindDirection_degrees;
    sprintf(oNewWindDirection.m_chaStatus, "%s", strStatus.c_str());

    m_voWindDirections_degrees.push_back(oNewWindDirection);
}

void cSpectrometerHDF5OutputFile::addTemperature(int64_t i64Timestamp_us, double dTemperature_degreesC, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewTemperature;
    oNewTemperature.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewTemperature.m_dValue = dTemperature_degreesC;
    sprintf(oNewTemperature.m_chaStatus, "%s", strStatus.c_str());

    m_voTemperatures_degreesC.push_back(oNewTemperature);
}

void cSpectrometerHDF5OutputFile::addAbsolutePressure(int64_t i64Timestamp_us, double dPressure_mbar, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewAbsolutePressure;
    oNewAbsolutePressure.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewAbsolutePressure.m_dValue = dPressure_mbar;
    sprintf(oNewAbsolutePressure.m_chaStatus, "%s", strStatus.c_str());

    m_voAbsolutePressures_mbar.push_back(oNewAbsolutePressure);
}

void cSpectrometerHDF5OutputFile::addRelativeHumidity(int64_t i64Timestamp_us, double dHumidity_percent, const string &strStatus)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    cTimestampedDouble oNewHumidity;
    oNewHumidity.m_dTimestamp_s = (double)i64Timestamp_us / 1e6;
    oNewHumidity.m_dValue = dHumidity_percent;
    sprintf(oNewHumidity.m_chaStatus, "%s", strStatus.c_str());

    m_voRelativeHumidities_percent.push_back(oNewHumidity);
}

// ROACH values

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
        oNewCoarseChannelSelection.m_dTimestamp_s = i64Timestamp_us / 1e6; //Convert from us to standard s
        oNewCoarseChannelSelection.m_u32Value = u32ChannelNo;
        sprintf(oNewCoarseChannelSelection.m_chaStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

        m_voROACHNBChannelSelects.push_back(oNewCoarseChannelSelection);
    }
}

void cSpectrometerHDF5OutputFile::setFrequencyFs(double dFrequencyFs_Hz)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    //This only a single value, not a history of values.
    //Update the current record whenever a new value is received

    m_dROACHFrequencyFs_Hz = dFrequencyFs_Hz;
}

void cSpectrometerHDF5OutputFile::setSizeOfCoarseFFT(uint32_t u32SizeOfCoarseFFT_nSamp)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_u32ROACHSizeOfCoarseFFT_nSamp = u32SizeOfCoarseFFT_nSamp;
    //cout << "cSpectrometerHDF5OutputFile::setSizeOfCoarseFFT(): wrote size of coarse FFT " << u32SizeOfCoarseFFT_nSamp << endl;
}

void cSpectrometerHDF5OutputFile::setSizeOfFineFFT(uint32_t u32SizeOfFineFFT_nSamp)
{
    boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

    m_u32ROACHSizeOfFineFFT_nSamp = u32SizeOfFineFFT_nSamp;
    //cout << "cSpectrometerHDF5OutputFile::setSizeOfFineFFT(): wrote size of fine FFT " << u32SizeOfFineFFT_nSamp << endl;
}

void cSpectrometerHDF5OutputFile::addCoarseFFTShiftMask(int64_t i64Timestamp_us, uint32_t u32ShiftMask)
{
    if (!m_voROACHCoarseFFTShiftMasks.size() || u32ShiftMask != m_voROACHCoarseFFTShiftMasks[m_voROACHCoarseFFTShiftMasks.size() - 1].m_u32Value)
    {
        cTimestampedUnsignedInt oNewCoarseFFTShift;
        oNewCoarseFFTShift.m_dTimestamp_s = i64Timestamp_us / 1e6;
        oNewCoarseFFTShift.m_u32Value = u32ShiftMask;
        sprintf(oNewCoarseFFTShift.m_chaStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

        boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

        m_voROACHCoarseFFTShiftMasks.push_back(oNewCoarseFFTShift);
    }
}

void cSpectrometerHDF5OutputFile::addAttenuationADCChan0(int64_t i64Timestamp_us, double dADCAttenuationChan0_dB)
{
    if (!m_voROACHADCAttenuationsLcp_dB.size() || dADCAttenuationChan0_dB != m_voROACHADCAttenuationsLcp_dB[m_voROACHADCAttenuationsLcp_dB.size() - 1].m_dValue)
    {
        cTimestampedDouble oNewAdcAttenuantion;
        oNewAdcAttenuantion.m_dTimestamp_s = i64Timestamp_us / 1e6;
        oNewAdcAttenuantion.m_dValue = dADCAttenuationChan0_dB;
        sprintf(oNewAdcAttenuantion.m_chaStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

        boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

        m_voROACHADCAttenuationsLcp_dB.push_back(oNewAdcAttenuantion);
    }
}

void cSpectrometerHDF5OutputFile::addAttenuationADCChan1(int64_t i64Timestamp_us, double dADCAttenuationChan1_dB)
{
    if (!m_voROACHADCAttenuationsRcp_dB.size() || dADCAttenuationChan1_dB != m_voROACHADCAttenuationsRcp_dB[m_voROACHADCAttenuationsRcp_dB.size() - 1].m_dValue)
    {
        cTimestampedDouble oNewAdcAttenuantion;
        oNewAdcAttenuantion.m_dTimestamp_s = i64Timestamp_us / 1e6;
        oNewAdcAttenuantion.m_dValue = dADCAttenuationChan1_dB;
        sprintf(oNewAdcAttenuantion.m_chaStatus, "nominal"); //This is hardcoded to always be nominal for now for compatability with KATDal. Adapt with available status data.

        boost::shared_lock<boost::shared_mutex> oLock(m_oAppendDataMutex);

        m_voROACHADCAttenuationsRcp_dB.push_back(oNewAdcAttenuantion);
    }
}
