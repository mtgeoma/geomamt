#include "atsheader.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <iomanip> // setw

void
read_header (std::fstream &ats, std::string parameter)
{
  using namespace std;
  using namespace boost;
  using namespace boost::gregorian;
  using namespace boost::posix_time;

  if ( parameter == "HeaderLength" )
  {
    ats.seekg(0x000,ios::beg);
    int16_t  HeaderLength;
    ats.read(reinterpret_cast<char *>(&HeaderLength), sizeof(HeaderLength));
    cout << " " << lexical_cast<string>(static_cast<int>(HeaderLength));
  }
  else if ( parameter == "HeaderVers" )
  {
    ats.seekg(0x002,ios::beg);
    int16_t  HeaderVers;
    ats.read(reinterpret_cast<char *>(&HeaderVers), sizeof(HeaderVers));
    cout << " " << lexical_cast<string>(static_cast<int>(HeaderVers));
  }
  else if ( parameter == "Samples" )
  {
    ats.seekg(0x004,ios::beg);
    uint32_t  Samples;
    ats.read(reinterpret_cast<char *>(&Samples), sizeof(Samples));
    cout << " " << lexical_cast<string>(static_cast<int>(Samples));
  }
  else if ( parameter == "SampleFreq" )
  {
    ats.seekg(0x008,ios::beg);
    float  SampleFreq;
    ats.read(reinterpret_cast<char *>(&SampleFreq), sizeof(SampleFreq));
    cout << " " << lexical_cast<string>(SampleFreq);
  }
  else if ( parameter == "StartDateTime" )
  {
    ats.seekg(0x00C,ios::beg);
    int32_t  StartDateTime;
    ats.read(reinterpret_cast<char *>(&StartDateTime), sizeof(StartDateTime));
    cout << " " << lexical_cast<string>(static_cast<int>(StartDateTime));
  }
  else if ( parameter == "StartDateTimeIso" )
  {
    ats.seekg(0x00C,ios::beg);
    int32_t  StartDateTime;
    ats.read(reinterpret_cast<char *>(&StartDateTime), sizeof(StartDateTime));
    ptime StartDateTimeIso(date(1970,1,1),
			   seconds(static_cast<int>(StartDateTime)));
    cout << " " << to_iso_extended_string(StartDateTimeIso);
  }
  else if ( parameter == "EndDateTime" ||  parameter == "EndDateTimeIso" )
  {
    ats.seekg(0x000,ios::beg);
    int16_t  HeaderLength;
    ats.read(reinterpret_cast<char *>(&HeaderLength), sizeof(HeaderLength));

    ats.seekg(0x008,ios::beg);
    float  SampleFreq;
    ats.read(reinterpret_cast<char *>(&SampleFreq), sizeof(SampleFreq));

    ats.seekg(0x00C,ios::beg);
    int32_t  StartDateTime;
    ats.read(reinterpret_cast<char *>(&StartDateTime), sizeof(StartDateTime));

    ats.seekg(0, std::ios::end);
    long ats_file_size = ats.tellg();

    long data_block_size = ats_file_size -
      static_cast<long>(HeaderLength);
    if( data_block_size % static_cast<long>(sizeof(float)) != 0 )
    {
      std::cerr << "data block size in ats file"
		<< " isn't a multiple of sizeof(float)";
      ats.close();
      exit(1);
    }

    data_block_size /= sizeof(float);
    long seconds_duration = data_block_size / static_cast<long>(SampleFreq);
    double microseconds_fraction =
      static_cast<double>(data_block_size) / static_cast<double>(SampleFreq);
    microseconds_fraction = (microseconds_fraction -
			     static_cast<double>(seconds_duration))*1e6;

    if ( parameter == "EndDateTime" )
    {
      seconds_duration += static_cast<int>(StartDateTime);
      cout << " " << lexical_cast<std::string>(seconds_duration) << "."
	   << lexical_cast<std::string>(microseconds_fraction);
    }
    else
    {
      ptime StartDateTimeIso(date(1970,1,1),
			     seconds(static_cast<int>(StartDateTime)));
      ptime EndDateTimeIso = StartDateTimeIso +
	seconds(seconds_duration)+microseconds(long(microseconds_fraction));
      cout << " " << to_iso_extended_string(EndDateTimeIso);
    }
  }
  else if ( parameter == "LSBMV" )
  {
    ats.seekg(0x010,ios::beg);
    double  LSBMV;
    ats.read(reinterpret_cast<char *>(&LSBMV), sizeof(LSBMV));
    cout << " " << lexical_cast<string>(LSBMV);
  }
  else if ( parameter == "GMTOffset" )
  {
    ats.seekg(0x018,ios::beg);
    int32_t  GMTOffset;
    ats.read(reinterpret_cast<char *>(&GMTOffset), sizeof(GMTOffset));
    cout << " " << lexical_cast<string>(static_cast<int>(GMTOffset));
  }
  else if ( parameter == "OrigSampleFreq" )
  {
    ats.seekg(0x01C,ios::beg);
    float  OrigSampleFreq;
    ats.read(reinterpret_cast<char *>(&OrigSampleFreq), sizeof(OrigSampleFreq));
    cout << " " << lexical_cast<string>(OrigSampleFreq);
  }
  else if ( parameter == "ADUSerNum" )
  {
    ats.seekg(0x020,ios::beg);
    int16_t  ADUSerNum;
    ats.read(reinterpret_cast<char *>(&ADUSerNum), sizeof(ADUSerNum));
    cout << " " << lexical_cast<string>(static_cast<int>(ADUSerNum));
  }
  else if ( parameter == "ADCSerNum" )
  {
    ats.seekg(0x022,ios::beg);
    int16_t  ADCSerNum;
    ats.read(reinterpret_cast<char *>(&ADCSerNum), sizeof(ADCSerNum));
    cout << " " << lexical_cast<string>(static_cast<int>(ADCSerNum));
  }
  else if ( parameter == "ChanNo" )
  {
    ats.seekg(0x024,ios::beg);
    char  ChanNo;
    ats.read(reinterpret_cast<char *>(&ChanNo), sizeof(ChanNo));
    cout << " " << lexical_cast<string>(static_cast<int>(ChanNo));
  }
  else if ( parameter == "Chopper" )
  {
    ats.seekg(0x025,ios::beg);
    char  Chopper;
    ats.read(reinterpret_cast<char *>(&Chopper), sizeof(Chopper));
    cout << " " << lexical_cast<string>(static_cast<int>(Chopper));
  }
  else if ( parameter == "ChanType" )
  {
    ats.seekg(0x026,ios::beg);
    char  ChanType[2];
    ats.read(reinterpret_cast<char *>(&ChanType), sizeof(ChanType));
    string StrChanType(ChanType,2);
    StrChanType=StrChanType.substr(0,StrChanType.find('\000'));
    cout << " " << StrChanType;
  }
  else if ( parameter == "SensorType" )
  {
    ats.seekg(0x028,ios::beg);
    char  SensorType[6];
    ats.read(reinterpret_cast<char *>(&SensorType), sizeof(SensorType));
    string StrSensorType(SensorType,6);
    StrSensorType=StrSensorType.substr(0,StrSensorType.find('\000'));
    cout << " " << StrSensorType;
  }
  else if ( parameter == "SensorSerNum" )
  {
    ats.seekg(0x02E,ios::beg);
    int16_t  SensorSerNum;
    ats.read(reinterpret_cast<char *>(&SensorSerNum), sizeof(SensorSerNum));
    cout << " " << lexical_cast<string>(static_cast<int>(SensorSerNum));
  }
  else if ( parameter == "PosX1" )
  {
    ats.seekg(0x030,ios::beg);
    float  PosX1;
    ats.read(reinterpret_cast<char *>(&PosX1), sizeof(PosX1));
    cout << " " << lexical_cast<string>(PosX1);
  }
  else if ( parameter == "PosY1" )
  {
    ats.seekg(0x034,ios::beg);
    float  PosY1;
    ats.read(reinterpret_cast<char *>(&PosY1), sizeof(PosY1));
    cout << " " << lexical_cast<string>(PosY1);
  }
  else if ( parameter == "PosZ1" )
  {
    ats.seekg(0x038,ios::beg);
    float  PosZ1;
    ats.read(reinterpret_cast<char *>(&PosZ1), sizeof(PosZ1));
    cout << " " << lexical_cast<string>(PosZ1);
  }
  else if ( parameter == "PosX2" )
  {
    ats.seekg(0x03C,ios::beg);
    float  PosX2;
    ats.read(reinterpret_cast<char *>(&PosX2), sizeof(PosX2));
    cout << " " << lexical_cast<string>(PosX2);
  }
  else if ( parameter == "PosY2" )
  {
    ats.seekg(0x040,ios::beg);
    float  PosY2;
    ats.read(reinterpret_cast<char *>(&PosY2), sizeof(PosY2));
    cout << " " << lexical_cast<string>(PosY2);
  }
  else if ( parameter == "PosZ2" )
  {
    ats.seekg(0x044,ios::beg);
    float  PosZ2;
    ats.read(reinterpret_cast<char *>(&PosZ2), sizeof(PosZ2));
    cout << " " << lexical_cast<string>(PosZ2);
  }
  else if ( parameter == "DipLength" )
  {
    ats.seekg(0x048,ios::beg);
    float  DipLength;
    ats.read(reinterpret_cast<char *>(&DipLength), sizeof(DipLength));
    cout << " " << lexical_cast<string>(DipLength);
  }
  else if ( parameter == "Angle" )
  {
    ats.seekg(0x04C,ios::beg);
    float  Angle;
    ats.read(reinterpret_cast<char *>(&Angle), sizeof(Angle));
    cout << " " << lexical_cast<string>(Angle);
  }
  else if ( parameter == "ProbeRes" )
  {
    ats.seekg(0x050,ios::beg);
    float  ProbeRes;
    ats.read(reinterpret_cast<char *>(&ProbeRes), sizeof(ProbeRes));
    cout << " " << lexical_cast<string>(ProbeRes);
  }
  else if ( parameter == "DCOffset" )
  {
    ats.seekg(0x054,ios::beg);
    float  DCOffset;
    ats.read(reinterpret_cast<char *>(&DCOffset), sizeof(DCOffset));
    cout << " " << lexical_cast<string>(DCOffset);
  }
  else if ( parameter == "PreGain" )
  {
    ats.seekg(0x058,ios::beg);
    float  PreGain;
    ats.read(reinterpret_cast<char *>(&PreGain), sizeof(PreGain));
    cout << " " << lexical_cast<string>(PreGain);
  }
  else if ( parameter == "PostGain" )
  {
    ats.seekg(0x05C,ios::beg);
    float  PostGain;
    ats.read(reinterpret_cast<char *>(&PostGain), sizeof(PostGain));
    cout << " " << lexical_cast<string>(PostGain);
  }
  else if ( parameter == "Latitude_msec" )
  {
    ats.seekg(0x060,ios::beg);
    int32_t  Latitude;
    ats.read(reinterpret_cast<char *>(&Latitude), sizeof(Latitude));
    cout << " " << lexical_cast<string>(static_cast<int>(Latitude));
  }
  else if ( parameter == "Latitude" )
  {
    ats.seekg(0x060,ios::beg);
    int32_t  Latitude;
    ats.read(reinterpret_cast<char *>(&Latitude), sizeof(Latitude));
    cout << " " << lexical_cast<string>(static_cast<int>(Latitude)/3600.e3);
  }
  else if ( parameter == "Longitude_msec" )
  {
    ats.seekg(0x064,ios::beg);
    int32_t  Longitude;
    ats.read(reinterpret_cast<char *>(&Longitude), sizeof(Longitude));
    cout << " " << lexical_cast<string>(static_cast<int>(Longitude));
  }
  else if ( parameter == "Longitude" )
  {
    ats.seekg(0x064,ios::beg);
    int32_t  Longitude;
    ats.read(reinterpret_cast<char *>(&Longitude), sizeof(Longitude));
    cout << " " << lexical_cast<string>(static_cast<int>(Longitude)/3600.e3);
  }
  else if ( parameter == "Elevation_cm" )
  {
    ats.seekg(0x068,ios::beg);
    int32_t  Elevation;
    ats.read(reinterpret_cast<char *>(&Elevation), sizeof(Elevation));
    cout << " " << lexical_cast<string>(static_cast<int>(Elevation));
  }
  else if ( parameter == "Elevation" )
  {
    ats.seekg(0x068,ios::beg);
    int32_t  Elevation;
    ats.read(reinterpret_cast<char *>(&Elevation), sizeof(Elevation));
    cout << " " << lexical_cast<string>(static_cast<int>(Elevation)/100.0);
  }
  else if ( parameter == "LatLongType" )
  {
    ats.seekg(0x06C,ios::beg);
    char  LatLongType;
    ats.read(reinterpret_cast<char *>(&LatLongType), sizeof(LatLongType));
    cout << " " << lexical_cast<string>(LatLongType);
  }
  else if ( parameter == "AddCoordType" )
  {
    ats.seekg(0x06D,ios::beg);
    char  AddCoordType;
    ats.read(reinterpret_cast<char *>(&AddCoordType), sizeof(AddCoordType));
    if( 0x00 == AddCoordType )
    {
      AddCoordType = ' ';
    }
    cout << " " << lexical_cast<string>(AddCoordType);
  }
  else if ( parameter == "RefMedian" )
  {
    ats.seekg(0x06E,ios::beg);
    int16_t  RefMedian;
    ats.read(reinterpret_cast<char *>(&RefMedian), sizeof(RefMedian));
    cout << " " << lexical_cast<string>(static_cast<int>(RefMedian));
  }
  else if ( parameter == "XCoord" )
  {
    ats.seekg(0x070,ios::beg);
    double  XCoord;
    ats.read(reinterpret_cast<char *>(&XCoord), sizeof(XCoord));
    cout << " " << lexical_cast<string>(XCoord);
  }
  else if ( parameter == "YCoord" )
  {
    ats.seekg(0x078,ios::beg);
    double  YCoord;
    ats.read(reinterpret_cast<char *>(&YCoord), sizeof(YCoord));
    cout << " " << lexical_cast<string>(YCoord);
  }
  else if ( parameter == "GPSStat" )
  {
    ats.seekg(0x080,ios::beg);
    char  GPSStat;
    ats.read(reinterpret_cast<char *>(&GPSStat), sizeof(GPSStat));
    cout << " " << lexical_cast<string>(GPSStat);
  }
  else if ( parameter == "GPSAccuracy" )
  {
    ats.seekg(0x081,ios::beg);
    char  GPSAccuracy;
    ats.read(reinterpret_cast<char *>(&GPSAccuracy), sizeof(GPSAccuracy));
    cout << " " << lexical_cast<string>(static_cast<int>(GPSAccuracy));
  }
  else if ( parameter == "UTCOffset" )
  {
    ats.seekg(0x082,ios::beg);
    int16_t  UTCOffset;
    ats.read(reinterpret_cast<char *>(&UTCOffset), sizeof(UTCOffset));
    cout << " " << lexical_cast<string>(static_cast<int>(UTCOffset));
  }
  else if ( parameter == "SystemType" )
  {
    ats.seekg(0x084,ios::beg);
    char  SystemType[12];
    ats.read(reinterpret_cast<char *>(&SystemType), sizeof(SystemType));
    string StrSystemType(SystemType,12);
    StrSystemType=StrSystemType.substr(0,StrSystemType.find('\000'));
    cout << " " << StrSystemType;
  }
  else if ( parameter == "SurveyHeaderName" )
  {
    ats.seekg(0x090,ios::beg);
    char  SurveyHeaderName[12];
    ats.read(reinterpret_cast<char *>(&SurveyHeaderName),
	     sizeof(SurveyHeaderName));
    string StrSurveyHeaderName(SurveyHeaderName,12);
    StrSurveyHeaderName =
      StrSurveyHeaderName.substr(0, StrSurveyHeaderName.find('\000'));
    cout << " " << StrSurveyHeaderName;
  }
  else if ( parameter == "MeasType" )
  {
    ats.seekg(0x09C,ios::beg);
    char  MeasType[4];
    ats.read(reinterpret_cast<char *>(&MeasType), sizeof(MeasType));
    string StrMeasType(MeasType,4);
    StrMeasType=StrMeasType.substr(0,StrMeasType.find('\000'));
    cout << " " << StrMeasType;
  }
  else if ( parameter == "LogFileName" )
  {
    ats.seekg(0x0A0,ios::beg);
    char  LogFileName[12];
    ats.read(reinterpret_cast<char *>(&LogFileName), sizeof(LogFileName));
    string StrLogFileName(LogFileName,12);
    StrLogFileName=StrLogFileName.substr(0,StrLogFileName.find('\000'));
    cout << " " << StrLogFileName;
  }
  else if ( parameter == "SelfTestResult" )
  {
    ats.seekg(0x0AC,ios::beg);
    char  SelfTestResult[2];
    ats.read(reinterpret_cast<char *>(&SelfTestResult), sizeof(SelfTestResult));
    string StrSelfTestResult(SelfTestResult,2);
    StrSelfTestResult =
      StrSelfTestResult.substr(0,StrSelfTestResult.find('\000'));
    cout << " " << StrSelfTestResult;
  }
  else if ( parameter == "Reserved5" )
  {
    ats.seekg(0x0AE,ios::beg);
    char  Reserved5[2];
    ats.read(reinterpret_cast<char *>(&Reserved5), sizeof(Reserved5));
    string StrReserved5(Reserved5,2);
    StrReserved5=StrReserved5.substr(0,StrReserved5.find('\000'));
    cout << " " << StrReserved5;
  }
  else if ( parameter == "CalFreqs" )
  {
    ats.seekg(0x0B0,ios::beg);
    int16_t  CalFreqs;
    ats.read(reinterpret_cast<char *>(&CalFreqs), sizeof(CalFreqs));
    cout << " " << lexical_cast<string>(static_cast<int>(CalFreqs));
  }
  else if ( parameter == "CalEntryLength" )
  {
    ats.seekg(0x0B2,ios::beg);
    int16_t  CalEntryLength;
    ats.read(reinterpret_cast<char *>(&CalEntryLength), sizeof(CalEntryLength));
    cout << " " << lexical_cast<string>(static_cast<int>(CalEntryLength));
  }
  else if ( parameter == "CalVersion" )
  {
    ats.seekg(0x0B4,ios::beg);
    int16_t  CalVersion;
    ats.read(reinterpret_cast<char *>(&CalVersion), sizeof(CalVersion));
    cout << " " << lexical_cast<string>(static_cast<int>(CalVersion));
  }
  else if ( parameter == "CalStartAddress" )
  {
    ats.seekg(0x0B6,ios::beg);
    int16_t  CalStartAddress;
    ats.read(reinterpret_cast<char *>(&CalStartAddress),
	     sizeof(CalStartAddress));
    cout << " " << lexical_cast<string>(static_cast<int>(CalStartAddress));
  }
  else if ( parameter == "LFFilters" )
  {
    ats.seekg(0x0B8,ios::beg);
    char  LFFilters[8];
    ats.read(reinterpret_cast<char *>(&LFFilters), sizeof(LFFilters));
    string StrLFFilters(LFFilters,8);
    StrLFFilters=StrLFFilters.substr(0,StrLFFilters.find('\000'));
    cout << " " << StrLFFilters;
  }
  else if ( parameter == "ADU06CalFilename" )
  {
    ats.seekg(0x0C0,ios::beg);
    char  ADU06CalFilename[12];
    ats.read(reinterpret_cast<char *>(&ADU06CalFilename),
	     sizeof(ADU06CalFilename));
    string StrADU06CalFilename(ADU06CalFilename,12);
    StrADU06CalFilename =
      StrADU06CalFilename.substr(0,StrADU06CalFilename.find('\000'));
    cout << " " << StrADU06CalFilename;
  }
  else if ( parameter == "ADUCalTimeDate" )
  {
    ats.seekg(0x0CC,ios::beg);
    int32_t  ADUCalTimeDate;
    ats.read(reinterpret_cast<char *>(&ADUCalTimeDate), sizeof(ADUCalTimeDate));
    cout << " " << lexical_cast<string>(static_cast<int>(ADUCalTimeDate));
  }
  else if ( parameter == "SensorCalFilename" )
  {
    ats.seekg(0x0D0,ios::beg);
    char  SensorCalFilename[12];
    ats.read(reinterpret_cast<char *>(&SensorCalFilename),
	     sizeof(SensorCalFilename));
    string StrSensorCalFilename(SensorCalFilename,12);
    StrSensorCalFilename =
      StrSensorCalFilename.substr(0,StrSensorCalFilename.find('\000'));
    cout << " " << StrSensorCalFilename;
  }
  else if ( parameter == "SensorCalTimeDate" )
  {
    ats.seekg(0x0DC,ios::beg);
    int32_t  SensorCalTimeDate;
    ats.read(reinterpret_cast<char *>(&SensorCalTimeDate),
	     sizeof(SensorCalTimeDate));
    cout << " " << lexical_cast<string>(static_cast<int>(SensorCalTimeDate));
  }
  else if ( parameter == "PowerlineFreq1" )
  {
    ats.seekg(0x0E0,ios::beg);
    float  PowerlineFreq1;
    ats.read(reinterpret_cast<char *>(&PowerlineFreq1), sizeof(PowerlineFreq1));
    cout << " " << lexical_cast<string>(PowerlineFreq1);
  }
  else if ( parameter == "PowerlineFreq2" )
  {
    ats.seekg(0x0E4,ios::beg);
    float  PowerlineFreq2;
    ats.read(reinterpret_cast<char *>(&PowerlineFreq2), sizeof(PowerlineFreq2));
    cout << " " << lexical_cast<string>(PowerlineFreq2);
  }
  else if ( parameter == "HFFilters" )
  {
    ats.seekg(0x0E8,ios::beg);
    char  HFFilters[8];
    ats.read(reinterpret_cast<char *>(&HFFilters), sizeof(HFFilters));
    string StrHFFilters(HFFilters,8);
    StrHFFilters=StrHFFilters.substr(0,StrHFFilters.find('\000'));
    cout << " " << StrHFFilters;
  }
  else if ( parameter == "CSAMTFreq" )
  {
    ats.seekg(0x0F0,ios::beg);
    float  CSAMTFreq;
    ats.read(reinterpret_cast<char *>(&CSAMTFreq), sizeof(CSAMTFreq));
    cout << " " << lexical_cast<string>(CSAMTFreq);
  }
  else if ( parameter == "CSAMTBlocks" )
  {
    ats.seekg(0x0F4,ios::beg);
    int16_t  CSAMTBlocks;
    ats.read(reinterpret_cast<char *>(&CSAMTBlocks), sizeof(CSAMTBlocks));
    cout << " " << lexical_cast<string>(static_cast<int>(CSAMTBlocks));
  }
  else if ( parameter == "CSAMTStacksPBlock" )
  {
    ats.seekg(0x0F6,ios::beg);
    int16_t  CSAMTStacksPBlock;
    ats.read(reinterpret_cast<char *>(&CSAMTStacksPBlock),
	     sizeof(CSAMTStacksPBlock));
    cout << " " << lexical_cast<string>(static_cast<int>(CSAMTStacksPBlock));
  }
  else if ( parameter == "CSAMTBlockLength" )
  {
    ats.seekg(0x0F8,ios::beg);
    int32_t  CSAMTBlockLength;
    ats.read(reinterpret_cast<char *>(&CSAMTBlockLength),
	     sizeof(CSAMTBlockLength));
    cout << " " << lexical_cast<string>(static_cast<int>(CSAMTBlockLength));
  }
  else if ( parameter == "ADBBoardType" )
  {
    ats.seekg(0x0FC,ios::beg);
    char  ADBBoardType[4];
    ats.read(reinterpret_cast<char *>(&ADBBoardType), sizeof(ADBBoardType));
    string StrADBBoardType(ADBBoardType,4);
    StrADBBoardType=StrADBBoardType.substr(0,StrADBBoardType.find('\000'));
    cout << " " << StrADBBoardType;
  }
  else if ( parameter == "Client" )
  {
    ats.seekg(0x100,ios::beg);
    char  Client[16];
    ats.read(reinterpret_cast<char *>(&Client), sizeof(Client));
    string StrClient(Client,16);
    StrClient=StrClient.substr(0,StrClient.find('\000'));
    cout << " " << StrClient;
  }
  else if ( parameter == "Contractor" )
  {
    ats.seekg(0x110,ios::beg);
    char  Contractor[16];
    ats.read(reinterpret_cast<char *>(&Contractor), sizeof(Contractor));
    string StrContractor(Contractor,16);
    StrContractor=StrContractor.substr(0,StrContractor.find('\000'));
    cout << " " << StrContractor;
  }
  else if ( parameter == "Area" )
  {
    ats.seekg(0x120,ios::beg);
    char  Area[16];
    ats.read(reinterpret_cast<char *>(&Area), sizeof(Area));
    string StrArea(Area,16);
    StrArea=StrArea.substr(0,StrArea.find('\000'));
    cout << " " << StrArea;
  }
  else if ( parameter == "SurveyID" )
  {
    ats.seekg(0x130,ios::beg);
    char  SurveyID[16];
    ats.read(reinterpret_cast<char *>(&SurveyID), sizeof(SurveyID));
    string StrSurveyID(SurveyID,16);
    StrSurveyID=StrSurveyID.substr(0,StrSurveyID.find('\000'));
    cout << " " << StrSurveyID;
  }
  else if ( parameter == "Operator" )
  {
    ats.seekg(0x140,ios::beg);
    char  Operator[16];
    ats.read(reinterpret_cast<char *>(&Operator), sizeof(Operator));
    string StrOperator(Operator,16);
    StrOperator=StrOperator.substr(0,StrOperator.find('\000'));
    cout << " " << StrOperator;
  }
  else if ( parameter == "Reserved" )
  {
    ats.seekg(0x150,ios::beg);
    char  Reserved[112];
    ats.read(reinterpret_cast<char *>(&Reserved), sizeof(Reserved));
    string StrReserved(Reserved,112);
    StrReserved=StrReserved.substr(0,StrReserved.find('\000'));
    cout << " " << StrReserved;
  }
  else if ( parameter == "XmlHeader" )
  {
    ats.seekg(0x1C0,ios::beg);
    char  XmlHeader[64];
    ats.read(reinterpret_cast<char *>(&XmlHeader), sizeof(XmlHeader));
    string StrXmlHeader(XmlHeader,64);
    StrXmlHeader=StrXmlHeader.substr(0,StrXmlHeader.find('\000'));
    cout << " " << StrXmlHeader;
  }
  else if ( parameter == "Comments" )
  {
    ats.seekg(0x200,ios::beg);
    char  Comments[512];
    ats.read(reinterpret_cast<char *>(&Comments), sizeof(Comments));
    string StrComments(Comments,512);
    StrComments=StrComments.substr(0,StrComments.find('\000'));
    cout << " " << StrComments;
  }
  else
  {
    cerr << "unknow read parameter " << parameter << "\n";
  }
}

void
read_header (std::fstream &ats)
{
  using namespace std;
  using namespace boost;
  using namespace boost::gregorian;
  using namespace boost::posix_time;

  ats.seekg(0x000,ios::beg);
  int16_t  HeaderLength;
  ats.read(reinterpret_cast<char *>(&HeaderLength), sizeof(HeaderLength));
  cout << left << setw(19) << "HeaderLength: "
       << lexical_cast<string>(static_cast<int>(HeaderLength)) << "\n";

  int16_t  HeaderVers;
  ats.read(reinterpret_cast<char *>(&HeaderVers), sizeof(HeaderVers));
  cout << left << setw(19) << "HeaderVers: "
       << lexical_cast<string>(static_cast<int>(HeaderVers)) << "\n";

  uint32_t  Samples;
  ats.read(reinterpret_cast<char *>(&Samples), sizeof(Samples));
  cout << left << setw(19) << "Samples: "
       << lexical_cast<string>(static_cast<int>(Samples)) << "\n";

  float  SampleFreq;
  ats.read(reinterpret_cast<char *>(&SampleFreq), sizeof(SampleFreq));
  cout << left << setw(19) << "SampleFreq: "
       << lexical_cast<string>(SampleFreq) << "\n";

  int32_t  StartDateTime;
  ats.read(reinterpret_cast<char *>(&StartDateTime), sizeof(StartDateTime));
  cout << left << setw(19) << "StartDateTime: "
       << lexical_cast<string>(static_cast<int>(StartDateTime)) << "\n";

  ptime StartDateTimeIso(date(1970,1,1),
			 seconds(static_cast<int>(StartDateTime)));
  cout << left << setw(19) << "StartDateTimeIso: "
       << to_iso_extended_string(StartDateTimeIso) << " [only read mode]\n";

  ats.seekg(0, std::ios::end);
  long ats_file_size = static_cast<long>(ats.tellg());
  long data_block_size = ats_file_size - static_cast<long>(HeaderLength);

  if( data_block_size % static_cast<long>(sizeof(float)) != 0 )
  {
    cerr << "data block size in ats file isn't a multiple of sizeof(float)\n";
    ats.close();
    exit(1);
  }

  data_block_size /= sizeof(float);
  long seconds_duration = data_block_size / long(SampleFreq);
  double microseconds_fraction =
    double(data_block_size) / double(SampleFreq);
  microseconds_fraction = (microseconds_fraction -
			   double(seconds_duration))*1e6;

  ptime EndDateTimeIso = StartDateTimeIso +
    seconds(seconds_duration) + microseconds(long(microseconds_fraction));

  seconds_duration += static_cast<int>(StartDateTime);
  cout << left << setw(19) << "EndDateTime: "
       << lexical_cast<std::string>(seconds_duration) << "."
       << lexical_cast<std::string>(int(microseconds_fraction))
       << " [only read mode]\n";

  cout << left << setw(19) << "EndDateTimeIso: "
       << to_iso_extended_string(EndDateTimeIso) << " [only read mode]\n";

  ats.seekg(0x010,ios::beg);
  double  LSBMV;
  ats.read(reinterpret_cast<char *>(&LSBMV), sizeof(LSBMV));
  cout << left << setw(19) << "LSBMV: " << lexical_cast<string>(LSBMV) << "\n";

  int32_t  GMTOffset;
  ats.read(reinterpret_cast<char *>(&GMTOffset), sizeof(GMTOffset));
  cout << left << setw(19) << "GMTOffset: "
       << lexical_cast<string>(static_cast<int>(GMTOffset)) << "\n";

  float  OrigSampleFreq;
  ats.read(reinterpret_cast<char *>(&OrigSampleFreq), sizeof(OrigSampleFreq));
  cout << left << setw(19) << "OrigSampleFreq: "
       << lexical_cast<string>(OrigSampleFreq) << "\n";

  int16_t  ADUSerNum;
  ats.read(reinterpret_cast<char *>(&ADUSerNum), sizeof(ADUSerNum));
  cout << left << setw(19) << "ADUSerNum: "
       << lexical_cast<string>(static_cast<int>(ADUSerNum)) << "\n";

  int16_t  ADCSerNum;
  ats.read(reinterpret_cast<char *>(&ADCSerNum), sizeof(ADCSerNum));
  cout << left << setw(19) << "ADCSerNum: "
       << lexical_cast<string>(static_cast<int>(ADCSerNum)) << "\n";

  char  ChanNo;
  ats.read(reinterpret_cast<char *>(&ChanNo), sizeof(ChanNo));
  cout << left << setw(19) << "ChanNo: "
       << lexical_cast<string>(static_cast<int>(ChanNo)) << "\n";

  char  Chopper;
  ats.read(reinterpret_cast<char *>(&Chopper), sizeof(Chopper));
  cout << left << setw(19) << "Chopper: "
       << lexical_cast<string>(static_cast<int>(Chopper)) << "\n";

  char  ChanType[2];
  ats.read(reinterpret_cast<char *>(&ChanType), sizeof(ChanType));
  string StrChanType(ChanType,2);
  StrChanType=StrChanType.substr(0,StrChanType.find('\000'));
  cout << left << setw(19) << "ChanType: " << StrChanType << "\n";

  char  SensorType[6];
  ats.read(reinterpret_cast<char *>(&SensorType), sizeof(SensorType));
  string StrSensorType(SensorType,6);
  StrSensorType=StrSensorType.substr(0,StrSensorType.find('\000'));
  cout << left << setw(19) << "SensorType: " << StrSensorType << "\n";

  int16_t  SensorSerNum;
  ats.read(reinterpret_cast<char *>(&SensorSerNum), sizeof(SensorSerNum));
  cout << left << setw(19) << "SensorSerNum: "
       << lexical_cast<string>(static_cast<int>(SensorSerNum)) << "\n";

  float  PosX1;
  ats.read(reinterpret_cast<char *>(&PosX1), sizeof(PosX1));
  cout << left << setw(19) << "PosX1: " << lexical_cast<string>(PosX1) << "\n";

  float  PosY1;
  ats.read(reinterpret_cast<char *>(&PosY1), sizeof(PosY1));
  cout << left << setw(19) << "PosY1: " << lexical_cast<string>(PosY1) << "\n";

  float  PosZ1;
  ats.read(reinterpret_cast<char *>(&PosZ1), sizeof(PosZ1));
  cout << left << setw(19) << "PosZ1: " << lexical_cast<string>(PosZ1) << "\n";

  float  PosX2;
  ats.read(reinterpret_cast<char *>(&PosX2), sizeof(PosX2));
  cout << left << setw(19) << "PosX2: " << lexical_cast<string>(PosX2) << "\n";

  float  PosY2;
  ats.read(reinterpret_cast<char *>(&PosY2), sizeof(PosY2));
  cout << left << setw(19) << "PosY2: " << lexical_cast<string>(PosY2) << "\n";

  float  PosZ2;
  ats.read(reinterpret_cast<char *>(&PosZ2), sizeof(PosZ2));
  cout << left << setw(19) << "PosZ2: " << lexical_cast<string>(PosZ2) << "\n";

  float  DipLength;
  ats.read(reinterpret_cast<char *>(&DipLength), sizeof(DipLength));
  cout << left << setw(19) << "DipLength: "
       << lexical_cast<string>(DipLength) << "\n";

  float  Angle;
  ats.read(reinterpret_cast<char *>(&Angle), sizeof(Angle));
  cout << left << setw(19) << "Angle: " << lexical_cast<string>(Angle) << "\n";

  float  ProbeRes;
  ats.read(reinterpret_cast<char *>(&ProbeRes), sizeof(ProbeRes));
  cout << left << setw(19) << "ProbeRes: "
       << lexical_cast<string>(ProbeRes) << "\n";

  float  DCOffset;
  ats.read(reinterpret_cast<char *>(&DCOffset), sizeof(DCOffset));
  cout << left << setw(19) << "DCOffset: "
       << lexical_cast<string>(DCOffset) << "\n";

  float  PreGain;
  ats.read(reinterpret_cast<char *>(&PreGain), sizeof(PreGain));
  cout << left << setw(19) << "PreGain: "
       << lexical_cast<string>(PreGain) << "\n";

  float  PostGain;
  ats.read(reinterpret_cast<char *>(&PostGain), sizeof(PostGain));
  cout << left << setw(19) << "PostGain: "
       << lexical_cast<string>(PostGain) << "\n";

  int32_t  Latitude_msec;
  ats.read(reinterpret_cast<char *>(&Latitude_msec), sizeof(Latitude_msec));
  cout << left << setw(19) << "Latitude_msec: "
       << lexical_cast<string>(static_cast<int>(Latitude_msec)) << "\n";

  cout << left << setw(19) << "Latitude: "
       << lexical_cast<string>(static_cast<int>(Latitude_msec)/3600.e3)
       << " [only read mode]\n";

  int32_t  Longitude_msec;
  ats.read(reinterpret_cast<char *>(&Longitude_msec), sizeof(Longitude_msec));
  cout << left << setw(19) << "Longitude_msec: "
       << lexical_cast<string>(static_cast<int>(Longitude_msec)) << "\n";

  cout << left << setw(19) << "Longitude: "
       << lexical_cast<string>(static_cast<int>(Longitude_msec)/3600.e3)
       << " [only read mode]\n";

  int32_t  Elevation_cm;
  ats.read(reinterpret_cast<char *>(&Elevation_cm), sizeof(Elevation_cm));
  cout << left << setw(19) << "Elevation_cm: "
       << lexical_cast<string>(static_cast<int>(Elevation_cm)) << "\n";

  cout << left << setw(19) << "Elevation: "
       << lexical_cast<string>(static_cast<int>(Elevation_cm)/100.0)
       << " [only read mode]\n";

  char  LatLongType;
  ats.read(reinterpret_cast<char *>(&LatLongType), sizeof(LatLongType));
  cout << left << setw(19) << "LatLongType: "
       << lexical_cast<string>(LatLongType) << "\n";

  char  AddCoordType;
  ats.read(reinterpret_cast<char *>(&AddCoordType), sizeof(AddCoordType));

  if( 0x00 == AddCoordType )
  {
    AddCoordType = ' ';
  }

  cout << left << setw(19) << "AddCoordType: "
       << lexical_cast<string>(AddCoordType) << "\n";

  int16_t  RefMedian;
  ats.read(reinterpret_cast<char *>(&RefMedian), sizeof(RefMedian));
  cout << left << setw(19) << "RefMedian: "
       << lexical_cast<string>(static_cast<int>(RefMedian)) << "\n";

  double  XCoord;
  ats.read(reinterpret_cast<char *>(&XCoord), sizeof(XCoord));
  cout << left << setw(19) << "XCoord: "
       << lexical_cast<string>(XCoord) << "\n";

  double  YCoord;
  ats.read(reinterpret_cast<char *>(&YCoord), sizeof(YCoord));
  cout << left << setw(19) << "YCoord: "
       << lexical_cast<string>(YCoord) << "\n";

  char  GPSStat;
  ats.read(reinterpret_cast<char *>(&GPSStat), sizeof(GPSStat));
  cout << left << setw(19) << "GPSStat: "
       << lexical_cast<string>(GPSStat) << "\n";

  char  GPSAccuracy;
  ats.read(reinterpret_cast<char *>(&GPSAccuracy), sizeof(GPSAccuracy));
  cout << left << setw(19) << "GPSAccuracy: "
       << lexical_cast<string>(static_cast<int>(GPSAccuracy)) << "\n";

  int16_t  UTCOffset;
  ats.read(reinterpret_cast<char *>(&UTCOffset), sizeof(UTCOffset));
  cout << left << setw(19) << "UTCOffset: "
       << lexical_cast<string>(static_cast<int>(UTCOffset)) << "\n";

  char  SystemType[12];
  ats.read(reinterpret_cast<char *>(&SystemType), sizeof(SystemType));
  string StrSystemType(SystemType,12);
  StrSystemType=StrSystemType.substr(0,StrSystemType.find('\000'));
  cout << left << setw(19) << "SystemType: " << StrSystemType << "\n";

  char  SurveyHeaderName[12];
  ats.read(reinterpret_cast<char *>(&SurveyHeaderName),
	   sizeof(SurveyHeaderName));
  string StrSurveyHeaderName(SurveyHeaderName,12);
  StrSurveyHeaderName =
    StrSurveyHeaderName.substr(0, StrSurveyHeaderName.find('\000'));
  cout << left << setw(19) << "SurveyHeaderName: "
       << StrSurveyHeaderName << "\n";

  char  MeasType[4];
  ats.read(reinterpret_cast<char *>(&MeasType), sizeof(MeasType));
  string StrMeasType(MeasType,4);
  StrMeasType=StrMeasType.substr(0,StrMeasType.find('\000'));
  cout << left << setw(19) << "MeasType: " << StrMeasType << "\n";

  char  LogFileName[12];
  ats.read(reinterpret_cast<char *>(&LogFileName), sizeof(LogFileName));
  string StrLogFileName(LogFileName,12);
  StrLogFileName=StrLogFileName.substr(0,StrLogFileName.find('\000'));
  cout << left << setw(19) << "LogFileName: " << StrLogFileName << "\n";

  char  SelfTestResult[2];
  ats.read(reinterpret_cast<char *>(&SelfTestResult), sizeof(SelfTestResult));
  string StrSelfTestResult(SelfTestResult,2);
  StrSelfTestResult =
    StrSelfTestResult.substr(0,StrSelfTestResult.find('\000'));
  cout << left << setw(19) << "SelfTestResult: " << StrSelfTestResult << "\n";

  char  Reserved5[2];
  ats.read(reinterpret_cast<char *>(&Reserved5), sizeof(Reserved5));
  string StrReserved5(Reserved5,2);
  StrReserved5=StrReserved5.substr(0,StrReserved5.find('\000'));
  cout << left << setw(19) << "Reserved5: " << StrReserved5 << "\n";

  int16_t  CalFreqs;
  ats.read(reinterpret_cast<char *>(&CalFreqs), sizeof(CalFreqs));
  cout << left << setw(19) << "CalFreqs: "
       << lexical_cast<string>(static_cast<int>(CalFreqs)) << "\n";

  int16_t  CalEntryLength;
  ats.read(reinterpret_cast<char *>(&CalEntryLength), sizeof(CalEntryLength));
  cout << left << setw(19) << "CalEntryLength: "
       << lexical_cast<string>(static_cast<int>(CalEntryLength)) << "\n";

  int16_t  CalVersion;
  ats.read(reinterpret_cast<char *>(&CalVersion), sizeof(CalVersion));
  cout << left << setw(19) << "CalVersion: "
       << lexical_cast<string>(static_cast<int>(CalVersion)) << "\n";

  int16_t  CalStartAddress;
  ats.read(reinterpret_cast<char *>(&CalStartAddress),
	   sizeof(CalStartAddress));
  cout << left << setw(19) << "CalStartAddress: "
       << lexical_cast<string>(static_cast<int>(CalStartAddress)) << "\n";

  char  LFFilters[8];
  ats.read(reinterpret_cast<char *>(&LFFilters), sizeof(LFFilters));
  string StrLFFilters(LFFilters,8);
  StrLFFilters=StrLFFilters.substr(0,StrLFFilters.find('\000'));
  cout << left << setw(19) << "LFFilters: " << StrLFFilters << "\n";

  char  ADU06CalFilename[12];
  ats.read(reinterpret_cast<char *>(&ADU06CalFilename),
	   sizeof(ADU06CalFilename));
  string StrADU06CalFilename(ADU06CalFilename,12);
  StrADU06CalFilename =
    StrADU06CalFilename.substr(0,StrADU06CalFilename.find('\000'));
  cout << left << setw(19) << "ADU06CalFilename: "
       << StrADU06CalFilename << "\n";

  int32_t  ADUCalTimeDate;
  ats.read(reinterpret_cast<char *>(&ADUCalTimeDate), sizeof(ADUCalTimeDate));
  cout << left << setw(19) << "ADUCalTimeDate: "
       << lexical_cast<string>(static_cast<int>(ADUCalTimeDate)) << "\n";

  char  SensorCalFilename[12];
  ats.read(reinterpret_cast<char *>(&SensorCalFilename),
	   sizeof(SensorCalFilename));
  string StrSensorCalFilename(SensorCalFilename,12);
  StrSensorCalFilename =
    StrSensorCalFilename.substr(0,StrSensorCalFilename.find('\000'));
  cout << left << setw(19) << "SensorCalFilename: "
       << StrSensorCalFilename << "\n";

  int32_t  SensorCalTimeDate;
  ats.read(reinterpret_cast<char *>(&SensorCalTimeDate),
	   sizeof(SensorCalTimeDate));
  cout << left << setw(19) << "SensorCalTimeDate: "
       << lexical_cast<string>(static_cast<int>(SensorCalTimeDate)) << "\n";

  float  PowerlineFreq1;
  ats.read(reinterpret_cast<char *>(&PowerlineFreq1), sizeof(PowerlineFreq1));
  cout << left << setw(19) << "PowerlineFreq1: "
       << lexical_cast<string>(PowerlineFreq1) << "\n";

  float  PowerlineFreq2;
  ats.read(reinterpret_cast<char *>(&PowerlineFreq2), sizeof(PowerlineFreq2));
  cout << left << setw(19) << "PowerlineFreq2: "
       << lexical_cast<string>(PowerlineFreq2) << "\n";

  char  HFFilters[8];
  ats.read(reinterpret_cast<char *>(&HFFilters), sizeof(HFFilters));
  string StrHFFilters(HFFilters,8);
  StrHFFilters=StrHFFilters.substr(0,StrHFFilters.find('\000'));
  cout << left << setw(19) << "HFFilters: " << StrHFFilters << "\n";

  float  CSAMTFreq;
  ats.read(reinterpret_cast<char *>(&CSAMTFreq), sizeof(CSAMTFreq));
  cout << left << setw(19) << "CSAMTFreq: "
       << lexical_cast<string>(CSAMTFreq) << "\n";

  int16_t  CSAMTBlocks;
  ats.read(reinterpret_cast<char *>(&CSAMTBlocks), sizeof(CSAMTBlocks));
  cout << left << setw(19) << "CSAMTBlocks: "
       << lexical_cast<string>(static_cast<int>(CSAMTBlocks)) << "\n";

  int16_t  CSAMTStacksPBlock;
  ats.read(reinterpret_cast<char *>(&CSAMTStacksPBlock),
	   sizeof(CSAMTStacksPBlock));
  cout << left << setw(19) << "CSAMTStacksPBlock: "
       << lexical_cast<string>(static_cast<int>(CSAMTStacksPBlock)) << "\n";

  int32_t  CSAMTBlockLength;
  ats.read(reinterpret_cast<char *>(&CSAMTBlockLength),
	   sizeof(CSAMTBlockLength));
  cout << left << setw(19) << "CSAMTBlockLength: "
       << lexical_cast<string>(static_cast<int>(CSAMTBlockLength)) << "\n";

  char  ADBBoardType[4];
  ats.read(reinterpret_cast<char *>(&ADBBoardType), sizeof(ADBBoardType));
  string StrADBBoardType(ADBBoardType,4);
  StrADBBoardType=StrADBBoardType.substr(0,StrADBBoardType.find('\000'));
  cout << left << setw(19) << "ADBBoardType: " << StrADBBoardType << "\n";

  char  Client[16];
  ats.read(reinterpret_cast<char *>(&Client), sizeof(Client));
  string StrClient(Client,16);
  StrClient=StrClient.substr(0,StrClient.find('\000'));
  cout << left << setw(19) << "Client: " << StrClient << "\n";

  char  Contractor[16];
  ats.read(reinterpret_cast<char *>(&Contractor), sizeof(Contractor));
  string StrContractor(Contractor,16);
  StrContractor=StrContractor.substr(0,StrContractor.find('\000'));
  cout << left << setw(19) << "Contractor: " << StrContractor << "\n";

  char  Area[16];
  ats.read(reinterpret_cast<char *>(&Area), sizeof(Area));
  string StrArea(Area,16);
  StrArea=StrArea.substr(0,StrArea.find('\000'));
  cout << left << setw(19) << "Area: " << StrArea << "\n";

  char  SurveyID[16];
  ats.read(reinterpret_cast<char *>(&SurveyID), sizeof(SurveyID));
  string StrSurveyID(SurveyID,16);
  StrSurveyID=StrSurveyID.substr(0,StrSurveyID.find('\000'));
  cout << left << setw(19) << "SurveyID: " << StrSurveyID << "\n";

  char  Operator[16];
  ats.read(reinterpret_cast<char *>(&Operator), sizeof(Operator));
  string StrOperator(Operator,16);
  StrOperator=StrOperator.substr(0,StrOperator.find('\000'));
  cout << left << setw(19) << "Operator: " << StrOperator << "\n";

  char  Reserved[112];
  ats.read(reinterpret_cast<char *>(&Reserved), sizeof(Reserved));
  string StrReserved(Reserved,112);
  StrReserved=StrReserved.substr(0,StrReserved.find('\000'));
  cout << left << setw(19) << "Reserved: " << StrReserved << "\n";

  char  XmlHeader[64];
  ats.read(reinterpret_cast<char *>(&XmlHeader), sizeof(XmlHeader));
  string StrXmlHeader(XmlHeader,64);
  StrXmlHeader=StrXmlHeader.substr(0,StrXmlHeader.find('\000'));
  cout << left << setw(19) << "XmlHeader: " << StrXmlHeader << "\n";

  char  Comments[512];
  ats.read(reinterpret_cast<char *>(&Comments), sizeof(Comments));
  string StrComments(Comments,512);
  StrComments=StrComments.substr(0,StrComments.find('\000'));
  cout << left << setw(19) << "Comments: " << StrComments << "\n";
}
