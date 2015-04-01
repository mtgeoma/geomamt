#include "atsheader.hpp"

void
write_header (std::fstream &ats, std::string parameter, std::string value)
{
  using namespace std;
  using namespace boost;
  if ( parameter == "HeaderLength" )
  {
    ats.seekp(0x000,ios::beg);
    int16_t HeaderLength = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&HeaderLength), sizeof(HeaderLength));
  }
  else if ( parameter == "HeaderVers" )
  {
    ats.seekp(0x002,ios::beg);
    int16_t HeaderVers = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&HeaderVers), sizeof(HeaderVers));
  }
  else if ( parameter == "Samples" )
  {
    ats.seekp(0x004,ios::beg);
    uint32_t Samples = lexical_cast<uint32_t>(value);
    ats.write(reinterpret_cast<char *>(&Samples), sizeof(Samples));
  }
  else if ( parameter == "SampleFreq" )
  {
    ats.seekp(0x008,ios::beg);
    float SampleFreq = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&SampleFreq), sizeof(SampleFreq));
  }
  else if ( parameter == "StartDateTime" )
  {
    ats.seekp(0x00C,ios::beg);
    int32_t StartDateTime = lexical_cast<int32_t>(value);
    ats.write(reinterpret_cast<char *>(&StartDateTime), sizeof(StartDateTime));
  }
  else if ( parameter == "LSBMV" )
  {
    ats.seekp(0x010,ios::beg);
    double LSBMV = lexical_cast<double>(value);
    ats.write(reinterpret_cast<char *>(&LSBMV), sizeof(LSBMV));
  }
  else if ( parameter == "GMTOffset" )
  {
    ats.seekp(0x018,ios::beg);
    int32_t GMTOffset = lexical_cast<int32_t>(value);
    ats.write(reinterpret_cast<char *>(&GMTOffset), sizeof(GMTOffset));
  }
  else if ( parameter == "OrigSampleFreq" )
  {
    ats.seekp(0x01C,ios::beg);
    float OrigSampleFreq = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&OrigSampleFreq),
	      sizeof(OrigSampleFreq));
  }
  else if ( parameter == "ADUSerNum" )
  {
    ats.seekp(0x020,ios::beg);
    int16_t ADUSerNum = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&ADUSerNum), sizeof(ADUSerNum));
  }
  else if ( parameter == "ADCSerNum" )
  {
    ats.seekp(0x022,ios::beg);
    int16_t ADCSerNum = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&ADCSerNum), sizeof(ADCSerNum));
  }
  else if ( parameter == "ChanNo" )
  {
    ats.seekp(0x024,ios::beg);
    int16_t ChanNo = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&ChanNo), 1);
  }
  else if ( parameter == "Chopper" )
  {
    ats.seekp(0x025,ios::beg);
    int16_t Chopper = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&Chopper), 1);
  }
  else if ( parameter == "ChanType" )
  {
    ats.seekp(0x026,ios::beg);
    if ( value.size() < 2 )
    {
      value.append(2-value.size(),'\000');
    }
    char ChanType[2];
    value.copy(ChanType,2);
    ats.write(reinterpret_cast<char *>(&ChanType), sizeof(ChanType));
  }
  else if ( parameter == "SensorType" )
  {
    ats.seekp(0x028,ios::beg);
    if ( value.size() < 6 )
    {
      value.append(6-value.size(),'\000');
    }
    char SensorType[6];
    value.copy(SensorType,6);
    ats.write(reinterpret_cast<char *>(&SensorType), sizeof(SensorType));
  }
  else if ( parameter == "SensorSerNum" )
  {
    ats.seekp(0x02E,ios::beg);
    int16_t SensorSerNum = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&SensorSerNum), sizeof(SensorSerNum));
  }
  else if ( parameter == "PosX1" )
  {
    ats.seekp(0x030,ios::beg);
    float PosX1 = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&PosX1), sizeof(PosX1));
  }
  else if ( parameter == "PosY1" )
  {
    ats.seekp(0x034,ios::beg);
    float PosY1 = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&PosY1), sizeof(PosY1));
  }
  else if ( parameter == "PosZ1" )
  {
    ats.seekp(0x038,ios::beg);
    float PosZ1 = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&PosZ1), sizeof(PosZ1));
  }
  else if ( parameter == "PosX2" )
  {
    ats.seekp(0x03C,ios::beg);
    float PosX2 = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&PosX2), sizeof(PosX2));
  }
  else if ( parameter == "PosY2" )
  {
    ats.seekp(0x040,ios::beg);
    float PosY2 = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&PosY2), sizeof(PosY2));
  }
  else if ( parameter == "PosZ2" )
  {
    ats.seekp(0x044,ios::beg);
    float PosZ2 = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&PosZ2), sizeof(PosZ2));
  }
  else if ( parameter == "DipLength" )
  {
    ats.seekp(0x048,ios::beg);
    float DipLength = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&DipLength), sizeof(DipLength));
  }
  else if ( parameter == "Angle" )
  {
    ats.seekp(0x04C,ios::beg);
    float Angle = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&Angle), sizeof(Angle));
  }
  else if ( parameter == "ProbeRes" )
  {
    ats.seekp(0x050,ios::beg);
    float ProbeRes = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&ProbeRes), sizeof(ProbeRes));
  }
  else if ( parameter == "DCOffset" )
  {
    ats.seekp(0x054,ios::beg);
    float DCOffset = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&DCOffset), sizeof(DCOffset));
  }
  else if ( parameter == "PreGain" )
  {
    ats.seekp(0x058,ios::beg);
    float PreGain = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&PreGain), sizeof(PreGain));
  }
  else if ( parameter == "PostGain" )
  {
    ats.seekp(0x05C,ios::beg);
    float PostGain = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&PostGain), sizeof(PostGain));
  }
  else if ( parameter == "Latitude_msec" )
  {
    ats.seekp(0x060,ios::beg);
    int32_t Latitude_msec = lexical_cast<int32_t>(value);
    ats.write(reinterpret_cast<char *>(&Latitude_msec), sizeof(Latitude_msec));
  }
  else if ( parameter == "Longitude_msec" )
  {
    ats.seekp(0x064,ios::beg);
    int32_t Longitude_msec = lexical_cast<int32_t>(value);
    ats.write(reinterpret_cast<char *>(&Longitude_msec),
	      sizeof(Longitude_msec));
  }
  else if ( parameter == "Elevation_cm" )
  {
    ats.seekp(0x068,ios::beg);
    int32_t Elevation_cm = lexical_cast<int32_t>(value);
    ats.write(reinterpret_cast<char *>(&Elevation_cm), sizeof(Elevation_cm));
  }
  else if ( parameter == "LatLongType" )
  {
    ats.seekp(0x06C,ios::beg);
    char LatLongType = lexical_cast<char>(value);
    ats.write(reinterpret_cast<char *>(&LatLongType), sizeof(LatLongType));
  }
  else if ( parameter == "AddCoordType" )
  {
    ats.seekp(0x06D,ios::beg);
    char AddCoordType = lexical_cast<char>(value);
    ats.write(reinterpret_cast<char *>(&AddCoordType), sizeof(AddCoordType));
  }
  else if ( parameter == "RefMedian" )
  {
    ats.seekp(0x06E,ios::beg);
    int16_t RefMedian = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&RefMedian), sizeof(RefMedian));
  }
  else if ( parameter == "XCoord" )
  {
    ats.seekp(0x070,ios::beg);
    double XCoord = lexical_cast<double>(value);
    ats.write(reinterpret_cast<char *>(&XCoord), sizeof(XCoord));
  }
  else if ( parameter == "YCoord" )
  {
    ats.seekp(0x078,ios::beg);
    double YCoord = lexical_cast<double>(value);
    ats.write(reinterpret_cast<char *>(&YCoord), sizeof(YCoord));
  }
  else if ( parameter == "GPSStat" )
  {
    ats.seekp(0x080,ios::beg);
    char GPSStat = lexical_cast<char>(value);
    ats.write(reinterpret_cast<char *>(&GPSStat), sizeof(GPSStat));
  }
  else if ( parameter == "GPSAccuracy" )
  {
    ats.seekp(0x081,ios::beg);
    int16_t GPSAccuracy = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&GPSAccuracy), 1);
  }
  else if ( parameter == "UTCOffset" )
  {
    ats.seekp(0x082,ios::beg);
    int16_t UTCOffset = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&UTCOffset), sizeof(UTCOffset));
  }
  else if ( parameter == "SystemType" )
  {
    ats.seekp(0x084,ios::beg);
    if ( value.size() < 12 )
    {
      value.append(12-value.size(),'\000');
    }
    char SystemType[12];
    value.copy(SystemType,12);
    ats.write(reinterpret_cast<char *>(&SystemType), sizeof(SystemType));
  }
  else if ( parameter == "SurveyHeaderName" )
  {
    ats.seekp(0x090,ios::beg);
    if ( value.size() < 12 )
    {
      value.append(12-value.size(),'\000');
    }
    char SurveyHeaderName[12];
    value.copy(SurveyHeaderName,12);
    ats.write(reinterpret_cast<char *>(&SurveyHeaderName),
	      sizeof(SurveyHeaderName));
  }
  else if ( parameter == "MeasType" )
  {
    ats.seekp(0x09C,ios::beg);
    if ( value.size() < 4 )
    {
      value.append(4-value.size(),'\000');
    }
    char MeasType[4];
    value.copy(MeasType,4);
    ats.write(reinterpret_cast<char *>(&MeasType), sizeof(MeasType));
  }
  else if ( parameter == "LogFileName" )
  {
    ats.seekp(0x0A0,ios::beg);
    if ( value.size() < 12 )
    {
      value.append(12-value.size(),'\000');
    }
    char LogFileName[12];
    value.copy(LogFileName,12);
    ats.write(reinterpret_cast<char *>(&LogFileName), sizeof(LogFileName));
  }
  else if ( parameter == "SelfTestResult" )
  {
    ats.seekp(0x0AC,ios::beg);
    if ( value.size() < 2 )
    {
      value.append(2-value.size(),'\000');
    }
    char SelfTestResult[2];
    value.copy(SelfTestResult,2);
    ats.write(reinterpret_cast<char *>(&SelfTestResult),
	      sizeof(SelfTestResult));
  }
  else if ( parameter == "Reserved5" )
  {
    ats.seekp(0x0AE,ios::beg);
    if ( value.size() < 2 )
    {
      value.append(2-value.size(),'\000');
    }
    char Reserved5[2];
    value.copy(Reserved5,2);
    ats.write(reinterpret_cast<char *>(&Reserved5), sizeof(Reserved5));
  }
  else if ( parameter == "CalFreqs" )
  {
    ats.seekp(0x0B0,ios::beg);
    int16_t CalFreqs = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&CalFreqs), sizeof(CalFreqs));
  }
  else if ( parameter == "CalEntryLength" )
  {
    ats.seekp(0x0B2,ios::beg);
    int16_t CalEntryLength = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&CalEntryLength),
	      sizeof(CalEntryLength));
  }
  else if ( parameter == "CalVersion" )
  {
    ats.seekp(0x0B4,ios::beg);
    int16_t CalVersion = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&CalVersion), sizeof(CalVersion));
  }
  else if ( parameter == "CalStartAddress" )
  {
    ats.seekp(0x0B6,ios::beg);
    int16_t CalStartAddress = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&CalStartAddress),
	      sizeof(CalStartAddress));
  }
  else if ( parameter == "LFFilters" )
  {
    ats.seekp(0x0B8,ios::beg);
    if ( value.size() < 8 )
    {
      value.append(8-value.size(),'\000');
    }
    char LFFilters[8];
    value.copy(LFFilters,8);
    ats.write(reinterpret_cast<char *>(&LFFilters), sizeof(LFFilters));
  }
  else if ( parameter == "ADU06CalFilename" )
  {
    ats.seekp(0x0C0,ios::beg);
    if ( value.size() < 12 )
    {
      value.append(12-value.size(),'\000');
    }
    char ADU06CalFilename[12];
    value.copy(ADU06CalFilename,12);
    ats.write(reinterpret_cast<char *>(&ADU06CalFilename),
	      sizeof(ADU06CalFilename));
  }
  else if ( parameter == "ADUCalTimeDate" )
  {
    ats.seekp(0x0CC,ios::beg);
    int32_t ADUCalTimeDate = lexical_cast<int32_t>(value);
    ats.write(reinterpret_cast<char *>(&ADUCalTimeDate),
	      sizeof(ADUCalTimeDate));
  }
  else if ( parameter == "SensorCalFilename" )
  {
    ats.seekp(0x0D0,ios::beg);
    if ( value.size() < 12 )
    {
      value.append(12-value.size(),'\000');
    }
    char SensorCalFilename[12];
    value.copy(SensorCalFilename,12);
    ats.write(reinterpret_cast<char *>(&SensorCalFilename),
	      sizeof(SensorCalFilename));
  }
  else if ( parameter == "SensorCalTimeDate" )
  {
    ats.seekp(0x0DC,ios::beg);
    int32_t SensorCalTimeDate = lexical_cast<int32_t>(value);
    ats.write(reinterpret_cast<char *>(&SensorCalTimeDate),
	      sizeof(SensorCalTimeDate));
  }
  else if ( parameter == "PowerlineFreq1" )
  {
    ats.seekp(0x0E0,ios::beg);
    float PowerlineFreq1 = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&PowerlineFreq1),
	      sizeof(PowerlineFreq1));
  }
  else if ( parameter == "PowerlineFreq2" )
  {
    ats.seekp(0x0E4,ios::beg);
    float PowerlineFreq2 = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&PowerlineFreq2),
	      sizeof(PowerlineFreq2));
  }
  else if ( parameter == "HFFilters" )
  {
    ats.seekp(0x0E8,ios::beg);
    if ( value.size() < 8 )
    {
      value.append(8-value.size(),'\000');
    }
    char HFFilters[8];
    value.copy(HFFilters,8);
    ats.write(reinterpret_cast<char *>(&HFFilters), sizeof(HFFilters));
  }
  else if ( parameter == "CSAMTFreq" )
  {
    ats.seekp(0x0F0,ios::beg);
    float CSAMTFreq = lexical_cast<float>(value);
    ats.write(reinterpret_cast<char *>(&CSAMTFreq), sizeof(CSAMTFreq));
  }
  else if ( parameter == "CSAMTBlocks" )
  {
    ats.seekp(0x0F4,ios::beg);
    int16_t CSAMTBlocks = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&CSAMTBlocks), sizeof(CSAMTBlocks));
  }
  else if ( parameter == "CSAMTStacksPBlock" )
  {
    ats.seekp(0x0F6,ios::beg);
    int16_t CSAMTStacksPBlock = lexical_cast<int16_t>(value);
    ats.write(reinterpret_cast<char *>(&CSAMTStacksPBlock),
	      sizeof(CSAMTStacksPBlock));
  }
  else if ( parameter == "CSAMTBlockLength" )
  {
    ats.seekp(0x0F8,ios::beg);
    int32_t CSAMTBlockLength = lexical_cast<int32_t>(value);
    ats.write(reinterpret_cast<char *>(&CSAMTBlockLength),
	      sizeof(CSAMTBlockLength));
  }
  else if ( parameter == "ADBBoardType" )
  {
    ats.seekp(0x0FC,ios::beg);
    if ( value.size() < 4 )
    {
      value.append(4-value.size(),'\000');
    }
    char ADBBoardType[4];
    value.copy(ADBBoardType,4);
    ats.write(reinterpret_cast<char *>(&ADBBoardType), sizeof(ADBBoardType));
  }
  else if ( parameter == "Client" )
  {
    ats.seekp(0x100,ios::beg);
    if ( value.size() < 16 )
    {
      value.append(16-value.size(),'\000');
    }
    char Client[16];
    value.copy(Client,16);
    ats.write(reinterpret_cast<char *>(&Client), sizeof(Client));
  }
  else if ( parameter == "Contractor" )
  {
    ats.seekp(0x110,ios::beg);
    if ( value.size() < 16 )
    {
      value.append(16-value.size(),'\000');
    }
    char Contractor[16];
    value.copy(Contractor,16);
    ats.write(reinterpret_cast<char *>(&Contractor), sizeof(Contractor));
  }
  else if ( parameter == "Area" )
  {
    ats.seekp(0x120,ios::beg);
    if ( value.size() < 16 )
    {
      value.append(16-value.size(),'\000');
    }
    char Area[16];
    value.copy(Area,16);
    ats.write(reinterpret_cast<char *>(&Area), sizeof(Area));
  }
  else if ( parameter == "SurveyID" )
  {
    ats.seekp(0x130,ios::beg);
    if ( value.size() < 16 )
    {
      value.append(16-value.size(),'\000');
    }
    char SurveyID[16];
    value.copy(SurveyID,16);
    ats.write(reinterpret_cast<char *>(&SurveyID), sizeof(SurveyID));
  }
  else if ( parameter == "Operator" )
  {
    ats.seekp(0x140,ios::beg);
    if ( value.size() < 16 )
    {
      value.append(16-value.size(),'\000');
    }
    char Operator[16];
    value.copy(Operator,16);
    ats.write(reinterpret_cast<char *>(&Operator), sizeof(Operator));
  }
  else if ( parameter == "Reserved" )
  {
    ats.seekp(0x150,ios::beg);
    if ( value.size() < 112 )
    {
      value.append(112-value.size(),'\000');
    }
    char Reserved[112];
    value.copy(Reserved,112);
    ats.write(reinterpret_cast<char *>(&Reserved), sizeof(Reserved));
  }
  else if ( parameter == "XmlHeader" )
  {
    ats.seekp(0x1C0,ios::beg);
    if ( value.size() < 64 )
    {
      value.append(64-value.size(),'\000');
    }
    char XmlHeader[64];
    value.copy(XmlHeader,64);
    ats.write(reinterpret_cast<char *>(&XmlHeader), sizeof(XmlHeader));
  }
  else if ( parameter == "Comments" )
  {
    ats.seekp(0x200,ios::beg);
    if ( value.size() < 512 )
    {
      value.append(512-value.size(),'\000');
    }
    char Comments[512];
    value.copy(Comments,512);
    ats.write(reinterpret_cast<char *>(&Comments), sizeof(Comments));
  }
  else
  {
    cerr << "unknow write parameter " << parameter << "\n";
  }
}
