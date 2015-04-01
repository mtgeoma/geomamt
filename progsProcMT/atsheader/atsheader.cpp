#include <cstdlib> // exit
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <iomanip> // setw
#include "atsheader.hpp"

int main(int argc, char* argv[]) {
  using namespace std;
  using namespace boost;
  string ats_file_name;

  if(argc<2)
  {
    cerr << argv[0] << " read an ats file and return selected header data\n"
	 << "  or change its value.\n\n"
	 << "  usage to read all parameters: " << argv[0] << " file.ats\n\n"
	 << "  usage to read some parameters: " << argv[0] << " file.ats"
	 << " parameter ...\n\n"
	 << "  usage to write: " << argv[0] << " file.ats parameter=value ...\n"
	 << "\n  First column in read all parameters shows which parameters"
	 << " could be used.\n"
	 << "  Note that some of them are only available in read mode.\n";
    exit(1);
  }

  ats_file_name=string(argv[1]);

  fstream ats;
  // before all, check ats size.
  // open ats file to read
  ats.open(ats_file_name.c_str(), ios_base::binary | ios_base::in);
  if(ats.fail())
  {
    cerr << "couldn't open '"<< ats_file_name << "' to read"<< endl;
    exit(1);
  }

  ats.seekg(0, std::ios::end);
  long ats_file_size = static_cast<long>(ats.tellg());
  if(ats_file_size == 0)
  {
    cerr << ats_file_name << ": empty file"<< endl;
    exit(1);
  }

  ats.seekg(0x000,ios::beg);
  int16_t  HeaderLength;
  ats.read(reinterpret_cast<char *>(&HeaderLength), sizeof(HeaderLength));
  if(ats_file_size < HeaderLength)
  {
    cerr << ats_file_name << ": file size is less than HeaderLength("
	 << HeaderLength << ")" << endl;
    exit(1);
  }
  ats.close();

  if (argc == 2)
  {
    // open ats file to read
    ats.open(ats_file_name.c_str(), ios_base::binary | ios_base::in);
    if(ats.fail())
    {
      cerr << "couldn't open '"<< ats_file_name << "' to read"<< endl;
      exit(1);
    }
    cout << left << setw(19) << "FileName: " << ats_file_name
	 << " [only read mode]\n";
    read_header(ats);
    ats.close();
  }
  else
  {
    bool only_read;
    if ( string(argv[2]).find("=") !=string::npos )
    {
      only_read = false;
      // open ats file to write
      ats.open(ats_file_name.c_str(), ios_base::binary |
	       ios_base::out | ios_base::in);
      if(ats.fail())
      {
	cerr << "couldn't open '"<< ats_file_name << "' to write"<< endl;
	exit(1);
      }
    }
    else
    {
      only_read = true;
      // open ats file to read
      ats.open(ats_file_name.c_str(), ios_base::binary | ios_base::in);
      if(ats.fail())
      {
	cerr << "couldn't open '"<< ats_file_name << "' to read"<< endl;
	exit(1);
      }
    }

    for(int i=2; i<argc; i++)
    {
      if(only_read)
      {
	string parameter = string(argv[i]);
	if ( parameter == "FileName" )
	{
	  cout << " " << ats_file_name;
	}
	else
	{
	  read_header(ats, parameter);
	}
      }
      else
      {
	string parameter = string(argv[i]);
	vector< string > v_str;
	split(v_str,parameter,is_any_of("="));
	if ( v_str.size() != 2 )
	{
	  cerr << " Error! Skiping writer parameter: " << parameter << "\n"
	       << " Writer parameter must have one \"=\" character\n";
	}
	cout << "will write " << v_str[0] << " with value " << v_str[1] << endl;
	write_header(ats, v_str[0], v_str[1] );
      }
    }

    if(only_read)
    {
      cout << endl;
    }

    ats.close();
  }
}
