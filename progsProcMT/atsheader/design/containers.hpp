
const int minimalHeaderLength=0x200+512; // highest addrs element plus its size

// make a class ats_header? Hide file access and keep different types members.
class header {
public:
  header();
  double get_double(std::string label);
  long get_long(std::string label);
  std::string get_string(std::string label);
private:
  // ats header access datails.

  enum MTXtype {String,Double,Float,qint16,qint32,quint32};

  struct headerElement{
    MTXtype type;
    int addr;
    int size;
    headerElement(MTXtype type_,int addr_,int size_)
      :type(type_),addr(addr_),size(size_){};
  };

  typedef  mapAddr std::map<std::string,headerElement>;

  mapATSHeader.insert(std::make_pair("HeaderLength",
				     headerElement(qint16,0x000,2)));
};

// don't appears to be a good choice:
class parameter {
public:
  parameter(long addr, long size, <T> value(?)):
  {};
  get_value(std::fstream ats);
private:
  long addr;
  long size;
  <T> value;
};

map<string,parameter> header;
header[label]=parameter(addr,size,<T>);

parameter::get_value (std::fstream ats)
{
  ats.seekg(addr,ios:beg);
  ats.read(reinterpret_cast<char *>(&value),
	   sizeof(size));
}
