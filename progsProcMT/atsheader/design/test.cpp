#include <iostream>
#include <string>
#include <map>

enum MTXtype {String,Double,Float,qint16,qint32,quint32};

struct headerElement{
  MTXtype type;
  int addr;
  int size;
  // headerElement(MTXtype type_=String,int addr_=0,int size_=0)
  headerElement(MTXtype type_,int addr_,int size_)
    :type(type_),addr(addr_),size(size_){};
};

int main()
{

  std::map<std::string,headerElement> metronixATSHeader;
  // headerElement HeaderLength(qint16,0x000,2);
  headerElement HeaderLength=headerElement(qint16,0x000,2);
  // qint16  HeaderLength      0x000   2
  // metronixATSHeader["HeaderLength"]=headerElement(qint16,0x000,2);
  metronixATSHeader.insert(std::make_pair("HeaderLength",
					  headerElement(qint16,0x000,2)));
  // std::cout << metronixATSHeader["HeaderLength"].size << "\n";
  std::cout << HeaderLength.type << "\n";
}
