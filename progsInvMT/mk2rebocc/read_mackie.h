#ifndef READ_MACKIE_H
#define READ_MACKIE_H
#include <string>
using namespace std;

typedef struct {
  double T;
  double rho;
  double phi;
  double err;
} data;

data read_line_mackie(string line, string mode_type);

#endif // READ_MACKIE_H
