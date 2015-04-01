#include "system_dependent.hpp"
#include <pwd.h> //passwd getpwuid
#include <unistd.h> // geteuid

std::string ats2asc::diretorio_do_usuario()
{
  struct passwd *p = getpwuid(geteuid());
  return std::string(p->pw_dir);
}
