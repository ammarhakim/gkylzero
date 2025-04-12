#include "atomic_data.h"
char* get_adas_file_loc(char* filename);
int file_isreg(const char *path);
atomic_data loadadf11(char* filename);
char* concat(const char *s1, const char *s2);
