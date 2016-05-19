#include "randombytes.h"
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace nfl {

static int fd = -1;

void randombytes(unsigned char *x, unsigned long long xlen) {
  int i;

  if (fd == -1) {
    for (;;) {
      fd = open("/dev/urandom", O_RDONLY);
      if (fd != -1) break;
      sleep(1);
    }
  }

  while (xlen > 0) {
    i = (xlen < 1048576) ? xlen : 1048576;
    i = read(fd, x, i);

    if (i < 1) {
      sleep(1);
      continue;
    }

    x += i;
    xlen -= i;
  }
}
}