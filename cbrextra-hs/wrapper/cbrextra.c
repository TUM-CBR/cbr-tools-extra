#include <stdlib.h>
#include "HsFFI.h"
#include "Bio\Cbr_stub.h"

HsBool cbrextra_init(void){
  int argc = 3;
  char *argv[] = { "libcbrextra", "+RTS", "-A32m", NULL };
  char **pargv = argv;

  // Initialize Haskell runtime
  hs_init_with_rtsopts(&argc, &pargv);

  // do any other initialization here and
  // return false if there was a problem
  return HS_BOOL_TRUE;
}

void cbrdemo() {
    demo();
}

void cbrextra_end(void){
  hs_exit();
}