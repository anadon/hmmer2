/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* sre_ctype.c
 *
 * For portability. Some systems have functions tolower, toupper
 * as macros (for instance, MIPS M-2000 RISC/os!)
 *
 * SVN $Id: sre_ctype.c 1530 2005-12-13 20:53:08Z eddy $
 */

#include "squidconf.h"

#include <ctype.h>
#include "squid.h"

int
sre_tolower(int c) {
  if (isupper(c)) return tolower(c);
  else return c;
}

int
sre_toupper(int c) {
  if (islower(c)) return toupper(c);
  else return c;
}

