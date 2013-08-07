/*
*  Copyright (C) 2007 David Chin, Julien Sylvestre
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#ifndef _LPC_H
#define _LPC_H

#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <lal/LALMalloc.h>

/*<lalLaTeX file="LPCH">
\section{Header \texttt{LPC.h}}
\label{s:LPC.h}

Functions for linear predictor filters.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LPC.h>
\end{verbatim}

\subsection*{Error Conditions}
\input{LPCHErrTab}
</lalLaTeX>*/

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID (LPCH, "$Id$");

  /******************************** <lalErrTable file="LPCHErrTab"> */
#define LPCH_EIN 1
#define LPCH_EMEM 2
#define LPCH_ENUM 3

#define LPCH_MSGEIN "invalid input"
#define LPCH_MSGEMEM "memory error"
#define LPCH_MSGENUM "numerical error"
/*************************************************** </lalErrTable> */

void LALPolystab(LALStatus *status,
		 REAL4Vector *a);

void LALLPC(LALStatus *status,
	    REAL4Vector *aout,    /* returned filter coefficients */
	    REAL4Vector *x,    /* training data */
	    UINT4 p            /* filter order */
	    );

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif