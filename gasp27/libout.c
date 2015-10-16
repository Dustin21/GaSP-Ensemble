/*****************************************************************/
/*   LOW-LEVEL OUTPUT ROUTINES                                   */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--2000.                  */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "define.h"
#include "implem.h"
#include "lib.h"

static FILE *LogFile = NULL;

static char    Buf[MAXTOK+1];
static Matrix  ErrorMat;

boolean        ErrorSave = NO;
string         ErrorVar  = NULL;
int            ErrorSeverityLevel = SEV_ERROR;
size_t         ErrorTry  = 0;

static string  SeverityStr[] = SEVERITY_STRS;

/******************************+++********************************/
/*                                                               */
/*   void      SetLogFile(FILE *fp)                              */
/*                                                               */
/*   Purpose:  Set pointer for the log file.                     */
/*                                                               */
/*   Version:  1990 September 30                                 */
/*                                                               */
/*****************************************************************/

void SetLogFile(FILE *fp)
{
     LogFile = fp;

     return;
}

/******************************+++********************************/
/*                                                               */
/*   FILE      *GetLogFile(void)                                 */
/*                                                               */
/*   Purpose:  Return a pointer to the log file.                 */
/*                                                               */
/*   Version:  1991 September 2                                  */
/*                                                               */
/*****************************************************************/

FILE *GetLogFile(void)
{
     return LogFile;
}

/******************************+++********************************/
void FileOutput(FILE *fp, const string Format, ...)
/*****************************************************************/
/*   Purpose:  Send output in variable argument list to fp.      */
/*             If fp == stdout, then also send output to the log */
/*             file.                                             */
/*                                                               */
/*   96.06.18: Put output in Buf to avoid calling vfprintf       */
/*             twice: First call changes Args with some          */
/*             compilers.                                        */
/*   97.07.04: fputs instead of fprintf - fprintf interprets "%" */
/*             in Buf as a conversion character.                 */
/*                                                               */
/*   Version:  1996.07.04                                        */
/*****************************************************************/
{
     va_list   Args;

     /* Put the output in Buf. */
     va_start(Args, Format);
     vsprintf(Buf, Format, Args);
     va_end(Args);

     /* Output to fp. */
     fputs(Buf, fp);

     if (fp == stdout && LogFile != NULL)
          /* Output to the log file. */
          fputs(Buf, LogFile);

     return;
}

/******************************+++********************************/
/*   void      Output(const string Format, ...)                  */
/*                                                               */
/*   Purpose:  Output variable argument list to stdout and to    */
/*             the log file.                                     */
/*             Error(), Fatal(), and Incompatibility() are       */
/*             similar except:                                   */
/*                                                               */
/*                                 precedes the output by        */
/*                                                               */
/*             Error()             "Error/Warning: "             */
/*             Fatal()             "Fatal error: "               */
/*             Incompatibility     "Incompatibility: "           */
/*                                                               */
/*             and Fatal() also calls exit(1).                   */
/*                                                               */
/* 1995.05.02:                                                   */
/* 2000.02.15: Output("\n") added to Fatal().                    */
/*****************************************************************/

void Output(const string Format, ...)
{
     va_list   Args;

     va_start(Args, Format);

     OutputVA(Format, Args);

     va_end(Args);

     return;
}

void Error(const string Format, ...)
{
     va_list   Args;

     va_start(Args, Format);

     if (ErrorSave)
          /* Save the message in ErrorMat. */
          ErrorToMat(SeverityStr[ErrorSeverityLevel], Format, Args);
     else
     {
          Output("%s: ", SeverityStr[ErrorSeverityLevel]);
          OutputVA(Format, Args);
     }

     va_end(Args);

     return;
}

void Fatal(const string Format, ...)
{
     va_list   Args;

     va_start(Args, Format);

     Output("%s", "Fatal error: ");
     OutputVA(Format, Args);
     Output("\n");

     va_end(Args);

     exit(1);
}

void Incompatibility(const string Format, ...)
{
     va_list   Args;

     va_start(Args, Format);

     if (ErrorSave)
          /* Save the message in ErrorMat. */
          ErrorToMat("Incompatibility", Format, Args);
     else
     {
          Output("%s", "Incompatibility: ");
          OutputVA(Format, Args);
     }

     va_end(Args);

     return;
}

/******************************+++********************************/
void OutputVA(const string Format, va_list Args)
/*****************************************************************/
/*   Purpose:  Send output in variable argument list to stdout   */
/*             and to the log file.                              */
/*                                                               */
/*   96.06.18: Put output in Buf to avoid calling vprintf and    */
/*             vfprintf: First call changes Args with some       */
/*             compilers.                                        */
/*   97.07.04: fputs instead of printf/fprintf as they interpret */
/*             "%" in Buf as a conversion character.             */
/*                                                               */
/*   Version:  1996.07.04                                        */
/*****************************************************************/
{
     /* Put the output in Buf. */
     vsprintf(Buf, Format, Args);

     /* Output to stdout (puts outputs newline character). */
     fputs(Buf, stdout);

     /* Output to the log file. */
     if (LogFile != NULL)
          fputs(Buf, LogFile);
}


/******************************+++********************************/
void ErrorToMat(const string Severity, const string Format,
          va_list Args)
/*****************************************************************/
/*   Purpose:  Save error, warning, etc. messages in a matrix.   */
/*                                                               */
/*   Version:  1995 March 10                                     */
/*****************************************************************/
{
     size_t    j, nRowsOld, LastTry;
     size_t    *Try;
     string    LastMess, LastVar, TermPtr;
     string    *Message, *Variable;

     if (!MatInitialized(&ErrorMat))
     {
          MatInit(RECT, MIXED, YES, &ErrorMat);
          MatPutText(&ErrorMat,
               "The following error messages were generated:\n");
     }

     Variable = MatStrColFind(&ErrorMat, VARIABLE, NO);
     Try      = MatSize_tColFind(&ErrorMat, "Try", NO);
     Message  = MatStrColFind(&ErrorMat, "Message", NO);

     nRowsOld = MatNumRows(&ErrorMat);

     LastVar = (Variable != NULL) ? Variable[nRowsOld - 1] : NULL;
     LastTry = (Try != NULL) ? Try[nRowsOld - 1] : 0;
     LastMess = (Message != NULL) ? Message[nRowsOld - 1] : NULL;

     /* Put the message in Buf. */
     vsprintf(Buf, Format, Args);

     /* Remove any terminating ".\n". */
     TermPtr = Buf + strlen(Buf) - 2;
     if (stricmp(TermPtr, ".\n") == 0)
          *TermPtr = NULL;

     if (stricmp(ErrorVar, LastVar) == 0 && ErrorTry == LastTry &&
               stricmp(Buf, LastMess) == 0)
          /* Do not repeat same message. */
          return;

     /* Allocate a new row for the new message. */
     MatReAlloc(nRowsOld + 1, MatNumCols(&ErrorMat), &ErrorMat);

     if (ErrorVar != NULL)
     {
          j = MatColumnAdd(VARIABLE, STRING, &ErrorMat);
          MatPutStrElem(&ErrorMat, nRowsOld, j, ErrorVar);
     }

     if (ErrorTry != 0)
     {
          j = MatColumnAdd("Try", SIZE_T, &ErrorMat);
          MatPutSize_tElem(&ErrorMat, nRowsOld, j, ErrorTry);
     }

     j = MatColumnAdd("Severity", STRING, &ErrorMat);
     MatPutStrElem(&ErrorMat, nRowsOld, j, Severity);

     j = MatColumnAdd("Message", STRING, &ErrorMat);
     MatPutStrElem(&ErrorMat, nRowsOld, j, Buf);

     return;
}

/******************************+++********************************/
void ErrorMatOut(void)
/*****************************************************************/
/*   Purpose:  Output matrix of error, warning, etc. messages.   */
/*                                                               */
/*   Version:  1995 May 11                                       */
/*****************************************************************/
{
     if (MatNumRows(&ErrorMat) > 0)
     {
          MatWriteBlock(&ErrorMat, NO, stdout);
          Output("\n");
          MatReAlloc(0, 0, &ErrorMat);
     }

     ErrorVar = NULL;
     ErrorTry = 0;
     ErrorSave = NO;
}

static size_t nTempCharsLast = 0;

/******************************+++********************************/
void OutputTemp(const string Format, ...)
/*****************************************************************/
/*   Purpose:  Output temporary message to stdout, which will be */
/*             overwritten by next temporary message.            */
/*                                                               */
/*   Version:  1996.04.05                                        */
/*****************************************************************/
{
     size_t    i, nTempChars;
     va_list   Args;

     va_start(Args, Format);

     /* Backspace previous message. */
     for (i = 0; i < nTempCharsLast; i++)
          printf("\b");

     /* Write new message. */
     nTempChars = vprintf(Format, Args);

     /* If new message is shorter, blank out remainder */
     /* of old message, then backspace blanks.         */
     for (i = nTempChars; i < nTempCharsLast; i++)
          printf(" ");
     for (i = nTempChars; i < nTempCharsLast; i++)
          printf("\b");

     va_end(Args);

     /* Needed for Sun. */
     fflush(stdout);

     nTempCharsLast = nTempChars;

     return;
}

/******************************+++********************************/
void OutTime(void)
/*****************************************************************/
/* Purpose:    Output time and date.                             */
/*                                                               */
/* 1991.06.11: Created.                                          */
/* 1999.06.17: "\n" added to Output;                             */
/*             Check if licence has expired.                     */
/* 2000.02.15: 1900 added to year for Y2K compatibility, and     */
/*             format for year changed from "%02d" to %4d".      */
/* 2011.01.19: "License" alternative spelling removed            */
/*****************************************************************/
{
     time_t    Tim;    /* Time.               */
     struct tm *TimSt;  /* Structure for time. */

     time(&Tim);
     TimSt = localtime(&Tim);
     Output("%02d:%02d:%02d on %4d/%02d/%02d\n",
               TimSt->tm_hour, TimSt->tm_min, TimSt->tm_sec,
               TimSt->tm_year + 1900, TimSt->tm_mon + 1,
               TimSt->tm_mday);

     if (TimSt->tm_year > YEAR_LIC)
          Fatal("Licence expired!");
}
