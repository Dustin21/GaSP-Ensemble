/*****************************************************************/
/*   CALLED BY ACED OR GASP MAIN TO INITIALIZE AND EXECUTE       */
/*   EVENT LOOP                                                  */
/*                                                               */
/*   Copyright (c) William J. Welch 1994--95.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "define.h"
#include "implem.h"
#include "lib.h"
#include "model.h"
#include "alex.h"

string         FuncName;

extern Matrix  DbStatus;

/*******************************+++*******************************/
void Run(int argc, char *argv[], size_t nFns, Function *ImpFn,
          const string Banner, const string Prompt)
/*****************************************************************/
/* Purpose:    Run an application.                               */
/*                                                               */
/* 1996.09.01: DEFAULT_OUT becomes LOG_FILE.                     */
/* 1999.06.16: Only one "\n" after OutTime (OutTime now outputs  */
/*             one "\n".                                         */
/*****************************************************************/
{
     char      *DotPtr;
     FILE      *JobFile, *LogFile;
     string    LogFileName;

     if (argc > 3)
          Fatal("Usage: %s [JobFile [LogFile]]\n", argv[0]);

     /* Get name of log file. */
     if (argc == 3)
          LogFileName = argv[2];
     else if (argc == 2)
     {
          /* If argv[1] has an extension, remove it. */
          if ( (DotPtr = strrchr(argv[1], '.')) != NULL)
               *DotPtr = NULL;

          LogFileName = StrPaste(2, argv[1], ".out");

          /* Restore extension to argv[1]. */
          if (DotPtr != NULL)
               *DotPtr = '.';
     }
     else
          LogFileName = LOG_FILE;

     /* Open log file. */
     LogFile = FileOpen(LogFileName, "a");

     if (argc == 2)
          AllocFree(LogFileName);

     SetLogFile(LogFile);

     Output(Banner);
     OutTime();
     Output("\n");

     /* Initialize database scalars and matrices. */
     DbScalarInit();
     DbMatInit();
     DbColTemplInit();

     /* Execute command line input file (if present). */
     if (argc >= 2 && (JobFile = FileOpen(argv[1], "r")) != NULL)
          EventLoop(argv[1], JobFile, nFns, ImpFn, Prompt);

     /* Return control to stdin (unless "Stop" encountered). */
     EventLoop("standard input", stdin, nFns, ImpFn, Prompt);

     exit(0);
}

/*******************************+++*******************************/
void EventLoop(string FileName, FILE *InpFile, size_t nFns,
          Function *ImpFn, const string Prompt)
/*****************************************************************/
/*   Purpose:  Execute events.                                   */
/*                                                               */
/*   Version:  1996.03.18                                        */
/*****************************************************************/
{
     char           Operator;
     DbMatrix       *D;
     int            ErrNum;
     size_t         FnIndex, Index;
     string         Buf, Token1, Token2, Option;
     const string   *Check;

     Output("Reading commands from %s.\n", FileName);

     /* Read commands until:                              */
     /* an error is encountered when reading from a file, */
     /* or EOF is encountered, or                         */
     /* "Stop" or "Quit" is encountered.                  */
     ErrNum = OK;
     while ( (ErrNum == OK || ErrNum == INCOMPAT_ERR
               || InpFile == stdin) &&
               (Buf = BufReadPrompt(InpFile, Prompt, YES)) != NULL)
     {
          /* Skip blank lines. */
          while (*Buf != NULL && isspace(*Buf))
               Buf++;
          if (*Buf == NULL)
               continue;

          Token1 = StrDup(BufToken("=<>", &Buf, &Operator));

          for (FnIndex = 0; FnIndex < nFns; FnIndex++)
               if (stricmp(Token1, ImpFn[FnIndex].FuncName) == 0)
                    break;

          if (FnIndex < nFns)
          {
               /* An implemented function (verb). */

               /* Check inputs. */

               FuncName = ImpFn[FnIndex].FuncName;

               Check = ImpFn[FnIndex].Check;
               ErrNum = OK;
               while (*Check != NULL)
                    if (DbCheck(*(Check++)) != OK)
                         ErrNum = !OK;
               Output("\n");
               MatWriteBlock(&DbStatus, NO, stdout);
               Output("\n");
               FuncName = NULL;

               if (ErrNum == OK)
               {
                    ErrNum = ExecuteFunc(*ImpFn[FnIndex].Func);
                    DbOutputMatStatus();
               }

               MatFree(&DbStatus);
          }

          else if (stricmp(Token1, "System") == 0)
               /* Pass the rest of the line to the */
               /* operating system.                */
               ErrNum = system(Buf);

          else if (stricmp(Token1, "Memory") == 0)
               Output("%u contiguous reals available.\n",
                         AllocMax(sizeof(real)));

          else if (stricmp(Token1, "Quit") == 0 ||
                    stricmp(Token1, "Stop") == 0)
          {
               Output("Goodbye.\n\n");
               exit(0);
          }

          else
          {
               /* Remaining commands need Token2. */
               Token2 = StrDup(BufTok(&Buf));

               if (*Token1 == NULL || *Token2 == NULL ||
                         strchr("=<>", Operator) == NULL)
               {
                    Error("Cannot parse this line.\n");
                    ErrNum = INPUT_ERR;
               }

               else if (Operator == '=')
               {
                    /* Scalar assignment. */
                    if ( (Index = DbScalIndex(Token1, YES))
                             == INDEX_ERR)
                         ErrNum = INPUT_ERR;
                    else
                         ErrNum = DbScalParse(Index, Token2);
               }

               else if (Operator == '<')
               {
                    /* Matrix input. */
                    if ( (D = DbMatFind(Token1, NO)) == NULL)
                         ErrNum = INPUT_ERR;
                    else
                         ErrNum = DbMatRead(D, Token2);
               }

               else if (Operator == '>')
               {
                    /* Matrix output. */
                    Option = BufTok(&Buf);   /* Blocked? */
                    if ( (D = DbMatFind(Token1, NO)) == NULL)
                         ErrNum = INPUT_ERR;
                    else
                         ErrNum = DbMatWrite(D, Token2, Option);
               }
               else
               {
                    Error("Cannot parse this line.\n");
                    ErrNum = INPUT_ERR;
               }

               AllocFree(Token2);
          }
          AllocFree(Token1);

     }

     if (Buf == NULL)
          Output("End of file %s encountered.\n", FileName);
     else
          Output("Terminating input from file %s.\n", FileName);

     return;
}

/*******************************+++*******************************/
int ExecuteFunc(int (*Func)(void))
/*****************************************************************/
/*   Purpose:  Execute a function (verb).                        */
/*                                                               */
/*   Version:  1996.03.18                                        */
/*****************************************************************/
{
     int       ErrNum;
     time_t    Finish, Start;

     time(&Start);
     ErrNum = Func();
     time(&Finish);

     Output("Seconds: %g\n\n",  difftime(Finish, Start));

     ErrorMatOut();

     return ErrNum;
}
