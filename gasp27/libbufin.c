/*****************************************************************/
/*   LOW-LEVEL ROUTINES FOR BUFFERED INPUT.                      */
/*                                                               */
/*   Copyright (c) William J. Welch 1991.                        */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "implem.h"
#include "lib.h"

static char Buffer[INPUT_COLS + 2];   /* Includes '\n' and NULL. */

/*******************************+++*******************************/
string BufReadPrompt(FILE *InpFile, const string Prompt,
          boolean Echo)
/*****************************************************************/
/*   Purpose:  Read a buffer (line) from a file with optional    */
/*             prompting and echoing.                            */
/*             Comments from '#' to the end of the line are      */
/*             deleted.                                          */
/*                                                               */
/*   Args:     InpFile   Input file.                             */
/*             Prompt    The prompt for input.                   */
/*             Echo      The line will be echoed [via Output()]  */
/*                       if Echo == YES.                         */
/*                                                               */
/*   Returns:  A pointer to the line or NULL if end of file.     */
/*                                                               */
/*   Version:  1995 May 18                                       */
/*****************************************************************/
{
     char *NewEnd;

     if (Prompt != NULL)
     {
          Output(Prompt);
          fflush(stdout);
     }

     if (fgets(Buffer, INPUT_COLS + 2, InpFile) == NULL)
           return NULL;

     else
     {
          if (Echo && InpFile != stdin)
               /* Echo to stdout and the log file. */
               Output(Buffer);
          else if (Echo)
               /* Only echo to the log file. */
               FileOutput(GetLogFile(), Buffer);

          /* Delete comments starting with '#'. */
          if ( (NewEnd = strchr(Buffer, '#')) != NULL)
               if (NewEnd == Buffer)
                    *NewEnd = NULL;
               else
               {
                    NewEnd[0] = '\n';
                    NewEnd[1] = NULL;
               }

          return Buffer;
     }
}

/*******************************+++*******************************/
/*                                                               */
/*   string    BufToken(const string Delim, string *BufPtr,      */
/*                  char *TermChar)                              */
/*                                                               */
/*   Purpose:  Separate the next token from a buffer (line).     */
/*                                                               */
/*   Args:     Delim     Token delimiters (in addition to white  */
/*                       space and comma).   NULL if none.       */
/*             BufPtr    On input, *BufPtr is the current buffer */
/*                       string.  On output, *BufPtr is the      */
/*                       string starting immediately after the   */
/*                       token found.                            */
/*             TermChar  On exit, if TermChar != NULL, *TermChar */
/*                       is the delimiter found.                 */
/*                                                               */
/*   Returns:  A pointer to the beginning of the token, which    */
/*             might be "".                                      */
/*                                                               */
/*   Version:  1991 June 12                                      */
/*                                                               */
/*****************************************************************/

string BufToken(const string Delim, string *BufPtr, char *TermChar)
{
     char      LocalTermChar;
     string    Buf;
     string    Token;

     Buf = *BufPtr;

     /* Skip white space. */
     while (isspace(*Buf))
          Buf++;

     /* Token now starts. */
     Token = Buf;

     /* Token is all characters until next */
     /* white space, comma, or delimiter.  */
     while (*Buf != NULL && !isspace(*Buf) && *Buf != ',' &&
               (Delim == NULL || strchr(Delim, *Buf) == NULL) )
          Buf++;

     LocalTermChar = *Buf;

     /* Terminate the token and reposition Buf to point */
     /* immediately after the terminating character.    */
     if (*Buf != NULL)
     {
          *Buf = NULL;
          Buf++;

          /* Skip any terminating white space,     */
          /* and update the terminating character. */
          if (isspace(LocalTermChar))
          {
               while(isspace(*Buf))
                    Buf++;
               if (*Buf == ',' || (*Buf != NULL && Delim != NULL &&
                         strchr(Delim, *Buf) != NULL) )
               {
                    LocalTermChar = *Buf;
                    Buf++;
               }
          }
     }

     *BufPtr = Buf;
     if (TermChar != NULL)
          *TermChar = LocalTermChar;

     return Token;
}

/*******************************+++*******************************/
/*                                                               */
/*   string    BufPromptToken(const string AppPrompt,            */
/*                  const string TokPrompt, const string Assign) */
/*                                                               */
/*   Purpose:  Prompt for a token.                               */
/*                                                               */
/*   Args:     AppPrompt ApplicationPrompt, e.g., "APP>".        */
/*             TokPrompt Token prompt, e.g., "n".                */
/*             Assign    Assignment string, e.g., " = ".         */
/*                                                               */
/*   Returns:  A pointer to the beginning of the token, which    */
/*             might be "".                                      */
/*                                                               */
/*   Version:  1992 April 16                                     */
/*                                                               */
/*****************************************************************/

string BufPromptToken(const string AppPrompt,
          const string TokPrompt, const string Assign)
{
     string    Buf, Prompt;

     Prompt = StrPaste(3, AppPrompt, TokPrompt, Assign);
     Buf = BufReadPrompt(stdin, Prompt, YES);
     AllocFree(Prompt);

     return BufTok(&Buf);
}

/*******************************+++*******************************/
/*                                                               */
/*   string    BufForceToken(FILE *InpFile, const string Prompt, */
/*                  boolean Echo, const string Delim,            */
/*                  string *BufPtr, char *TermChar)              */
/*                                                               */
/*   Purpose:  Separate the next token from a buffer (line).     */
/*             If the buffer is empty, then the buffer is        */
/*             refilled until a token is found (i.e. a token is  */
/*             forced).                                          */
/*             Note that a buffer consisting of, for example,    */
/*             just white space and a comma will lead to the     */
/*             return of a token of length zero.                 */
/*                                                               */
/*   Args:     InpFile   Input file.                             */
/*             Prompt    The prompt for input.                   */
/*             Echo      The line will be echoed [via Output()]  */
/*                       if Echo == YES.                         */
/*             Delim     Token delimiters (in addition to white  */
/*                       space and comma).  NULL for none.       */
/*             BufPtr    On input, *BufPtr is the current buffer */
/*                       string.  On output, *BufPtr is the      */
/*                       string starting immediately after the   */
/*                       token found.                            */
/*             TermChar  On exit, if TermChar != NULL, *TermChar */
/*                       is the delimiter found.                 */
/*                                                               */
/*   Returns:  A pointer to the beginning of the token or NULL   */
/*             if end of file.                                   */
/*                                                               */
/*   Version:  1991 June 12                                      */
/*                                                               */
/*****************************************************************/

string BufForceToken(FILE *InpFile, const string Prompt,
          boolean Echo, const string Delim, string *BufPtr,
          char *TermChar)
{
     char     LocalTermChar;
     string   Buf, Token;

     Buf = *BufPtr;

     if (Buf == NULL)
          Buf = BufReadPrompt(InpFile, Prompt, Echo);

     if (Buf != NULL)
          Token = BufToken(Delim, &Buf, &LocalTermChar);

     while (Buf != NULL && *Token == NULL &&
               (LocalTermChar == NULL || isspace(LocalTermChar)))
     {
          Buf = BufReadPrompt(InpFile, Prompt, Echo);
          if (Buf != NULL)
               Token = BufToken(Delim, &Buf, &LocalTermChar);
     }

     *BufPtr = Buf;

     if (Buf != NULL)
     {
          if (TermChar != NULL)
               *TermChar = LocalTermChar;
          return Token;
     }
     else
          return NULL;
}
