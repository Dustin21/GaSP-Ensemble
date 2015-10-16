/*****************************************************************/
/*                                                               */
/*   FILE-MANIPULATION ROUTINES                                  */
/*                                                               */
/*   Copyright (c) William J. Welch 1991.                        */
/*   All rights reserved.                                        */
/*                                                               */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "implem.h"
#include "lib.h"

/*******************************+++*******************************/
/*                                                               */
/*   FILE      *FileOpen(const string FileName,                  */
/*                  const string Mode)                           */
/*                                                               */
/*   Purpose:  Version of fopen() that:                          */
/*             (1) automatically backs up the file if it already */
/*             exists and is being opening for writing;          */
/*             (2) issues an error message if the file cannot be */
/*             opened.                                           */
/*                                                               */
/*   Args:     FileName  File name.                              */
/*             Mode      Input-output mode.                      */
/*                                                               */
/*   Returns:  Pointer to the file or NULL if unsuccessful.      */
/*                                                               */
/*   Version:  1991 July 8                                       */
/*                                                               */
/*****************************************************************/

FILE *FileOpen(const string FileName, const string Mode)
{
     FILE *FilePtr;

     if (stricmp(Mode, "w") == 0 &&
              (FilePtr = fopen(FileName, "r")) != NULL)
     {
          /* Writing to an existing file: back it up. */
          FileBackup(FileName, FilePtr);
     }

     if ( (FilePtr = fopen(FileName, Mode)) == NULL)
          Error("File %s cannot be opened.\n", FileName);

     return FilePtr;
}

/*******************************+++*******************************/
void FileBackup(const string FileName, FILE *Source)
/*****************************************************************/
/*   Purpose:  Copy a file to a backup file.                     */
/*             The backup name is FileName with the extension    */
/*             replaced by ".bak" or, if there is no extension,  */
/*             with ".bak" appended.                             */
/*                                                               */
/*   Version:  1995 May 2                                        */
/*****************************************************************/
{
     char BackupName[MAXTOK+1];
     char *Dot;
     FILE *Target;
     int  c;

     strcpy(BackupName, FileName);

     if ( (Dot = strrchr(BackupName, '.')) != NULL)
          strcpy(Dot + 1, "bak");
     else
          strcat(BackupName, ".bak");

     Output("Backing up %s to %s.\n", FileName, BackupName);

     if ( (Target = fopen(BackupName, "w")) == NULL)
          Error("Backup file %s cannot be opened.\n", BackupName);
     else
          while ( (c = getc(Source)) != EOF)
               putc(c, Target);

     fclose(Source);
     fclose(Target);
}
