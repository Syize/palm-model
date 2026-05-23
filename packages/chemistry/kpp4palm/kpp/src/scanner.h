#ifndef __SCANNER_H
#define __SCANNER_H

void CmdStoicmat( char *cmd );
void CheckAll();
void LookAtAll();
void TransportAll();
void DefineInitializeNbr( char *cmd );
void DefineXGrid( char *cmd );
void DefineYGrid( char *cmd );
void DefineZGrid( char *cmd );
void SparseData( char *cmd );
void AddUseFile( char *fname );
int ParseEquationFile( char * filename );
#endif