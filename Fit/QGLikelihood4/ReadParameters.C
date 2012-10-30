#include <stdio.h>
#include <string>

#ifndef READ_
#define READ_

#define MAX_STR_LENGTH 1023

class Read
{
public:
	//costruttore
	Read(){};
	//distruttore
	~Read(){};
	//funzioni
	char* ReadParameterFromFile(const char*fileName,const char*parName);
	float ConvertToFloat(const char*Par);
private:
	//strutture dati
};

//==============================
char* Read::ReadParameterFromFile(const char*fileName,const char * parName)
{
FILE *fr=fopen(fileName,"r");
if(fr==NULL) return NULL;
char *R=new char[MAX_STR_LENGTH]; //must be a pointer - should survive outside
char P[MAX_STR_LENGTH];
char S[MAX_STR_LENGTH];

//leggi una linea su S
int STATUS=1;
while(STATUS!=EOF)
{
	char c='\0';
	int i=0;
	while ( (STATUS=fscanf(fr,"%c",&c))!=EOF ){
		if(c=='\0') break;
		if(c=='\n') break;
		S[i]=c;i++;
		}
	S[i]='\0';
	i=0;
	while(S[i]!='=' && S[i]!='\0'){P[i]=S[i];i++;}
	P[i]='\0';int j=0;
	if(S[i]=='=')i++;
	while(S[i]!='\0'){R[j]=S[i];i++;j++;}
	R[j]='\0';
	bool isPar=true;
	for(i=0;;i++)
		{
		if(P[i]=='\0' && parName[i]=='\0')break;
		if(P[i]=='\0' || parName[i]=='\0'){isPar=false; break;} // but not together
		if(parName[i]!=P[i]){isPar=false;break;}
		}
	
	if(isPar){fclose(fr);return R;}
}
fclose(fr);
return NULL;
}

float Read::ConvertToFloat(const char*Par)
{
}
#endif
