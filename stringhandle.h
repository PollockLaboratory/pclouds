#ifndef stringhandle_h__
#define stringhandle_h__

#include <string>

void strright(int i, char *source, char *target);

int straddright(char *source, char *target, char *add);

void strleft(int i, char* source, char* target);

int strcompare(int i, char *source, char *target);

int strreplaceright(char *source, char *target, char *replace);

void getsubstring(char* pchSource, char *pchTarget, int nStart, int nEnd);

int stringcompare(char *pchSource, char *pchTarget);

int stringcompare(const char *pchSource, const char *pchTarget);

int stringinitial(char *pchSource, int iLength);

int stringcopy(char *pchSource, char *pchTarget, int iLength);

void getreversecomplement(const char* pchSource, char* pchTarget);

int stringtonumber(const std::string& m_strTemp);

unsigned int stringtolargenumber(const std::string& m_strTemp);

void stringtoarray(const std::string& src, char* target);

bool isonessr(char *pchPattern);

bool istwossr(char *pchPattern);

bool isthreessr(char *pchPattern);

bool isfourssr(char *pchPattern);

bool issegmentvalid(char *pchPattern);

#endif // stringhandle_h__
