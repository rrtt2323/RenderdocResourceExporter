#include "CSVFile.h"
#include "FileStream.h"
#include "checked.h"

inline void  StringSplit(std::string& s, const char* delim, std::vector< std::string >& ret)
{
	size_t last = 0;
	size_t index = s.find_first_of(delim, last);
	while (index != std::string::npos)
	{
		ret.push_back(s.substr(last, index - last));
		last = index + 1;
		index = s.find_first_of(delim, last);
	}
	if (index - last > 0)
	{
		ret.push_back(s.substr(last, index - last));
	}
}

void TrimSpace(std::string& strParam)
{
	if (strParam.empty())
		return;

	const char* start = strParam.c_str();
	int eraseCount = 0;
	while (*start && isspace(*start))
	{
		start++;
		eraseCount++;
	}
	if (eraseCount)
		strParam = strParam.erase(0, eraseCount);

	eraseCount = 0;
	const char* end = strParam.c_str() + strParam.size() - 1;
	while (end > start && isspace(*end))
	{
		end--;
		eraseCount++;
	}
	if (eraseCount)
		strParam = strParam.erase(strParam.size() - eraseCount, eraseCount);
}

int	 ConvertUTF8ToUnicode(const char* szUtf8, std::wstring& strOut)
{
	std::string utf8Str = szUtf8;
//#ifdef _DEBUG
//	std::string::iterator end_it = core_utf8::find_invalid(utf8Str.begin(), utf8Str.end());
//	if (end_it != utf8Str.end()) {
//		assert(false);
//	}
//#endif 
	core_utf8::utf8to16(utf8Str.begin(), utf8Str.end(), std::back_inserter(strOut));
	return strOut.size(); 
}

int  ConvertUnicodeToUTF8(const wchar_t* srcUniStr, std::string& strOut)
{
	std::wstring utf16Str = srcUniStr;

	core_utf8::utf16to8(utf16Str.begin(), utf16Str.end(), std::back_inserter(strOut));
	return strOut.size(); 
}
 
void TrimQuotes(std::string& strParam)
{
	if(strParam.size()>2)
	{
		if(strParam[0] == '"' && strParam[strParam.size()-1]=='"')
		{
			strParam = strParam.substr(1,strParam.size()-2);
		}
	}
}

void CCSVFileRow::ParseRow(std::string& strLine)
{
	StringSplit(strLine,",",m_vecData);
	for(int iCol = 0;iCol<m_vecData.size();iCol++)
	{
		TrimQuotes(m_vecData[iCol]);
	}
}

void CCSVFileRow::Reserve(int iElemCount)
{
	m_vecData.reserve(iElemCount);
	m_vecData.resize(iElemCount);
}


CCSVFile::CCSVFile(const char* szCSVFile)
{
	OpenCSVFile(szCSVFile);
}

CCSVFile::~CCSVFile()
{
	m_mapFields.clear();
	m_vecRows.clear(); 
}

bool	CCSVFile::OpenCSVFile(const char* file)
{
	bool res = false;
	//
	m_mapFields.clear();
	m_vecRows.clear();

 
	CFileStream* pFileStream = CreateFileStream();
 
	if(pFileStream)
	{
		if(pFileStream->Open(file,CFileStream::EFILEACCESS_READ) == CFileStream::EFILEOP_OK)
		{
			std::string strBuffer;
			
			pFileStream->ReadLine(strBuffer);

			std::vector<std::string> vecFiledStrings;
			StringSplit(strBuffer,",",vecFiledStrings);

			for(int iField = 0;iField<vecFiledStrings.size();iField++)
			{
				TrimQuotes(vecFiledStrings[iField]);
				TrimSpace(vecFiledStrings[iField]);
				m_mapFields[vecFiledStrings[iField]] = iField;
			}
			while(1)
			{
				if(pFileStream->GetStatus() == CFileStream::EFILEOP_EOS) 
					break;
				
				pFileStream->ReadLine(strBuffer);
				if(strBuffer.empty())
					continue;
				CCSVFileRow row;
				m_vecRows.push_back(row);
				m_vecRows.back().ParseRow(strBuffer);
			}
			pFileStream->Close();

			res = true;
		}
 
 
		DestroyFileStream(pFileStream);
 
	}
	else
		res = false;
	return res;
}

void	CCSVFile::WriteCSVFile(const char* szFileName)
{
#ifndef USE_PACK_FILE
	CFileStream* pFileStream = CreateFileStream();

	if(pFileStream)
	{
		if(pFileStream->Open(szFileName,CFileStream::EFILEACCESS_WRITE) == CFileStream::EFILEOP_OK)
		{
			char buf[2048]={0};
			unsigned long dwWrotten = 0;
			unsigned long dwLineSize = 0;

			std::vector<std::string>  vecFields;
			vecFields.reserve(m_mapFields.size());
			vecFields.resize(m_mapFields.size());
			for(std::map<std::string,int>::iterator itr = m_mapFields.begin();
				itr != m_mapFields.end();itr++)
			{
				vecFields[itr->second] = itr->first;
			}


			for(int iField = 0;iField < vecFields.size() - 1;iField++)
			{
				strcat_s(buf + dwLineSize ,2048 - dwLineSize,vecFields[iField].c_str());
				dwLineSize += vecFields[iField].size();
				strcat_s(buf + dwLineSize ,2048 - dwLineSize,",");
				dwLineSize += 1;
			}

			strcat_s(buf + dwLineSize ,2048 - dwLineSize,vecFields.back().c_str());
			dwLineSize += vecFields.back().size();
			strcat_s(buf + dwLineSize ,2048 - dwLineSize,"\n");
			dwLineSize += 1;

			pFileStream->Write(dwLineSize, buf,&dwWrotten);


			for(int iRow =0;iRow <m_vecRows.size();iRow++)
			{
				CCSVFileRow& row = m_vecRows[iRow]; 
				dwWrotten = 0;
				dwLineSize = 0;
				memset(buf,0,sizeof(buf));
				for(int iCol = 0;iCol < row.Size()-1;iCol++)
				{
					strcat_s(buf + dwLineSize ,2048 - dwLineSize,row[iCol].c_str());
					dwLineSize += row[iCol].size();
					strcat_s(buf + dwLineSize ,2048 - dwLineSize,",");
					dwLineSize += 1;
				}
				strcat_s(buf + dwLineSize ,2048 - dwLineSize,row[row.Size() -1].c_str());
				dwLineSize += row[row.Size() -1].size();
				strcat_s(buf + dwLineSize ,2048 - dwLineSize,"\n");
				dwLineSize += 1;

				pFileStream->Write(dwLineSize, buf,&dwWrotten);
			}
		}

		DestroyFileStream(pFileStream);
	}
#endif
}

int	CCSVFile::GetRowNum()
{
	return m_vecRows.size();
}

bool	CCSVFile::GetCellValue(const char* fieldName,int rowIndex,std::wstring& resValue)
{ 
	int colIndex = m_mapFields[fieldName];
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	if(m_vecRows[rowIndex][colIndex] == "null")
		resValue.clear();
	else
		ConvertUTF8ToUnicode(m_vecRows[rowIndex][colIndex].c_str(),resValue);

	return true;
}
bool	CCSVFile::GetCellValue(int colIndex,int rowIndex,std::wstring& resValue)
{
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false; 
	if(m_vecRows[rowIndex][colIndex] == "null")
		resValue.clear();
	else
		ConvertUTF8ToUnicode(m_vecRows[rowIndex][colIndex].c_str(),resValue);

	return true;
}

bool	CCSVFile::GetCellValue(const char* fieldName,int rowIndex,std::string& resValue)
{
 
	int colIndex = m_mapFields[fieldName];
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	resValue = m_vecRows[rowIndex][colIndex];
	if(resValue == "null")
		resValue.clear();
	return true;
}
bool	CCSVFile::GetCellValue(int colIndex,int rowIndex,std::string& resValue)
{
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	resValue = m_vecRows[rowIndex][colIndex];
	if(resValue == "null")
		resValue.clear();
	return true;
}
bool	CCSVFile::GetCellValue(const char* fieldName,int rowIndex,int& resValue)
{	
 
	int colIndex = m_mapFields[fieldName];
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	if(m_vecRows[rowIndex][colIndex] == "null")
		resValue = 0;
	else
		resValue = atoi(m_vecRows[rowIndex][colIndex].c_str());
	 
	return true;
}
bool	CCSVFile::GetCellValue(int colIndex,int rowIndex,int& resValue)
{	
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	if(m_vecRows[rowIndex][colIndex] == "null")
		resValue = 0;
	else
		resValue = atoi(m_vecRows[rowIndex][colIndex].c_str());
		return true;
}

bool	CCSVFile::GetCellValue(const char* fieldName,int rowIndex,unsigned char& resValue)
{	 
	int colIndex = m_mapFields[fieldName];
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	if(m_vecRows[rowIndex][colIndex] == "null")
		resValue = 0;
	else
		resValue = atoi(m_vecRows[rowIndex][colIndex].c_str());

	return true;
}
bool	CCSVFile::GetCellValue(int colIndex,int rowIndex,unsigned char& resValue)
{	
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	if(m_vecRows[rowIndex][colIndex] == "null")
		resValue = 0;
	else
		resValue = atoi(m_vecRows[rowIndex][colIndex].c_str());
	return true;
}



bool	CCSVFile::GetCellValue(const char* fieldName,int rowIndex,float& resValue)
{ 
	int colIndex = m_mapFields[fieldName];
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	if(m_vecRows[rowIndex][colIndex] == "null")
		resValue = 0.0f;
	else
		resValue = atof(m_vecRows[rowIndex][colIndex].c_str());
	return true;
}
bool	CCSVFile::GetCellValue(int colIndex,int rowIndex,float& resValue)
{
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	if(m_vecRows[rowIndex][colIndex] == "null")
		resValue = 0.0f;
	else
		resValue = atof(m_vecRows[rowIndex][colIndex].c_str());
	return true;
}

bool	CCSVFile::GetCellValue(const char* fieldName,int rowIndex,bool& resValue)
{ 
	int colIndex = m_mapFields[fieldName];
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	if(m_vecRows[rowIndex][colIndex] == "null" || m_vecRows[rowIndex][colIndex] == "FALSE" ||m_vecRows[rowIndex][colIndex] == "false")
		resValue = false;
	else
		resValue = true;
	return true;
}
bool	CCSVFile::GetCellValue(int colIndex,int rowIndex,bool& resValue)
{
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	if(m_vecRows[rowIndex][colIndex] == "null" || m_vecRows[rowIndex][colIndex] == "FALSE" ||m_vecRows[rowIndex][colIndex] == "false")
		resValue = false;
	else
		resValue = true;
	return true;
}

bool	CCSVFile::SetCellValue(const char* fieldName,int rowIndex,std::wstring& resValue)
{ 
	int colIndex = m_mapFields[fieldName];
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	if(resValue.empty())
		m_vecRows[rowIndex][colIndex] = "null";
	else
	{
		m_vecRows[rowIndex][colIndex].clear();
		ConvertUnicodeToUTF8(resValue.c_str(),m_vecRows[rowIndex][colIndex]);
	}

	return true;
}
bool	CCSVFile::SetCellValue(int colIndex,int rowIndex,std::wstring& resValue)
{
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false; 
	if(resValue.empty())
		m_vecRows[rowIndex][colIndex] = "null";
	else
	{
		m_vecRows[rowIndex][colIndex].clear();
		ConvertUnicodeToUTF8(resValue.c_str(),m_vecRows[rowIndex][colIndex]);
	}
	return true;
}

bool	CCSVFile::SetCellValue(const char* fieldName,int rowIndex,std::string& resValue)
{ 
	int colIndex = m_mapFields[fieldName];
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	if(resValue.empty())
		m_vecRows[rowIndex][colIndex] = "null";
	else
		m_vecRows[rowIndex][colIndex] = resValue; 
	return true;
}
bool	CCSVFile::SetCellValue(int colIndex,int rowIndex,std::string& resValue)
{
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	if(resValue.empty())
		m_vecRows[rowIndex][colIndex] = "null";
	else
		m_vecRows[rowIndex][colIndex] = resValue; 
	return true;
}
bool	CCSVFile::SetCellValue(const char* fieldName,int rowIndex,int resValue)
{	 
	int colIndex = m_mapFields[fieldName];
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false; 
	char str[32] = {0};
	sprintf_s(str,"%d",resValue);
	m_vecRows[rowIndex][colIndex] = str; 

	return true;
}
bool	CCSVFile::SetCellValue(int colIndex,int rowIndex,int resValue)
{	
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	char str[32] = {0};
	sprintf_s(str,"%d",resValue);
	m_vecRows[rowIndex][colIndex] = str; 
	return true;
}

bool	CCSVFile::SetCellValue(const char* fieldName,int rowIndex,unsigned char resValue)
{	 
	int colIndex = m_mapFields[fieldName];
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	char str[32] = {0};
	sprintf_s(str,"%d",resValue);
	m_vecRows[rowIndex][colIndex] = str; 

	return true;
}
bool	CCSVFile::SetCellValue(int colIndex,int rowIndex,unsigned char resValue)
{	
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	char str[32] = {0};
	sprintf_s(str,"%d",resValue);
	m_vecRows[rowIndex][colIndex] = str; 
	return true;
}



bool	CCSVFile::SetCellValue(const char* fieldName,int rowIndex,float resValue)
{ 
	int colIndex = m_mapFields[fieldName];
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	char str[32] = {0};
	sprintf_s(str,"%f",resValue);
	m_vecRows[rowIndex][colIndex] = str; 
	return true;
}
bool	CCSVFile::SetCellValue(int colIndex,int rowIndex,float resValue)
{
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	char str[32] = {0};
	sprintf_s(str,"%f",resValue);
	m_vecRows[rowIndex][colIndex] = str; 
	return true;
}

bool	CCSVFile::SetCellValue(const char* fieldName,int rowIndex,bool resValue)
{ 
	int colIndex = m_mapFields[fieldName];
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false; 
	m_vecRows[rowIndex][colIndex] = resValue?"TRUE":"FALSE"; 
	return true;
}
bool	CCSVFile::SetCellValue(int colIndex,int rowIndex,bool resValue)
{
	if(rowIndex >= m_vecRows.size())
		return false;
	if(colIndex >= m_mapFields.size())
		return false;
	m_vecRows[rowIndex][colIndex] = resValue?"TRUE":"FALSE"; 
	return true;
}