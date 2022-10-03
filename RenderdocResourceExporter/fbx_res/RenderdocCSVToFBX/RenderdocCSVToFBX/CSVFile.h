#ifndef CSVFILE_H_INCLUDE
#define CSVFILE_H_INCLUDE

#include <vector>
#include <string>
#include <map>
 
class   CCSVFileRow
{
public:
	CCSVFileRow()
	{

	}

	~CCSVFileRow()
	{
		m_vecData.clear();
	}
	std::string const& operator[](std::size_t index) const
	{
		return m_vecData[index];
	}
	std::string & operator[](std::size_t index) 
	{
		return m_vecData[index];
	}

	std::size_t Size() const
	{
		return m_vecData.size();
	}
	void ParseRow(std::string& strLine);

	void Reserve(int iElemCount);
private:
	std::vector<std::string>    m_vecData;
};

class   CCSVFile
{
public:
	CCSVFile(const char* szCSVFile);
	~CCSVFile();

	bool	OpenCSVFile(const char* file);
	void	WriteCSVFile(const char* szFileName);

	int		GetRowNum();
 
	bool	GetCellValue(const char* fieldName,int rowIndex,std::wstring& resValue);
	bool	GetCellValue(int colIndex,int rowIndex,std::wstring& resValue);
	bool	GetCellValue(const char* fieldName,int rowIndex,std::string& resValue);
	bool	GetCellValue(int colIndex,int rowIndex,std::string& resValue);
	bool	GetCellValue(const char* fieldName,int rowIndex,int& resValue);
	bool	GetCellValue(int colIndex,int rowIndex,int& resValue);
	bool	GetCellValue(const char* fieldName,int rowIndex,unsigned char& resValue);
	bool	GetCellValue(int colIndex,int rowIndex,unsigned char& resValue);
	bool	GetCellValue(const char* fieldName,int rowIndex,float& resValue);
	bool	GetCellValue(int colIndex,int rowIndex,float& resValue);
	bool	GetCellValue(const char* fieldName,int rowIndex,bool& resValue);
	bool	GetCellValue(int colIndex,int rowIndex,bool& resValue);


	bool	SetCellValue(const char* fieldName,int rowIndex,std::wstring& resValue);
	bool	SetCellValue(int colIndex,int rowIndex,std::wstring& resValue);
	bool	SetCellValue(const char* fieldName,int rowIndex,std::string& resValue);
	bool	SetCellValue(int colIndex,int rowIndex,std::string& resValue);
	bool	SetCellValue(const char* fieldName,int rowIndex,int resValue);
	bool	SetCellValue(int colIndex,int rowIndex,int resValue);
	bool	SetCellValue(const char* fieldName,int rowIndex,unsigned char resValue);
	bool	SetCellValue(int colIndex,int rowIndex,unsigned char resValue);
	bool	SetCellValue(const char* fieldName,int rowIndex,float resValue);
	bool	SetCellValue(int colIndex,int rowIndex,float resValue);
	bool	SetCellValue(const char* fieldName,int rowIndex,bool resValue);
	bool	SetCellValue(int colIndex,int rowIndex,bool resValue);

	std::map<std::string,int>&	GetFieldsMap()
	{
		return m_mapFields;
	}

	std::vector<CCSVFileRow>&	GetRowsVector()
	{
		return m_vecRows; 
	}
private:
	std::map<std::string,int>	m_mapFields;
	std::vector<CCSVFileRow> m_vecRows; 
};
 


#endif