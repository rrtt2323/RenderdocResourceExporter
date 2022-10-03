#include "FileStream.h"
#include "assert.h"

std::string FindExtensionFromFileName(const std::string& name)
{
	int pos = name.rfind('.');
	if(pos != std::string::npos)
	{
		pos++;
		return name.substr(pos,name.size()- pos);
	}
	else
		return std::string("");	 
}


std::string FindShortNameFromFileName(const std::string& name,bool bIncludeExt)
{ 
	int pos = name.rfind('\\');
 
 
	if(pos != std::string::npos)
	{
		pos++;
		std::string fileName = name.substr(pos,name.size()- pos);
		if(bIncludeExt)
			return fileName;
		else
		{
			pos = fileName.rfind('.');
			if(pos != std::string::npos)
			{
				return fileName.substr(0,pos);
			}
			else
				return std::string("");	 
		}
	}
	else
		return std::string("");	 
}


 

//-----------------------------------------------------------------------------
bool FileDelete(const char * name)
{
   if(!name || (strlen(name) >= 260))
      return(false);
   return(remove(name)== 0);
}

 
//-----------------------------------------------------------------------------
// Self-explanatory.
//-----------------------------------------------------------------------------
bool CFileStream::HasCapability(EFILEACCESS_MODE_FLAG cap) const
{
    return (0 != (unsigned long(cap) & m_eFileAccesssMode));
}

//-----------------------------------------------------------------------------
// Constructors & Destructor
//-----------------------------------------------------------------------------
CFileStream::CFileStream()
: m_currentStatus(EFILEOP_CLOSED), m_eFileAccesssMode(EFILEACCESS_READ) 
{
	m_pMemBuffer = NULL;
	m_iMemPosition = m_dwMemSize = 0;
	m_bStream = false;
}


//-----------------------------------------------------------------------------
// Sets and returns the currentStatus to status.
//-----------------------------------------------------------------------------
CFileStream::EFILEOP_STATUS CFileStream::SetStatus(CFileStream::EFILEOP_STATUS status)
{
    return m_currentStatus = status;
}

bool	CFileStream::ReadInt(int&	value)
{
	int tmp = 0;
	if(Read(sizeof(int),(char*)&tmp) == EFILEOP_OK)
	{
		value = tmp;
		return true;
	}
	else
		return false;
}
bool	CFileStream::WriteInt(int	value)
{
	if(Write(sizeof(int),(const char*)&value) == EFILEOP_OK)
		return true;
	else
		return false;
}
bool	CFileStream::ReadDWORD(unsigned long&  value)
{
	unsigned long tmp = 0;
	if(Read(sizeof(unsigned long),(char*)&tmp) == EFILEOP_OK)
	{
		value = tmp;
		return true;
	}
	else
		return false;
}
bool	CFileStream::WriteDWORD(unsigned long value)
{
	if(Write(sizeof(unsigned long),(const char*)&value) == EFILEOP_OK)
		return true;
	else
		return false;
}

bool	CFileStream::WriteWORD(unsigned short value)
{
	if(Write(sizeof(unsigned short),(const char*)&value) == EFILEOP_OK)
		return true;
	else
		return false;
}
bool	CFileStream::ReadWORD(unsigned short& value)
{
	unsigned short tmp = 0;
	if(Read(sizeof(unsigned short),(char*)&tmp) == EFILEOP_OK)
	{
		value = tmp;
		return true;
	}
	else
		return false;
}

bool	CFileStream::WriteShort(short value)
{
	if(Write(sizeof(short),(const char*)&value) == EFILEOP_OK)
		return true;
	else
		return false;
}
bool	CFileStream::ReadShort(short& value)
{
	short tmp = 0;
	if(Read(sizeof(short),(char*)&tmp) == EFILEOP_OK)
	{
		value = tmp;
		return true;
	}
	else
		return false;
}

bool	CFileStream::WriteBYTE(unsigned char value)
{
	if(Write(sizeof(unsigned char),(const char*)&value) == EFILEOP_OK)
		return true;
	else
		return false;
}
bool	CFileStream::ReadBYTE(unsigned char& value)
{
	unsigned char tmp = 0;
	if(Read(sizeof(unsigned char),(char*)&tmp) == EFILEOP_OK)
	{
		value = tmp;
		return true;
	}
	else
		return false;
}

bool	CFileStream::Readbool(bool& value)
{
	char tmp = 0;
	if(Read(sizeof(char),(char*)&tmp) == EFILEOP_OK)
	{
		value = tmp;
		return true;
	}
	else
		return false;
}
bool	CFileStream::Writebool(bool value)
{
	char bValue = value;
	if(Write(sizeof(char),(const char*)&bValue) == EFILEOP_OK)
		return true;
	else
		return false;
}

bool	CFileStream::ReadBool(bool& value)
{
	char tmp = 0;
	if(Read(sizeof(char),(char*)&tmp) == EFILEOP_OK)
	{
		value = tmp;
		return true;
	}
	else
		return false;
}
bool	CFileStream::WriteBool(bool value)
{
	char bValue = value;
	if(Write(sizeof(char),(const char*)&bValue) == EFILEOP_OK)
		return true;
	else
		return false;
} 

bool	CFileStream::ReadFloat(float& value)
{
	float tmp = 0;
	if(Read(sizeof(float),(char*)&tmp) == EFILEOP_OK)
	{
		value = tmp;
		return true;
	}
	else
		return false;
}
bool	CFileStream::WriteFloat(float value)
{
	if(Write(sizeof(float),(const char*)&value) == EFILEOP_OK)
		return true;
	else
		return false;
}

bool CFileStream::WriteString(const char *string, int maxLen)
{
   int len = string ? unsigned long(strlen(string)) : 0;
   if(len > maxLen)
		return false;

   if(!WriteWORD(unsigned short(len)))
	   return false;
   if(len)
      return Write(len, string);

   return true;
}

bool CFileStream::ReadString(char* buf,int bufSize)
{
   unsigned short len;
   if(!ReadWORD(len))
	   return false;
	if(len >= bufSize)
		return false;
	return Read(int(len), buf) == CFileStream::EFILEOP_OK;
}

bool	CFileStream::ReadString(std::string& str,bool sizeWord)
{
	char buf[256]={0};
	if(sizeWord)
	{
		if(ReadString(buf,256))
		{
			str = buf;
			return true;
		}
	}
	else
	{
	   int len;
	   if(!ReadInt(len))
		   return false;
		if(len >= 256)
			return false;
		if(Read(int(len), buf) == CFileStream::EFILEOP_OK)
		{
			str = buf;
			return true;
		}
	}
	return false;
}

bool	CFileStream::ReadLine(std::string& outStr)
{
    if(m_bStream) 
	{
		if(EFILEOP_OK != m_currentStatus)
			return false;
		int linesize = 0;
		char buf[1024]={0};
		unsigned long dwReaded = 0;
		Read(1,&buf[linesize],&dwReaded);
		while(buf[linesize]!= '\n' && m_currentStatus == EFILEOP_OK)
		{
			linesize++;
			assert(linesize <= 1024);
			Read(1,&buf[linesize],&dwReaded);
		}
		outStr = buf;
		return true;
	}
    else
    {
		int linesize = 0;
		while(m_pMemBuffer[m_iMemPosition + linesize]!= '\n' && m_pMemBuffer[m_iMemPosition + linesize]!= '\r' && m_iMemPosition + linesize < m_dwMemSize)
		{
			linesize++;
		}

		assert(linesize <= 2048);
		char buf[2048]={0};
		memcpy(buf,m_pMemBuffer + m_iMemPosition,linesize);
		m_iMemPosition += linesize;
		outStr = buf;
		if(m_iMemPosition == m_dwMemSize)
		{
			m_currentStatus = EFILEOP_EOS; 
		}
		else
		{
			m_iMemPosition++;  //move to /n next pos
		}
		return true;
    }
}
//-----------------------------------------------------------------------------
// After construction, the currentStatus will be EFILEOP_CLOSED and the capabilities
// will be 0.
//-----------------------------------------------------------------------------
CSynFileStream::CSynFileStream()
:CFileStream()
{
    m_handle = NULL;
	m_Type = EFILESTREAM_TYPE_SYN;
}

//-----------------------------------------------------------------------------
// insert a copy constructor here... (currently disabled)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------------
CSynFileStream::~CSynFileStream()
{
    Close();
    m_handle = NULL;
	m_Type = EFILESTREAM_TYPE_SYN;
}


//-----------------------------------------------------------------------------
// Open a file in the mode specified by openMode (Read, Write, or ReadWrite).
// Truncate the file if the mode is either Write or ReadWrite and truncate is
// true.
//
// Sets capability appropriate to the openMode.
// Returns the currentStatus of the file.
//-----------------------------------------------------------------------------
CFileStream::EFILEOP_STATUS CSynFileStream::Open(const char *filename, const EFILEACCESS_MODE openMode,bool bStream)
{
 
	m_strFileName = filename; 


	m_bStream = bStream;

	assert(NULL != filename);
	 
    // Close the file if it was already open...
    if (EFILEOP_CLOSED != m_currentStatus)
        Close();

    // create the appropriate type of file...
    switch (openMode)
    {
    case EFILEACCESS_READ:
	{
		int err = fopen_s(&m_handle, m_strFileName.c_str(), "rb");
		printf("openfailed:%d", err);
	}
        break;
    case EFILEACCESS_WRITE:
		fopen_s(&m_handle, m_strFileName.c_str(),"wb");
		m_bStream = true;
        break;
    case EFILEACCESS_READWRITE:
		fopen_s(&m_handle, m_strFileName.c_str(),"a+b");//如果文件不存在，自动创建
		if(m_handle)
		{
			fclose(m_handle);
		}
		fopen_s(&m_handle, m_strFileName.c_str(),"r+b");
		m_bStream = true;
        break;
    case EFILEACCESS_WRITEAPPEND:
		fopen_s(&m_handle, m_strFileName.c_str(),"a+b");
		m_bStream = true;
        break;

    default:
		return m_currentStatus = EFILEOP_INVALIDCALL;

	}
    
    if (m_handle == NULL)                // handle not created successfully
	{
		 return ResetStatus();
	}
    else
    {
        // successfully created file, so set the file capabilities...
		m_eFileAccesssMode = openMode;
        switch (openMode)
        {
        case EFILEACCESS_READ:
        case EFILEACCESS_WRITE:
        case EFILEACCESS_WRITEAPPEND:
        case EFILEACCESS_READWRITE:
            break;
        default:
			return m_currentStatus = EFILEOP_INVALIDCALL;	
        }

		if(!bStream && m_eFileAccesssMode & unsigned long(EFILEACCESS_MODE_FLAG_READ))
		{
			delete []m_pMemBuffer;
			long  oldCursor = ftell(m_handle);
			fseek(m_handle, 0L, SEEK_END);
			m_dwMemSize = ftell(m_handle);
			fseek(m_handle,oldCursor,SEEK_SET);
			m_pMemBuffer = new unsigned char[m_dwMemSize];
			unsigned long	dwRead = 0;
			unsigned long dwReadTotal = 0;
			while(dwReadTotal < m_dwMemSize)
			{
				dwRead = fread(m_pMemBuffer,sizeof(unsigned char),m_dwMemSize - dwReadTotal,m_handle);
				if(!dwRead && feof(m_handle) == 0)
				{
					delete [] m_pMemBuffer;
					return ResetStatus();
				}
				dwReadTotal += dwRead;
			}
		}
		
                                   // success!
    }
	return m_currentStatus = EFILEOP_OK;     
}

//-----------------------------------------------------------------------------
// Get the current position of the file pointer.
//-----------------------------------------------------------------------------
unsigned long CSynFileStream::GetPosition()
{
	if(m_bStream)
	{
		unsigned long pos = ftell(m_handle);
		if(pos == 0xffffffff)
			m_currentStatus = EFILEOP_IOERR;
		return pos;
	}
	else
		return m_iMemPosition;
}

//-----------------------------------------------------------------------------
// Set the position of the file pointer.
// Absolute and relative positioning is supported via the absolutePos
// parameter.
//
// If positioning absolutely, position MUST be positive - an EFILEOP_IOERR results if
// position is negative.
// Position can be negative if positioning relatively, however positioning
// before the start of the file is an EFILEOP_IOERR.
//
// Returns the currentStatus of the file.
//-----------------------------------------------------------------------------
CFileStream::EFILEOP_STATUS CSynFileStream::SetPosition(int position, bool absolutePos)
{

    if (EFILEOP_OK != m_currentStatus && EFILEOP_EOS != m_currentStatus)
        return m_currentStatus;

	unsigned long finalPos;
	if(m_bStream)
	{
		if(absolutePos)
		{                                                  // absolute position
			assert(0 <= position);
	        
			// position beyond EFILEOP_EOS is OK
			if(fseek(m_handle,position,SEEK_SET)!=0)
				return ResetStatus();
			finalPos = ftell(m_handle);
	 
		}
		else
		{	
			// relative position
			assert( (GetPosition() >= (unsigned long)abs(position) && 0 > position) || 0 <= position);
 
			if(fseek(m_handle,position,SEEK_CUR)!=0)
				return ResetStatus();
			finalPos = ftell(m_handle);
		}

		if (0xffffffff == finalPos)
			return ResetStatus();                                        // unsuccessful
		else if (finalPos >= GetSize())
			return m_currentStatus = EFILEOP_EOS;                                // success, at end of file
		else
			return m_currentStatus = EFILEOP_OK;                                // success!
	}
	else
	{
		if(absolutePos)
			finalPos = position;
		else
			finalPos  = m_iMemPosition + position;

		if (finalPos >= GetSize())
		{
			m_iMemPosition = GetSize();
			return m_currentStatus = EFILEOP_EOS;                                // success, at end of file
		}
		else
		{
			if(finalPos < 0)
				m_iMemPosition = 0;
			m_iMemPosition = finalPos;
			return m_currentStatus = EFILEOP_OK;                                // success!
		}
	}

}

//-----------------------------------------------------------------------------
// Get the size of the file in bytes.
// It is an error to query the file size for a EFILEOP_CLOSED file, or for one with an
// error status.
//-----------------------------------------------------------------------------
unsigned long CSynFileStream::GetSize() const
{
 
    if (EFILEOP_OK == m_currentStatus || EFILEOP_EOS == m_currentStatus)
    {
		if(m_bStream)
		{
			long  oldCursor = ftell(m_handle);
			fseek(m_handle, 0L, SEEK_END);
			unsigned long high = ftell(m_handle);
			fseek(m_handle,oldCursor,SEEK_SET);
			return high;                // success!
		}
		else
			return m_dwMemSize;
    }
    else
        return 0;                                                // unsuccessful
}

//-----------------------------------------------------------------------------
// Flush the file.
// It is an error to flush a read-only file.
// Returns the currentStatus of the file.
//-----------------------------------------------------------------------------
CFileStream::EFILEOP_STATUS CSynFileStream::Flush()
{
 
    if (0 != fflush(m_handle))
        return ResetStatus();                                        // unsuccessful
    else
        return m_currentStatus = EFILEOP_OK;                                // success!
}

//-----------------------------------------------------------------------------
// Close the File.
//
// Returns the currentStatus
//-----------------------------------------------------------------------------
CFileStream::EFILEOP_STATUS CSynFileStream::Close()
{
    // check if it's already closed...
    if (EFILEOP_CLOSED == m_currentStatus)
        return m_currentStatus;

    // it's not, so close it...
    if (m_handle != NULL)
    {
        if (fclose(m_handle)!= 0)
            return ResetStatus();                                    // unsuccessful
    }
    m_handle = NULL;

	delete []m_pMemBuffer;
	m_dwMemSize = m_iMemPosition = 0;
    return m_currentStatus = EFILEOP_CLOSED;
}

//-----------------------------------------------------------------------------
// Self-explanatory.
//-----------------------------------------------------------------------------
CFileStream::EFILEOP_STATUS CSynFileStream::GetStatus() const
{
    return m_currentStatus;
}

//-----------------------------------------------------------------------------
// Sets and returns the currentStatus when an error has been encountered.
//-----------------------------------------------------------------------------
CFileStream::EFILEOP_STATUS CSynFileStream::ResetStatus()
{
	return m_currentStatus = EFILEOP_IOERR;
}



//-----------------------------------------------------------------------------
// Read from a file.
// The number of bytes to read is passed in size, the data is returned in src.
// The number of bytes read is available in bytesRead if a non-Null pointer is
// provided.
//-----------------------------------------------------------------------------
CFileStream::EFILEOP_STATUS CSynFileStream::Read(unsigned long size, char *dst, unsigned long *bytesRead)
{
	if(m_bStream)
	{
		assert(EFILEOP_CLOSED != m_currentStatus);
		assert(m_handle != NULL);
		assert(NULL != dst);
		//assert(true == HasCapability(EFILEACCESS_MODE_FLAG_READ));
		assert(0 != size);
	    
		if (EFILEOP_OK != m_currentStatus || 0 == size)
			return m_currentStatus;
		else
		{
			unsigned long lastBytes;
			unsigned long *bytes = (NULL == bytesRead) ? &lastBytes : (unsigned long *)bytesRead;
			*bytes = fread(dst,sizeof(unsigned char),size,m_handle);
			if (size != *bytes)    //这里应该是循环调用，还没有完全更改，以后看到再改
			{
				if(feof(m_handle)!= 0)
					return m_currentStatus = EFILEOP_EOS;                        // end of stream
				else
					return ResetStatus(); 
			} 
				                                   // unsuccessful
		}
		return m_currentStatus = EFILEOP_OK;                                    // successfully read size bytes
	}
	else
	{
		if (EFILEOP_OK != m_currentStatus || 0 == size)
			return m_currentStatus;
		else
		{
			unsigned long lastBytes;
			unsigned long *bytes = (NULL == bytesRead) ? &lastBytes : (unsigned long *)bytesRead;
			if(m_iMemPosition + size > m_dwMemSize)
				size = m_dwMemSize - m_iMemPosition;// unsuccessful
			memcpy(dst,m_pMemBuffer + m_iMemPosition,size);
			*bytes = size;
			m_iMemPosition += size;

			if(m_iMemPosition == m_dwMemSize)
			{
				m_currentStatus = EFILEOP_EOS; 
				return EFILEOP_OK;
			}
		}
		return m_currentStatus = EFILEOP_OK;                                    // successfully read size bytes
	}
}

//-----------------------------------------------------------------------------
// Write to a file.
// The number of bytes to write is passed in size, the data is passed in src.
// The number of bytes written is available in bytesWritten if a non-Null
// pointer is provided.
//-----------------------------------------------------------------------------
CFileStream::EFILEOP_STATUS CSynFileStream::Write(unsigned long size, const char *src, unsigned long *bytesWritten)
{
    assert(EFILEOP_CLOSED != m_currentStatus);
    assert(m_handle!= NULL);
    //assert(true == HasCapability(EFILEACCESS_MODE_FLAG_WRITE));
    assert(0 != size);
    
    if ((EFILEOP_OK != m_currentStatus && EFILEOP_EOS != m_currentStatus) || 0 == size)
        return m_currentStatus;
    else
    {
        unsigned long lastBytes;
        unsigned long *bytes = (NULL == bytesWritten) ? &lastBytes : (unsigned long *)bytesWritten;
		
		*bytes =  fwrite(src, sizeof(unsigned char),size, m_handle);
        if (*bytes  == size)
		{   
			return m_currentStatus = EFILEOP_OK;                            // success!
		}
        else
		{
            return ResetStatus();                                    // unsuccessful
		}
    }
}


CFileStream*	CreateFileStream(CFileStream::EFILESTREAM_TYPE type)
{
	CFileStream* pFileStream = NULL;
	switch (type)
	{
	case CFileStream::EFILESTREAM_TYPE_SYN:
		pFileStream = new CSynFileStream();
		break; 
	}
	return pFileStream;
}
void	DestroyFileStream(CFileStream* pFileStream)
{ 
	delete pFileStream; 
}