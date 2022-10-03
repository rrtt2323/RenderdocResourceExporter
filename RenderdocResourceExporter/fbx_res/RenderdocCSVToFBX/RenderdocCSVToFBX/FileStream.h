#ifndef _FILESTREAM_H_INCLUDE
#define _FILESTREAM_H_INCLUDE

#include <string>

//class CFileStream;
//class  CFileMaster
//{ 
//public:
//	CFileMaster()
//	{
//		m_bFileLoaded = false;
//	}
//public:
//	virtual void	OnFileLoaded(CFileStream* pFileStream);
//	virtual void	ProcessFileLoaded(CFileStream* pFileStream)= 0;
//	bool	IsFileLoaded()
//	{
//		return m_bFileLoaded;
//	}
//private:
//	bool	m_bFileLoaded;
//};

class  CFileStream
{
public:
	enum  EFILESTREAM_TYPE
	{
		EFILESTREAM_TYPE_SYN,
		EFILESTREAM_TYPE_ASYN,
		EFILESTREAM_TYPE_SYNPACK,
		EFILESTREAM_TYPE_ASYNPACK,
	};
   /// What is the status of our file handle?
   enum EFILEOP_STATUS
   {
      EFILEOP_OK = 0,           ///< EFILEOP_OK!
      EFILEOP_IOERR,          ///< Read or Write error
      EFILEOP_EOS,              ///< End of Stream reached (mostly for reads)
      EFILEOP_INVALIDCALL,      ///< An unsupported operation used.  Always accompanied by AssertWarn
      EFILEOP_CLOSED,           ///< Tried to operate on a closed stream (or detached filter)
	  EFILEOP_LOADING,			///< Asynchronous  loading now
	  EFILEOP_PRELOAD,			///< Asynchronous  pre for load
	  EFILEOP_LOADED,
	  EFILEOP_WAITINGDEL,
	  EFILEOP_WAITWRITE,			///<  
      EFILEOP_UNKNOWNERR      ///< Catchall
   };

   /// How are we accessing the file?
   enum EFILEACCESS_MODE
   {
      EFILEACCESS_READ         = 1,  ///< Open for read only, starting at beginning of file.
      EFILEACCESS_WRITE        = 2,  ///< Open for write only, starting at beginning of file; will blast old contents of file.
      EFILEACCESS_READWRITE    = 3,  ///< Open for read-write.
      EFILEACCESS_WRITEAPPEND  = 6   ///< Write-only, starting at end of file.
   };

	/// Flags used to indicate what we can do to the file.
	enum EFILEACCESS_MODE_FLAG
	{
		EFILEACCESS_MODE_FLAG_READ       = 0,
		EFILEACCESS_MODE_FLAG_WRITE		 = 1
	};


	CFileStream();
	virtual ~CFileStream(){};
protected:
	EFILEOP_STATUS				m_currentStatus;   ///< Current status of the file (EFILEOP_OK, EFILEOP_IOERR, etc.).
	EFILEACCESS_MODE			m_eFileAccesssMode;         ///< Keeps track of file capabilities.
	//CFileMaster*				m_pFileMaster;
	EFILESTREAM_TYPE			m_Type;
 
	unsigned char*						m_pMemBuffer;
	unsigned long						m_dwMemSize;
	std::string					m_strFileName;
	int							m_iMemPosition;

	bool						m_bStream;			//aysn file stream not support stream mode,only synFile support stream mode,
										//when a filestream is open as stream,caller read data when need,not read all file into the memory
public:
	EFILEACCESS_MODE	GetAccessMode()
	{
		return m_eFileAccesssMode;
	}

	EFILESTREAM_TYPE	GetType()
	{
		return m_Type;
	}

	unsigned char*	GetContent()
	{
		return (unsigned char*)m_pMemBuffer;
	}

   //! name	Open.
   //! param	const char * filename - 
   //! param	const EFILEACCESS_MODE openMode - 
   //! param	bool bStream -bStream if this param set true ,the file is open as a stream if this param set false,the file will be read to the memory
   //! note		Opens a file for access using the specified EFILEACCESS_MODE
   //! return 	Status The status of the file
   virtual EFILEOP_STATUS Open(const char *filename, const EFILEACCESS_MODE openMode,bool bStream= true) = 0;
 
   //! name	GetPosition.
   //! note	 Gets the current position in the file
   //! return 	unsigned long	This is in bytes from the beginning of the file.
   virtual unsigned long			GetPosition()  = 0;

   /// Sets the current cursor position in the file.
   virtual EFILEOP_STATUS	SetPosition(int position, bool absolutePos = true) = 0;
   
   /// Returns the size of the file
   virtual unsigned long GetSize()const 
   {
	   return m_dwMemSize;
   }
 
   virtual const std::string &GetFileName() const
   {
	   return m_strFileName;
   }
   
   /// Make sure everything that's supposed to be written to the file gets written.
   ///
   /// @returns The status of the file.
   virtual EFILEOP_STATUS Flush() = 0;
   
   /// Closes the file
   ///
   /// @returns The status of the file.
   virtual EFILEOP_STATUS Close() = 0;
   
   /// Gets the status of the file
   virtual EFILEOP_STATUS GetStatus() const = 0;
   
   /// Reads "size" bytes from the file, and dumps data into "dst".
   /// The number of actual bytes read is returned in bytesRead
   /// @returns The status of the file
   virtual  EFILEOP_STATUS Read(unsigned long size, char *dst, unsigned long *bytesRead = NULL) = 0;
   
   /// Writes "size" bytes into the file from the pointer "src".
   /// The number of actual bytes written is returned in bytesWritten
   /// @returns The status of the file
   virtual EFILEOP_STATUS Write(unsigned long size, const char *src, unsigned long *bytesWritten = NULL) = 0;

   /// Returns whether or not this file is capable of the given function.
   bool HasCapability(EFILEACCESS_MODE_FLAG cap) const;
	
   EFILEOP_STATUS SetStatus(EFILEOP_STATUS status);    ///< Setter for the current status.

   //syn load function used in syn load
   virtual    void			ReadFileSize(){};
   virtual    void			AllocFileBuffer(){};
	virtual EFILEOP_STATUS	LoadFile(){return EFILEOP_OK;};  //load file to stream memory
	virtual EFILEOP_STATUS	SaveFile(){return EFILEOP_OK;};  //load file to stream memory
 
    virtual EFILEOP_STATUS	ResetStatus()= 0;                 ///< Called after error encountered.
	
	virtual void	NotifyMasterLoaded() = 0;

public:
	bool	ReadInt(int&	value);
	bool	WriteInt(int	value);

	bool	ReadDWORD(unsigned long&  value);
	bool	WriteDWORD(unsigned long value);

	bool	WriteWORD(unsigned short value);
	bool	ReadWORD(unsigned short& value);
 
	bool	WriteShort(short value);
	bool	ReadShort(short& value);

	bool	WriteBYTE(unsigned char value);
	bool	ReadBYTE(unsigned char& value);

	bool	Readbool(bool& value);
	bool	Writebool(bool value);


	bool	ReadBool(bool& value);
	bool	WriteBool(bool value);

	bool	ReadFloat(float& value) ;
	bool	WriteFloat(float value) ;
 
	bool	WriteString(const char *string, int maxLen = 255);
	bool	ReadString(char* buf,int bufSize);
	bool	ReadString(std::string& str,bool sizeWord = true);

	bool	ReadLine(std::string& outStr);
};


class CSynFileStream: public CFileStream
{
public:
   CSynFileStream();                     ///< Default constructor
   virtual ~CSynFileStream();            ///< Destructor

   /// Opens a file for access using the specified EFILEACCESS_MODE
   ///
   /// @returns The status of the file
   EFILEOP_STATUS	Open(const char *filename, const EFILEACCESS_MODE openMode,bool bStream = true);


   /// Gets the current position in the file
   ///
   /// This is in bytes from the beginning of the file.
   unsigned long			GetPosition() ;

   /// Sets the current position in the file.
   EFILEOP_STATUS	SetPosition(int position, bool absolutePos = true);
   
   /// Returns the size of the file
   unsigned long GetSize() const;
   
   /// Make sure everything that's supposed to be written to the file gets written.
   ///
   /// @returns The status of the file.
   EFILEOP_STATUS	Flush();
   
   /// Closes the file
   ///
   /// @returns The status of the file.
   EFILEOP_STATUS	Close();
   
   /// Gets the status of the file
   EFILEOP_STATUS	GetStatus() const;
   
   /// Reads "size" bytes from the file, and dumps data into "dst".
   /// The number of actual bytes read is returned in bytesRead
   /// @returns The status of the file
   EFILEOP_STATUS	Read(unsigned long size, char *dst, unsigned long *bytesRead = NULL);
   
   /// Writes "size" bytes into the file from the pointer "src".
   /// The number of actual bytes written is returned in bytesWritten
   /// @returns The status of the file
   EFILEOP_STATUS	Write(unsigned long size, const char *src, unsigned long *bytesWritten = NULL);

 
protected:
   EFILEOP_STATUS	ResetStatus();                 ///< Called after error encountered.
   virtual void		NotifyMasterLoaded() {};
protected:
	FILE*		m_handle;           ///< Pointer to the file handle.
};

//extent function
 std::string FindExtensionFromFileName(const std::string& name); 
 std::string FindShortNameFromFileName(const std::string& name,bool bIncludeExt= false);

//-------------------------------------- Helper Functions
//将"/"转换为"/"
static void forwardslash(char *str)
{
   while(*str)
   {
      if(*str == '\\')
         *str = '/';
      str++;
   }
}

//将"/"转换为"/"
static void backslash(char *str)
{
   while(*str)
   {
      if(*str == '/')
         *str = '\\';
      str++;
   }
}

 CFileStream*	CreateFileStream(CFileStream::EFILESTREAM_TYPE type = CFileStream::EFILESTREAM_TYPE_SYN);
 void			DestroyFileStream(CFileStream* pFileStream);


 
#endif // _NBFILE_H_INCLUDE
