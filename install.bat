set src=RenderdocResourceExporter

set dst=%APPDATA%\qrenderdoc\extensions\%src%
if not exist "%dst%" mkdir "%dst%"
xcopy "%~dp0%src%\*" "%dst%" /i /e /Y /C

::pause