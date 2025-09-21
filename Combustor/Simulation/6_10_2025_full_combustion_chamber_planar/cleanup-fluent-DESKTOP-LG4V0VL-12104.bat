echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v252\fluent/ntbin/win64/winkill.exe"

start "tell.exe" /B "C:\PROGRA~1\ANSYSI~1\v252\fluent\ntbin\win64\tell.exe" DESKTOP-LG4V0VL 64248 CLEANUP_EXITING
timeout /t 1
"C:\PROGRA~1\ANSYSI~1\v252\fluent\ntbin\win64\kill.exe" tell.exe
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 29052) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 13848) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 8048) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 27252) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 12104) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 16132)
del "C:\Users\rjmcm\Documents\GitHub\Turbojet\Combustor\Simulation\6_10_2025_full_combustion_chamber_planar\cleanup-fluent-DESKTOP-LG4V0VL-12104.bat"
