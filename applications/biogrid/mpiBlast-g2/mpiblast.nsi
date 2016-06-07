;NSIS Modern User Interface version 1.63
;Basic Example Script
;Written by Joost Verburg

!define MUI_PRODUCT "mpiBLAST" ;Define your own software name here
!define MUI_VERSION "1.2.1" ;Define your own software version here
!define MPIBLAST_DIR "C:\Development\mpiblast"

!include "MUI.nsh"

;--------------------------------
;Configuration

  ;General
  OutFile "${MUI_PRODUCT}-${MUI_VERSION}.exe"

  ;Folder selection page
  InstallDir "$PROGRAMFILES\${MUI_PRODUCT}"
  
  ;Remember install folder
  InstallDirRegKey HKCU "Software\${MUI_PRODUCT}" ""

;--------------------------------
;Modern UI Configuration

  !define MUI_LICENSEPAGE
  !define MUI_COMPONENTSPAGE
  !define MUI_DIRECTORYPAGE
  
  !define MUI_ABORTWARNING
  
  !define MUI_UNINSTALLER
  !define MUI_UNCONFIRMPAGE
  
;--------------------------------
;Languages
 
  !insertmacro MUI_LANGUAGE "English"
  
;--------------------------------
;Language Strings

  ;Description
  LangString DESC_SecCopyUI ${LANG_ENGLISH} "Install the mpiblast executables"

;--------------------------------
;Data
  
  LicenseData "${MPIBLAST_DIR}\LICENSE"

;--------------------------------
;Installer Sections

Section "mpiBLAST executables" SecCopyUI
  ;ADD YOUR OWN STUFF HERE!

  SetOutPath "$INSTDIR"
  File "${MPIBLAST_DIR}\AUTHORS"
  File "${MPIBLAST_DIR}\blast.cgi"
  File "${MPIBLAST_DIR}\ChangeLog"
  File "${MPIBLAST_DIR}\COPYING"
  File "${MPIBLAST_DIR}\LICENSE"
  File "${MPIBLAST_DIR}\mpiblast.conf"
  File "${MPIBLAST_DIR}\msvc_projects\Release\mpiblast.exe"
  File "${MPIBLAST_DIR}\msvc_projects\Release\mpiformatdb.exe"
  File "${MPIBLAST_DIR}\NEWS"
  File "${MPIBLAST_DIR}\README"
  File "${MPIBLAST_DIR}\TODO"
  File "${MPIBLAST_DIR}\WWWBlastWrap.pl"

  ;Store install folder
  WriteRegStr HKCU "Software\${MUI_PRODUCT}" "" $INSTDIR
  
  ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"

SectionEnd

Section "mpiBLAST source code" "Source code to used build mpiBLAST"

  SetOutPath "$INSTDIR"
  File "${MPIBLAST_DIR}\configure.in"
  File "${MPIBLAST_DIR}\depcomp"
  File "${MPIBLAST_DIR}\INSTALL"
  File "${MPIBLAST_DIR}\install-sh"
  File "${MPIBLAST_DIR}\Makefile.am"
  File "${MPIBLAST_DIR}\missing"
  File "${MPIBLAST_DIR}\mkinstalldirs"
  File "${MPIBLAST_DIR}\mpiblast.nsi"
  File "${MPIBLAST_DIR}\readdb-patch-2002-08"
  File "${MPIBLAST_DIR}\readdb.patch"
  File "${MPIBLAST_DIR}\readdb.patch-2002-12"
  File "${MPIBLAST_DIR}\patch-NCBIToolbox_Nov14_2003"

  SetOutPath "$INSTDIR\src"
  File "${MPIBLAST_DIR}\src\blastjob.cpp"
  File "${MPIBLAST_DIR}\src\blastjob.hpp"
  File "${MPIBLAST_DIR}\src\blast_hooks.c"
  File "${MPIBLAST_DIR}\src\blast_hooks.h"
  File "${MPIBLAST_DIR}\src\db_spec.cpp"
  File "${MPIBLAST_DIR}\src\db_spec.hpp"
  File "${MPIBLAST_DIR}\src\distributed_bioseq.c"
  File "${MPIBLAST_DIR}\src\distributed_bioseq.h"
  File "${MPIBLAST_DIR}\src\embed_rank.cpp"
  File "${MPIBLAST_DIR}\src\embed_rank.hpp"
  File "${MPIBLAST_DIR}\src\file_util.cpp"
  File "${MPIBLAST_DIR}\src\file_util.hpp"
  File "${MPIBLAST_DIR}\src\formatdb_hooks.h"
  File "${MPIBLAST_DIR}\src\fragment_list.cpp"
  File "${MPIBLAST_DIR}\src\fragment_list.hpp"
  File "${MPIBLAST_DIR}\src\getopt.c"
  File "${MPIBLAST_DIR}\src\getopt.h"
  File "${MPIBLAST_DIR}\src\getopt1.c"
  File "${MPIBLAST_DIR}\src\Makefile.am"
  File "${MPIBLAST_DIR}\src\mpiblast.cpp"
  File "${MPIBLAST_DIR}\src\mpiblast.hpp"
  File "${MPIBLAST_DIR}\src\mpiblast_config.cpp"
  File "${MPIBLAST_DIR}\src\mpiblast_config.hpp"
  File "${MPIBLAST_DIR}\src\mpiblast_types.h"
  File "${MPIBLAST_DIR}\src\mpiblast_util.cpp"
  File "${MPIBLAST_DIR}\src\mpiblast_util.hpp"
  File "${MPIBLAST_DIR}\src\mpiformatdb.cpp"
  File "${MPIBLAST_DIR}\src\ncbi_sizeof.c"
  File "${MPIBLAST_DIR}\src\ncbi_sizeof.h"

  SetOutPath "$INSTDIR\msvc_projects"
  File "${MPIBLAST_DIR}\msvc_projects\mpiblast.sln"
  File "${MPIBLAST_DIR}\msvc_projects\mpiblast.vcproj"
  File "${MPIBLAST_DIR}\msvc_projects\mpiformatdb.vcproj"
SectionEnd

;Display the Finish header
;Insert this macro after the sections if you are not using a finish page
!insertmacro MUI_SECTIONS_FINISHHEADER

;--------------------------------
;Descriptions

!insertmacro MUI_FUNCTIONS_DESCRIPTION_BEGIN
  !insertmacro MUI_DESCRIPTION_TEXT ${SecCopyUI} $(DESC_SecCopyUI)
!insertmacro MUI_FUNCTIONS_DESCRIPTION_END
 
;--------------------------------
;Uninstaller Section

Section "Uninstall"

  ;ADD YOUR OWN STUFF HERE!

  Delete "$INSTDIR\Uninstall.exe"
  Delete /REBOOTOK "$INSTDIR\AUTHORS"
  Delete /REBOOTOK "$INSTDIR\blast.cgi"
  Delete /REBOOTOK "$INSTDIR\ChangeLog"
  Delete /REBOOTOK "$INSTDIR\COPYING"
  Delete /REBOOTOK "$INSTDIR\LICENSE"
  Delete /REBOOTOK "$INSTDIR\mpiblast.conf"
  Delete /REBOOTOK "$INSTDIR\mpiblast.exe"
  Delete /REBOOTOK "$INSTDIR\mpiformatdb.exe"
  Delete /REBOOTOK "$INSTDIR\NEWS"
  Delete /REBOOTOK "$INSTDIR\readdb-patch-2002-08"
  Delete /REBOOTOK "$INSTDIR\readdb.patch"
  Delete /REBOOTOK "$INSTDIR\readdb.patch-2002-12"
  Delete /REBOOTOK "$INSTDIR\README"
  Delete /REBOOTOK "$INSTDIR\TODO"
  Delete /REBOOTOK "$INSTDIR\WWWBlastWrap.pl"

  Delete /REBOOTOK "$INSTDIR\configure.in"
  Delete /REBOOTOK "$INSTDIR\depcomp"
  Delete /REBOOTOK "$INSTDIR\INSTALL"
  Delete /REBOOTOK "$INSTDIR\install-sh"
  Delete /REBOOTOK "$INSTDIR\Makefile.am"
  Delete /REBOOTOK "$INSTDIR\missing"
  Delete /REBOOTOK "$INSTDIR\mkinstalldirs"
  Delete /REBOOTOK "$INSTDIR\mpiblast.nsi"

  Delete /REBOOTOK "$INSTDIR\src\blastjob.cpp"
  Delete /REBOOTOK "$INSTDIR\src\blastjob.hpp"
  Delete /REBOOTOK "$INSTDIR\src\blast_hooks.c"
  Delete /REBOOTOK "$INSTDIR\src\blast_hooks.h"
  Delete /REBOOTOK "$INSTDIR\src\db_spec.cpp"
  Delete /REBOOTOK "$INSTDIR\src\db_spec.hpp"
  Delete /REBOOTOK "$INSTDIR\src\distributed_bioseq.c"
  Delete /REBOOTOK "$INSTDIR\src\distributed_bioseq.h"
  Delete /REBOOTOK "$INSTDIR\src\embed_rank.cpp"
  Delete /REBOOTOK "$INSTDIR\src\embed_rank.hpp"
  Delete /REBOOTOK "$INSTDIR\src\file_util.cpp"
  Delete /REBOOTOK "$INSTDIR\src\file_util.hpp"
  Delete /REBOOTOK "$INSTDIR\src\formatdb_hooks.h"
  Delete /REBOOTOK "$INSTDIR\src\fragment_list.cpp"
  Delete /REBOOTOK "$INSTDIR\src\fragment_list.hpp"
  Delete /REBOOTOK "$INSTDIR\src\getopt.c"
  Delete /REBOOTOK "$INSTDIR\src\getopt.h"
  Delete /REBOOTOK "$INSTDIR\src\getopt1.c"
  Delete /REBOOTOK "$INSTDIR\src\Makefile.am"
  Delete /REBOOTOK "$INSTDIR\src\mpiblast.cpp"
  Delete /REBOOTOK "$INSTDIR\src\mpiblast.hpp"
  Delete /REBOOTOK "$INSTDIR\src\mpiblast_config.cpp"
  Delete /REBOOTOK "$INSTDIR\src\mpiblast_config.hpp"
  Delete /REBOOTOK "$INSTDIR\src\mpiblast_types.h"
  Delete /REBOOTOK "$INSTDIR\src\mpiblast_util.cpp"
  Delete /REBOOTOK "$INSTDIR\src\mpiblast_util.hpp"
  Delete /REBOOTOK "$INSTDIR\src\mpiformatdb.cpp"
  Delete /REBOOTOK "$INSTDIR\src\ncbi_sizeof.c"
  Delete /REBOOTOK "$INSTDIR\src\ncbi_sizeof.h"

  Delete /REBOOTOK "$INSTDIR\msvc_projects\mpiblast.sln"
  Delete /REBOOTOK "$INSTDIR\msvc_projects\mpiblast.vcproj"
  Delete /REBOOTOK "$INSTDIR\msvc_projects\mpiformatdb.vcproj"

  RMDir "$INSTDIR\msvc_projects"
  RMDir "$INSTDIR\src"
  RMDir "$INSTDIR"

  DeleteRegKey /ifempty HKCU "Software\${MUI_PRODUCT}"
  
  ;Display the Finish header
  !insertmacro MUI_UNFINISHHEADER

SectionEnd