# $Id$

if (MAINTAINER)
  add_custom_target (distrib-man)
  add_custom_target (man ALL)
endif (MAINTAINER)

set (MANPAGES)
set (USAGE)
set (HTMLMAN)

if (MAINTAINER)
  set (MANDIR ${CMAKE_CURRENT_BINARY_DIR})
else (MAINTAINER)
  set (MANDIR ${CMAKE_CURRENT_SOURCE_DIR})
endif (MAINTAINER)

foreach (TOOL ${TOOLS})
  set (MANPAGES ${MANPAGES} ${MANDIR}/${TOOL}.1)
  set (USAGE ${USAGE} ${MANDIR}/${TOOL}.usage)
  set (HTMLMAN ${HTMLMAN} ${MANDIR}/${TOOL}.1.html)
  if (MAINTAINER)
    add_custom_command (OUTPUT ${TOOL}.1
      COMMAND pod2man --center=\"GeographicLib Utilities\"
	--release=\"GeographicLib ${GeographicLib_VERSION_STRING}\"
        ${CMAKE_CURRENT_SOURCE_DIR}/${TOOL}.pod > ${TOOL}.1
      COMMENT "Building man page for ${TOOL}"
      MAIN_DEPENDENCY ${TOOL}.pod)
    add_custom_command (OUTPUT ${TOOL}.1.html
      COMMAND
        pod2html --noindex ${CMAKE_CURRENT_SOURCE_DIR}/${TOOL}.pod |
	sed -e 's%<head>%<head><link href="http://search.cpan.org/s/style.css" rel="stylesheet" type="text/css">%'
	  -e 's%<code>\\\([^<>]*\\\)\(\\\(.\\\)\)</code>%<a href="\\1.\\2.html">&</a>%'g > ${TOOL}.1.html
      COMMENT "Building html version of man page for ${TOOL}"
      MAIN_DEPENDENCY ${TOOL}.pod)
    add_custom_command (OUTPUT ${TOOL}.usage
      COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/makeusage.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/${TOOL}.pod > ${TOOL}.usage
      COMMENT "Building usage code for ${TOOL}"
      MAIN_DEPENDENCY ${TOOL}.pod)
  endif (MAINTAINER)
endforeach (TOOL)

if (MAINTAINER)
  add_custom_target (manpages ALL DEPENDS ${MANPAGES}
    COMMENT "Building all the man pages")
  add_custom_target (usage ALL DEPENDS ${USAGE}
    COMMENT "Converting the man pages to online usage")
  add_custom_target (htmlman ALL DEPENDS ${HTMLMAN}
    COMMENT "Building html versions of the man pages")
  add_dependencies (man manpages usage htmlman)
  add_dependencies (distrib-man man)
  add_custom_command (TARGET distrib-man
    COMMAND
      for f in ${MANPAGES} ${USAGE} ${HTMLMAN}\; do
        cmp "$$f" ${CMAKE_CURRENT_SOURCE_DIR}/`basename "$$f"` > /dev/null 2>&1 ||
        install -m 644 "$$f" ${CMAKE_CURRENT_SOURCE_DIR}\; done
    COMMENT "Installing man documentation page in source tree")
else (MAINTAINER)
  install (FILES ${MANPAGES} DESTINATION man/man1)
endif (MAINTAINER)
