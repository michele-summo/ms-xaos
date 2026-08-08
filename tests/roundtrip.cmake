# Saving a position and loading it back has to land on the same image.
#
# Run as: cmake -DXAOS=<binary> -DPOSITIONS=<dir> -DWORK=<dir> -P roundtrip.cmake
#
# For each position: render it, replay it with a savepos appended, render what
# that wrote, and compare the two images. Any difference means the file does
# not describe the state it was saved from -- which is invisible until someone
# reopens a picture months later and finds it changed.
#
# The comparison is on the rendered pixels rather than on the file text,
# because the text legitimately differs: the saver writes defaults the original
# left implicit, more digits than an older version did, and its own version
# number in the header.

file(GLOB positions ${POSITIONS}/*.xpf)
list(LENGTH positions count)
if(count EQUAL 0)
   message(FATAL_ERROR "no positions found in ${POSITIONS}")
endif()

set(failures "")

foreach(position ${positions})
   get_filename_component(name ${position} NAME_WE)
   set(dir ${WORK}/${name})
   file(REMOVE_RECURSE ${dir})
   file(MAKE_DIRECTORY ${dir})
   file(COPY ${position} DESTINATION ${dir})
   get_filename_component(basename ${position} NAME)

   # XaoS resolves the savepos path against the working directory, so
   # everything happens inside dir with plain relative names.
   file(READ ${position} text)
   file(WRITE ${dir}/replay.xaf "${text}\n(savepos \"out.xpf\")\n")

   execute_process(COMMAND ${XAOS} -render ${basename} -size 80x60 -basename a
                   WORKING_DIRECTORY ${dir} OUTPUT_QUIET ERROR_QUIET)
   execute_process(COMMAND ${XAOS} -render replay.xaf -size 80x60 -basename t
                   WORKING_DIRECTORY ${dir} OUTPUT_QUIET ERROR_QUIET)

   if(NOT EXISTS ${dir}/out.xpf)
      list(APPEND failures "${name}: nothing was saved")
      continue()
   endif()

   execute_process(COMMAND ${XAOS} -render out.xpf -size 80x60 -basename b
                   WORKING_DIRECTORY ${dir} OUTPUT_QUIET ERROR_QUIET)

   if(NOT EXISTS ${dir}/a000000.png OR NOT EXISTS ${dir}/b000000.png)
      list(APPEND failures "${name}: one of the two renders produced no image")
      continue()
   endif()

   file(MD5 ${dir}/a000000.png before)
   file(MD5 ${dir}/b000000.png after)
   if(NOT before STREQUAL after)
      list(APPEND failures "${name}: reloading it renders differently")
      message(STATUS "  ${name}: kept ${dir} for inspection")
   else()
      message(STATUS "ok  ${name}")
      file(REMOVE_RECURSE ${dir})
   endif()
endforeach()

if(failures)
   foreach(f ${failures})
      message(SEND_ERROR "FAIL ${f}")
   endforeach()
   message(FATAL_ERROR "${count} positions checked, save/load is not faithful")
endif()

message(STATUS "${count} positions survive a save and load unchanged")
