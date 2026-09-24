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
#
# Every render also has to exit cleanly. The images alone do not show it: from
# 1.6 to 1.7.3 each run crashed on its way out, after the PNG was written, and
# this test went on passing because it looked only at the pictures.

# Renders input inside dir as prefix000000.png. A non-zero exit -- or the text
# CMake reports in place of a number when the process crashes -- fails the
# position, and whatever it said on stderr is shown.
macro(xaos_render input prefix)
   execute_process(COMMAND ${XAOS} -render ${input} -size 80x60 -basename ${prefix}
                   WORKING_DIRECTORY ${dir} OUTPUT_QUIET
                   RESULT_VARIABLE render_rc ERROR_VARIABLE render_err)
   if(NOT render_rc STREQUAL "0")
      list(APPEND failures "${name}: rendering ${input} exited with ${render_rc}")
      string(STRIP "${render_err}" render_err)
      if(render_err)
         message(STATUS "  ${name}: ${input} said: ${render_err}")
      endif()
      set(clean FALSE)
   endif()
endmacro()

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

   set(clean TRUE)
   xaos_render(${basename} a)
   xaos_render(replay.xaf t)

   if(NOT EXISTS ${dir}/out.xpf)
      list(APPEND failures "${name}: nothing was saved")
      continue()
   endif()

   xaos_render(out.xpf b)

   if(NOT EXISTS ${dir}/a000000.png OR NOT EXISTS ${dir}/b000000.png)
      list(APPEND failures "${name}: one of the two renders produced no image")
      continue()
   endif()

   # Compared even after a bad exit, so that one run reports both.
   file(MD5 ${dir}/a000000.png before)
   file(MD5 ${dir}/b000000.png after)
   if(NOT before STREQUAL after)
      list(APPEND failures "${name}: reloading it renders differently")
      message(STATUS "  ${name}: kept ${dir} for inspection")
   elseif(NOT clean)
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
   message(FATAL_ERROR "${count} positions checked, not all of them render, "
                       "save and reload cleanly")
endif()

message(STATUS "${count} positions survive a save and load unchanged, "
               "every render exiting cleanly")
