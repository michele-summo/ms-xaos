How to create (pre)releases by using the scripts in this folder?
================================================================

In this folder we provide scripts that can create installation bundles
for various platforms.

We assume that you are already familiar with the
[steps for compiling XaoS](https://github.com/xaos-project/XaoS/wiki/Developer's-Guide#build-instructions).

In a nutshell, you can build XaoS via qmake or cmake. When using qmake,
it puts the binary into `../bin`. By contrast, cmake will use a folder
you set at the configuration step; most users choose the folder `../build`
(but it's up to you). When using cmake, you need to set (see below) explicitly
where the binary has been created.

Debian Linux and variants
-------------------------

Simply run the script `create-deb` from command line. See the comments
in the script in the first few rows to learn which packages must be
installed in advance. Also, if you built XaoS via cmake, you need to
set the build path given as an argument.

Microsoft Windows
-----------------

First, make sure that the required DLLs are collected in the folder `../bin`.
This can be done via the script `create-windows-bin-folder`. When
using cmake, enter the build directory as first (and only) argument.

Before running the script `deploy-win.bat`, you need to add the path of
tools `windeployqt.exe` and `binarycreator.exe` to the system
environment variable PATH. These tools are usually located in Qt's
`Tools\QtInstallerFramework\<VERSION>\bin` and
`Tools\QtDesignStudio\qt6_design_studio_reduced_version\bin`.
When installing Qt, make sure that the required packages are present.

To avoid problems with finding certain files, you should make sure that
there is no special character in the full path of the `XaoS` folder.
Otherwise some files may be missing from the installation bundle (for
example, the .cat files).

The steps above will create the file `..\xaossetup.exe` which can be
used to install XaoS via a simple Qt-based installer.

Alternatively, we provide a .zip file as well. To create it, use the
script `create-windows-zip`. It assumes that the Qt-based installer
has already been created.

MS XaoS: the whole release in one command
-----------------------------------------

The scripts below this one are upstream's, and are built around MSYS2 and the
Qt Installer Framework. `ms-release` is built around what this fork actually
compiles with -- cmake, MinGW, and the `deploy` target, which fills `../bin`
with the Qt and compiler runtime through windeployqt.

```
cd tools
./ms-release              # build the zip and say what to do with it
./ms-release --tag        # and make the tag and push it
./ms-release --publish    # and create the release on GitHub
```

It never asks which version is being released: it reads that from
`src/include/config.h`, so the tag, the name of the zip and the About box
cannot disagree. The tag is `v` and the version -- `v1.6`.

What it does, in order:

1. **Refuses to go on if the working tree has uncommitted changes.** A release
   built from a tree nobody can check out again is the one thing a release must
   not be.
2. Finds the configured build folder under `build/`, builds, and **runs the
   tests**. Failing tests stop it.
3. Runs the `deploy` target, which is windeployqt plus libquadmath -- the one
   thing windeployqt does not know about, since only `XaoS-quad` links it and it
   comes from the toolchain rather than from Qt.
4. Puts `bin/` together with `catalogs/`, `tutorial/` and `examples/` -- the
   three the program looks for at run time -- plus `COPYING`, the guide and the
   release notes, and zips the lot.
5. With `--tag`, tags and pushes. With `--publish`, calls `gh release create`.

`--publish` wants the GitHub CLI: `winget install GitHub.cli`, then
`gh auth login` once. Without it the script prints the command to run, and the
release can equally be made by hand on github.com -- Releases, Draft a new
release, choose the tag, paste `doc/RELEASE-NOTES.md`, attach the zip.

The notes live in `doc/RELEASE-NOTES.md` and are written by hand. Update them
before tagging, not after.

### Doing it without the script

```
cmake --build build/<your-build-folder> --target deploy
git tag -a v1.6 -m "MS XaoS 1.6"
git push origin v1.6
```

then zip `bin catalogs tutorial examples` and attach that to the release.

### What is not automated

`.github/workflows/build.yml` compiles on Linux and macOS at every push. It has
no Windows job and produces no artefacts, so nothing about the release happens
by itself. Making it happen on a tag push would mean adding a Windows job that
installs Qt, runs `deploy`, and uploads -- worth doing if releases become
frequent, and not before.

MacOS
-----

Before running the script `deploy-mac` from command line you need to add
the path of tool `macdeployqt` to the system environment variable PATH.

Depending on which build system you used (qmake or cmake) you may need
to add an argument that points to the build folder.

WebAssembly
-----------

Since Qt 6.8 the WebAssembly build is straightforward. It can be
performed via both build systems (qmake or cmake).

After a successful build, you need the files qtloader.js, qtlogo.svg,
xaos.data, xaos.html, xaos.js, xaos.wasm and xaos.worker.js. Please
follow the steps in the wiki if you want to install them on your web
server, similarly to the publicly available website, https://xaos.app.
