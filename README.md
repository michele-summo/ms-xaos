MS XaoS
=======

A fork of [XaoS](https://github.com/xaos-project/XaoS) 4.3.3 by
[Michele Summo](https://www.michelesummo.it). Everything the original does, it
still does; what has been added is described in
**[the guide](doc/ms-xaos-guide.md)**, and in short it is:

- **A second binary at 128-bit precision**, so a zoom can go some seven orders
  of magnitude deeper before the arithmetic runs out.
- **Coherent noise in the user formula** — `randsc` and its family, which lay a
  random field over the plane in blobs, squares, irregular polygons, hexagons
  or triangles, with a kaleidoscope to fold the plane first.
- **Orbit traps and stripe averaging**, as ordinary formula functions.
- **Eighteen more colouring modes**, inside and out, and smooth colouring made
  to work for the user formula, for every bailout shape and for the Newton
  emulation.
- **A user formula language that says what it means**: every function
  documented with what it takes and what it defaults to, `poly`, `ifiter` and
  its relatives, arguments that may be left out by leaving their place empty,
  and `p1` to `p9999` for reaching back into the orbit.
- **A reference window** listing every function, variable and notation the
  formula language accepts, beside the formula being written.

335 regression tests, run at both precisions.

Supporting this fork
--------------------

MS XaoS is free software and is written in spare time. If it is useful to you,
you can support the work at **[Buy Me a Coffee](https://buymeacoffee.com/michelesummo)**.
There is more about the project, and about what else I do, at
[www.michelesummo.it](https://www.michelesummo.it).

About XaoS
----------

XaoS (pronounced *chaos*) is a realtime interactive fractal zoomer. This
means that it lets you zoom smoothly into any place in the fractal you
choose without the long calculation required by other fractal generators.
It has many other features too, like different fractal types, autopilot,
special coloring modes, random palette generation, color cycling, etc.

XaoS (since version 4.3, September 2023) is available as a full-featured web application as well.
You can try it at [xaos.app](https://xaos.app).
Also, you can try a simplified version of the XaoS zooming algorithm online using
[XaoS.js](https://xaos-project.github.io/XaoSjs/).

XaoS is based on [Qt](http://www.qt.io), and has been tested on Windows, Mac and Linux. It
should also work on any other platform supporting Qt Widgets, like the BSDs.

- Original authors: [Jan Hubička](https://www.ucw.cz/~hubicka/) and [Thomas Marsh](https://www.linkedin.com/in/thomasmarsh).
- Current maintainers: [J.B. Langston](https://www.linkedin.com/in/jblangston/) and [Zoltán Kovács](https://matek.hu/zoltan).
- See [CREDITS](CREDITS.md) for a complete list of contributors.

Project Resources
-----------------

This fork:

- [GitHub Repo](https://github.com/michele-summo/ms-xaos)
- [What this fork adds](doc/ms-xaos-guide.md)
- [Support the work](https://buymeacoffee.com/michelesummo)

The original project:

- [GitHub Repo](https://github.com/xaos-project/XaoS)
- [Binary Releases](https://github.com/xaos-project/XaoS/releases)
- [Documentation](https://github.com/xaos-project/XaoS/wiki)
- [Issue Tracker](https://github.com/xaos-project/XaoS/issues)
- [XaoS website](http://xaos-project.github.io/)
- [Users' Google Group](http://groups.google.com/group/xaos-users)
- [Developer's Guide](https://github.com/xaos-project/XaoS/wiki/Developer's-Guide)

Bugs and Support
----------------

XaoS is a community-supported free software project. The Google Group provides
a place for XaoS users to help each other get the most out of XaoS. XaoS
developers also monitor this discussion and answer questions from time to time.
You can browse the archives freely but to prevent spam, you must join the group
in order to post. Google Groups provides options for participation from a
traditional mailing list to a completely web-based forum, so you don’t have to
get emails if you don’t want to.

If you think you have found a bug in XaoS or have an idea for a new feature,
please let us know about it -- for this fork, on
[its own issue tracker](https://github.com/michele-summo/ms-xaos/issues), and
for anything the original does too, upstream. XaoS is developed on a volunteer basis and the
developers work on it in their spare time. Therefore, we can’t guarantee that
issues will be addressed in a certain timeframe. If you are able to fix a bug
or implement a new feature yourself, pull requests are very welcome.

License
-------

Copyright © 1996-2024 XaoS Contributors

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
