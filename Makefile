# Copyright (C) 2015 Pierre de Buyl

# This file is part of vmf90

# vmf90 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# vmf90 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with vmf90.  If not, see <http://www.gnu.org/licenses/>.

build:
	@mkdir build

hmf: | build 
	make -C build -f ../scripts/Makefile $@ 

fel: | build
	make -C build -f ../scripts/Makefile $@

clean:
	make -C build -f ../scripts/Makefile $@

doc:
	(cd doc ; doxygen doxy_conf)

.PHONY: doc
