#!/bin/bash
set -e
cd "$(dirname "$0")"

# only if template.zip actually exists (README says the library ships unpacked)
[ -f template.zip ] && unzip -o -q template.zip

cd external_tools
for t in fiberdock multiprot naccess; do unzip -o -q "$t.zip"; done

cd naccess && csh install.scr && cd ..
cd BeEM-master && g++ -O3 BeEM.cpp -o BeEM && cd ..

chmod 755 multiprot/multiprot.Linux
chmod 755 fiberdock/FiberDock fiberdock/nma