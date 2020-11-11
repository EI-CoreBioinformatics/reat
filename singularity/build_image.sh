#!/usr/bin/env bash
set -euxo
version=0.0.1
rundir=${PWD}
cd $(mktemp -d)
sudo singularity build reat.img ${rundir}/reat_singularity.def
mkdir -p /ei/software/testing/reat/${version}/x86_64/bin
cp reat.img /ei/software/testing/reat/${version}/x86_64/reat-${version}.img
cat > singularity.exec <<SINGULARITY_EXEC
#!/bin/bash
DIR=\`dirname \$(readlink -f \$0)\`
singularity exec \$DIR/../reat-${version}.img \$(basename "\$0") \$@
SINGULARITY_EXEC