set -euxo
version=0.0.4
rundir=${PWD}
cd $(mktemp -d)
cp ${rundir}/reat_singularity.def reat.def
sudo singularity build reat.img reat.def
mkdir -p /ei/software/testing/reat/${version}/x86_64/bin
cp reat.img /ei/software/testing/reat/${version}/x86_64/reat-${version}.img
cat > singularity.exec <<SINGULARITY_EXEC
#!/bin/bash
DIR=\`dirname \$(readlink -f \$0)\`
singularity exec \$DIR/../reat-${version}.img \$(basename "\$0") \$@
SINGULARITY_EXEC
