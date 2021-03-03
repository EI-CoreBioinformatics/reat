set -euxo
version=0.0.9
rundir=${PWD}
cd $(mktemp -d)
cp ${rundir}/reat_singularity.def reat.def
sudo singularity build reat.img reat.def
mkdir -p /ei/software/testing/reat/${version}/x86_64/bin
cp reat.img /ei/software/testing/reat/${version}/x86_64/reat-${version}.img