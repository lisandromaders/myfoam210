#! /bin/sh

OLD=$1
NEW=$2

pushd $NEW
 rename $OLD $NEW *
 sed -i s/$OLD/$NEW/g *
 sed -i s/$OLD/$NEW/g ./Make/*

 sed -i s/FOAM_LIBBIN/FOAM_USER_LIBBIN/g Make/files
popd

# wclean
# wmake
