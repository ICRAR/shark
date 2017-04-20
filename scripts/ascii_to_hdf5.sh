#!/bin/bash
#
# Ascii to HDF5 converter for descendants.txt file
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2017
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307  USA
#

nhalos=$(head --lines 1 descendants.txt | cut -f1 -d" ")
echo "Converting $nhalos rows"

# Separate each column into its own file
# We 'sed 1d' to ignore the first line
cat descendants.txt | sed 1d | tr -s ' ' | cut -d' ' -f1 > halo_ids
cat descendants.txt | sed 1d | tr -s ' ' | cut -d' ' -f2 > halo_snapshots
cat descendants.txt | sed 1d | tr -s ' ' | cut -d' ' -f3 > descendant_ids
cat descendants.txt | sed 1d | tr -s ' ' | cut -d' ' -f4 > descendant_snapshots

# And now import
rm -f descendants.hdf5
h5import \
        halo_ids             -d $nhalos -p Halo_IDs             -t TEXTIN -s 64 \
		  halo_snapshots       -d $nhalos -p Halo_Snapshots       -t TEXTIN \
		  descendant_ids       -d $nhalos -p Descendant_IDs       -t TEXTIN -s 64 \
		  descendant_snapshots -d $nhalos -p Descendant_Snapshots -t TEXTIN \
		  -o descendants.hdf5

rm halo_ids halo_snapshots descendant_ids descendant_snapshots
