#!/bin/bash 
# Roby Joehanes
# Copyright 2007 Roby Joehanes
# This file is distributed under the GNU General Public License version 3.0.
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#DEFAULT_MEM=`free | awk '/Mem/ {print $4}'`
DEFAULT_MEM=`free | awk '/-/ {print $4}'`
DEFAULT_MEM=$(( DEFAULT_MEM * 2 / 1024 / 3 ))
DEFAULT_MEM="-Xmx${DEFAULT_MEM}M"

ALL_JARS="*.jar"
ALL_LIBS=`echo $ALL_JARS | sed 's/.jar\ /.jar\:/g'`
#ALL_LIBS=gpgpu.jar:qgenerics.jar:javacl-1.0.0-RC1-shaded.jar:jdistlib-0.3.5-bin.jar:javacsv.jar:commons-compress-1.3.jar:jgrapht.jar

echo $ALL_LIBS

/usr/java/latest/bin/java -server $DEFAULT_MEM -cp $ALL_LIBS gov.nih.eqtl.QeQTLAnalysis $*
#java -server $DEFAULT_MEM -cp $ALL_LIBS gov.nih.eqtl.QeQTLAnalysis $*

