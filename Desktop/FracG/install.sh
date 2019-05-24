#! /bin/bash
clear
echo "FracG is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FracG is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FracG.  If not, see <https://www.gnu.org/licenses/>.
"
read -rsn1 -p"Press any key to continue";echo
clear
echo "***************************
       FracG - V.0.1 
Fault and Fracture Analysis 
***************************"

if `dpkg -s  xutils-dev | grep -q Status;`
   then
      echo "Found XUTILS"
   else
      echo "
      ERROR: XUTILS is not installed!
      Please refer to manual to install it."
      return
fi

if `dpkg -s libboost-dev | grep -q Status;`
   then
      echo "Found BOOST "
   else
      echo "
      ERROR: BOOST is not installed!
      Please refer to manual to install it."
      exit
fi

if `dpkg -s  gdal-bin | grep -q Status;`
   then
      echo "Found GDAL "
   else
      echo "
      ERROR: GDAL is not installed!
      Please refer to manual to install it."
      exit
fi

if `dpkg -s  libarmadillo-dev | grep -q Status;`
   then
      echo "Found ARMADILLO
      "
   else
      echo "
      ERROR: ARMADILLO is not installed
      Please refer to manual to install it."
      exit
fi

xmkmf -a
echo "Excuting make file"
make 
echo "Finished buidling"
chmod -x FracG 
cp FracG /usr/local/bin/
chmod -u+x /usr/local/bin/FracG
echo "Run with FracG"
