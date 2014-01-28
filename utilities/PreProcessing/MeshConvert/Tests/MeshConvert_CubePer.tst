<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Meshconvert with Periodic Boundary condition and Boundary Layer </description>
    <executable>MeshConvert</executable>
    <parameters> -m peralign:dir=y:surf1=3:surf2=5 -m bl:surf=4,6:layers=4:r=3:nq=7 cube.dat cube_nek.xml </parameters>
    <files>
        <file description="Input File">cube.dat</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="cube_nek.xml">
                <sha1>a92cce580691c38d1b7a33f4d1ced3a3ad3f4c2d</sha1>
             </file>
         </metric>
    </metrics>
</test>
