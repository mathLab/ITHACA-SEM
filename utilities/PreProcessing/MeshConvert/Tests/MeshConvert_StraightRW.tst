<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Meshconvert with Spherigons and variable Boundary Layer </description>
    <executable>MeshConvert</executable>
    <parameters> -m spherigon:surf=10:surf=13 -m spherigon:surf=8:surf=9 -m bl:surf=3,10,13:layers=4:r="1.7*( 1-x/0.3 )+1":nq=7 -m bl:surf=2,8,9:layers=4:r="1.7*(1-(x-0.27)/0.078)+1":nq=7 SL_NEK.dat StraightRWGeom.xml </parameters>
    <files>
        <file description="Input File">SL_NEK.dat</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="StraightRWGeom.xml">
                <sha1>72a34f2f0895058700c64f671ed67006bd2cb40f</sha1>
             </file>
         </metric>
    </metrics>
</test>
