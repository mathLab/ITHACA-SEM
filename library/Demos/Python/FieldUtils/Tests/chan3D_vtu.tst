<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D vtu output </description>
    <executable python="true">chan3D_vtu.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3D.vtu">
                <sha1>ac74dd20883482c3b529e4aa82fbf1291c0d17d0</sha1>
             </file>
         </metric>
    </metrics>
</test>

