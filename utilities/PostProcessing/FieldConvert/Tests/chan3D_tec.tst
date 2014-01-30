<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D tecplot output </description>
    <executable>FieldConvert</executable>
    <parameters> -e chan3D.xml chan3D.fld chan3D.dat</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3D.dat">
                <sha1>5775c16cc4c8f1a2e059170c6fbbbcb4739c8659</sha1>
             </file>
         </metric>
    </metrics>
</test>

