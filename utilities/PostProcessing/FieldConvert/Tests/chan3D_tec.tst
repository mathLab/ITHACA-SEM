<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D tecplot output </description>
    <executable>FieldConvert</executable>
    <parameters> chan3D.xml chan3D.fld chan3D.dat</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3D.dat">
                <sha1>c2c71cd8e94a7489a446e81c4a8017e65d1692b8</sha1>
             </file>
         </metric>
    </metrics>
</test>

