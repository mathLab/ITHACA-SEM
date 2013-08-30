<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D vtu output </description>
    <executable>FieldConvert</executable>
    <parameters> chan3D.xml chan3D.fld chan3D.vtu</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3D.vtu">
                <sha1>6d7edc7dac35ae0b7626e90564d31c98eaaa1d9d</sha1>
             </file>
         </metric>
    </metrics>
</test>

