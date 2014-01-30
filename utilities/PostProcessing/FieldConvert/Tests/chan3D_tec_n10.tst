<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D tecplot output with n=10 physical points</description>
    <executable>FieldConvert</executable>
    <parameters> -e -n 10 chan3D.xml chan3D.fld chan3D.dat</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3D.dat">
                <sha1>bafe913a59cdd16b5bc298d3379b14bbd6a8d387</sha1>
             </file>
         </metric>
    </metrics>
</test>

