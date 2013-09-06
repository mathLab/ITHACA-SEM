<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D tecplot output with n=10 physical points</description>
    <executable>FieldConvert</executable>
    <parameters> -n 10 chan3D.xml chan3D.fld chan3D.dat</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3D.dat">
                <sha1>65bfbc766df0551d63e243d8dfbaad8ebcd14747</sha1>
             </file>
         </metric>
    </metrics>
</test>

