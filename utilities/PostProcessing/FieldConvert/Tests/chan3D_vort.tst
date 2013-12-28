<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D vorticity output </description>
    <executable>FieldConvert</executable>
    <parameters> -m vorticity chan3D.xml chan3D.fld chan3D_vort.fld</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3D_vort.fld">
                <sha1>4c38e8eebe5433d6d74f7f053a25297352f3503f</sha1>
             </file>
         </metric>
    </metrics>
</test>

