<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D vorticity output </description>
    <executable>FieldConvert</executable>
    <parameters> -m vorticity chan3D.xml chan3D.fld chan3D_vort.fld</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3D_vort.fld">
                <sha1>532fd3d5b4c4b0560936f09da2db1d1063191de3</sha1>
             </file>
         </metric>
    </metrics>
</test>

