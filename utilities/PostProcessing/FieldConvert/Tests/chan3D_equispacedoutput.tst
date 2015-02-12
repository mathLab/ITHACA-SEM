<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a field into a equi-spaced tecplot file</description>
    <executable>FieldConvert</executable>
    <parameters>-m equispacedoutput chan3D.xml chan3D.fld equispacedoutput.dat</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
    </files>
    <metrics>
        <metric type="file" id="1">
            <file filename="equispacedoutput.dat">
                <sha1>7f37aca6d15da8a50ef2feb736a47036d4de175b</sha1>
             </file>
         </metric>
    </metrics>
</test>

