<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a .pts file to .dat </description>
    <executable>FieldConvert</executable>
    <parameters>-m interppoints:fromxml=chan3D.xml:fromfld=chan3D.fld chan3D_probe.pts chan3D_probe.dat</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
        <file description="Session File">chan3D_probe.pts</file>
    </files>
    <metrics>
        <metric type="file" id="1">
            <file filename="chan3D_probe.dat">
                <sha1>26a1c31c159b798db40df148d29355e9c378d97b</sha1>
             </file>
         </metric>
    </metrics>
</test>

