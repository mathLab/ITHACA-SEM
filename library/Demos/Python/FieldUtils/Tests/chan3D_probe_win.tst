<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a .pts file to .dat </description>
    <executable python="true">chan3D_probe_win.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
        <file description="Session File">chan3D_probe.pts</file>
    </files>
    <metrics>
        <metric type="file" id="1">
            <file filename="chan3D_probe.dat">
                <sha1>fbc0a3682b31a300334961bfd910364b3a229d5c</sha1>
             </file>
         </metric>
    </metrics>
</test>

