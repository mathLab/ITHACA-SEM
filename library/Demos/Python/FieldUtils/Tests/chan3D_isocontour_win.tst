<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Extract a isocontour</description>
    <executable python="true">chan3D_isocontour_win.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
    </files>
    <metrics>
        <metric type="file" id="1">
            <file filename="isocontour.dat">
                <sha1>ed9c242cd1a9797e7dc7cc950500a3f16f554c07</sha1>
             </file>
         </metric>
    </metrics>
</test>

